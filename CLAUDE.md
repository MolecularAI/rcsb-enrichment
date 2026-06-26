# CLAUDE.md — rcsb-enrichment

Design decisions, API quirks, and non-obvious constraints for this codebase.
Update this file whenever a decision is made that future contributors (or Claude) would otherwise have to rediscover.

---

## Project overview

Single-file CLI script (`src/rcsb_enrichment/rcsb_enrichment.py`) that reads a CSV of PDB IDs, queries several REST APIs, and writes an enriched CSV. The script is intentionally monolithic for now; refactoring into modules is planned.

---

## API endpoints used

### RCSB Data API (`https://data.rcsb.org/rest/v1/core`)
- `/entry/{pdb_id}` — entry-level metadata, quality summary blocks, entity IDs
- `/polymer_entity/{pdb_id}/{entity_id}` — sequence, UniProt IDs, target cofactors
- `/polymer_entity_instance/{pdb_id}/{asym_id}` — per-chain validation features (RSRZ, outlier flags)
- `/nonpolymer_entity/{pdb_id}/{entity_id}` — ligand entity metadata, instance IDs
- `/nonpolymer_entity_instance/{pdb_id}/{asym_id}` — per-ligand validation scores, contact neighbours

### RCSB Search API (`https://search.rcsb.org/rcsbsearch/v2/query`)
- Used for UniProt-based related entry search, sequence similarity search, InChIKey→CCD lookup, and holo structure search.

### PDBe API (`https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites`)
- SIFTS-curated binding site descriptions.

---

## Key API field names (verified against live data)

These were confirmed by live API inspection; do not change without re-verifying.

### Entry-level quality
- `pdbx_vrpt_summary_geometry[0]` — clashscore, percent_ramachandran_outliers, percent_rotamer_outliers, bonds_RMSZ, angles_RMSZ
- `pdbx_vrpt_summary_diffraction[0].percent_RSRZ_outliers` — **NOT** in geometry block; only present for X-ray structures with EDS processing
- `rcsb_entry_container_identifiers.non_polymer_entity_ids` — nonpolymer entity IDs (note: `non_polymer_entity_ids`, not `nonpolymer_entity_ids`)

### Ligand validation
- `rcsb_nonpolymer_instance_validation_score[0].RSCC` — uppercase; `.rscc` returns None
- `rcsb_nonpolymer_instance_validation_score[0].RSR` — uppercase
- `rcsb_nonpolymer_instance_validation_score[0].mogul_bonds_RMSZ` — bond geometry Z-score
- `rcsb_nonpolymer_instance_validation_score[0].mogul_angles_RMSZ` — angle geometry Z-score
- `rcsb_nonpolymer_instance_validation_score[0].intermolecular_clashes` — integer count; field `intermolecular_clashscore` does not exist
- `rcsb_nonpolymer_instance_validation_score[0].is_subject_of_investigation` — True/False/None

### Contact neighbours
- `rcsb_target_neighbors` on nonpolymer_entity_instance — array of {atom_id, comp_id, distance, target_asym_id, target_atom_id, target_auth_seq_id, target_comp_id, target_entity_id, target_model_id, target_seq_id}
- Uses internal asym_id, not auth chain ID, for the polymer instance endpoint

### Per-residue validation
- `rcsb_polymer_instance_feature` — array of feature objects; two encoding layouts:
  - **Dense**: RSRZ, RSR, RSCC, OWAB, Q_SCORE, ASA — `values[]` is a continuous array, `beg_seq_id` is the start residue
  - **Sparse**: ANGLE_OUTLIERS, BOND_OUTLIERS, CLASHES, RAMACHANDRAN_OUTLIER, ROTAMER_OUTLIER, etc. — one element per outlier residue, `values[0]` = count
- `pdbx_vrpt_summary_entity_geometry` — per-chain aggregate only (bonds_RMSZ, angles_RMSZ), not per-residue

### Search API
- `rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id` — correct attribute for ligand CCD code in search; `rcsb_nonpolymer_entity.comp_id` does not exist
- InChIKey→CCD: service=`text_chem`, attribute=`rcsb_chem_comp_descriptor.InChIKey`, return_type=`mol_definition`; value must be bare InChIKey (no `InChIKey=` prefix); `rcsb_chem_comp_descriptor.descriptor` is not searchable
- `rcsb_target_cofactors[].neighbor_flag`: `"Y"` = binder targets a neighbouring protein, NOT this one; `"N"` = direct binder. Only `"N"` entries are used.

---

## Design decisions

### Related entry split: siblings vs full-length (1.4× chain length threshold)
**Decision:** entries whose deposited chain length exceeds the query sequence length × 1.4 are classified as "full-length" (`fulllength_pdb_ids`); shorter entries are "siblings" (`related_pdb_ids`).  
**Why:** a common PDB pattern is a fragment crystallised in isolation plus the same protein as a subunit in a larger complex. Mixing them misleads SAR analysis. The 1.4× factor gives a ~40% buffer for expression tags, disordered termini, and minor construct differences.  
**Implementation:** two RCSB search requests with `entity_poly.rcsb_sample_sequence_length` `less_or_equal`/`greater` filters — no per-hit extra fetches.  
**Rule when only one class populated:** use `related_pdb_ids` (column 1).

### Binder lookup: three tiers
**Decision:** known binders are collected from (1) the query structure's own entities, (2) sibling structures unconditionally, (3) full-length structures only if the binder's `cofactor_chem_comp_id` matches a ligand co-crystallised in the query fragment.  
**Why for tier 3:** `rcsb_target_cofactors` has no binding-residue fields, and `reference_sequence_identifiers` gives only fractional coverage (no start/end residue positions). The only available evidence that a binder contacts the fragment domain is if its CCD code appears in the fragment's `ligands_present`.  
**Deduplication:** by InChIKey, falling back to name. Capped at `_MAX_RELATED_BINDER_ENTRIES = 5` sibling/full-length lookups each.

### neighbor_flag filtering
**Decision:** only `neighbor_flag == "N"` cofactors are kept (direct binders to this protein).  
**Why:** `neighbor_flag == "Y"` means the molecule is annotated against a *neighbouring* protein in DrugBank's interaction network — not this target. Including them overstates binding evidence.

### Holo structure lookup for apo/fragment inputs
**Decision:** if the query structure has no co-crystallised ligand (most fragment/apo structures), `ligand_binding_quality` is derived from holo structures found by searching UniProt + CCD code.  
**Why:** the traffic-light score was always empty for apo inputs before this change.  
**Flow:** known binder → CCD (from `cofactor_chem_comp_id` or InChIKey→CCD search) → UniProt+CCD holo search → `get_ligand_quality()` on holo entries.  
**Limit:** `_MAX_HOLO_ENTRIES = 3` holo structures per binder.

### Ligand interest classification
**Decision:** two output columns — `ligands_interesting` (drug-like) and `ligands_noninteresting` (ions/solvents/cofactors/detergents).  
**Why:** co-crystallised non-polymer entities include many biologically relevant but not drug-targeted molecules (Mg²⁺, ATP, HEPES, PEG) that would corrupt the traffic-light score if mixed with actual drug candidates.  
**Primary signal:** RCSB `is_subject_of_investigation` flag from `rcsb_nonpolymer_instance_validation_score`. Three-way logic: `True` → interesting, `False` → not interesting, `None` (absent in older entries) → fall back to `_NON_INTERESTING_CCD` exclusion set (~100 CCD codes).  
**Traffic-light scope:** `ligand_binding_quality` aggregates only from interesting ligands.

### Traffic-light scoring
**Decision:** six metrics each scored 0/1/2 (good/fair/bad), mean threshold < 0.67 → "good", < 1.33 → "fair", ≥ 1.33 → "bad". Hard override: RSCC < 0.50 → "bad" regardless.  
**Thresholds rationale:**
- RSCC: 0.80/0.60 — standard crystallographic acceptability thresholds (RCSB uses 0.80 as "good")
- RSR: 0.20/0.35 — standard RCSB thresholds; RSR > 0.35 is flagged as poor in wwPDB validation reports
- RMSZ bonds/angles: 1.5/2.5 — RCSB uses 2.0 as the outlier threshold; 1.5 is conservative
- Intermolecular clashes: 0/2 — any clash is a concern; ≥3 is clearly problematic
- Contact outlier fraction: 0.10/0.25 — >10% of contact residues with outliers is notable; >25% is a quality concern
**Missing metrics** (e.g. RSCC absent for cryo-EM or old entries without EDS) are excluded from the mean rather than treated as zero.

### Contact residue cutoff: 4.0 Å
**Decision:** residues with any atom within 4.0 Å of any ligand atom are counted as contact residues.  
**Why:** 4 Å captures hydrogen-bond donors/acceptors and van der Waals contacts; tighter (3.5 Å) misses many relevant interactions, wider (5 Å) includes too many peripheral residues that may have outliers unrelated to binding.

### Excel PDB ID repair
**Decision:** `_normalise_pdb_id()` strips commas and spaces before taking the first 4 characters.  
**Why:** Excel reformats PDB IDs beginning with a digit (e.g. `6BHD`) using thousands separators, producing `6,000 BHD`. Removing commas and spaces recovers the original ID. A warning is logged whenever normalisation changes the raw value.

### TLS / CA bundle (macOS only)
**Decision:** `_build_ca_bundle()` merges macOS system keychain certs with certifi's bundle into a temp PEM file at process start.  
**Why:** corporate TLS-inspecting proxies install a root CA in the system keychain that certifi does not include, causing SSL errors on all HTTPS requests. The merge happens once per process; the temp file is not cleaned up (it is small and deterministic).

### RSRZ is not in `pdbx_vrpt_summary_geometry`
**Why recorded:** `percent_RSRZ_outliers` appears to belong in the geometry validation block but is actually in `pdbx_vrpt_summary_diffraction`, a separate block only present for X-ray structures processed by DCC/EDS. This was a live-API-verified correction; do not move it back.

---

## Planned refactoring

The script is intentionally kept as a single file during the exploratory phase. When refactoring into modules, the natural split is:

```
src/rcsb_enrichment/
    client.py          # RCSBClient, _build_ca_bundle
    quality.py         # get_entry_quality, get_ligand_quality, traffic_light
    entities.py        # get_polymer_entities, _extract_direct_binders
    related.py         # get_related_by_uniprot*, get_related_by_sequence
    holo.py            # get_holo_ligand_quality, _inchikey_to_ccd, _find_holo_entries
    binding_sites.py   # get_pdbe_binding_sites
    ligand_filter.py   # _NON_INTERESTING_CCD, _is_interesting_ligand
    enrich.py          # enrich_row (orchestration)
    cli.py             # main, detect_pdb_col, detect_uniprot_col
```
