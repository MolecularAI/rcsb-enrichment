# CLAUDE.md — rcsb-enrichment

Design decisions, API quirks, and non-obvious constraints for this codebase.
Update this file whenever a decision is made that future contributors (or Claude) would otherwise have to rediscover.

---

## Project overview

CLI tool that reads a CSV of PDB IDs, queries several REST APIs, and writes an enriched CSV.
Package layout:

```
src/rcsb_enrichment/
    client.py          # RCSBClient, _build_ca_bundle, URL constants
    quality.py         # get_entry_quality, get_ligand_quality, traffic_light
    entities.py        # get_polymer_entities, extract_direct_binders
    related.py         # get_related_by_uniprot*, get_related_by_sequence
    holo.py            # get_holo_ligand_quality, _inchikey_to_ccd, _find_holo_entries
    binding_sites.py   # get_pdbe_binding_sites
    ligand_filter.py   # _NON_INTERESTING_CCD, is_interesting_ligand
    enrich.py          # enrich_row, _normalise_pdb_id, _is_valid_pdb_id,
                       #   _entity_matches_names, _collect_entity_names
    cli.py             # main, detect_pdb_col, detect_uniprot_col
```

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

### Polymer entity (`/polymer_entity/{pdb_id}/{entity_id}`)
- `rcsb_polymer_entity.pdbx_description` — human-readable molecule name (e.g. "Tubulin alpha-1B chain"); stored as `entity["description"]` by `get_polymer_entities`

### Nonpolymer entity (`/nonpolymer_entity/{pdb_id}/{entity_id}`)
- `rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id` — CCD code (NOT `comp_id`)
- `rcsb_nonpolymer_entity_container_identifiers.asym_ids` — list of instance asym_ids (NOT `nonpolymer_entity_instance_ids`)
- `rcsb_nonpolymer_entity.pdbx_description` — human-readable ligand name (e.g. "GUANOSINE-5'-TRIPHOSPHATE"); stored as `description` in each ligand metric dict by `get_ligand_quality`

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

### Related entry three-way split: fragment / sibling / full-length
**Decision:** same-UniProt entries are classified into three groups by deposited chain length relative to the query sequence:
- **Fragment** (`fragment_pdb_ids`): chain length < query × 0.8 — a subset/domain of the query.
- **Sibling** (`sibling_pdb_ids`): chain length in [query × 0.8, query × 1.4] — same construct or minor variant.
- **Full-length** (`fulllength_pdb_ids`): chain length > query × 1.4 — query is a domain/fragment of this entry.

**Why:** a query that is itself a fragment (e.g. a domain construct) may have PDB depositions of the isolated domain (sibling), the full protein (full-length), and even shorter sub-domain constructs (fragment). Mixing all three misleads SAR analysis; the three-way split makes the relationship explicit.  
**Thresholds:** upper 1.4× gives ~40% buffer for expression tags and disordered termini. Lower 0.8× gives 20% tolerance for minor N/C-terminal truncations before calling something a subset.  
**Implementation:** three RCSB search requests using `return_type: "polymer_entity"` with `entity_poly.rcsb_sample_sequence_length` `less`/`less_or_equal`/`greater` filters. The sibling set is derived as the `≤upper` result minus the `<lower` fragment set (avoids a fourth request). `return_type` must be `polymer_entity`, not `entry` — with `entry`, the length filter matches any chain in the entry, not necessarily the UniProt-annotated one, causing false fragment/full-length classifications.  
**Binder lookup tier for fragments:** fragments are treated like siblings (unconditional binder inclusion) since they cover a subset of the query domain and their cofactors are directly relevant.

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

### Peptide ligand detection
**Decision:** polymer chains classified as peptide binders are written to a separate `peptide_ligands` column, not mixed with `ligands_interesting`/`ligands_noninteresting` (which cover small-molecule non-polymer entities only).  
**Classification logic (in `get_polymer_entities`):**  
1. BIRD-annotated chains (`rcsb_polymer_entity_container_identifiers.bird_id` is non-null) — always peptide ligands; BIRD is the PDB's formal oligopeptide ligand registry.  
2. Protein-type chains with no UniProt mapping **and** sequence length ≤ `_PEPTIDE_LEN_THRESHOLD` (30 residues) — treated as short uncharacterised peptide binders.  
**Column value:** comma-separated labels — BIRD ID when available, otherwise the one-letter sequence.  
**Downstream effect:** peptide ligand chains are excluded from the receptor entity list; they do not contribute to UniProt resolution, related-entry search, cofactor lookup, or traffic-light scoring.

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
**Decision:** `_build_ca_bundle()` merges macOS system keychain certs with certifi's bundle into a temp PEM file at process start. Each keychain cert is individually probed with `ssl.SSLContext.load_verify_locations` and dropped if OpenSSL rejects it.  
**Why:** corporate TLS-inspecting proxies install a root CA in the system keychain that certifi does not include, causing SSL errors on all HTTPS requests. The per-cert probe is needed because Python 3.14's OpenSSL enforces strict RFC 5280 compliance: legacy macOS certs whose `Basic Constraints` extension is not marked critical are now rejected with `SSLCertVerificationError`. Filtering them out before writing the bundle prevents the error without breaking the proxy CA inclusion. The merge happens once per process; the temp file is not cleaned up (it is small and deterministic).

### RSRZ is not in `pdbx_vrpt_summary_geometry`
**Why recorded:** `percent_RSRZ_outliers` appears to belong in the geometry validation block but is actually in `pdbx_vrpt_summary_diffraction`, a separate block only present for X-ray structures processed by DCC/EDS. This was a live-API-verified correction; do not move it back.

### `iridium_score`: design and known gaps vs. OpenEye Iridium

**Source:** Warren et al., "Essential considerations for using protein–ligand structures in drug discovery", *Drug Discovery Today* 17(23–24), 1270–1281, December 2012.

**What we implemented (`quality.iridium_score`):**  
A weighted-mean composite grade ("good" / "fair" / "bad") over six global structure criteria and one binding-site criterion:

| Criterion | Good | Fair | Bad | Weight |
|-----------|------|------|-----|--------|
| Resolution (Å) | ≤ 2.5 | ≤ 3.0 | > 3.0 | 1 (X-ray only) |
| R-free | ≤ 0.25 | ≤ 0.30 | > 0.30 | 1 (X-ray only) |
| Clashscore | ≤ 10 | ≤ 25 | > 25 | 1 |
| Ramachandran outliers % | ≤ 0.5 | ≤ 2.0 | > 2.0 | 1 |
| Rotamer outliers % | ≤ 1.0 | ≤ 5.0 | > 5.0 | 1 |
| RSRZ outliers % | ≤ 5.0 | ≤ 10.0 | > 10.0 | 1 (X-ray + EDS only) |
| Best ligand `binding_quality` (traffic-light) | good | fair | bad | 2 |

Weighted mean < 0.67 → "good", < 1.33 → "fair", ≥ 1.33 → "bad". Missing metrics are excluded (not penalised); the score degrades gracefully for cryo-EM and NMR entries.

**What the actual Iridium paper does (verified from paper):**  
Iridium classifies X-ray structures into three tiers — **HT** (Highly Trustworthy), **MT** (Mildly Trustworthy), **NT** (Not Trustworthy) — using a sequential decision tree, not a weighted mean. The tree (Chart 1 in the paper) filters in this order:

1. **Electron density available** (public EDS data) — structures without it cannot enter the pipeline at all. This alone cut 728 → 233 structures (68% rejected).
2. **R-factor < 0.40** — structures at ≥ 0.40 go to NT; at this value a model is statistically indistinguishable from random.
3. **Density for ligands** — no density → removed; partial density (RSR 0.10–0.30, RSCC 0.50–0.90) → MT; complete density (RSR < 0.10, RSCC > 0.90) required for HT.
4. **Non-covalent ligands only** — covalent ligands → MT (pose/affinity prediction requires knowledge of the covalent bond, which most docking algorithms lack).
5. **No alternate conformations for active-site side chains** (within 5 Å of any ligand atom) → MT if present.
6. **No alternate conformations for ligand atoms** → MT if present.
7. **No crystal packing contacts** (no atom from a neighbouring protein in the lattice within 6 Å of any ligand atom) — electrostatic surface residues can stabilise ligand conformations not present in solution → MT if violated.
8. **Complete density in active site** — all active-site residue and co-factor atoms have density.
9. **Occupancy = 1.0 in active site** — ligand and active-site atoms must be fully occupied → MT if < 1.0.

Only structures passing all nine steps reach the HT set (121 of 728, 17%).

**Key corrections to the prior CLAUDE.md entry (which was recalled from training data):**

- **Tier names:** HT / MT / **NT** (Not Trustworthy), not "HT / MT / LT".
- **Clashscore and Ramachandran outlier % are NOT Iridium criteria.** They are wwPDB/MolProbity validation metrics. The original Iridium paper does not use them. We added them to our score because they are available from RCSB and are reasonable proxies, but they are our addition, not Iridium's.
- **B-factor is NOT a filter criterion in the Iridium decision tree.** The paper discusses the occupancy-weighted B-factor (owB-factor) as a useful local quality indicator alongside RSR and RSCC, but it does not appear as a step in Chart 1. It was incorrect to call it "the single biggest gap".
- **RSCC thresholds in the paper:** complete density = RSCC > 0.90 / RSR < 0.10; partial density = RSCC 0.50–0.90 / RSR 0.10–0.30; no density = RSCC < 0.50 / RSR > 0.30. Our traffic-light uses 0.80/0.60 for RSCC, which is more lenient than the paper's 0.90 threshold for "complete density".
- **R-factor (not R-free) < 0.40** is the NT hard cutoff. The paper also recommends R-free < 0.45 and R − R-free difference ≤ 0.05 (Table 1). We currently use R-free ≤ 0.25/0.30 as good/fair thresholds — a stricter criterion on R-free but we do not check R − R-free divergence.

**Gaps in our implementation vs. the actual Iridium paper (in priority order):**

1. **Ligand occupancy < 1.0** — a hard HT criterion; structures with any active-site atom at < 1.0 occupancy go to MT. Available from RCSB: `rcsb_nonpolymer_instance_validation_score[0].average_occupancy`. Threshold: < 1.0 → downgrade.
2. **Crystal packing contacts within 6 Å of ligand** — not tracked. Would require fetching the full unit cell or symmetry-mate coordinates; not directly available via the RCSB REST API.
3. **Alternate conformations for ligand and active-site side chains** — not tracked. Would require inspecting `pdbx_struct_conf` or the full PDB/mmCIF file.
4. **R − R-free divergence > 0.05** (overfitting indicator) — R-work is available from RCSB; we could add the difference check.
5. **RSCC "good" threshold should arguably be 0.90, not 0.80** — the paper equates 0.90 with "complete density". Our more lenient 0.80 threshold is closer to RCSB's own "acceptable" flag.

**What our score does that Iridium does not:**  
- Accepts cryo-EM and NMR (Iridium is X-ray only; we skip X-ray-only sub-criteria gracefully).
- Penalises intermolecular clashes and contact-residue outlier fraction via the `binding_quality` traffic-light (weight 2 in our composite).
- Uses Clashscore and Ramachandran / rotamer outlier % from wwPDB validation as protein geometry proxies.

### Multi-UniProt related-entry search

**Decision:** when a structure has multiple receptor chains (e.g. a heterodimer), `get_related_by_uniprot_split` / `get_related_by_uniprot` is called once per resolved UniProt ID, and results are merged with deduplication.  
**Why:** the original code only searched `all_uniprot[0]`. For a heteromeric complex, this silently skipped all chains beyond the first, returning an incomplete set of related structures.  
**Implementation:** a `seen` set tracks already-added PDB IDs within the sibling and full-length lists to prevent duplicates from overlapping search results.

### Entity name column (`entity_names`)

**Decision:** a comma-separated `entity_names` column lists the `pdbx_description` of every polymer and non-polymer entity for that structure. Polymer descriptions appear first, then non-polymer. Duplicates (e.g. homodimer chains with the same description) are deduplicated. The column is placed immediately after the PDB ID column in the output CSV.  
**Scope:** populated for both primary rows (from `enrich_row`) and related-entry ligand sub-rows (from `_fetch_related_ligand_data`). The column is included in `tags` passed to `build_ligand_rows`, so it propagates onto sibling and full-length sub-rows.  
**API source:** `rcsb_polymer_entity.pdbx_description` (polymer) and `rcsb_nonpolymer_entity.pdbx_description` (non-polymer).

### Entity name filtering (`--entity-names`)

**Decision:** an optional `--entity-names` CLI argument accepts comma-separated substrings. After `get_polymer_entities` returns, receptor entities are filtered to those whose `pdbx_description` contains at least one term as a whole whitespace-delimited token (case-insensitive). Filtered-out entities do not contribute to UniProt resolution, related-entry search, cofactor lookup, or binding-site scoring.  
**Matching rule:** split description on whitespace, check if any filter term equals a token. Hyphenated compounds are single tokens: `"Tubulin"` matches `"Tubulin alpha-1B chain"` but NOT `"Tubulin-Tyrosine Ligase"`.  
**Why token-based not substring:** substring matching would make `"Tubulin"` spuriously match `"Tubulin-Tyrosine Ligase"`, returning four entities for 5S5V instead of the expected two.

### Cross-row deduplication cache

**Decision:** `cli.main()` creates two dicts — `related_data_cache` and `binders_cache` — that are shared across all `enrich_row` calls in a run. They are keyed by PDB ID and cache the results of `_fetch_related_ligand_data` and `extract_direct_binders` respectively. When multiple input PDB IDs resolve to the same set of related structures (e.g. 5S5V and 4X1Y are both tubulin structures), each related entry is fetched only once.  
**Why dicts not sets:** the cache must store the full fetched data, not just presence; `_fetch_related_ligand_data` is expensive (multiple API calls per related entry).  
**Mutation safety:** `extract_direct_binders` results are cached as returned; callers must `dict(b)` before adding `binder_source_type` to avoid mutating the cached original. This is enforced at the call sites in `enrich_row`.  
**New output columns:** `all_related_pdb_ids` and `all_fulllength_pdb_ids` on primary rows list the complete search-result sets (all found siblings/full-lengths, after self and other-input exclusion), independently of whether those entries have meaningful ligands.

### UniProt search result cache (`uniprot_search_cache`)

**Decision:** `cli.main()` creates a `uniprot_search_cache` dict shared across all `enrich_row` calls. It is keyed by `(uniprot_id, seq_len, max_rows)` and stores the `(siblings, full_length)` tuple returned by `get_related_by_uniprot_split` (or the flat list from `get_related_by_uniprot` when seq_len is 0).  
**Why:** when multiple input PDB IDs resolve to the same set of UniProt IDs (e.g. two tubulin structures both filtered to P81947 + Q6B856), each search fires 2 × 2 = 4 API calls. The cache collapses these to 4 calls total for the entire run regardless of how many input rows share those UniProt IDs.  
**Cache key includes `max_rows`** so a future call with a different `--max-related` value gets a fresh search rather than a truncated cached result.

---

