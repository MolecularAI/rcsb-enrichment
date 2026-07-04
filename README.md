# rcsb-enrichment

Enrich a CSV of PDB IDs with structure quality metrics, related entries, known binders, and ligand binding quality scores — all sourced from the RCSB PDB and PDBe REST APIs.

## Features

- **Structure quality score** — composite Iridium-like grade (`good` / `fair` / `bad`) combining resolution, R-free, clashscore, Ramachandran/rotamer outliers, RSRZ outliers %, and ligand binding quality
- **Crystal structure quality** — resolution, R-work/R-free, clashscore, Ramachandran/rotamer outliers %, RSRZ outliers %
- **Entity names** — comma-separated list of all polymer and non-polymer entity descriptions per structure, placed immediately after the PDB ID column
- **Related structures** — split into same-region siblings vs full-length proteins via a chain-length threshold (1.4×), searched by UniProt accession or sequence similarity; all resolved UniProt IDs are searched when a structure has multiple receptor chains
- **Known binders** — direct pharmacological binders (ChEMBL/DrugBank, `neighbor_flag=N` only) from the query structure, its siblings, and full-length proteins
- **Holo structure lookup** — for each known binder, finds co-crystallised PDB entries via InChIKey→CCD→UniProt+CCD search, even when the query structure is apo or a fragment
- **Ligand binding quality** — per-ligand traffic-light score (`good` / `fair` / `bad`) combining RSCC, RSR, Mogul RMSZ (bonds + angles), intermolecular clashes, and fraction of contact residues (≤4 Å) carrying geometry/density outliers
- **Ligand classification** — splits co-crystallised ligands into drug-like (`ligands_interesting`) and non-interesting (`ligands_noninteresting`: ions, solvents, cofactors, detergents) using the RCSB `is_subject_of_investigation` flag with CCD exclusion list fallback
- **Binding site annotations** — UniProt/CSA site features and PDBe-KB/SIFTS binding site descriptions
- **Entity name filtering** — restrict processing to specific receptor chains by molecule name substring (`--entity-names`)
- **Cross-row deduplication** — related structures that appear in search results for multiple input PDB IDs are fetched only once per run
- **Excel PDB ID repair** — recovers IDs mangled by Excel's thousands-separator formatting (e.g. `6,000 BHD` → `6BHD`)

## Installation

```bash
git clone https://github.com/MolecularAI/rcsb-enrichment.git
cd rcsb-enrichment
pip install -e .
```

**Requirements:** Python ≥ 3.9, `requests`, `certifi`, `pandas`

## Usage

```bash
rcsb-enrich --input proteins.csv --output enriched.csv

# Full options
rcsb-enrich \
    --input proteins.csv \
    --output enriched.csv \
    --pdb-col PDBID \
    --uniprot-col Uniprot \
    --seq-identity 0.9 \
    --max-related 25 \
    --delay 0.1 \
    --entity-names "Tubulin,Kinase"
```

The `--pdb-col` and `--uniprot-col` arguments are auto-detected from common column name aliases if omitted.

`--entity-names` accepts a comma-separated list of molecule name substrings. When set, only receptor chains whose `pdbx_description` contains at least one term as a whole whitespace-delimited word are processed; other chains are ignored. This is useful for heteromeric complexes where only specific subunits are of interest. For example, `--entity-names Tubulin` applied to PDB entry 5S5V (which contains alpha-tubulin, beta-tubulin, Stathmin-4, and Tubulin-Tyrosine Ligase) retains only the two tubulin chains. Note that hyphenated compound words are treated as single tokens: `Tubulin` does **not** match `Tubulin-Tyrosine`.

## Output columns

### Row metadata

| Column | Description |
|---|---|
| `row_type` | `primary` for the input PDB entry; `ligand` for per-ligand sub-rows from related structures |
| `parent_pdb_id` | PDB ID of the input structure (set on ligand sub-rows) |
| `entity_names` | Comma-separated descriptions of all polymer and non-polymer entities (placed immediately after the PDB ID column) |
| `related_pdb_ids` | PDB ID of the sibling structure this ligand sub-row comes from |
| `fulllength_pdb_ids` | PDB ID of the full-length structure this ligand sub-row comes from |

### Structure quality

| Column | Description |
|---|---|
| `species` | Source organism(s), comma-separated |
| `structure_quality` | Composite grade: `good` / `fair` / `bad` (Iridium-like) |
| `exp_method` | Experimental method (X-RAY DIFFRACTION, ELECTRON MICROSCOPY, …) |
| `resolution_A` | Diffraction resolution in Å |
| `r_work` / `r_free` | Crystallographic R-factors |
| `clashscore` | MolProbity all-atom clashscore |
| `ramachandran_outliers_pct` | % Ramachandran outliers |
| `rotamer_outliers_pct` | % rotamer outliers |
| `rsrz_outliers_pct` | % RSRZ outliers (X-ray + EDS only) |
| `bonds_rmsz` / `angles_rmsz` | Bond/angle RMSZ from wwPDB validation |

### Per-ligand detail (ligand sub-rows)

| Column | Description |
|---|---|
| `ligand_type` | `small_molecule` or `peptide` |
| `ligand_id` | CCD code or BIRD ID / one-letter sequence (peptides) |
| `chain_id` | Asymmetric unit chain ID |
| `ligand_rscc` | Real-space correlation coefficient |
| `ligand_rsr` | Real-space R-value |
| `ligand_rmsz_bonds` / `ligand_rmsz_angles` | Mogul bond/angle geometry Z-scores |
| `ligand_intermolecular_clashes` | Count of intermolecular clashes |
| `contact_residue_count` | Number of protein residues within 4 Å |
| `contact_outlier_fraction` | Fraction of contact residues carrying validation outliers |
| `contact_residues` | Semicolon-separated contact residue labels |
| `binding_quality` | Per-ligand traffic-light: `good` / `fair` / `bad` |

### Ligands summary (primary rows)

| Column | Description |
|---|---|
| `ligands_present` | All CCD codes of co-crystallised non-polymer entities |
| `ligands_interesting` | Drug-like ligands (ISI flag / not in exclusion list) |
| `ligands_noninteresting` | Ions, solvents, cofactors, detergents |
| `ligand_quality` | JSON — per-instance metrics for interesting ligands |
| `holo_quality` | JSON — ligand quality from holo structures (populated for apo/fragment inputs) |

### Related entries

| Column | Description |
|---|---|
| `all_related_pdb_ids` | All sibling PDB IDs found by search (up to `--max-related`), regardless of ligand content |
| `all_fulllength_pdb_ids` | All full-length PDB IDs found by search (up to `--max-related`), regardless of ligand content |
| `related_pdb_ids_no_ligand` | Sibling PDB IDs with no meaningful ligand |
| `related_pdb_ids_no_ligand_count` | Count of the above |
| `fulllength_pdb_ids_no_ligand` | Full-length PDB IDs with no meaningful ligand |
| `fulllength_pdb_ids_no_ligand_count` | Count of the above |
| `related_search_method` | `uniprot_split`, `uniprot`, `sequence_id_<threshold>`, or `none` |

The `related_pdb_ids` and `fulllength_pdb_ids` tag columns (one PDB ID per ligand sub-row) identify which related structure a ligand sub-row's data came from.

### Binding sites and binders

| Column | Description |
|---|---|
| `binding_site_sources` | Contributing databases: `UniProt/CSA`, `ChEMBL/DrugBank`, `PDBe-KB/SIFTS` |
| `binding_site_notes` | Human-readable binding site descriptions |
| `known_binders` | Names of direct pharmacological binders |
| `known_binder_smiles` | SMILES of known binders (same order) |

## TLS / proxy note

On corporate networks with TLS-inspecting proxies, the script merges macOS system keychain certificates with certifi's CA bundle at startup. This is transparent on non-macOS systems.

## API rate limits

RCSB asks for polite access. The default `--delay 0.1` (100 ms between requests) is conservative. Increase it if you see 429 responses.

## Development

```bash
pip install -e .
pytest
```
