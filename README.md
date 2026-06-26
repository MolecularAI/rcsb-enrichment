# rcsb-enrichment

Enrich a CSV of PDB IDs with structure quality metrics, related entries, known binders, and ligand binding quality scores — all sourced from the RCSB PDB and PDBe REST APIs.

## Features

- **Crystal structure quality** — resolution, R-work/R-free, clashscore, Ramachandran/rotamer outliers %, RSRZ outliers %
- **Related structures** — split into same-region siblings vs full-length proteins via a chain-length threshold (1.4×), searched by UniProt accession or sequence similarity
- **Known binders** — direct pharmacological binders (ChEMBL/DrugBank, `neighbor_flag=N` only) from the query structure, its siblings, and full-length proteins
- **Holo structure lookup** — for each known binder, finds co-crystallised PDB entries via InChIKey→CCD→UniProt+CCD search, even when the query structure is apo or a fragment
- **Ligand binding quality** — per-ligand traffic-light score (`good` / `fair` / `bad`) combining RSCC, RSR, Mogul RMSZ (bonds + angles), intermolecular clashes, and fraction of contact residues (≤4 Å) carrying geometry/density outliers
- **Ligand classification** — splits co-crystallised ligands into drug-like (`ligands_interesting`) and non-interesting (`ligands_noninteresting`: ions, solvents, cofactors, detergents) using the RCSB `is_subject_of_investigation` flag with CCD exclusion list fallback
- **Binding site annotations** — UniProt/CSA site features and PDBe-KB/SIFTS binding site descriptions
- **Excel PDB ID repair** — recovers IDs mangled by Excel's thousands-separator formatting (e.g. `6,000 BHD` → `6BHD`)

## Installation

```bash
git clone https://github.com/<org>/rcsb-enrichment.git
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
    --delay 0.1
```

The `--pdb-col` and `--uniprot-col` arguments are auto-detected from common column name aliases if omitted.

## Output columns

| Column | Description |
|---|---|
| `exp_method` | Experimental method (X-RAY DIFFRACTION, ELECTRON MICROSCOPY, …) |
| `resolution_A` | Diffraction resolution in Å |
| `r_work` / `r_free` | Crystallographic R-factors |
| `clashscore` | MolProbity all-atom clashscore |
| `ramachandran_outliers_pct` | % Ramachandran outliers |
| `rotamer_outliers_pct` | % rotamer outliers |
| `rsrz_outliers_pct` | % RSRZ outliers (X-ray + EDS only) |
| `bonds_rmsz` / `angles_rmsz` | Bond/angle RMSZ from wwPDB validation |
| `ligands_present` | All CCD codes of co-crystallised non-polymer entities |
| `ligands_interesting` | Drug-like ligands (ISI flag / not in exclusion list) |
| `ligands_noninteresting` | Ions, solvents, cofactors, detergents |
| `ligand_quality` | JSON — per-instance metrics for interesting ligands |
| `cofactor_ion_quality` | JSON — per-instance metrics for non-interesting ligands |
| `ligand_binding_quality` | Worst traffic-light across interesting ligands: `good`/`fair`/`bad` |
| `holo_quality` | JSON — per-binder holo structure lookup results |
| `related_pdb_ids` | Same-region sibling PDB entries (same UniProt, chain length ≤ 1.4×) |
| `related_pdb_count` | Count of siblings |
| `fulllength_pdb_ids` | Full-length protein PDB entries (chain length > 1.4× query) |
| `fulllength_pdb_count` | Count of full-length entries |
| `related_search_method` | How related entries were found (`uniprot_split`, `uniprot`, `sequence_id_*`, `none`) |
| `known_binders` | Names of direct pharmacological binders (ChEMBL/DrugBank) |
| `known_binder_smiles` | SMILES of known binders |
| `has_binding_site` | Boolean — any binding site evidence found |
| `binding_site_sources` | Which databases contributed (`UniProt/CSA`, `ChEMBL/DrugBank`, `PDBe-KB/SIFTS`) |
| `binding_site_notes` | Human-readable binding site descriptions |

## TLS / proxy note

On corporate networks with TLS-inspecting proxies, the script merges macOS system keychain certificates with certifi's CA bundle at startup. This is transparent on non-macOS systems.

## API rate limits

RCSB asks for polite access. The default `--delay 0.1` (100 ms between requests) is conservative. Increase it if you see 429 responses.

## Development

```bash
pip install -e ".[dev]"
pytest
```
