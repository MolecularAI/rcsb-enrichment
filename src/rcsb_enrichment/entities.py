"""Polymer entity data: UniProt IDs, sequence, cofactors, peptide ligand detection."""

from .client import RCSBClient, RCSB_DATA

# Protein chains shorter than this with no UniProt mapping are treated as
# peptide ligands.  BIRD-annotated chains are always peptide ligands regardless.
_PEPTIDE_LEN_THRESHOLD = 30


def get_polymer_entities(client: RCSBClient, pdb_id: str, entity_ids: list) -> list:
    entities = []
    for eid in entity_ids:
        data = client.get(f"{RCSB_DATA}/polymer_entity/{pdb_id}/{eid}")
        if not data:
            continue

        entity: dict = {}

        entity_poly = data.get("entity_poly") or {}
        if entity_poly.get("rcsb_entity_polymer_type") != "Protein":
            continue

        container = data.get("rcsb_polymer_entity_container_identifiers") or {}
        entity["uniprot_ids"] = container.get("uniprot_ids") or []
        entity["sequence"] = entity_poly.get("pdbx_seq_one_letter_code_can") or ""
        entity["description"] = (data.get("rcsb_polymer_entity") or {}).get("pdbx_description") or ""

        src_organisms = data.get("rcsb_entity_source_organism") or []
        entity["species"] = src_organisms[0].get("ncbi_scientific_name") if src_organisms else None

        bird_id = container.get("bird_id")
        seq_len = len(entity["sequence"])
        is_peptide = bool(bird_id) or (
            not entity["uniprot_ids"] and 0 < seq_len <= _PEPTIDE_LEN_THRESHOLD
        )
        entity["is_peptide_ligand"] = is_peptide
        entity["bird_id"] = bird_id

        if is_peptide:
            entities.append(entity)
            continue

        cofactors = []
        for c in data.get("rcsb_target_cofactors") or []:
            if c.get("neighbor_flag") == "Y":
                continue
            cofactors.append({
                "name": c.get("cofactor_name"),
                "smiles": c.get("cofactor_SMILES"),
                "inchikey": c.get("cofactor_InChIKey"),
                "chem_comp_id": c.get("cofactor_chem_comp_id"),
                "source": c.get("resource_name", ""),
            })
        entity["cofactors"] = cofactors

        site_types = {"BINDING_SITE", "ACTIVE_SITE", "METAL_COORDINATION", "SITE"}
        site_features = []
        for f in data.get("rcsb_polymer_entity_feature") or []:
            if f.get("type") in site_types:
                site_features.append({
                    "name": f.get("name"),
                    "feature_id": f.get("feature_id"),
                    "source": f.get("provenance_source"),
                    "type": f.get("type"),
                })
        entity["site_features"] = site_features

        entities.append(entity)

    return entities


def extract_direct_binders(client: RCSBClient, pdb_id: str) -> list:
    """Return all direct (neighbor_flag=='N') binders from all protein entities of pdb_id."""
    entry = client.get(f"{RCSB_DATA}/entry/{pdb_id}")
    if not entry:
        return []
    entity_ids = (
        (entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids") or []
    )
    binders = []
    for eid in entity_ids:
        data = client.get(f"{RCSB_DATA}/polymer_entity/{pdb_id}/{eid}")
        if not data:
            continue
        if (data.get("entity_poly") or {}).get("rcsb_entity_polymer_type") != "Protein":
            continue
        for c in data.get("rcsb_target_cofactors") or []:
            if c.get("neighbor_flag") == "Y":
                continue
            binders.append({
                "name": c.get("cofactor_name"),
                "smiles": c.get("cofactor_SMILES"),
                "inchikey": c.get("cofactor_InChIKey"),
                "chem_comp_id": c.get("cofactor_chem_comp_id"),
                "source": c.get("resource_name", ""),
                "from_pdb": pdb_id,
            })
    return binders
