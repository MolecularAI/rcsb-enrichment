"""Entry-level and per-ligand quality metrics, plus traffic-light scoring."""

from .client import RCSBClient, RCSB_DATA
from .ligand_filter import is_interesting_ligand

CONTACT_DIST_CUTOFF = 4.0  # Å — ligand–residue contact threshold

# Sparse feature types: one element per outlier residue, values[0] = count
_SPARSE_FEATURE_TYPES = {
    "ANGLE_OUTLIERS", "BOND_OUTLIERS", "CLASHES", "SYMM_CLASHES",
    "RAMACHANDRAN_OUTLIER", "ROTAMER_OUTLIER", "CHIRAL_OUTLIERS",
    "PLANE_OUTLIERS", "STEREO_OUTLIER", "CIS-PEPTIDE",
    "MOGUL_ANGLE_OUTLIERS", "MOGUL_BOND_OUTLIERS",
}

# Feature types where values[] is a dense per-residue array starting at beg_seq_id
_DENSE_FEATURE_TYPES = {"RSRZ", "RSR", "RSCC", "OWAB", "Q_SCORE", "ASA"}


def get_entry_quality(client: RCSBClient, pdb_id: str) -> dict:
    data = client.get(f"{RCSB_DATA}/entry/{pdb_id}")
    if not data:
        return {}

    out = {}

    exptl = data.get("exptl") or []
    out["exp_method"] = exptl[0].get("method") if exptl else None

    entry_info = data.get("rcsb_entry_info") or {}
    dres = entry_info.get("diffrn_resolution_high") or {}
    out["resolution_A"] = dres.get("value")

    refine = data.get("refine") or []
    if refine:
        out["r_work"] = refine[0].get("ls_R_factor_R_work")
        out["r_free"] = refine[0].get("ls_R_factor_R_free")
    else:
        out["r_work"] = None
        out["r_free"] = None

    geom = data.get("pdbx_vrpt_summary_geometry") or []
    if geom:
        g = geom[0]
        out["clashscore"] = g.get("clashscore")
        out["ramachandran_outliers_pct"] = g.get("percent_ramachandran_outliers")
        out["rotamer_outliers_pct"] = g.get("percent_rotamer_outliers")
        out["bonds_rmsz"] = g.get("bonds_RMSZ")
        out["angles_rmsz"] = g.get("angles_RMSZ")
    else:
        out["clashscore"] = None
        out["ramachandran_outliers_pct"] = None
        out["rotamer_outliers_pct"] = None
        out["bonds_rmsz"] = None
        out["angles_rmsz"] = None

    # percent_RSRZ_outliers lives in the diffraction block, not geometry (X-ray + EDS only)
    diff = data.get("pdbx_vrpt_summary_diffraction") or []
    out["rsrz_outliers_pct"] = diff[0].get("percent_RSRZ_outliers") if diff else None

    ligands = entry_info.get("nonpolymer_bound_components") or []
    out["ligands_present"] = ",".join(ligands) if ligands else ""

    identifiers = data.get("rcsb_entry_container_identifiers") or {}
    out["_polymer_entity_ids"] = identifiers.get("polymer_entity_ids") or []
    out["_nonpolymer_entity_ids"] = identifiers.get("non_polymer_entity_ids") or []

    return out


def _parse_residue_features(features: list) -> dict:
    """Return {seq_id: {feature_type: value}} from rcsb_polymer_instance_feature."""
    residue_data: dict = {}
    for feat in features:
        ftype = feat.get("type", "")
        for pos in feat.get("feature_positions") or []:
            values = pos.get("values") or []
            beg = pos.get("beg_seq_id")
            if beg is None or not values:
                continue
            if ftype in _DENSE_FEATURE_TYPES:
                for i, val in enumerate(values):
                    rid = beg + i
                    residue_data.setdefault(rid, {})[ftype] = val
            elif ftype in _SPARSE_FEATURE_TYPES:
                residue_data.setdefault(beg, {})[ftype] = values[0]
    return residue_data


def traffic_light(
    rscc,
    rsr,
    rmsz_bonds,
    rmsz_angles,
    intermolecular_clashes,
    contact_outlier_fraction: float,
) -> str:
    """Aggregate per-ligand and per-contact-residue metrics into a traffic-light score.

    Scoring table (each criterion contributes 0=good, 1=fair, 2=bad):
      RSCC       ≥0.80→0   ≥0.60→1   <0.60→2   (hard override: <0.50 → "bad")
      RSR        ≤0.20→0   ≤0.35→1   >0.35→2
      RMSZ bonds ≤1.5→0    ≤2.5→1    >2.5→2
      RMSZ angles≤1.5→0    ≤2.5→1    >2.5→2
      i-mol clashes 0→0    1-2→1     ≥3→2
      contact outlier fraction ≤0.10→0 ≤0.25→1 >0.25→2

    Mean score: <0.67→"good"  <1.33→"fair"  ≥1.33→"bad"
    Missing metrics are excluded from the mean rather than treated as zero.
    """
    if rscc is not None and rscc < 0.50:
        return "bad"

    scores = []
    if rscc is not None:
        scores.append(0 if rscc >= 0.80 else (1 if rscc >= 0.60 else 2))
    if rsr is not None:
        scores.append(0 if rsr <= 0.20 else (1 if rsr <= 0.35 else 2))
    if rmsz_bonds is not None:
        scores.append(0 if rmsz_bonds <= 1.5 else (1 if rmsz_bonds <= 2.5 else 2))
    if rmsz_angles is not None:
        scores.append(0 if rmsz_angles <= 1.5 else (1 if rmsz_angles <= 2.5 else 2))
    if intermolecular_clashes is not None:
        scores.append(0 if intermolecular_clashes == 0 else (1 if intermolecular_clashes <= 2 else 2))
    scores.append(0 if contact_outlier_fraction <= 0.10 else (1 if contact_outlier_fraction <= 0.25 else 2))

    if not scores:
        return ""
    mean = sum(scores) / len(scores)
    return "good" if mean < 0.67 else ("fair" if mean < 1.33 else "bad")


def iridium_score(quality: dict, best_ligand_traffic: str | None) -> str:
    """Composite structure quality grade analogous to OpenEye Iridium (HT/MT/LT).

    Combines global crystallographic metrics (resolution, R-free, clashscore,
    Ramachandran/rotamer outliers, RSRZ) with ligand binding-site quality.
    Returns "good", "fair", or "bad".

    Scoring table (0=good, 1=fair, 2=bad):
      resolution_A          ≤2.5→0  ≤3.0→1  >3.0→2   (X-ray only)
      r_free                ≤0.25→0 ≤0.30→1 >0.30→2  (X-ray only)
      clashscore            ≤10→0   ≤25→1   >25→2
      ramachandran_outliers ≤0.5%→0 ≤2.0%→1 >2.0%→2
      rotamer_outliers      ≤1.0%→0 ≤5.0%→1 >5.0%→2
      rsrz_outliers         ≤5%→0   ≤10%→1  >10%→2   (X-ray + EDS only)
      best ligand traffic   good→0  fair→1  bad→2    (weight 2×)

    Missing metrics are excluded from the weighted mean.
    Thresholds: mean <0.67→"good", <1.33→"fair", ≥1.33→"bad".
    """
    exp = (quality.get("exp_method") or "").upper()
    is_xray = "X-RAY" in exp

    weighted_scores: list[tuple[float, float]] = []  # (score, weight)

    if is_xray:
        res = quality.get("resolution_A")
        if res is not None:
            weighted_scores.append((0 if res <= 2.5 else (1 if res <= 3.0 else 2), 1))
        rfree = quality.get("r_free")
        if rfree is not None:
            weighted_scores.append((0 if rfree <= 0.25 else (1 if rfree <= 0.30 else 2), 1))

    clash = quality.get("clashscore")
    if clash is not None:
        weighted_scores.append((0 if clash <= 10 else (1 if clash <= 25 else 2), 1))

    rama = quality.get("ramachandran_outliers_pct")
    if rama is not None:
        weighted_scores.append((0 if rama <= 0.5 else (1 if rama <= 2.0 else 2), 1))

    rota = quality.get("rotamer_outliers_pct")
    if rota is not None:
        weighted_scores.append((0 if rota <= 1.0 else (1 if rota <= 5.0 else 2), 1))

    rsrz = quality.get("rsrz_outliers_pct")
    if rsrz is not None:
        weighted_scores.append((0 if rsrz <= 5.0 else (1 if rsrz <= 10.0 else 2), 1))

    if best_ligand_traffic in ("good", "fair", "bad"):
        lig_score = {"good": 0, "fair": 1, "bad": 2}[best_ligand_traffic]
        weighted_scores.append((lig_score, 2))

    if not weighted_scores:
        return ""
    total_weight = sum(w for _, w in weighted_scores)
    mean = sum(s * w for s, w in weighted_scores) / total_weight
    return "good" if mean < 0.67 else ("fair" if mean < 1.33 else "bad")


def _get_contact_residue_outlier_fraction(
    client: RCSBClient,
    pdb_id: str,
    neighbors: list,
    cache: dict,
) -> float:
    """Fraction of contact residues (≤CONTACT_DIST_CUTOFF Å) that carry any validation outlier."""
    contact_seq_ids: dict = {}
    for n in neighbors:
        if (n.get("distance") or 999) > CONTACT_DIST_CUTOFF:
            continue
        asym = n.get("target_asym_id")
        seq_id = n.get("target_seq_id")
        if asym and seq_id is not None:
            contact_seq_ids.setdefault(asym, set()).add(seq_id)

    if not contact_seq_ids:
        return 0.0

    total = 0
    outliers = 0
    for asym, seq_ids in contact_seq_ids.items():
        if asym not in cache:
            inst = client.get(f"{RCSB_DATA}/polymer_entity_instance/{pdb_id}/{asym}")
            features = (inst.get("rcsb_polymer_instance_feature") or []) if inst else []
            cache[asym] = _parse_residue_features(features)
        res_data = cache[asym]
        for sid in seq_ids:
            total += 1
            rdata = res_data.get(sid, {})
            if any(k in _SPARSE_FEATURE_TYPES for k in rdata):
                outliers += 1

    return outliers / total if total else 0.0


def get_ligand_quality(client: RCSBClient, pdb_id: str, nonpolymer_entity_ids: list) -> list:
    """Return per-ligand quality metrics including binding-site residue validation."""
    results = []
    chain_feature_cache: dict = {}

    for eid in nonpolymer_entity_ids:
        entity_data = client.get(f"{RCSB_DATA}/nonpolymer_entity/{pdb_id}/{eid}")
        if not entity_data:
            continue
        container = entity_data.get("rcsb_nonpolymer_entity_container_identifiers") or {}
        comp_id = container.get("nonpolymer_comp_id") or eid
        instance_ids = container.get("asym_ids") or []
        np_description = (entity_data.get("rcsb_nonpolymer_entity") or {}).get("pdbx_description") or ""

        for asym_id in instance_ids:
            inst = client.get(f"{RCSB_DATA}/nonpolymer_entity_instance/{pdb_id}/{asym_id}")
            if not inst:
                continue

            scores = (inst.get("rcsb_nonpolymer_instance_validation_score") or [None])[0] or {}
            rscc = scores.get("RSCC")
            rsr = scores.get("RSR")
            rmsz_bonds = scores.get("mogul_bonds_RMSZ")
            rmsz_angles = scores.get("mogul_angles_RMSZ")
            intermolecular_clashes = scores.get("intermolecular_clashes")
            isi = scores.get("is_subject_of_investigation")
            interesting = is_interesting_ligand(comp_id, isi)

            neighbors = inst.get("rcsb_target_neighbors") or []
            outlier_frac = _get_contact_residue_outlier_fraction(
                client, pdb_id, neighbors, chain_feature_cache
            )

            n_contacts = sum(
                1 for n in neighbors
                if (n.get("distance") or 999) <= CONTACT_DIST_CUTOFF
            )
            contact_residues = sorted({
                f"{n.get('target_comp_id','?')}{n.get('target_auth_seq_id','?')}"
                f"({n.get('target_asym_id','')})"
                for n in neighbors
                if (n.get("distance") or 999) <= CONTACT_DIST_CUTOFF
            })

            tl = traffic_light(
                rscc, rsr, rmsz_bonds, rmsz_angles,
                intermolecular_clashes, outlier_frac,
            )

            results.append({
                "ligand_id": comp_id,
                "chain_id": asym_id,
                "description": np_description,
                "is_interesting": interesting,
                "ligand_rscc": rscc,
                "ligand_rsr": rsr,
                "ligand_rmsz_bonds": rmsz_bonds,
                "ligand_rmsz_angles": rmsz_angles,
                "ligand_intermolecular_clashes": intermolecular_clashes,
                "contact_residue_count": n_contacts,
                "contact_outlier_fraction": round(outlier_frac, 3),
                "contact_residues": ";".join(contact_residues),
                "binding_quality": tl,
            })

    return results
