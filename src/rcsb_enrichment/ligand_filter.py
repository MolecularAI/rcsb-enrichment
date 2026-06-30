"""Ligand interest classification: drug-like vs ions/solvents/cofactors/detergents."""

# CCD codes that are never drug-like ligands of interest.
_NON_INTERESTING_CCD: frozenset = frozenset({
    # Water / solvents / cryoprotectants
    "HOH", "DOD", "EDO", "EGL", "PEG", "PE3", "PE4", "PE5", "PE7", "PE8",
    "GOL", "DMS", "ACT", "ACE", "MPD", "IPA", "EOH", "PGE", "PGO", "TFP",
    "TRS", "BU3", "BU2", "BU1", "1BO", "BME", "2ME", "MLI",
    # Inorganic ions
    "MG", "ZN", "CA", "NA", "CL", "FE", "MN", "CO", "NI", "CU", "K",
    "BR", "IOD", "F", "CD", "HG", "PT", "AU", "AG", "LI", "RB", "CS",
    "BA", "SR", "AL", "GA", "IN", "TL", "PB", "BI", "SB", "AS", "SE",
    "TE", "PO4", "SO4", "SO3", "NO3", "NO2", "CO3", "CO2", "OH", "OXY",
    "FE2", "FE3", "FES", "SF4", "F3S", "CLF", "ZNO", "CUA", "CUB", "MO",
    "MO3", "MO5", "MO6", "W", "WO4", "VO4", "REO",
    # Detergents / lipids / amphiphiles
    "OLC", "OLA", "LMT", "BOG", "DDM", "DM", "OG", "NG", "HG", "C8E",
    "PLC", "LPC", "LPG", "LPE", "LPS", "PGP", "DAG", "TAG", "MAG",
    "SDS", "DPC", "LDAO", "CHAPSO", "CHAPS",
    # Buffers / additives
    "MES", "HEPES", "TRIS", "PIPES", "MOPS", "BICINE", "CAPS", "CHES",
    "EPPS", "BIS", "TCEP", "DTT", "BME", "EDO", "PEG", "MPD", "PG4",
    # Common non-drug cofactors / metals in enzyme active sites
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "CDP", "UTP", "UDP",
    "NAD", "NAP", "NDP", "FAD", "FMN", "HEM", "HEC", "HEA", "COA", "COB",
    "SAM", "SAH", "THF", "TPP", "PLP", "PLR", "BIO", "FES", "SF4",
    "MTE", "F43", "CLA", "BCL", "PLM", "RET", "SRT",
})


# CCD codes that are always uninteresting regardless of the ISI flag.
# UNX = "unknown atom or ion" — identity not established, no drug-like value.
_ALWAYS_UNINTERESTING: frozenset = frozenset({"UNX"})


def is_interesting_ligand(comp_id: str, is_subject_of_investigation) -> bool:
    """Return True if the ligand is likely a drug-like molecule of interest.

    Uses the RCSB ISI flag as primary signal; falls back to the CCD exclusion
    list for older entries where the flag is absent.
    Codes in _ALWAYS_UNINTERESTING are excluded regardless of the ISI flag.
    """
    if comp_id.upper() in _ALWAYS_UNINTERESTING:
        return False
    if is_subject_of_investigation is True:
        return True
    if is_subject_of_investigation is False:
        return False
    return comp_id.upper() not in _NON_INTERESTING_CCD
