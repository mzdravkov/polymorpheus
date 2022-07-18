from hgvs.parser import Parser

THREE_LETTER_CODE_TO_ONE_LETTER = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Asx": "B",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Glx": "Z",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}


def __get_valid_hgvs(hgvs):
    if hgvs.startswith('p.'):
        # Missing transcript id.
        # This is the expected format from snpEff
        return 'MOCK:' + hgvs
    return hgvs


def parse_hgvs(hgvs):
    """
    Expects an HGVS string and returns
    a tuple of start_pos, end_pos, alt_protein
    """
    parser = Parser()
    variant = parser.parse(hgvs)
    posedit = variant.posedit
    start = posedit.pos.start.base - 1
    end = posedit.pos.end.base - 1
    alt = posedit.edit.alt
    return start, end, alt


def get_protein_variant(protein_seq, hgvs):
    valid_hgvs = __get_valid_hgvs(hgvs)
    start, end, alt = parse_hgvs(valid_hgvs)

    if alt.endswith('*'):
        return protein_seq[:start] + alt
    
    return protein_seq[:start] + alt + protein_seq[end + 1:]
    

def get_range_from_hgvs(hgvs):
    valid_hgvs = __get_valid_hgvs(hgvs)
    start, end, _ = parse_hgvs(valid_hgvs)
    return start, end
