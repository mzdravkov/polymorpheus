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


def get_protein_variant(protein_seq, hgvs):
    valid_hgvs = hgvs
    # Missing transcript id.
    # This is the expected format from snpEff
    if hgvs.startswith('p.'):
        valid_hgvs = 'MOCK:' + hgvs

    parser = Parser()
    variant = parser.parse(valid_hgvs)
    posedit = variant.posedit
    start = posedit.pos.start.base - 1
    end = posedit.pos.end.base - 1
    alt = posedit.edit.alt

    if alt.endswith('*'):
        return protein_seq[:start] + alt
    
    return protein_seq[:start] + alt + protein_seq[end + 1:]
    