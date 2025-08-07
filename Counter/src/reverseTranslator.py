import random

AA2NA = {
    "A": list("GCT,GCC,GCA,GCG".split(",")),
    "R": list("CGT,CGC,CGA,CGG,AGA,AGG".split(",")),
    "N": list("AAT,AAC".split(",")),
    "D": list("GAT,GAC".split(",")),
    "C": list("TGT,TGC".split(",")),
    "Q": list("CAA,CAG".split(",")),
    "E": list("GAA,GAG".split(",")),
    "G": list("GGT,GGC,GGA,GGG".split(",")),
    "H": list("CAT,CAC".split(",")),
    "I": list("ATT,ATC,ATA".split(",")),
    "L": list("TTA,TTG,CTT,CTC,CTA,CTG".split(",")),
    "K": list("AAA,AAG".split(",")),
    "M": list("ATG".split(",")),
    "F": list("TTT,TTC".split(",")),
    "P": list("CCT,CCC,CCA,CCG".split(",")),
    "S": list("TCT,TCC,TCA,TCG,AGT,AGC".split(",")),
    "T": list("ACT,ACC,ACA,ACG".split(",")),
    "W": list("TGG".split(",")),
    "Y": list("TAT,TAC".split(",")),
    "V": list("GTT,GTC,GTA,GTG".split(",")),
    "*": list("TAA,TGA,TAG".split(","))
}

def aa2na(seq):
    na_seq = [random.choice(AA2NA.get(c, ["---"])) for c in seq]
    return "".join(na_seq)
