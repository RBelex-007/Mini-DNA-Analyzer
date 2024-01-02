from collections import Counter
Nucleotides = ["A", "C", "G", "T"]
DNA_ReverseCompl = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}


#Check the sequence to make sure it is a string
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

def CountNUCfreq(seq):
    tmpFreqDict = {"A": 0, "T": 0, "C": 0, "G": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

def transcription(seq):
    return seq.replace("T", "U")

def Complement(seq):
    return ''.join([DNA_ReverseCompl[nuc] for nuc in seq])

def Reverse_Complement(seq):
    return ''.join([DNA_ReverseCompl[nuc] for nuc in seq])[::-1]

def GC_content(seq):
    return round((seq.count('C') + seq.count('G')/len(seq) * 100))

def GC_content_subsec(seq, k=20):
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subsec = seq[i:i + k]
        res.append(GC_content(subsec))
    return res

def translate_seq(seq, init_pos=0):
    return [DNA_Codons[seq[pos:pos +3]]for pos in range(init_pos, len(seq)-2, 3)]

def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict

