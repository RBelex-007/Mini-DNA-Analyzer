from DNA_Toolkit import *
import random


#Creating a random DNA sequence
randDNAstr = ''.join([random.choice(Nucleotides)
                      for nuc in range(50)])

rndDNAdstr = "actcgatgca"
DNAstr = validateSeq(randDNAstr)
print(CountNUCfreq(DNAstr))

print("[1]. " + DNAstr)

print("[2]. " + transcription(DNAstr))

print(f"[3]  + DNA String +Reverse Complement:\n5'  {DNAstr}  3'\n")

print(f"    {''.join(['|' for c in range(len(DNAstr))])}")

print(f"3'  {Complement(DNAstr)}  5'\n")

print(f"[4] Reverse complement: 5'  {Reverse_Complement(DNAstr)}  3'\n")

print(f"[5] GC Content:  {GC_content(DNAstr)}%\n")

print(f"[6] GC Content Subsection 5:  {GC_content_subsec(DNAstr, k=5)}\n")

print(f"[7] Amino-acid sequence from DNA:  {translate_seq(DNAstr, 0)}\n")

print(f'[8] Codon frequency (L):  {codon_usage(DNAstr, "L")}\n')