The script functions.py contains basic functions used for the detection of Non-B DNA motifs in the MPRA data.

The functions accept DNA sequence strings as input. They also constrain by the spacer and arm lengths for IRs, DRs and MRs. G-quadruplexes were identiified in both orienttations by reverse complementing the DNA sequence string.

Detection of G4 motifs was performed with:
G4s=re.finditer("([gG]{3,}\w{1,7}){3,}[gG]{3,}",seq)

Detection of IRs, MRs and DRs was performed with:
IRs=IR_function(seq,0,5,10)
DRs=DR_function(seq,0,5,10)
MRs=MR_function(seq,0,5,10)

STRs were identified with:
STRs=master_function_STR(seq)

Z-DNA were identified with:
Z_DNAs=Z_DNA_find(seq)

GC content correction was performed for the MPRA expression sccores with:
coefs=np.polyfit(GC_contentL,ScoresL,1)



