import sys

def get_protein_seq(seq):
    bases= "tcag"
    bases= bases.upper()
    amino_acids= 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codons= []
    for a in bases:
        for b in bases:
            for c in bases:
               codons.append(a + b + c)
    codon_table = dict(zip(codons, amino_acids))
    
    protein_seq= ""
    for i in range(0, len(seq), 3):
        codon= seq[i:i+3]
        amino_acid = codon_table[codon]
        protein_seq += amino_acid

    return(protein_seq)

try:
    with open('output3.txt', 'w') as output3:
        with open("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/PDB/Toxo_Ser_Thr_kinase.fasta") as ser_kinase_file:
            for line in ser_kinase_file:
                line= line.rstrip()
                if line.startswith(">"):
                    ID= line
                else:
                    seq=line
                    protein_seq= get_protein_seq(seq)
                    output3.write(f'{ID}\n{protein_seq}\n')
except FileNotFoundError:
    sys.exit(f'Cannot open C:/Users/Aakanksha Choudhary/OneDrive/Desktop/PDB/Toxo_Ser_Thr_kinase.fasta. Check the filename is correct')


with open("output4.txt", 'w') as output4:

    with open("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/PDB/GenesWithSignalPeptide_SignalP.txt") as gene_signal_file:
            header= next(gene_signal_file)
            header = header.rstrip() 
            header= header.split('\t')
            header= 'ID\t'+ 'Signal_start\t' + 'Signal_End\t' + 'Signal_seq\n'
            output4.write(header)
            for line1 in gene_signal_file:
                line1= line1.rstrip()
                line1= line1.split('\t')
                ID= line1[0]
                signal_start= int(line1[1])
                signal_end= int(line1[2])

                with open("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/PDB codes/output3.txt") as signal_pep:
                    for line2 in signal_pep:
                            line2= line2.rstrip()
                            if line2.startswith(">"):
                                gene_id= line2.lstrip('>')
                            else:
                                prot_seq= line2
                                signal_pep_seq= prot_seq[signal_start -1:signal_end]
                                output4.write(f'{gene_id}\t{signal_start}\t{signal_end}\t{signal_pep_seq}\n')


                                

                
        