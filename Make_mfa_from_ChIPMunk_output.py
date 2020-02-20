###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes output from ChIP-Munk (table with seqs)
#and prepares multy-fasta file ready for logo visualization (e.g. with WebLogo online service).
###############################################

#######
#Packages to be imported.
#######

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#######
#Variables to be defined.
#######

#Path to ChIP-Munk output.
Munk_input="C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Peaks_motif\DRIP_Seq_CTD_minus_Rif_minus_forward_peaks_threshold_300_sequences_longer_20bp_no_rRNA_LacI_peaks_motif.txt"

#Path to multi-fasta output.
Munk_output="C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Peaks_motif\DRIP_Seq_CTD_minus_Rif_minus_forward_peaks_threshold_300_sequences_longer_20bp_no_rRNA_LacI_peaks_motif.fasta"

def read_convert(input_path, output_path, reversecomplement):
    filein=open(input_path, 'r')
    fileout=open(output_path, 'w')
    for line in filein:
        line=line.rstrip().split('|')
        if line[0]=='WORD':
            data=line[1].split('\t')
            seq=data[2]
            seq_bio=Seq(seq, generic_dna)
            seq_bio_rc=seq_bio.reverse_complement()
            if reversecomplement==False:
                fileout.write(f'>{data[0]}\n{seq}\n')
            elif reversecomplement==True:
                fileout.write(f'>{data[0]}\n{str(seq_bio_rc)}\n')
    filein.close()
    fileout.close()
    print('Done!')
    return

read_convert(Munk_input, Munk_output, False)