# ORF_length_evolution
Script to recreate the Data analysis 

#Step1: Detect ORF in de novo transcripts

This part of the analysis takes a set of de novo transcripts from 7 populations of D.melanogaster and detects all ORFs in forward direction. Following this the getORF output is filtered using a custom python script to get the final set of ORFs to use for the follow up analysis.

Input:  
-De novo transcripts (fasta file) from seven populations, genome assembly, transcriptome assembly and annotation

Output: 
-A file with information on the genomic position of the ORFs (Chromosome, Start, End, Strand, Stop)
- The corrected nucleotide (with added STOP) and protein ORF files where the following was removed: Duplicate ORFs (different spliceform but same genomic positions) and ORFs that do not have a stop codon at their end

Run like this:

#Run Get ORF (run similarly for each population)
#DNA
getorf -sequence AK5DeNovoTranscripts.fa  -outseq AK5ORFs_DNA.fa -minsize 90 -reverse No -find 3
#Protein
getorf -sequence AK5DeNovoTranscripts.fa  -outseq AK5ORFs_Prot.fa -minsize 90 -reverse No -find 1
