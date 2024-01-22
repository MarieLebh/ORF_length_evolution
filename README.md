# ORF_length_evolution
Script to recreate the Data analysis 

# Step1: Detect ORF in de novo transcripts

This part of the analysis takes a set of de novo transcripts from 7 populations of D.melanogaster and detects all ORFs in forward direction. Following this the getORF output is filtered using a custom python script to get the final set of ORFs to use for the follow up analysis.

Input:  
- De novo transcripts (fasta file) from seven populations, genome assembly, transcriptome assembly and annotation

Output: 
- A file with information on the genomic position of the ORFs (TranscriptName,StartORFinGenome,EndORFinGenome,StopCodon,Chrom,Population,Strand)
- The corrected nucleotide (with added STOP) and protein ORF files where the following was removed: Duplicate ORFs (different spliceform but same genomic positions) and ORFs that do not have a stop codon at their end

Run like this:
1) Run Get ORF (run similarly for each population)
``` 
#DNA
getorf -sequence AK5DeNovoTranscripts.fa  -outseq AK5ORFs_DNA.fa -minsize 90 -reverse No -find 3
#Protein
getorf -sequence AK5DeNovoTranscripts.fa  -outseq AK5ORFs_Prot.fa -minsize 90 -reverse No -find 1
```
2) Run the scripts Part1_Pg1 - Pg7 (from https://github.com/AnnaGrBio/Proto-gene_emergence). You will need to make the following changes:
- The output file with the information of the stop codon of each ORF needs to have these columns (Chrom, Population and Strand need to be added): TranscriptName,StartORFinGenome,EndORFinGenome,StopCodon,Chrom,Population,Strand
- In cases two alternate spliceforms of a transcript contain the same ORF (= same genomic position, population + strand) only choose one.
- Adapt the script so it works for multiple ORFs per transcript and not just one.

# Step2: Build ORF orthogroups and detect length changes
In this part orthogroups are built from the de novo ORFs. ORFs are grouped into an orthogroup when they (i) share the same neighbouring genes and (ii) fullfill certain blast criteria (100 % query coverage, 90 % identity, E-value < 0.0001). We only considered orthogroups without paralogs containing at least two populations for further analysis. Following this length changes are assessed in suitable orthogroups. 

Input:
- De novo ORFs fasta files (protein)

Output:
- Orthogroups
- Info File with length changes
- File with longest sequence mapping to shorter ORFs

Run like this:
``` 
python3 The_long_survive_Part1.py
``` 

# Step3: Compare codons in long and short ORFs 
This part compares the codons in long and short ORFs.

Input:
- File with longest sequence mapping to shorter ORFs (from Step2)
- De novo ORF fasta files (nucleotide) (from Step1)
- Info File with ORF position in the genome (from Step1)
- De novo transcripts OR transcriptome assembly fasta

Output:
- File with info on whether long ORFs have a 1 nt neighbour of a stop codon at the end of alignment with the short ORF (or reverse with a start codon)
- File checking if the short ORFs have a stop codon in the spliced transcript (or reverse with a start codon)

Run like this:

1) Merge all de novo ORFs together
``` 
cat *DeNovoORFsDNA.fa > Merged_DeNovoORFsDNA.fa
``` 
2) Run python script (change all input paths at the begining)
``` 
The_long_survive_Part2.py
``` 

# Step4: Detect noncoding ORFs in syntenic regions
This part detects noncoding homologs of the transcribed ORFs in the syntenic regions of the target populations (without the transcribed ORF). It then adds checks if the non transcribed ORFs fullfill the previously set orthogroup criteria and (if yes) adds them to the orthogroups and assesses length changes again.

Input:
- Genome + Transcriptome assemblies (and corresponding gtf files)
- File where longest seq is mapped to the shorter ORFs (from Step2)
- File with genomic positions of each ORF (from Step1)
- ORF fasta files (with introns)

Output:
- BLASTN results
- Enabling mutations in homologs
- File with length change info of transcribed and nontranscribed ORFs
- File for R Figure 5

Run like this:

1) Prepare the files
``` 
python3 The_long_survive_Part3.py
``` 
2) Donwload the code from Github (https://github.com/AnnaGrBio/Proto-gene_emergence) and run Part2_Pg7-12
3) Combine transcribed and nontranscribed ORFs into one dataframe with length change information
```
python3 The_long_survive_Part4.py #Restructure file with transcribed ORFs
python3 The_long_survive_Part5.py #Add nontranscribed ORFs
``` 
4) Reformat the file for R (Figure 5)
``` 
python3 The_long_survive_Part6.py
```

# Step5: Plot the results
This part shows how to recreate the figures 4 and 5.

Input:
- Combined file with non transcribed and transcribed ORFs (from Step 4)
- File to create Figure 5 (from Step 4)
- File where each short ORF is aligned to the longest (transcribed only, from Step 2)
- The two output files from Step3 (for transcribed) plus one for non transcribed

Output:
- Figure 4 and 5

How to run:
1) Run the script The_long_survive_Part7.R in R 
