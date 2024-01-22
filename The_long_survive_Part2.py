#!/usr/bin/env python

"""
#################################
Created: Oct, 22 2023
Last edited: Jan 22
#################################
Analyse Codons (only in transcribed ORFS!)

This program takes as an input:

1. File1: A tsv file from Finalpipeline2
2. ORF fasta files (nucleotide, spliced) but merged
3. File2: Info on positions of ORFs in the genome (STOP codon positions are already reversed) 
4. De novo transcripts fasta OR transcriptome assembly fasta


It returns two output files:

1. A tsv file with info on whether long sequences have a 1 nt neighbour of a STOP codon at the end of alignment with the short seq (or reverse with Start)
2. A tsv file checking if the short seqs have a STOP codon in the spliced transcript (or reverse with Start)

"""

#########################################################################################
#Add input files and paths here
#########################################################################################
#PART1
File1 ="/home/m_lebh01/Documents/The_long_survive_13_01/Check_alignment_longest_vs_short_indels.txt"
File2 = "/global/students/research/m_lebh01/Protein_length_project/DetectORF/DeNovoORF_ALL_InfoFile"
Fasta ="/global/students/research/m_lebh01/Protein_length_project/DetectORF/DeNovoORF_ALL.fa"

#Path to de novo transcripts OR transcriptome assembly
AK5_TA = "/global/students/research/m_lebh01/Protein_length_project/De_novo_transcript_fasta/AK5_de_novo_transcript.fa"
DK5_TA = "/global/students/research/m_lebh01/Protein_length_project/De_novo_transcript_fasta/DK5_de_novo_transcript.fa"
GI5_TA = "/global/students/research/m_lebh01/Protein_length_project/De_novo_transcript_fasta/GI5_de_novo_transcript.fa"
SW5_TA = "/global/students/research/m_lebh01/Protein_length_project/De_novo_transcript_fasta/SW5_de_novo_transcript.fa"
UM_TA = "/global/students/research/m_lebh01/Protein_length_project/De_novo_transcript_fasta/UM_de_novo_transcript.fa"
YE_TA = "/global/students/research/m_lebh01/Protein_length_project/De_novo_transcript_fasta/YE_de_novo_transcript.fa"
Zamb_TA = "/global/students/research/m_lebh01/Protein_length_project/De_novo_transcript_fasta/Zamb_de_novo_transcript.fa"

#import the necessary packages and modules
import re
from Bio import SeqIO
import pandas as pd

#########################################################################################
#Examine the STOP codons in long ORF and spliced transcript
#########################################################################################

#This function identifies all 1nt neighbours of a codon
def Make_neighbour(seq):
    nuc1 = ["A", "T", "C", "G"]
    nuc2 = ["A", "T", "C", "G"]
    nuc3 = ["A", "T", "C", "G"]
    if len(seq) == 3:
        nuc1.remove(seq[0])
        nb1 = nuc1[0] + seq[1:]
        nb2 = nuc1[1] + seq[1:]
        nb3 = nuc1[2] + seq[1:]

        nuc2.remove(seq[1])
        nb4 = seq[0] + nuc2[0] + seq[2]
        nb5 = seq[0] + nuc2[1] + seq[2]
        nb6 = seq[0] + nuc2[2] + seq[2]

        nuc3.remove(seq[2])
        nb7 = seq[0:2] + nuc3[0]
        nb8 = seq[0:2] + nuc3[1]
        nb9 = seq[0:2] + nuc3[2]

        Neighbour_list = [nb1, nb2, nb3, nb4, nb5, nb6, nb7, nb8, nb9]
    else:
        print("Your input is not a Codon")
        Neighbour_list = 0
    return Neighbour_list

#This function finds all 1nt neighbours of the STOP of a shorter ORF and searches if they are present in the long ORF
def Detect_stop_neighbour_in_long(File1, File2, Fasta):
    file = open(File1, "r") #File with long <-> short comparison
    file2 = open(File2, "r").readlines() #File with genomic pos + stop codon of each ORF
    Outfile = open("Stop_neighbour_info_long_short_indels.txt", "w")
    Outfile.write("seq_id" + "\t" + "Orthogroup_ID" + "\t" + "Neighbour_present" + "\t" + "Searched_codon_in_short_seq" + "\t" + "Codon_at_same_pos_in_long_seq"  + "\t" + "Neighbour_present2" + "\t" + "Searched_codon_in_short_seq2" + "\t" + "Codon_at_same_pos_in_long_seq2" + "\t" + "Type" + "\n")
    test = 0
    info_dict = {}
    seq_dict = SeqIO.to_dict(SeqIO.parse(Fasta, 'fasta'))
    for line in file:
        for line2 in file2:
            l, l2 = line.strip().split("\t"), line2.strip().split(",")

            if l2[0] == l[3].split("::")[0]: #search for short seq in second file
                seq_id = l[2]
                STOP = l2[3]
                NSTOP = Make_neighbour(STOP) #Get all 1nt neighbours of that STOP codon
                record = seq_dict[l[2].split("::")[0]]
                seq = record.seq

                #Extract position of STOP neighbour
                start = int(l[5]) * 3
                end = start + 3
                fasta = seq[start:end]

                #Extract position of potential start
                START = "ATG"
                NBSTART = Make_neighbour(START)
                start2 = int(l[4]) * 3 - 3
                end2 = start2 + 3
                fasta2 = seq[start2:end2]

                if fasta in NSTOP:
                    test1 = "Stop_neighbour"
                    codon1 =  STOP
                    cod1 = fasta
                elif fasta == STOP:
                    test1 = "Stop_present"
                    codon1 = STOP
                    cod1 = fasta
                elif fasta not in NSTOP:
                    test1 = "No_Stop_or_neighbour"
                    codon1 = STOP
                    cod1 = fasta
                if fasta2 in NBSTART:
                    test2 = "Start_neighbour_present"
                    codon2 = "ATG"
                    cod2 = fasta2
                elif fasta2 == "ATG":
                    test2 = "Start_codon_present"
                    codon2 = "ATG"
                    cod2 = fasta2
                elif fasta2 not in NBSTART:
                    test2 = "No_start_neighbour_present"
                    codon2 = "ATG"
                    cod2 = fasta2
                Outfile.write(seq_id + "\t" + l[0] + "\t" + test1 + "\t" + codon1 + "\t" + str(cod1)  +  "\t" + test2 + "\t" + codon2 + "\t" + str(cod2)  +  "\t" + str(l[8]) + "\n")
                test1 = 0
                test2  = 0
                codon1 = 0
                cod1 = 0
                codon2 = 0
                cod2 = 0

#This function searches for the START/STOP codon of the long sequence in the spliced transcript of the short one
def check_for_stop_in_transcript(File1, File2):
    Outfile = open("Check_for_Stop_Codon_transcript_version2.csv", "w")
    Outfile.write("Short_seq,Stop_present,Stop_transcript,Stop_long,Start_present,Start_transcript,Start_long,Type"+ "\n")
    file = open(File1, "r")
    file2 = open(File2, "r").readlines()
    test1 = 0
    test2 = 0
    result_dict = {}
    Stop = {}
    STARTNB = Make_neighbour("ATG")

    #Sequence dictionary for each transcript fasta file
    AK5_dict = SeqIO.to_dict(SeqIO.parse(AK5_TA, "fasta"))
    DK5_dict = SeqIO.to_dict(SeqIO.parse(DK5_TA, "fasta"))
    GI5_dict = SeqIO.to_dict(SeqIO.parse(GI5_TA, "fasta"))
    SW5_dict = SeqIO.to_dict(SeqIO.parse(SW5_TA, "fasta"))
    UM_dict = SeqIO.to_dict(SeqIO.parse(UM_TA, "fasta"))
    YE_dict = SeqIO.to_dict(SeqIO.parse(YE_TA, "fasta"))
    Zamb_dict = SeqIO.to_dict(SeqIO.parse(Zamb_TA, "fasta"))
    #Save the Stop codon of each sequence id in a dictionary
    for line in file2:
        l = line.strip().split(",")
        Stop[l[0]] = l[3]
    print("Finished with the dictionarys")

    for line in file:
        for line2 in file2:
            l, l2 = line.strip().split("\t"), line2.strip().split(",")
            if l2[0] == l[3].split("::")[0]:

                ID = l[3]
                #Do this for STOP codon
                split = l2[0].split("_")
                Start = int(split[4]) - (int(l[7]) - int(l[5])) *3 #Length long seq - End alignment long seq
                End = Start + 3
                Population = split[0]
                x = Population.rstrip() + "_dict"
                seq = eval(x)[split[1]][Start:End]
                seq = str(seq.seq)
                seq = seq.upper()

		#Do this for START codon
                Start2 = int(split[3]) + 1 - (int(l[4]) * 3 -1) #Start ORF in transcript (+1 cause python 0 based) #Start alignment long (-1 to include aligned start pos), * 3 cause prot -> nuc
                End2 = Start2 + 3
                Population = split[0]
                x = Population.rstrip() + "_dict"
                seq2 = eval(x)[split[1]][Start2:End2]

                seq2 = str(seq2.seq)
                seq2 = seq2.upper()
                STOP = Make_neighbour(Stop[l[2].split("::")[0]])

                if seq in STOP:
                    test1 = "Stop_neighbour"
                    codon1 =  Stop[l[2].split("::")[0]]
                    cod1 = seq
                elif seq == Stop[l[2].split("::")[0]]:
                    test1 = "Stop_present"
                    codon1 = Stop[l[2].split("::")[0]]
                    cod1 = seq
                elif len(seq) < 3:
                    test1 = "Search_disrupted_by_TTS"
                    codon1 = Stop[l[2].split("::")[0]]
                    cod1 = seq
                elif seq not in STOP:
                    test1 = "No_Stop_or_neighbour"
                    codon1 = Stop[l[2].split("::")[0]]
                    cod1 = seq
                if End2 < 0:
                    test2 = "Search_disrupted_by_TSS"
                    codon2 = "ATG"
                    cod2 = " "
                elif seq2 in STARTNB:
                    test2 = "Start_neighbour_present"
                    codon2 = "ATG"
                    cod2 = seq2
                elif seq2 == "ATG":
                    test2 = "Start_codon_present"
                    codon2 = "ATG"
                    cod2 = seq2
                elif seq2 not in STARTNB:
                    test2 = "No_start_neighbour_present"
                    codon2 = "ATG"
                    cod2 = seq2
                print(ID, test1, Start, End, test2, Start2, End2)
                result_dict[ID] = [test1, cod1, codon1,test2, cod2, codon2, l[8]]
                Outfile.write(ID + "," + test1 + "," + str(cod1) + "," +  str(codon1) + ","  + test2 + "," + str(cod2) + "," +  str(codon2) + "," + str(l[8]) + "\n")
                test1 = 0
                test2  = 0
                codon1 = 0
                cod1 = 0
                codon2 = 0
                cod2 = 0


#Run in python
#COMPARE STOP CODONS
check_for_stop_in_transcript(File1, File2)
Detect_stop_neighbour_in_long(File1, File2, Fasta)




