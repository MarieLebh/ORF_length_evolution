#!/usr/bin/env python


"""
#################################
Created: Nov, 2023
Last edited: Jan 22, 2024
#################################
PREPARE DATA FOR SYNTENY ANALYSIS

This script creates the necessary input file to run the synteny analysis found at:

https://github.com/AnnaGrBio/Proto-gene_emergence

After running this code, you can download the necessary scripts from Github(Part2_Pg7-12) and run them.

"""

#Import all modules and packages
from Bio import Seq
from Bio import SeqIO, SearchIO
from collections import Counter
import random
import re
import pandas as pd
import gffutils
import json
from bisect import bisect_left, bisect_right
import subprocess
import os
import numpy as np
import Levenshtein


#############################################################################################
#Add input files and path to them here
#############################################################################################
#ORF info files
File1 = "/home/m_lebh01/Documents/The_long_survive_13_01/Check_alignment_longest_vs_short_indels.txt"
File2 = "/global/students/research/m_lebh01/Protein_length_project/DetectORF/DeNovoORF_ALL_InfoFile"
Fasta = "/home/m_lebh01/Documents/The_long_survive_synteny/ORF_Fasta_files/All_ORFs_Introns.fa"

#############################################################################################
#Part 1: Prepare the input files to match the follow up scripts

#This part required three input files:

#1. File 1 (From "The_long_survive_Part1.py"): Longest seq mapped to each shortest
#2. File 2: Contains genomic positions of each ORF
#3. File 3: ORF fasta files (with introns)

#It will generate all necessary output files to run the code from AnnaGrandchBio/Proto-gene_emergence.
#Contrarily it always chooses the longest sequence for the blast search.

#############################################################################################

def get_longest_ref_seq(File1, File2, Fasta):
    IDs = []
    Pops = []
    dict_x = {}
    Populations = ["AK5", "DK5", "GI5", "SW5", "UM", "YE", "Zamb"]
    #Make a dictionary with all populations that are not present per orthogroup
    with open(File1, "r") as File:
        next(File)
        for line in File:
            l = line.strip().split("\t")
            Orthogroup = l[0]
            Pop1 = l[2].split("::")[1]
            Pop2 = l[3].split("::")[1]
            try:
                dict_x[Orthogroup].append(Pop1)
                dict_x[Orthogroup].append(Pop2)
            except KeyError:
                dict_x[Orthogroup] = [Pop1]
                dict_x[Orthogroup].append(Pop2)
        for lst in dict_x.values():
            lst[:] = list(set(Populations) - set(lst))
    #Create the necessary outfiles and save the info there
    Outfile = open("InfoFile.csv", "w")
    Outfile.write("Orthogroup,NameDeNovo,NamePop,Exons,Chrom,Start,End" + "\n")
    Outfasta = open("OrthogroupsRefSeq.fa", "w")
    Outpops = open("popstosearchin.csv", "w")
    F2 = open(File2, "r").readlines()
    ORF_dict = SeqIO.to_dict(SeqIO.parse(Fasta, "fasta"))
    with open(File1, "r") as File:
        next(File)
        for line in File:
            l = line.strip().split("\t")
            ID = l[2]
            if ID not in IDs:
                IDs.append(l[2])
                for line2 in F2:
                    l2 = line2.split(",")
                    if l2[0] == ID.split("::")[0]:
                        Start =l2[1]
                        End = l2[2]
                        Chromosome = l2[4]
                        Orthogroup = l[0]
                        ID = l2[0]
                        NamePop = l2[5]
                        Search_pop = dict_x[l[0]]
                        if len(Search_pop) > 0:
                            Outfile.write(Orthogroup + "," + ID + "," + NamePop + "," + "NA" + "," + Chromosome + "," + str(Start) + "," + str(End) + "\n")
                            Outfasta.write(">" + ID + "\n")
                            Outfasta.write(str(ORF_dict[ID].seq) + "\n")
                            Outpops.write(",".join(Search_pop) + "\n")


#Run like this in python
get_longest_ref_seq(File1, File2, Fasta)

















































































































































































