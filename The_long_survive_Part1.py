#!/usr/bin/env python

"""
Created: Sept, 2023
Last edited: Jan 22, 2023

ANALYZE LENGTH CHANGES IN ORTHOGROUPS

This program runs the following analysis:
1. Build orthogroups using protein BLAST from ORFs from different populations
2. Selects all orthogroups for length analysis
3. Tests that the longest sequence in each suitable orthogroup matches all others and removed the ones where thats not true.
4. Maps the shorter sequences to the longest sequence to see where they are different and returns this as a result.

This program takes as an input:
(1) The fasta files (with protein sequences) of all seven populations of Drosophila melanogaster
(--> Add them to the function: change_header_all_files())
BEFORE: Make sure that the ORFs contain a correct STOP codon and also remove duplicate ORFs (from alternatively spliced transcripts).

The rest will run based on that.

It returns:

1. 7 fasta files with a renamed header and a merged fasta file with all
2. 7 lists with the IDs of all ORFS matching the fasta header
3. The individual Blast results for each population and the blast database
4. The filtered blast hits used to build the orthogroups
5. A file with the Orthogroups according to the Blast parameters
6. A file with all ambigous Orthogroups removed
7. A file with suitable orthogroups for further analysis of length changes
8. A file with information of length changes
9. A file with information whether the longest sequences within an orthogroup match with the rest
10. A file with information on how the shorter sequences map to the longest sequence
"""

#########################################################################################
#Add input files and path to them here
#########################################################################################
#Path to fasta files 7 populations:
Path = "/global/students/research/m_lebh01/Protein_length_project/DetectORF/"

#Filenames of de novo ORFs (protein sequences)
AK5fasta = "AK5_DeNovoORFwithStop_Protein.fa"
DK5fasta = "DK5_DeNovoORFwithStop_Protein.fa"
GI5fasta = "GI5_DeNovoORFwithStop_Protein.fa"
SW5fasta = "SW5_DeNovoORFwithStop_Protein.fa"
UMfasta = "UM_DeNovoORFwithStop_Protein.fa"
YEfasta = "YE_DeNovoORFwithStop_Protein.fa"
Zambfasta = "Zamb_DeNovoORFwithStop_Protein.fa" 


#import the necessary packages and modules
import csv
import pandas as pd
import os
from Bio import SeqIO
import subprocess
from itertools import permutations
from itertools import combinations
import numpy as np
from bisect import bisect_left, bisect_right

#########################################################################################
#Part1:Prepare files and run protein blast
#########################################################################################

#CRITERIA FOR ORF ORTHOGROUPS
#Change header of ORF fasta files in the fasta file and return a textfile with the header of each population
def change_fasta(fasta, outfile, population, header_file):
    records = []
    count = 0
    file = open(header_file, "w")
    for record in SeqIO.parse(fasta, "fasta"):
        count = count + 1
        record.id = record.id + "::" + population
        record.description = record.id
        records.append(record)
        x = record.id
        file.write(x + "\n")
    SeqIO.write(records, outfile, "fasta")

#Change header of ORF for all fasta files
#This is the only function that requires input! Add the correct path to the fasta files here.
def change_header_all_files():
      path = Path
      change_fasta(path + AK5fasta, "AK5_ORFs.fa", "AK5", "ORFak5")
      change_fasta(path + DK5fasta, "DK5_ORFs.fa", "DK5", "ORFdk5")
      change_fasta(path + GI5fasta, "GI5_ORFs.fa", "GI5", "ORFgi5")
      change_fasta(path + SW5fasta, "SW5_ORFs.fa", "SW5", "ORFsw5")
      change_fasta(path + UMfasta, "UM_ORFs.fa", "UM",  "ORFum")
      change_fasta(path + YEfasta, "YE_ORFs.fa", "YE",  "ORFye")
      change_fasta(path + Zambfasta, "Zamb_ORFs.fa", "Zamb",  "ORFzamb")

#Build Blast database and Blast the ORFs from each population against the other populations
def run_blast():
      print("Starting to build the blast database.")
      subprocess.call("cat *_ORFs.fa > merged_ORFs.fa", shell=True)
      subprocess.call("makeblastdb -in merged_ORFs.fa -dbtype prot", shell=True)
      print("Starting the blast.")
      subprocess.call("blastp -db merged_ORFs.fa -query AK5_ORFs.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 -num_threads 5  > AK5_blast_ORF.txt", shell = True)
      subprocess.call("blastp -db merged_ORFs.fa -query DK5_ORFs.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 -num_threads 5  > DK5_blast_ORF.txt", shell = True)
      subprocess.call("blastp -db merged_ORFs.fa -query GI5_ORFs.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 -num_threads 5  > GI5_blast_ORF.txt", shell = True)
      subprocess.call("blastp -db merged_ORFs.fa -query SW5_ORFs.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 -num_threads 5  > SW5_blast_ORF.txt", shell = True)
      subprocess.call("blastp -db merged_ORFs.fa -query UM_ORFs.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 -num_threads 5   > UM_blast_ORF.txt", shell = True)
      subprocess.call("blastp -db merged_ORFs.fa -query YE_ORFs.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 -num_threads 5  > YE_blast_ORF.txt", shell = True)
      subprocess.call("blastp -db merged_ORFs.fa -query Zamb_ORFs.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 -num_threads 5 > Zamb_blast_ORF.txt", shell = True)
      print("Successfully finished with the blast for all populations.")

#########################################################################################
#Part2: Determine synteny for each ORF and filter blast output with blast and synteny criteria.
#########################################################################################

def create_file():
    df2 = pd.read_csv(Path + "DeNovoORF_ALL_InfoFile")
    info_list = []

    for index, row in df2.iterrows():
        #Create the necessary input file for the following function
        info_dict = {}
        info_dict['ID'] = row['TranscriptName']
        info_dict['Seq'] = "NA" #This is just to fit the format of the next function
        info_dict['Chrom'] = row['Chrom']
        info_dict['NamePop'] = row['Population']
        info_dict['Start'] = row['StartORFinGenome']
        info_dict['End'] = row['EndORFinGenome']
        info_dict['Orthogroup'] = "NA"  #This is just to fit the format of the next function
        info_dict['StopCodon'] = row['StopCodon']
        info_dict['SearchPop'] = "NA" #This is just to fit the format of the next function

        info_list.append(info_dict)

    return info_list

#Detect the neighboring gene for each ORF (taken from: https://github.com/AnnaGrBio/Proto-gene_emergence)
def neighboring_gene(info_list, neighboring_degree=1, query=None):
    info_list_neigh = []
    if query is None:
        query = info_list
    if type(query) == dict:
        pop = query.get('NamePop')
        #print(pop)
        chrom = query.get('Chrom')
        start = query.get('Start')
        end = query.get('End')

        # Path to genome annotation directory
        annotated = pd.read_csv(
            "../The_long_survive_synteny/Genome_annotations/" + pop + "_genes.csv", delimiter=",")
        chrom_filter = annotated[annotated["Chrom"] == str(chrom)]
        array_stop = pd.array(chrom_filter["Stop"])
        array_start = pd.array(chrom_filter["Start"])
        # If start position is greater end position 3'-5' orientation is assumed, switches start and end information
        if start >= end:
            start = query.get('End')
            end = query.get('Start')
        # Neighboring gene is determined by bisection of the arrays
        try:
            lower = array_stop[bisect_left(array_stop, start) - neighboring_degree]
            above = array_start[bisect_right(array_start, end) + (neighboring_degree - 1)]

            lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
            above_neigh = chrom_filter.loc[chrom_filter["Start"] == above]

            # No gene left, RoiStart = Chromosome start
            if lower > above:
                query['GeneLeft'] = 'NaN'

            query['GeneRight'] = above_neigh['Gene'].to_string(index=False)
            query['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)

            # update info_list_neigh if entry found
            if query in info_list_neigh:
                idx = info_list_neigh.index(query)
                info_list_neigh.remove(query)
                info_list_neigh.insert(idx, query.copy())

        except IndexError:
            lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
            query['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)
            query['GeneRight'] = 'NaN'
            # update info_list_neigh if entry found
            if query in info_list_neigh:
                idx = info_list_neigh.index(query)
                info_list_neigh.remove(query)
                info_list_neigh.insert(idx, query.copy())
        return query

    else:
        for dictionary in info_list:
            pop = dictionary.get('NamePop')
            chrom = dictionary.get('Chrom')
            start = dictionary.get('Start')
            end = dictionary.get('End')
            annotated = pd.read_csv("../The_long_survive_synteny/Genome_annotations/" + pop + "_genes.csv", delimiter=",")
            chrom_filter = annotated[annotated["Chrom"] == str(chrom)]
            array_stop = pd.array(chrom_filter["Stop"])
            array_start = pd.array(chrom_filter["Start"])

            # DeNovo candiate is 3'-5', swap start and end to get Roi
            if start >= end:
                start = dictionary.get('End')
                end = dictionary.get('Start')

            try:
                lower = array_stop[bisect_left(array_stop, start) - neighboring_degree]
                above = array_start[bisect_right(array_start, end) + (neighboring_degree - 1)]

                lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
                above_neigh = chrom_filter.loc[chrom_filter["Start"] == above]

                # No gene left, RoiStart = Chromosome start
                if lower > above:
                    dictionary['GeneLeft'] = 'NaN'

                dictionary['GeneRight'] = above_neigh['Gene'].to_string(index=False)
                dictionary['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)
                info_list_neigh.append(dictionary.copy())

            # No gene right, RoiEnd should be chromosome end
            except IndexError:
                lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
                dictionary['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)
                dictionary['GeneRight'] = 'NaN'
                info_list_neigh.append(dictionary.copy())
                continue
    return info_list_neigh

#Filter the blast matches by criteria and only return the results where this matches as a df
def filter_matches(blast,  left_gene_dict, right_gene_dict)-> "df":
      with open(blast) as file_in:
              qseqid = []
              sseqid = []
              pident = []
              length = []
              qlen = []
              slen = []
              qstart = []
              qend = []
              sstart = []
              send = []
              evalue= []
              bitscore = []
              qcovhsp = []
              population_query = []
              population_hit = []
              for line in file_in:
                      l = line.split()
		      #ADD SYNTENY HERE(!!!!!)
                      x = l[0].split("::")[0]
                      y = l[1].split("::")[0]
                      if l[0] != l[1] and float(l[8]) < 0.0001 and float(l[2]) > 90 and left_gene_dict[x] == left_gene_dict[y] and right_gene_dict[x] == right_gene_dict[y]:
                            qseqid.append(l[0])
                            sseqid.append(l[1])
                            pident.append(l[2])
                            length.append(l[3])
                            qstart.append(l[4])
                            qend.append(l[5])
                            sstart.append(l[6])
                            send.append(l[7])
                            evalue.append(l[8])
                            bitscore.append(l[9])
                            qcovhsp.append(l[10])
                            qlen.append(l[11])
                            slen.append(l[12])
                            population_query.append(l[0].split('::')[-1])
                            population_hit.append(l[1].split('::')[-1])
               #All query-subject matches that fullfill the conditions are appended to a df
              df = pd.DataFrame()
              df["qseqid"] = qseqid
              df["sseqid"] = sseqid
              df["pident"] = pident
              df["alignment_length"] = length
              df["query_length"] = qlen
              df["subject_length"] = slen
              df["qstart"] = qstart
              df["qend"] = qend
              df["sstart"] = sstart
              df["send"] = send
              df["evalue"] = evalue
              df["bitscore"] = bitscore
              df["qdovhsp"] = qcovhsp
              df["population_query"] = population_query
              df["population_hit"] = population_hit
              df["query_hit"] = df["qseqid"].astype(str) + " " + df["sseqid"].astype(str)
              all_hits = qseqid + sseqid
      return all_hits, df #returns a list of all hits as well as a df with all filtered hits

#read in a textile and save the ids saved there as a list
def list_de_novo(de_novo)-> "list":
        with open(de_novo) as file_in:
                lines = [line.rstrip() for line in file_in]
        return lines

#Add all ORFs that do not have a blast hit
def add_no_matches(a, b, df):
    population = []
    a = set(a) #list of all ORFs
    print("Results for ", population)
    print("The number of all detected ORFs is: ", len(a))
    b = set(b) # list of all blast hit ORFs
    print("The number of ORFs with a blast hit matching the filtering criteria is: ", len(b))
    c = list(a - b) #list of all ORFs  without blast hits
    print("The number of specific ORFs without a sufficient blast hit is: ", len(c))
    count = 0
    for i in c:
        population.append(c[count].split("::")[1])
        count = count + 1
    df2 = pd.DataFrame()
    df2["qseqid"] = c
    df2["population_query"] = population
    final_df = pd.concat([df, df2])
    return final_df

#Filter the blast hits for one population
def filter_hits_population(blast, de_novo,  left_gene_dict, right_gene_dict):
      b, df  = filter_matches(blast,  left_gene_dict, right_gene_dict)
      a = list_de_novo(de_novo)
      final_df = add_no_matches(a, b, df)
      return final_df

#Final function for part 1
def final_function_part1():
    subprocess.call("cat *_blast_ORF.txt > All_Blast_Hits.txt", shell = True)
    subprocess.call("cat ORF* > AllORF", shell = True)
    #Do the synteny check
    info_list = create_file()
    print(info_list[:10])
    info_list_neigh = neighboring_gene(info_list, neighboring_degree=1, query=None)
    info_df = pd.DataFrame(info_list_neigh)
    left_gene_dict = info_df.set_index('ID')['GeneLeft'].to_dict()
    right_gene_dict =  info_df.set_index('ID')['GeneRight'].to_dict()
    df = filter_hits_population("All_Blast_Hits.txt", "AllORF", left_gene_dict, right_gene_dict)
    df.to_csv("Filtered_Blast_Hits_temp.csv", index = False, sep = "\t")


#########################################################################################
#Part3: Build orthogroups
#########################################################################################

#Build orthogroups based on the filtered Blast file (adapted from: https://github.com/AnnaGrBio/Transcripts-gain-and-loss)
def build_orthogroups():
    Dict_Groups = {}
    c = 1
    F = open("Filtered_Blast_Hits_temp.csv", "r")
    L=F.readlines()
    for i in L[1:]:
        line = i.split("\t")
        if len(Dict_Groups) == 0:
            Group = "Orthogroup"+str(c)
            c+=1
            Dict_Groups[Group] = [line[0], line[1]]
        else:
            seq1 = line[0]
            seq2 = line[1]
            if not seq2:
                seq2 = seq1
            Attribute = False
            for NameGroups in Dict_Groups.keys():
                MemberOfGroup = Dict_Groups[NameGroups]
                if seq1 in MemberOfGroup and seq2 in MemberOfGroup:

                    Attribute = True
                    break
                elif seq1 in MemberOfGroup and seq2 not in MemberOfGroup:

                    Dict_Groups[NameGroups].append(seq2)
                    Attribute = True
                    break
                elif seq2 in MemberOfGroup and seq1 not in MemberOfGroup:

                    Dict_Groups[NameGroups].append(seq1)
                    Attribute = True
                    break
            if Attribute == False:
                Group = "Orthogroup"+str(c)
                c+=1
                if seq1 == seq2:
                    Dict_Groups[Group] = [seq1]
                else:

                    Dict_Groups[Group] = [seq1, seq2]
    return Dict_Groups

#This function merges the dictionaries of orthogropus transcripts (adapted from: https://github.com/AnnaGrBio/Transcripts-gain-and-loss)
def merge_dics(Dict):
    for NameGroups in list(Dict.keys()):
        if NameGroups in Dict.keys():
            l = Dict[NameGroups]
            for Name in l:
                for NewNames in list(Dict.keys()):
                    if NewNames != NameGroups:
                        Newl = Dict[NewNames]
                        if Name in Newl:
                            for k in Newl:
                                if k not in l:
                                    Dict[NameGroups].append(k)
                            Dict.pop(NewNames)
    return Dict
    
#Build a final file for those orthogroups. Also add which ORFs are considered to be present in the ancestor.
#Here we assume that an ORFs present in Zambia and at least one additional population is present in the ancestor
#AK5 count =  x.count("AK5") > 1 -> column: count AK5
def build_final_file(x):
            F = open("Orthogroups_protBlast.txt", "w")
            F.write("Number_Orthogroup;IDs;Number_seqs_in_orthogroup;Number_populations_in_orthogroup;Ancestral")
            F.write("\n")
            count = 0
            count_pops = 0
            Zamb_Pop = "False"
            for key, value in x.items():
                  count = count + 1
                  if any("AK5" in s for s in value): #Count if each population occurs
                        count_pops = count_pops + 1
                  if any("DK5" in s for s in value):
                        count_pops = count_pops + 1
                  if any("GI5" in s for s in value):
                        count_pops = count_pops + 1
                  if any("SW5" in s for s in value):
                        count_pops = count_pops + 1
                  if any("UM" in s for s in value):
                        count_pops = count_pops + 1
                  if any("YE" in s for s in value):
                        count_pops = count_pops + 1
                  if any("Zamb" in s for s in value):
                        count_pops = count_pops + 1
                        if count_pops >1:
                             Zamb_Pop = "True"
                  F.write(key  + ";" + ",".join(value) + ";" + str(len(value))+ ";" + str(count_pops) + ";" + Zamb_Pop)
                  F.write("\n")
                  count_pops = 0
                  Zamb_Pop = "False"
            F.close()
            print("The number of all detected orthogroups is: ", str(count))

#Final function for building the orthogroups
def final_function_part2():
      y = build_orthogroups()
      y= merge_dics(y)
      print("Finished building the orthogroups")
      final =  build_final_file(y)
      print("Finished building the final file")

#########################################################################################
#Part4:Filter orthogroups and select suitable ones for analysis
#########################################################################################
    
#Generate all possible combinations of ids within an orthogroup
#This is only necessary when running the function ambigous orthogroups
def id_combinations(group_ids):
    ids = group_ids.split(",")
    return list(permutations(ids,2))

#Check whether there are ambigous orthogroups
#By default this is NOT included in the run but you can add it if necessary
def ambigous_orthogroups():
    df1 = pd.read_csv("Orthogroups_protBlast.txt", sep = ";")
    df2 = pd.read_csv("Filtered_Blast_Hits_temp.csv", sep = "\t")
    query_hit_combinations = set(tuple(sorted(str(row).split())) for row in df2["query_hit"] if row)
    removed_groups = []
    for index, row in df1.iterrows():
        group = row["Number_Orthogroup"]
        group_ids = row["IDs"]
        possible_combinations = id_combinations(group_ids)
        all_exists = all(tuple(sorted(combination)) in query_hit_combinations for combination in possible_combinations)
        if not all_exists:
            removed_groups.append(group)
    df1 = df1[~df1["Number_Orthogroup"].isin(removed_groups)]
    print("This number of orthogroups contained ambigous content:", len(removed_groups))
    print("These orthogroups were ambiguous and were thus removed from the dataset:", removed_groups)
    df1.to_csv("Orthogroups_ORF_prot_no_ambigous.txt", sep = ";", index = False)

#This function selects all orthogroups suitable for further analysis (i.e. no duplicated ORFs, 1 < Nr Pops < 7)
def select_suitable_OGs(Orthogroup_file):
    OG_df2 = pd.read_csv(Orthogroup_file, sep = ";")
    OG_df3 = OG_df2.loc[OG_df2["Number_populations_in_orthogroup"] > 1]
    OG_df4 = OG_df3.loc[OG_df3["Number_seqs_in_orthogroup"] == OG_df3["Number_populations_in_orthogroup"]]
    OG_df5 = OG_df4.loc[OG_df4["Number_populations_in_orthogroup"] <= 7]
    OG_df5.to_csv("Suitable_Orthogroups.csv", sep = ";", index = False)
    print("The number of all orthogroups suitable for follow up analyses (specific, 2-6 pops): ", len(OG_df5))
    print("The number of orthogroups that are specific: ", OG_df2[(OG_df2.Number_seqs_in_orthogroup == 1)].count())
    print("The number of orthogroups with duplicates: ", OG_df2[(OG_df2.Number_populations_in_orthogroup != OG_df2.Number_seqs_in_orthogroup)].count())
    print("The number of orthogroups with 7 specific seqs: ", OG_df2[(OG_df2.Number_populations_in_orthogroup == 7) & (OG_df2.Number_seqs_in_orthogroup == 7)].count())
    return OG_df2

#########################################################################################
#Part5: Detect length changes and stats in suitable orthogroups
#########################################################################################

#This function compares the longest sequence in the orthogroup to the shortest one and returns the max length difference per OG
def get_length_changes(OG_file, fasta):
    OG_fasta_dict = {}
    len_dict = {}
    len_count_dict = {}
    population = 0
    with open(OG_file, "r") as file:
        for line in file:
            line = line.split(";")
            OG_ID = line[0]
            orthogroups = line[1].split(",")
            for ID in orthogroups:
                for record in SeqIO.parse(fasta, "fasta"):
                    if record.id == ID:
                        try:
                            OG_fasta_dict[OG_ID].append(record.seq)
                        except KeyError:
                            OG_fasta_dict[OG_ID] = [record.seq]
    for key, value in OG_fasta_dict.items():
        lengths = [len(i) for i in value]
        seqs_in_OG = len(value)
        MAD = np.median(np.absolute(lengths - np.median(lengths)))
        SD = np. std(lengths)
        Max_length = max(lengths)
        Min_length = min(lengths)
        Mean_length = np.mean(lengths)
        Biggest_difference = Max_length - Min_length
        Median_length = np.median(lengths)
        lengths = ",".join(str(x) for x in lengths)
        len_dict[key] = [lengths, seqs_in_OG, MAD, SD, Max_length, Min_length, Biggest_difference, Mean_length, Median_length]
    headers = ["Number_Orthogroup", "Lengths", "Seqs_in_OG", "MAD", "SD", "Max_length", "Min_length", "Biggest_difference", "Mean_length", "Median_length"]
    df = pd.DataFrame([[k] + v for k,v in len_dict.items()], columns=headers)
    df["Change"] = "Length_variation"
    df.loc[df["Biggest_difference"] == 0, "Change"] = "No_change"
    #Print which orthogroups have length changes
    print("Number of Orthogroups with ORF length change and without length change:")
    print(df["Change"].value_counts())
    return df

#This function does the same but only returns the length and no other feature
def get_length_changes_length_only(OG_file, fasta):
    OG_fasta_dict = {}
    len_dict = {}
    len_count_dict = {}
    population = 0
    with open(OG_file, "r") as file:
        for line in file:
            line = line.split(";")
            OG_ID = line[0]
            orthogroups = line[1].split(",")
            for ID in orthogroups:
                for record in SeqIO.parse(fasta, "fasta"):
                    if record.id == ID:
                        #population = ID.split("::")[1]
                        try:
                            OG_fasta_dict[OG_ID].append(record.seq)
                        except KeyError:
                            OG_fasta_dict[OG_ID] = [record.seq]
    for key, value in OG_fasta_dict.items():
        lengths = [len(i) for i in value]
        lengths = ",".join(str(x) for x in lengths)
        len_dict[key] = [lengths]
    headers = ["Number_Orthogroup", "Lengths"]
    df = pd.DataFrame([[k] + v for k,v in len_dict.items()], columns=headers)
    return df

#This function merges the last part of the analysis
def final_function_part3():
    select_suitable_OGs("Orthogroups_protBlast.txt")
    df = get_length_changes("Suitable_Orthogroups.csv", "merged_ORFs.fa")
    df2 = pd.read_csv("Orthogroups_protBlast.txt", sep = ";")
    merged = pd.merge(df, df2, on = "Number_Orthogroup", how = "left")
    merged.to_csv("Length_variation_within_orthogroup_all.csv", sep = "\t", index = False)


#########################################################################################
#Part6: Select all orthogroups where the first longest seq matches all other seqs
#########################################################################################

#Identify orthogroups where there are ORF length changes.
def process_orthogroups(orthogroup_file, fasta_dir, blast_file2):
    # Read the Orthogroup file and store the sequences in each group
    orthogroups = {}
    with open(orthogroup_file, "r") as og_file:
        next(og_file)  # Skip the header
        for line in og_file:
            line = line.strip().split("\t")
            og_id, seq_ids = line[0], line[11].split(",")
            orthogroups[og_id] = seq_ids
    # Find the longest sequences in each orthogroup
    longest_sequences = {}
    for og_id, seq_ids in orthogroups.items():
        max_length = 0
        longest_seqs = []
        for seq_id in seq_ids:
            for record in SeqIO.parse(fasta_dir, "fasta"):
                if record.id == seq_id:
                    seq_length = len(record.seq)
                    if seq_length > max_length:
                        max_length = seq_length
                        longest_seqs = [seq_id]
                    elif seq_length == max_length:
                        longest_seqs.append(seq_id)
        longest_sequences[og_id] = longest_seqs
    # Read the Blast file and check if the longest sequences have matches with every other sequence in the orthogroup
    orthogroups_with_match = set()
    for og_id, seq_ids in orthogroups.items():
        longest_seqs = longest_sequences[og_id]
        query_hit_combinations = {(query, hit) for query in longest_seqs for hit in seq_ids if hit != query}
        found_combinations = set()
        with open(blast_file2, "r") as blast_file:
            next(blast_file)  # Skip the header
            for line in blast_file:
                fields = line.strip().split("\t")
                query, hit = fields[0], fields[1]
                if (query, hit) in query_hit_combinations:
                    found_combinations.add((query, hit))
                if (hit, query) in query_hit_combinations:
                    found_combinations.add((hit, query))
        if query_hit_combinations == found_combinations:
            orthogroups_with_match.add(og_id)
    # Filter orthogroups where the longest sequences have matches with every other sequence
    filtered_orthogroups = {og_id: seq_ids for og_id, seq_ids in orthogroups.items() if og_id in orthogroups_with_match}
    others = {og_id: seq_ids for og_id, seq_ids in orthogroups.items() if og_id not in orthogroups_with_match}
    print("In these orthogroups the longest sequence/sequences do not (all) match:", len(others))
    print("These are the orthogroups where the longest sequence(s) not all match the other sequences")
    print(others)
    print("In these orthogroup the longest sequence(s) have a blast hit with all the others:", len(filtered_orthogroups))
    # Write the filtered Orthogroups to the output file
    with open("Longest_seq_match_all.txt", "w") as out_file:
        out_file.write("Orthogroup IDs\tSequences\n")
        for og_id, seq_ids in filtered_orthogroups.items():
            out_file.write(f"{og_id}\t{','.join(seq_ids)}\n")

#########################################################################################
#Part7: Map the short sequences to the longest one and see if/where they differ in length
#########################################################################################

#This function maps the shorter seqs to the longest seq of the orthogroup
def map_short_to_long(orthogroup_file, fasta_dir, blast_file2):
    # Read the Orthogroup file and store the sequences in each group
    orthogroups = {}
    with open(orthogroup_file, 'r') as og_file:
        next(og_file)  # Skip the header
        for line in og_file:
            line = line.strip().split('\t')
            og_id, seq_ids = line[0], line[1].split(',')
            orthogroups[og_id] = seq_ids
    # Find the longest sequences in each orthogroup
    longest_sequences = {}
    for og_id, seq_ids in orthogroups.items():
        max_length = 0
        longest_seqs = []
        for seq_id in seq_ids:
            for record in SeqIO.parse(fasta_dir, "fasta"):
                if record.id == seq_id:
                    seq_length = len(record.seq)
                    if seq_length > max_length:
                        max_length = seq_length
                        longest_seqs = [seq_id]
        longest_sequences[og_id] = longest_seqs
    # Read the Blast file and check if the longest sequences have matches with every other sequence in the orthogroup
    orthogroups_with_match = set()
    Outfile = open("Mapped_short_seqs_to_first_longest_out_with_indels.txt", "w")
    Outfile.write("Number_Orthogroup" + "\t" + "Number_Seqs_Orthogroup" + "\t" + "Longest_Seq" + "\t" + "Other_Seq" + "\t" + "Start_alignmentLong"  + "\t" + "End_alignmentLong"  +"\t" +  "LenShort"  + "\t" + "LenLong"  +"\t" "Comparison" +  "Indel_present" + "\t" + "Number_indels_long" + "\t" + "Number_indels_short" + "\t" +  "Difference_indels" + "\t" + "Length_diff_to_start" + "\t" + "Length_diff_to_stop" + "\t" +"Difference_long_short_aa" + "\n")
    for og_id, seq_ids in orthogroups.items():
        Comparison = 0
        indel = 0
        Num_q = 0
        Num_s = 0
        Diff = 0
        indel_abs = 0
        diff_start_dir = 0
        diff_end_dir = 0
        longest_seqs = longest_sequences[og_id]
        num_seqs = len(longest_seqs)
        num2 = len(seq_ids)
        for query in longest_seqs:
            for hit in seq_ids:
                if hit != query:
                    with open(blast_file2, 'r') as blast_file:
                        next(blast_file)  # Skip the header
                        for line in blast_file:
                        #Map each shorter sequence to the long ones and compare start and end
                            fields = line.strip().split('\t')
                            query1, hit1 = fields[0], fields[1]


                            #Pident = float(fields[2])

                            if (query, hit) == (query1, hit1):

                                AL_q = int(fields[7]) - int(fields[6]) + 1 #Alignment length query (no indels)
                                AL_s = int(fields[9]) - int(fields[8])  + 1 #Alignment length subject (no indels)
                                Len_query = int(fields[4])
                                Len_subject = int(fields[5])
                                AL = int(fields[3])

                                if AL > AL_q and AL == AL_s:
                                    indel = "indel_in_long_seq"
                                    Num_q = AL - AL_q
                                    Num_s = 0
                                elif AL > AL_s and AL_q == AL:
                                    indel = "indel_in_short_seq"
                                    Num_s = AL - AL_s
                                    Num_q = 0
                                elif AL > AL_q and AL > AL_s:
                                    indel = "indel_in_short_and_long"
                                    Num_q = AL - AL_q
                                    Num_s = AL - AL_s
                                else:
                                    indel = "no_indel"
                                    Num_q = 0
                                    Num_s = 0

                                indel_abs = abs(Num_q - Num_s)

                                if int(fields[6]) ==  1 and int(fields[4]) == int(fields[7]):
                                    Comparison = "Same_start_and_end"
                                    Diff = int(fields[4]) - int(fields[5])
                                    diff_start_dir = 0
                                    diff_end_dir = 0

                                elif int(fields[6]) == 1 and int(fields[4]) > int(fields[7]):
                                    Comparison = "Same_start_end_earlier"
                                    Diff = int(fields[4]) - int(fields[5])
                                    diff_start_dir = 0
                                    diff_end_dir = Diff - Num_s + Num_q - diff_start_dir

                                elif int(fields[6]) > 1  and int(fields[4]) == int(fields[7]):
                                    Comparison = "Same_end_and_start_later"
                                    Diff = int(fields[4]) - int(fields[5])
                                    diff_start_dir = int(fields[6]) -1 - Num_s + Num_q
                                    diff_end_dir = 0


                                elif int(fields[6]) > 1 and int(fields[4]) > int(fields[7]):
                                    Comparison = "Short_seq_in_long_seq"
                                    Diff = int(fields[4]) - int(fields[5])
                                    diff_start_dir = int(fields[6]) -1 - Num_s + Num_q
                                    diff_end_dir = Diff - Num_s + Num_q - diff_start_dir

                                else:
                                    Comparison = "There must be a mistake"  #if this occurs check the data > there is a mistake then




                                Outfile.write(og_id + "\t" + str(num2) + "\t" + query + "\t" + hit + "\t" + str(fields[6])  + "\t" + str(fields[7]) +"\t" +  str(fields[5])  + "\t" + str(fields[4]) + "\t" + Comparison  +"\t" + indel + "\t" + str(Num_q) + "\t" + str(Num_s) + "\t" +  str(indel_abs) + "\t" + str(diff_start_dir) + "\t" + str(diff_end_dir) +  "\t" + str(Diff) + "\n")

                                Comparison = 0
                                indel = 0
                                Diff = 0
                                Num = 0
                                indel_abs = 0
                                diff_start_dir = 0
                                diff_end_dir = 0

                            elif (query, hit) == (hit1, query1):

                                AL_q = int(fields[7]) - int(fields[6]) + 1 #Alignment length query (no indels)
                                AL_s = int(fields[9]) - int(fields[8])  + 1 #Alignment length subject (no indels)
                                Len_query = int(fields[4])
                                Len_subject = int(fields[5])
                                AL = int(fields[3])
                                indel_abs = abs(Num_s - Num_q)

                                if AL > AL_q and AL == AL_s:
                                    indel = "indel_in_short_seq"
                                    Num_q = AL - AL_q
                                    Num_s = 0
                                elif AL > AL_s and AL_q == AL:
                                    indel = "indel_in_long_seq"
                                    Num_s = AL - AL_s
                                    Num_q = 0
                                elif AL > AL_q and AL > AL_s:
                                    indel = "indel_in_short_and_long"
                                    Num_q = AL - AL_q
                                    Num_s = AL - AL_s
                                else:
                                    indel = "no_indel"
                                    Num_q = 0
                                    Num_s = 0

                                indel_abs = abs(Num_s - Num_q)

                                if int(fields[8]) ==  1 and int(fields[5]) == int(fields[9]):
                                    Comparison = "Same_start_and_end"
                                    Diff = int(fields[5]) - int(fields[4])
                                    diff_start_dir = 0
                                    diff_end_dir = 0

                                elif int(fields[8]) == 1 and int(fields[5]) > int(fields[9]):
                                    Comparison = "Same_start_end_earlier"
                                    Diff = int(fields[5]) - int(fields[4])
                                    diff_start_dir = 0
                                    diff_end_dir = Diff + Num_s - Num_q - diff_start_dir

                                elif int(fields[8]) > 1  and int(fields[5]) == int(fields[9]):
                                    Comparison = "Same_end_and_start_later"
                                    Diff = int(fields[5]) - int(fields[4])
                                    diff_start_dir = int(fields[8]) -1 + Num_s - Num_q
                                    diff_end_dir = 0

                                elif int(fields[8]) > 1 and int(fields[5]) > int(fields[9]):
                                    Comparison = "Short_seq_in_long_seq"
                                    Diff = int(fields[5]) - int(fields[4])
                                    diff_start_dir = int(fields[8]) -1 + Num_s - Num_q
                                    diff_end_dir = Diff + Num_s - Num_q - diff_start_dir
                                else:
                                    Comparison = "There must be a mistake" #if this occurs check the data > there must be a mistake then


                                Outfile.write(og_id + "\t" + str(num2) + "\t" + query + "\t" + hit + "\t" + str(fields[8])  + "\t" + str(fields[9]) +"\t" +  str(fields[4])  + "\t" + str(fields[5]) + "\t" + Comparison  +"\t" + indel + "\t" + str(Num_s) + "\t" + str(Num_q) + "\t" +  str(indel_abs) + "\t" + str(diff_start_dir) + "\t" + str(diff_end_dir) +  "\t" + str(Diff) + "\n")

                                Comparison = 0
                                indel = 0
                                Num = 0
                                Diff = 0
                                indel_abs = 0
                                diff_start_dir = 0
                                diff_end_dir = 0

#########################################################################################
#Part8: Run everything at once
#########################################################################################

#This function can be used to run the complete analysis
def run_complete_script():
    change_header_all_files()
    run_blast()
    final_function_part1()
    final_function_part2()
    final_function_part3()
    process_orthogroups("Length_variation_within_orthogroup_all.csv", "merged_ORFs.fa", "Filtered_Blast_Hits_temp.csv")
    map_short_to_long("Longest_seq_match_all.txt", "merged_ORFs.fa", "Filtered_Blast_Hits_temp.csv")
    subprocess.call("awk '!NF || !seen[$0]++' Mapped_short_seqs_to_first_longest_out_with_indels.txt > Check_alignment_longest_vs_short_indels.txt", shell = True)

#########################################################################################
#Run like this in python (Default is to run the complete analysis)
#########################################################################################

run_complete_script()










































