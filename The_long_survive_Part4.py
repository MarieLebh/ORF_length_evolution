#!/usr/bin/env python
"""
#################################
Created: Dec 12, 2023
Last edited: Dec 14, 2023
#################################

This program takes as an input:

1. Short ORFs mapped to long (File)

It will then restructure the file like this:

Output:
1. Merged file with the following columns:

	Column 1 : Orthogroup ID
	Column 2 : Lines containing the ORF <comma separated>
	Column 3 : Lengths <comma separated, same order as col 2>
	Column 4 : Start coordinate wrt longest <comma separated, same order as col 2>
	Column 5 : End coordinate wrt longest <comma separated, same order as col 2>
	Column 6: Change or No-change
	
"""

#import the necessary packages and modules
import pandas as pd

#############################################################################################
#Add input files and path to them here
#############################################################################################

infile = "../The_long_survive_13_01/Check_alignment_longest_vs_short_indels.txt"

#############################################################################################
#Restructure the input file
#############################################################################################

#Reformat file and get final outfile
def reformat_file(infile):
    Pops = []
    Lens = []
    Change = 0
    Orthogroup = "Remove"
    Starts = 0
    Ends = 0
    Start = []
    End = []
    dict_x = {}
    #Populations = ["AK5", "DK5", "GI5", "SW5", "UM", "YE", "Zamb"]
    Orthogroup = 0
    with open(infile, "r") as File:
        next(File)
        for line in File:
            l = line.strip().split("\t")
                   
            if Orthogroup != l[0]:
                if len(set(Lens)) > 1:
                    Change = "Change"
                if len(set(Lens)) == 1:
                    Change = "No_Change"
                dict_x[Orthogroup] = [",".join(Pops), ",".join(Lens), ",".join(Start), ",".join(End), Change]
                
                Orthogroup = 0
                Pops = []
                Lens = []
                Start =  []
                End = []
                Change = 0
            
                Pop1 = l[2].split("::")[1]
                Pop2 = l[3].split("::")[1] 
                Starts = "1"
                Ends = l[7]
                Start.append(Starts)
                End.append(Ends)
                Start.append(l[4])
                End.append(l[5])
                Len_long = l[7]  
                Len_short = l[6]         
                Orthogroup = l[0] 
                Pops.append(Pop1)
                Pops.append(Pop2)
                Lens.append(Len_long)
                Lens.append(Len_short)
                
            elif Orthogroup == l[0]:
                Pop2 = l[3].split("::")[1] 
                Len_short = l[6]         
                Orthogroup = l[0] 
                Start.append(l[4])
                End.append(l[5])
                Pops.append(Pop2)
                Lens.append(Len_short) 

        if len(set(Lens)) > 1:
            Change = "Change"
        if len(set(Lens)) == 1:
            Change = "No_Change"
        dict_x[Orthogroup] = [",".join(Pops), ",".join(Lens), ",".join(Start), ",".join(End), Change]     
    #print(dict_x)   
    removed_value = dict_x.pop(0)
    headers=["Orthogroup_ID", "Populations", "Lengths", "Start", "End", "Change"]
    df = pd.DataFrame([[k] + v for k,v in dict_x.items()], columns=headers)                  
    #df = pd.DataFrame(dict_x.items(), columns=["Orthogroup_ID", "Populations", "Lengths", "Start_transcript_long", "End_transcript_long", "Change"])
    df.to_csv("Restructured_Outfile_ORF_lengths_v2.csv", sep = "\t", index = False)


#############################################################################################
#Run this in python
#############################################################################################
         
reformat_file(infile)
