#!/usr/bin/env python

"""
#################################
Creatd: Jan 10, 2024
Last edited: Jan 11, 2024
#################################
This program reformats the outputfile created before to  create Fig.5

As an input it takes the file from the prior program

You get one output file that can be used to create the final R plot.
"""

import pandas as pd

#############################################################################################
#Add input files and path to them here
#############################################################################################

Infile = "Merged_File_HOM_ORF.txt"

#############################################################################################
#Add input files and path to them here
#############################################################################################

def Restructure_df():
    df = pd.read_csv(Infile ,sep='\t')
    df["Start"] = df["Start"].apply(lambda x: list(map(int, x.split(","))))
    df["End"] = df["End"].apply(lambda x: list(map(int, x.split(","))))
    df["Lengths"] = df["Lengths"].apply(lambda x: list(map(int, x.split(","))))
    #Create the outputfile
    with open("Results.txt", "w") as outfile:
        outfile.write("Orthogroup,Population1,Population2,Type,Difference"+ "\n") #header (aka categories we want)
        for index, row in df.iterrows():
            start_values = row["Start"]
            end_values = row["End"]
            lengths = row["Lengths"]
            orthogroup_id = row["Orthogroup_ID"]
            populations = row["Populations"].split(",")
            # Compare the first ORF (= longest ORF) with all following ORFs (= short ORFs)
            for i in range(1, len(start_values)):
                #print(i)
                length = lengths[0] - lengths[i] #Calculate difference long ORF - short ORF
                length = str(length)
                #Check category
                if start_values[i] == 1 and end_values[0] == end_values[i]:
                    outfile.write(orthogroup_id + "," + populations[0] + ","+  populations[i] + "," + "Same_start_and_end" + "," + length + "\n")
                elif start_values[i] > 1 and end_values[0] == end_values[i]:
                    outfile.write(orthogroup_id + "," + populations[0] + ","+  populations[i] + "," + "5'truncation" + "," + length +"\n")
                elif start_values[i] == 1 and end_values[0] > end_values[i]:
                    outfile.write(orthogroup_id + "," + populations[0] + ","+  populations[i] +"," + "3'truncation" + "," + length +"\n")
                else:
                    outfile.write(orthogroup_id + "," + populations[0] + ","+  populations[i] + "," + "short_seq_in_long" + "," + length + "\n")

#############################################################################################
#Run the script here
#############################################################################################

Restructure_df()

