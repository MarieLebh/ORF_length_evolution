#!/usr/bin/env python
"""
#################################
Created: Dec 14, 2023
Last edited: Dec 14, 2023
#################################

Prior to running this script you need to run Part4 Reformat ORFs!

This program takes as an input:

1. Blast output file
2. ORF (protein) fasta files
3. Info file (from synteny analysis)
4. Merged info file with all ORFs 

It will then do the following analysis:

1. Select the best hit and translate it
2. Select the longest sequence from each orthogroup 
3. Blast them and see if the best hit has a sufficient match
4. Merge the sufficient hits with the transcribed ORFs

Output:
1. Merged file with transcribed and nontranscribed ORFs
"""

import subprocess
from Bio import Seq
import pandas as pd
from Bio import SeqIO
import os

#############################################################################################
#Add input files and path to them here
#############################################################################################
#ORF info file
Info = "InfoFile.csv"
#ORF fasta file
FastaORF = "/home/m_lebh01/Documents/The_long_survive/merged_ORFs.fa"
#Blast result file
Bfile = "FinalBLAST.rst"
#Merged transcript file
OFile = "Restructured_Outfile_ORF_lengths_v2.csv"


#############################################################################################
#Extract the best hits and blast them against their respective reference ORF (transcribed)
#############################################################################################

def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L


#Get the best blast hit per sequence, from: github.com/AnnaGrBio/Proto-gene_emergence) 
def access_best_hit_rsts(File):
    Liste = []
    print ("All hits : "+str(len(File)))
    OldQueryName = ""
    OldQueryPop = ""
    for i in File[1:]:
        #print(i)
        ligne1 = i.split("\n")[0]
        ligne2 = ligne1.split(",")
        NameQueryDeNovo = ligne2[0]
        #print(NameQueryDeNovo)
        NameQueryPop = ligne2[11]
        if NameQueryDeNovo!=OldQueryName:
            OldQueryName = NameQueryDeNovo
            OldQueryPop = NameQueryPop
            Liste.append(ligne1)
        elif NameQueryDeNovo == OldQueryName and NameQueryPop!=OldQueryPop:
            Liste.append(ligne1)
            OldQueryPop = NameQueryPop
    print ("Final best hits : "+str(len(Liste)))
    return Liste

#Translate the best hit
#Check that it has an in frame Start and Stop codon. If not return an empty string as seq.    
def translate_best_hit(Liste):
    newListe = []
    Stops = ["TAA", "TAG", "TGA"]
    for i in Liste:
        item = i.split(",")
        protein = Seq.translate(item[12], to_stop = True)
        #Check if there is an inframe Start codon (aka M)
        n = protein.find("M")
        t = len(protein) 
        #Get positions of the stop codon       
        s = n*3 + t*3
        e = s+3
        #This means no start at all        
        if n < 0:
            protein = ""
        #No hit with synteny    
        elif item[1] == "0":
            protein = ""
        #Checks that the end is a stop codon        
        elif item[12][s:e] not in Stops:
            protein = ""       
            stop = "NA"
        #Takes the protein from the start codon              
        else:
            protein = protein[n:]
        if item[12][s:e] in Stops:
            stop = item[12][s:e]
        seq_id = item[0]
        population = item[11]
        newListe.append([seq_id, population, protein, stop])
    return newListe

#Run blastp against the longest sequence  
#Run blastp against the longest sequence  
def Run_blast(Liste, Fasta):
    ORF_dict = SeqIO.to_dict(SeqIO.parse(Fasta, "fasta"))
    sstart = "no_hit"
    send = "no_hit"
    qlen = "no_hit"
    count = 0
    count_empty = 0
    Result_Liste = []
    #Get the necessary sequences
    for item in Liste:
        count = count+1
        Id = item[0]
        QueryPop = item[0].split("_")
        Id = Id + "::" + QueryPop[0]
        Pop = item[1]
        Hom_prot = item[2]
        Ref = str(ORF_dict[Id].seq)
        stop = item[3]
        #Homolog(=query)
        Fasta1 = open("Temp1.fa", "w")
        Fasta1.write(">" + Id + "\n")
        Fasta1.write(Hom_prot + "\n")
        Fasta1.close()
        #LongORF(=target)
        Fasta1 = open("Temp2.fa", "w")
        Fasta1.write(">" + Id + "\n")
        Fasta1.write(Ref + "\n")
        Fasta1.close()

        subprocess.call("blastp -query Temp1.fa -subject Temp2.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen slen' -qcov_hsp_perc 100 > Out_Temp.txt", shell = True)
        
        print("Start blasting homolog from ", Pop, " of ORF ", Id)
        print("Query: ", Hom_prot)
        print("Subject: ", Ref)

        with open("Out_Temp.txt", "r") as Blast:
            if not Blast.read(1):
                count_empty = count_empty + 1
                match = False
                Id = Id
                Pop = Pop
                sstart = "no_hit"
                send = "no_hit"
                qlen = "no_hit"
                print("No hit at all/ or no ORF", "\n")
                Result_Liste.append([Id, Pop, sstart, send, qlen])
        #Implement to first check if there even was a hit
            for line in Blast:
                match = False
                Id = Id
                Pop = Pop
                sstart = "no_hit"
                send = "no_hit"
                qlen = "no_hit"

                if float(line.split("\t")[8]) < 0.0001 and float(line.split("\t")[2]) > 90:
                    l = line.split("\t")
                    match = True
                    Id = Id
                    Pop = Pop
                    sstart = str(l[6])
                    send = str(l[7])
                    qlen = str(len(Hom_prot))
                    Result_Liste.append([Id, Pop, sstart, send, qlen,stop])
                    print("Sufficient hit", "\n")
                    break
                else:
                    match = False
                    Id = Id
                    Pop = Pop
                    sstart = "no_hit"
                    send = "no_hit"
                    qlen = "no_hit"
                    Result_Liste.append([Id, Pop, sstart, send, qlen, stop])
                    print("No sufficient hit", "\n")
                    break
                    
             #Save this as result liste. This should have all you need for the output file.          
    print(len(Result_Liste), "= Length Liste")        
    with open("TempOutfileBlast2.txt","w") as Outfile:
        Outfile.write("Seq_id,Pop,Start,End,lengthORF,stop" + "\n")
        for item in Result_Liste:
            item2 = ",".join(item)
            Outfile.write(item2  + "\n")
    print("Num of all homologs:", count, "Number of seqs with no hit:", count_empty)
    return Result_Liste
                   
#############################################################################################
#Part 2: Merge the transcribed ORFs with the suitable non transcribed ORFs
#############################################################################################
                  
#Add the orthogroup info to the file
def Add_orthogroup_Info(InfoFile):
    Info = openFile(InfoFile)
    ORFs = openFile("TempOutfileBlast2.txt")
    with open("TempOutfileBlastFinal.txt","w") as Outfile:
        Outfile.write("Orthogroup_Id,Seq_id,Pop,Start,End,lengthORF,stop" + "\n")         
        for item in ORFs:
            for i in Info:
                if item.split(",")[0].split("::")[0] == i.split(",")[1]:
                    Outfile.write(i.split(",")[0] + "," + item)

#Add the sufficient hits to the file
def Merge_ORF_and_homologs(HomFile, ORFFile):
    ORF = openFile(ORFFile)
    Hom = openFile(HomFile)
    Save = [] 
    count = 0
    with open("Merged_File_HOM_ORF.txt", "w") as Outfile:
        Outfile.write("Orthogroup_ID" + "\t" + "Populations" + "\t"  + "Lengths" + "\t" +  "Start" + "\t" +  "End" + "\t" +  "Change" + "\t" + "Transcription_status" + "\n" )
        for line in ORF[1:]:
            l = line.strip().split("\t")
            Orthogroup_ID = l[0]
            Populations = l[1]
            Lengths = l[2]
            Start = l[3]
            End = l[4]
            n = len(Populations.split(","))
            #print(Populations)
            lol = ["Transcript"] * n
            #print(n)
            Transcription_status = ",".join(lol)
            for line2 in Hom:
                l2 = line2.strip().split(",") 
                ID = l2[0]
                Pop = l2[2]
                Len = l2[5]
                Star = l2[3]
                En = l2[4]
                if En == "no_hit":
                    continue
                if Orthogroup_ID == ID:
                    Lengths = Lengths + "," + Len
                    Populations = Populations + "," + Pop
                    Start = Start + "," + Star
                    End = End + "," + En
                    Transcription_status = Transcription_status + "," + "Not_transcribed"
            if len(set(Lengths.split(","))) == 1:
                Change = "No_change"
            else:
                Change = "Change"
            Save.append([Orthogroup_ID, Populations, Lengths, Start, End, Change, Transcription_status])
            Outfile.write(Orthogroup_ID + "\t" + Populations + "\t"  + Lengths + "\t" +  Start + "\t" +  End + "\t" +  Change + "\t" + Transcription_status + "\n" )

#############################################################################################
#Run in python
#############################################################################################
File = openFile(Bfile)    
ListeBestHits = access_best_hit_rsts(File)
Liste = translate_best_hit(ListeBestHits)
Result = Run_blast(Liste, FastaORF)
Add_orthogroup_Info(Info)
Merge_ORF_and_homologs("TempOutfileBlastFinal.txt", OFile)



