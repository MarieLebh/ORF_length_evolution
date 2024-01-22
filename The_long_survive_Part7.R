############################################################################
#R script to generate Figures 4 and 5:
############################################################################
library(readr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(rstatix)

windowsFonts(sans = windowsFont("Helvetica"))

############################################################################
#Load and restructure the dataframe

df <- read_delim("Merged_File_HOM_ORF.txt", 
                                  delim = "\t", quote = "'", escape_backslash = TRUE, 
                                  locale = locale(grouping_mark = "."), 
                                  trim_ws = TRUE)
df$New = strsplit(df$Lengths, ",")
df$New = lapply(df$New, as.numeric)
df$MAD = lapply(df$New, mad)
df$median_length = lapply(df$New, median)
df$New = strsplit(df$Populations, ",")
df = df %>% 
  mutate(New2 = lengths(New))

#Get summary stats of MAD and median length
df$median_length = as.numeric(df$median_length)
df %>%
  group_by(New2, Change) %>%
  get_summary_stats(median_length, show = c("mean", "sd", "median"))
df$MAD = as.numeric(df$MAD)
df %>%
  group_by(New2) %>%
  get_summary_stats(MAD, show = c("mean", "sd", "median"))

#Count the number of popularions per orthogroup
table(df$New2)

############################################################################
#Fig 4A  (ORFs in orthogroup)
a = ggplot(data = df, aes(x = as.character(New2), fill = Change))+
  geom_histogram(position = "stack", stat = "count")+
  ylab("Number of \northogropus")+
  xlab("Number of lines in \n orthogroup")+
  theme_bw()+
  theme(axis.title = element_text(size = 19), 
        axis.text = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_blank())+
  guides(fill=guide_legend(ncol=2))

a = a + theme(text = element_text(family = "sans"))

#Fig4B (Median length)
b = ggplot(data = df, aes(x = as.character(New2), y = as.numeric(median_length), fill = Change))+
  geom_boxplot(alpha = 0.3, outlier.alpha = 1)+
  ylab("Median length")+
  xlab("Number of lines in \n orthogroup")+ 
  theme_bw()+
  theme(axis.title = element_text(size = 19), 
        axis.text = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_blank())+
  stat_summary(fun=mean, geom="point", 
               shape=18, size=4, show.legend=TRUE, aes(group=Change, col = Change), position=position_dodge(.78))+
  scale_y_continuous(trans = "log2")

b = b + theme(text = element_text(family = "sans"))

#Merge figures  a and b and save them
x = ggarrange(a,b, labels = c("(A)", "(B)"),  font.label= list(size = 19), common.legend = TRUE, legend = "bottom") 
ggsave("Fig4.svg", width = 20, height = 11, unit = "cm")

#GLM median length
model1 = glm(data = df, as.numeric(median_length) ~ New2*Change , family = "inverse.gaussian")
summary(model1)

#Create subsets with and without length variation
change = subset(df, Change == "Change")
no_change = subset(df, Change == "No_change")

#Variation stats
model1 = glm(data = change, as.numeric(median_length) ~ New2 , family = "inverse.gaussian")
summary(model1)

#No variation stats
model1 = glm(data = no_change, as.numeric(median_length) ~ New2, family = "inverse.gaussian")
summary(model1)

#Supplementary plot MAD
df2 = subset(df, Change != "No_change") #subset with only length variation
df2$MAD = as.numeric(df2$MAD)

#Get summary stats
df2 %>%
  group_by(New2) %>%
  get_summary_stats(MAD, show = c("mean", "sd", "median"))

a = ggplot(data = df2, aes(x = as.character(New2), y = as.numeric(MAD)))+
  geom_boxplot()+
  ylab("MAD \n (ORF length aa)")+
  xlab("Number of lines in \n orthogroup")+ 
  stat_summary(fun=mean, colour = "darkred", geom="point", 
               shape=18, size=4, show.legend=FALSE)+
  theme_bw()+
  theme(axis.title = element_text(size = 19), 
        axis.text = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 19),
        legend.position = "bottom",
        legend.title = element_blank())

a + theme(text = element_text(family = "sans"))

ggsave("Supp_figure_MAD.pdf", width = 10, height = 11, units = "cm")

#GLM with MAD
model = glm(data = df2, as.numeric(MAD) ~ New2)
summary(model)

############################################################################
#Now plot the length change trough 3' and 5' truncation

#Non transcried included
df <- read_csv("Results.txt") 
df <- subset(df, Type == "3'truncation" | Type == "5'truncation")
df$Dataset = "Nontranscribed \nincluded"

#Summary stats
df %>%
  group_by(Type) %>%
  get_summary_stats(Difference, show = c("mean", "sd", "median"))

#Only transcribed
df2 <- read_delim("Fig_9_transcribedonly.txt", delim = "\t") #Its the same df from Step1 (Check alignment ...)
df2$Comparison[grep("Same_start_end_earlier", df2$Comparison)] <- "3'truncation"
df2$Comparison[grep("Same_end_and_start_later", df2$Comparison)] <- "5'truncation"
df2$Type = df2$Comparison
df2 <- subset(df2, Type == "3'truncation" | Type == "5'truncation")
df2$Dataset = "Transcribed \nonly"

#Summary stats
df2 %>%
  group_by(Type) %>%
  get_summary_stats(Difference_long_short_aa, show = c("mean", "sd", "median"))

#Plot Fig.5C
lengths = ggplot()+
  geom_boxplot(data = df, aes(y = Difference, x = Dataset, fill = Type))+
  geom_boxplot(data = df2, aes(y = Difference_long_short_aa, x = Dataset, fill = Type))+
  scale_fill_manual(values = c("red", "green"))+
  ylab("Difference \nlong vs short ORF")+
  xlab("Category")+
  theme_bw()+
  guides(fill=guide_legend(nrow = 2))+
  theme(axis.title = element_text(size = 19),
        strip.text = element_text(size = 19),
        axis.text = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 19),
        legend.position = "bottom",
        legend.title = element_blank())

lengths + theme(text = element_text(family = "sans"))

ggsave("Fig5C.pdf", width = 10, height = 11, units = "cm")


############################################################################
#Now plot Figure 5 A and B

#Search for short STOP/START in long ORF
shortORF_df <- read_delim("Stop_neighbour_info_long_short_indels.txt", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)

shortORF_df$DF = "Transcribed only"

#Get counts per category
shortORF_df %>%
  group_by(Type,Neighbour_present2) %>%
  dplyr::summarise(count=n())

shortORF_df %>%
  group_by(Type,Neighbour_present) %>%
  dplyr::summarise(count=n())


df1 <- subset(shortORF_df, Type == "Same_start_end_earlier")
df1$Type[grep("Same_start_end_earlier", df1$Type)] <- "3'truncation"
df1$Neighbour_present[grep("No_Stop_or_neighbour", df1$Neighbour_present)] <- "Codon_missing"
df1$Codon_present = df1$Neighbour_present

df2 <- subset(shortORF_df, Type == "Same_end_and_start_later")
df2$Type[grep("Same_end_and_start_later", df2$Type)] <- "5'truncation"
df2$Neighbour_present2[grep("No_Start_present", df2$Neighbour_present2)] <- "Codon_missing"
df2$Neighbour_present2[grep("Start_codon_present", df2$Neighbour_present2)] <- "Start"
df2$Neighbour_present2[grep("Start_neighbour_present", df2$Neighbour_present2)] <- "Start_neighbour"
df2$Codon_present = df2$Neighbour_present2

#Transcribed and non transcribed
shortORF_df2 <- read_delim("CodonInfo.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
#Get counts per category
shortORF_df2 %>%
  group_by(Type,Neighbour_present2) %>%
  dplyr::summarise(count=n())

shortORF_df2 %>%
  group_by(Type,Neighbour_present) %>%
  dplyr::summarise(count=n())

shortORF_df2$DF = "Nontranscribed included"

df5 <- subset(shortORF_df2, Type == "Same_start_end_earlier" | Type == "3'truncation")
df5$Type[grep("Same_start_end_earlier", df5$Type)] <- "3'truncation"
df5$Neighbour_present[grep("No_Stop_or_neighbour", df5$Neighbour_present)] <- "Codon_missing"
df5$Codon_present = df5$Neighbour_present

df6 <- subset(shortORF_df2, Type == "Same_end_and_start_later" | Type == "5'truncation")
df6$Type[grep("Same_end_and_start_later", df6$Type)] <- "5'truncation"
df6$Neighbour_present2[grep("No_Start_present", df6$Neighbour_present2)] <- "Codon_missing"
df6$Neighbour_present2[grep("Start_codon_present", df6$Neighbour_present2)] <- "Start"
df6$Neighbour_present2[grep("Start_neighbour_present", df6$Neighbour_present2)] <- "Start_neighbour"
df6$Codon_present = df6$Neighbour_present2

#Figure 5A)
merge = (ggplot(data = df1, aes(x = Type, fill = Codon_present))+
           geom_histogram(position = "stack", stat = "count")+
           geom_histogram(data = df2, position = "stack", stat = "count")+
           geom_histogram(data = df5, position = "stack", stat = "count")+
           geom_histogram(data = df6, position = "stack", stat = "count")+
           facet_wrap(~ DF, nrow = 2)+
           scale_fill_manual(values = c("darkgrey","darkgreen", "lightgreen", "darkred", "#4d6b53"))+
           ylab("Number long vs short \nORF comparisons")+
           xlab("Category")+
           theme_bw()+
           guides(fill=guide_legend(ncol=3, nrow = 2))+
           theme(text = element_text(size = 19),
                 axis.title = element_text(size = 19),
                 axis.text = element_text(size = 15, colour = "black"),
                 legend.text = element_text(size = 12),
                 legend.position = "bottom",
                 legend.title = element_blank()))

ggsave("Fig5A.pdf", width = 10, height = 11, units = "cm")

############################################################################
#Figure 5B)

#Search for long STOP/START in transcript of short
transcr_df <- read_delim("Check_for_Stop_Codon_transcript_version2.csv", 
                  delim = ",", escape_double = FALSE, 
                  trim_ws = TRUE)

#Get counts per category
transcr_df %>%
  group_by(Type,Start_present) %>%
  dplyr::summarise(count=n())

transcr_df %>%
  group_by(Type,Stop_present) %>%
  dplyr::summarise(count=n())


df3 <- subset(transcr_df, Type == "Same_start_end_earlier")
df3$Type[grep("Same_start_end_earlier", df3$Type)] <- "3' trunctuation 2"
df3$Stop_present[grep("No_Stop_or_neighbour", df3$Stop_present)] <- "Codon_missing"
df3$Stop_present[grep("Stop_neighbour", df3$Stop_present)] <- "Stop \nneighbour"
df3$Codon_present = df3$Stop_present

df4 <- subset(transcr_df, Type == "Same_end_and_start_later")
df4$Type[grep("Same_end_and_start_later", df4$Type)] <- "5' trunctuation 2"
df4$Start_present[grep("No_start_neighbour_present", df4$Start_present)] <- "Codon_missing"
df4$Start_present[grep("Search_disrupted_by_TSS", df4$Start_present)] <- "Transcript \ntruncation"
df4$Start_present[grep("Start_codon_present", df4$Start_present)] <- "Start"
df4$Start_present[grep("Start_neighbour_present", df4$Start_present)] <- "Start_neighbour"
df4$Codon_present = df4$Start_present

#Plot Figure 5B
merge1 = (ggplot(data = df3, aes(x = Type, fill = Codon_present))+
    geom_histogram(position = "stack", stat = "count")+
    geom_histogram(data = df4, position = "stack", stat = "count")+
    scale_fill_manual(values = c("darkgrey","darkgreen", "lightgreen", "darkred", "#4d6b53"))+
    ylab("Number long vs short \nORF comparisons")+
    xlab("Category")+
    theme_bw()+
    guides(fill=guide_legend(ncol=2, nrow = 3))+
    theme(text = element_text(size = 19),
          axis.title = element_text(size = 19),
          axis.text = element_text(size = 15, colour = "black"),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_blank(),)+
      scale_x_discrete(limits=c("5' trunctuation 2", "3' trunctuation 2"),
                       labels= c("5' truncation", "3' truncation"))
    
    )

merge1 + theme(text = element_text(family = "sans"))

ggsave("5B.pdf", width = 10, height = 11, units = "cm")
