#load Zeisel data
#get cell types with highest SST expression
#look up p-values for GO intersection - if tested

library(readr)
library(reshape2)
library(homologene)
library(ROCR)
library(optparse)
library(dplyr)

level <- "TranscriptomicName"
source("./R/ZeiselProcessingCode.R") #load data

rpkmByType <- group_by(rpkmDataCore, geneName, CellClassID) %>%summarise( log1Expression = mean(log1Expression))

#Cdk6
print.data.frame(filter(rpkmByType, geneName == "Cdk6") %>% arrange(desc(log1Expression)))
filter(rpkmByType, geneName == "Sst") %>% arrange(desc(log1Expression))


#Vip Mybpc1

print.data.frame(filter(rpkmByType, geneName == "Mybpc1") %>% arrange(desc(log1Expression)))
filter(rpkmByType, geneName == "Vip") %>% arrange(desc(log1Expression))

#oligos
print.data.frame(filter(rpkmByType, geneName == "Opalin") %>% arrange(desc(log1Expression)))

print.data.frame(filter(rpkmByType, geneName == "9630013A20Rik") %>% arrange(desc(log1Expression)))
