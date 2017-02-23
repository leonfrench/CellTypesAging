load("/Users/lfrench/Desktop/results/CellTypesAging/results/AW.Tasic.emperical.iter.1.TranscriptomicName.1479920613/aging_genes_aw_fisher_data.RData",v=T)

aw_result_mouseFilter$rank <- rank(-1*aw_result_mouseFilter$agingPValuesWithDirection, ties.method = "random")

#get middle percentile for p=0
1-(aw_result_mouseFilter %>% mutate(absRank = abs(agingPValuesWithDirection)) %>% arrange(absRank))[1, "rank"]/nrow(aw_result_mouseFilter)
aw_result_mouseFilter %>% filter(round(nrow(aw_result_mouseFilter)/2) == rank)


load("/Users/lfrench/Desktop/results/CellTypesAging/results/Mistry.Zeisel.emperical.iter.1.TranscriptomicName.1479938797/aging_genes_aw_fisher_data.RData",v=T)

aw_result_mouseFilter$rank <- rank(-1*aw_result_mouseFilter$agingPValuesWithDirection, ties.method = "random")

#get middle percentile for p=0
1-(aw_result_mouseFilter %>% mutate(absRank = abs(agingPValuesWithDirection)) %>% arrange(absRank))[1, "rank"]/nrow(aw_result_mouseFilter)
aw_result_mouseFilter %>% filter(round(nrow(aw_result_mouseFilter)/2) == rank)


load("/Users/lfrench/Desktop/results/CellTypesAging/server results/Blood.Tasic.emperical.iter.10000.TranscriptomicName.1479930030/aging_genes_aw_fisher_data.RData",v=T)

aw_result_mouseFilter$rank <- rank(-1*aw_result_mouseFilter$agingPValuesWithDirection, ties.method = "random")

#get middle percentile for p=0
1-(aw_result_mouseFilter %>% mutate(absRank = abs(agingPValuesWithDirection)) %>% arrange(absRank))[1, "rank"]/nrow(aw_result_mouseFilter)
aw_result_mouseFilter %>% filter(round(nrow(aw_result_mouseFilter)/2) == rank)
