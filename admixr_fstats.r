#!/usr/bin/env Rscript

pops1and2=$1

library(admixr)
library(readr)


dataset_eigs55=eigenstrat(geno = "data_1240k.geno", snp = "data_1240k.snp", ind = "data_1240k.ind")

popsmerged=data.frame(read.table(pops1and2))


f4_v0=f4(W=popsmerged[,1] , X=popsmerged[,1] , Y=popsmerged[,1] , Z= "YRI",data=dataset_eigs55)

write.table(f4_v0,"f4_v0_results", quote=F)


