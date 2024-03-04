library(ape)
library(arules)
library(Biostrings)
#library(btools)
library(cba)
library(cowplot)
library(dada2)
library(DECIPHER)
library(dplyr)
library(dunn.test)
library(flextable)
#library(ggordiplots)
library(ggplot2)
library(ggplotify)
library(gmodels)
library(grid)
library(gtools)
library(indicspecies)
library(iNEXT)
library(kmer)
library(metagMisc)
library(officer)
library(philentropy)
library(phyloseq)
library(phyloseqCompanion)
library(plyr)
library(qiime2R)
library(RColorBrewer)
library(Rcpp)
library(readr)
library(Rtsne)
library(sangeranalyseR)
library(seqinr)
library(stringr)
library(tidyverse)
library(vegan)





#remotes::install_github("jbisanz/qiime2R")
library(qiime2R)

#remotes::install_github("twbattaglia/btools")
#library(btools)

#remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)

#used packageVersion of sangeranalyseR 1.2.0
if(!require(sangeranalyseR)){
  install.packages("sangeranalseR")
  library(sangeranalyseR)
}

knitr::write_bib(c(.packages(),"bookdown"),"packages.bib")
