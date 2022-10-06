library(vegan)
library(plyr)
library(Rtsne)
library(Rcpp)
library(tidyverse)
library(philentropy)
library(metagMisc)
library(phyloseq)
library(cowplot)
library(ggordiplots)
library(grid)
library(ggplotify)
library(indicspecies)
library(iNEXT)
library(tidyverse)
library(plyr)
library(cowplot)
library(dunn.test)
library(ggplot2)
library(ape)



library(gtools)





#remotes::install_github("jbisanz/qiime2R")
library(qiime2R)

#remotes::install_github("twbattaglia/btools")
library(btools)

#remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)

#used packageVersion of sangeranalyseR 1.2.0
if(!require(sangeranalyseR)){
  install.packages("sangeranalseR")
  library(sangeranalyseR)
}

library(Biostrings)
library(dada2)
library(stringr)
