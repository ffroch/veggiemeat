#ONLY-FIRST-TIME - data structure setup
#source("A_MakeDataStructure.R")

#load libraries
source("B_LoadLibraries.R")

#import and prepare data
source("C_PrepareData.R")



#OPTIONAL-doing some general read distribution stats
#source("D0_SequencingStats.R")


#OPTIONAL-rarefaction plots with iNEXT
#source("E0_AlphaDiv_iNEXT.R")



#we decided to rarefy based on coverage (option 1) - the second option was size based (option 2)
#for option 1 we set a coverage of 99.5%. Since all 28 samples had a high coverage, we used all of them
#for option 2 we excluded sample C4 with 439 reads and rarefied all samples to 2463 reads

#OPTIONAL - alpha diversities with both options in comparison
#source("E1_AlphaDiv_Comparison.R)

#alpha diversity analysis with coverage based rarefaction
source("E2_AlphaDiv.R")


#classic alpha diversity with coverage rarefied data - Shannon, Simpson, Faith's PD 
source("E3_AlphaDiv_classic.R")


#beta-diversity is done with 100 iterations of coverage based rarefaction
source("F0_RarefactionAndDistanceMatrices.R")


#OPTIONAL-NMDS plots
#source("F1_BetaDiv_NMDS.R")



#tSNE plots
source("F2_BetaDiv_tSNE.R")

#PERMANOVA----
source("F3_BetaDiv_PERMANOVA.R")

#LEfse with relative abundance of coverage based rarified data rarec
#export for LEfse tool with bioconda see ./07-lefse/01_lefse_pipe.sh
source("G1_LEfseExport.R")

#use LEfse prepared data for group comparison including correction for muliltple comparisons
source("G2_LEfse_multicomp.R")

#OPTIONAL: hand made LEfse plots
#source("G3_LEfse_plots.R")


#ONLY THE FIRST TIME OR WHEN CHANGING QUALITY PARAMETERS processing the sanger sequenced isolates
#source("H1_IsolateSeq_Import_Trimming.R")

#ONLY THE FIRST TIME OR WHEN CHANGING QUALITY PARAMETERS call trimmed Sanger sequences
#source("H2_IsolateSeq_FastaImport.R")

#ONLY THE FIRST TIME OR WHEN CHANGING QUALITY PARAMETERS assign taxonomy to SILVA 138.1 database
#source("H3_IsolateSeq_TaxonomyR.R")

#blastcorrection 
source("I1_Blastcorrection.R")

source("I2_Isolates.R")

#same sample isolate sequences clustering
source("J1_ClusteringIsolates.R")

#plot figure 2
source("K0_Fig2_3percent_grouped_by_group.R")
