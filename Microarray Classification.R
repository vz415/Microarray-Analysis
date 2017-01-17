# These are the instructions for analyzing the microarray samples. Followed step-by-step
library(Biobase)
library("ALL")
data("ALL")
# A lot of base packages are masked after loading this package. Watch out and see if it makes a difference in your coding.

# Data consists of microarrray samples from 128 individuals with acute lymphoblastic leukemia. We'll focus on the B-cell ALL
# which consists of 95 samples.
# Data is divided into several groups. First, have the assay data with the gene expression levels matrix. For this, have about
# 12,625 genes and 128 samples. The object also contains a lot of meta-data about the samples of the experiment. This
# includes the phenoData part with information on the sample names and several associated co-variates. Also includes info on the features
# (i.e. genes) as well as annotations of the genes from biomedical databases. Finally, object also contaions information that 
# describes the experiment.

# Will start by obtaining some info on co-variates associated to each sample:
pD <- phenoData(ALL)
varMetadata(pD)
# Playing/viewing the data
table(ALL$BT)
table(ALL$mol.biol)
table(ALL$BT, ALL$mol.bio)

# First two statements obtain the names an descriptions of the existing co-variates. We then obtain some information on the 
# distribution of the samples across two main co-variates: the BT variable that determines the type of acute lympoblastic leukemia
# and the mol.bio variable that describes the cytogenetic abnormality found on each sample (NEG represents no abnormalities)

# Obtaining some info on the genes and samples:
featureNames(ALL)[1:10]
sampleNames(ALL)[1:5]
# Shows the names of the first 10 genes and the names of the first 5 samples.

# Obtaining the subset of data that we will use:
tgt.cases <- which(ALL$BT %in% levels(ALL$BT)[1:5] & # OH SHIT! LEVELS CAN ISOLATE THOSE THINGS IN A LIST!
                     ALL$mol.bio %in% levels(ALL$mol.bio)[1:4]) # HOW DID I NOT KNOW THIS!
ALLb <- ALL[,tgt.cases]
# Checks out. Has 94 samples and everything else looks about right.

