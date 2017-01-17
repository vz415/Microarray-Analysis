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

# Updating the BT and mol.bio variables to factors
ALLb$BT <- factor(ALLb$BT)
ALLb$mol.bio <- factor(ALLb$mol.bio)

#Pre-processing steps done.
# ------------------------------------------------
# Exploring the dataset
es <- exprs(ALLb)
dim(es)
# Matrix has 12,625 rows (genes/features) and 94 columns (samples/cases)

# In terms of dimensionality, there are way too many variables for the number of cases. First goal will be to eliminate some of the
# genes from our analysis to reduce dimensionality. To do this, we will explore the expression levels in the data.
summary(as.vector(es))

# Better overview can be graphically achieved with 'genefilter' package. Install first.
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
library("genefilter")
hist(as.vector(es), breaks = 80, prob=T,
     xlab = 'Expression Levels',
     main = 'Histogram of Overall Expression Levels')
abline(v=c(median(as.vector(es)),
           shorth(as.vector(es)), # wut does this mean.
           quantile(as.vector(es), c(0.25, 0.75))),
       lty=2, col=c(2,3,4,4))
legend('topright', c('Median','Shorth','1stQ','3rdQ'),
       lty=2,col=c(2,3,4,4))
# Dotted lines show the median, first and third quartiles, and the shorth. A shorth is the mean of the values in the shortest half.
# Shorth is the midpoint of that half, which is the least median of squares estimate of location, and length of the shortest half.
# Shorth is a robust estimator of the centrality of a continuous distribution that is implemented by the function shorth().

# Now want to know, are teh distributions of the gene expression levels of the samples eith a particular mutation different from
# each other? 
sapply(levels(ALLb$mol.bio), function(x) summary(as.vector(es[,
                                                              which(ALLb$mol.bio == x)])))
# Things, whatever that means, are fairly similar across these subsets of samples and, moreover, they are similar to the global 
# distribution of expression levels... Tryng to say this sample is a fair representation of the population?

#----------------------------------------------------------
# Gene (Feature) Selection
