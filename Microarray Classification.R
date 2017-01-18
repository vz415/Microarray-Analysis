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
# Generally, there are two types of approaches to feature selection: 1. Filters 2. Wrappers
# Filter approaches are carried out in a single step and use statistical properties of the features to select the final set.
# Wrappers typically involve a search process where we iterativel look for the set of featuresthat is more adequate for the
# data mining tools that are applied.==> Feature wrappers have a clear overhead in terms of computational resources....
# Feature Wrappers Steps: Running the full filter -> model -> evaluate cycle several times until some CONVERGENCE CRITERIA are met.

# First gene filtering methods to described based on info concerning distribution of expression levels.
# This will show which genes are not expressed at all or show very small variability. The latter means that this cannot be used
# to differentiate among samples. 

# We can get an overall idea of the dist. of the expression levels of each gene across all individuals with the following graph.
# Graph will use the median and IQR as representatives of these distributions.
rowIQRs <- function(em){
  rowQ(em,ceiling(0.75*ncol(em))) - rowQ(em, floor(0.25*ncol(em)))
}    # Created a custom function for IQRs
# rowMedians() function obtains a vector of medians per row. This gets the median expression level of each gene.
plot(rowMedians(es), rowIQRs(es),
     xlab='Median expression level',
     ylab='IQR expression level',
     main='Main Characteristics of Gene Expression Levels')
# We can see that a large proportion of the genes have very low variability (IQRs approaching 0). 
# If hene has a very low variability across all samples, then it is reasonably safe to conclude that it won't be useful 
# in discriminating the B-Cells. Caveat is, they are dependent, and, in this case, use the RELIEF method.

# Removing any genes with variability < 1/5 --> Just a heuristic threshold. nsFilter() from library(genefilter) can do this filtering.
library(genefilter)
ALLb <- nsFilter(ALLb,
                 var.func=IQR,
                 var.cutoff=IQR(as.vector(es))/5,
                 feature.exclude="^AFFX")    # What dis do?
# Outputs a cool filter log that says: 1. how many duplicates were removed. 2. number of low variance genes removed. 3. features removed
# 4. and numRemoved ENTREZID, whatever that means. 
# Were able to reduce to 3,942 genes from initial ~12k. Still have a ways to go, but a good start.
ALLb

# Update original dataset
ALLb <- ALLb$eset
es <- exprs(ALLb)
dim(es)

# ANOVA Filters
# If a gene has a distribution of expression values (say, the mean) that is similar across all possible values of the target variable,
# then it won't be useful to discriminate among these variables. i.e. Genes which have high statistical conf of having the same mean
# expression level across the roups of samples belonging to each mutation will be discarded from further analysis.

# We have four groups of cases, one for each of the gene mutations of B-cell ALL we are considering. Carry out as follows:
f <- Anova(ALLb$mol.bio, p = 0.01) # dat p value. mol.bio is the subgroup used in the analysis.
ff <- filterfun(f) # Creates the filter from the ANOVA
selGenes <- genefilter(exprs(ALLb), ff)
# Gets down to 746!
ALLb <- ALLb[selGenes,]
ALLb

es <- exprs(ALLb)
plot(rowMedians(es), rowIQRs(es),
     xlab='Median expression level',
     ylab='IQR expression level',
     main='Distribution properties of the Selected Genes')    # Muuuuuuch better.
# At this point, could normalize, but some modeling techniques are affected by this result, so abstain for the moment.

#---------------------------------------------
# Filtering Using Random Forests
# Random forests can be used to obtain a ranking of the features in terms of their usefulness for a classification task.