### This notebook will cover gene ontology, accessing the GO databased, annotations, enrichments analysis, enrichment package in R
### We consider functional annotation, an atuomatical interpretation of gene lists producted by the analyses
# as of ontologu version 1.2854, dated 06.04.2012, there were 
# 22506 biological process terms
# 2980 cellular component terms
#9341 molecular function terms
#1630 obsolete terms

### GO OBJECT IN R
### ACCESSING the GO databased in R
#### GO databared is available in R though package GO.db
# this is an examples of a more general database object, which can store various type of information

library(GO.db)
GO.db # to see the version of GO

### GO OBJECT KEYS
### to access the information in the databased, we can use one of several keys. For some databased, there may be several key types available,
#the list of key types can be obtained with keytypes() functions.
length(keys(GO.db)) # the number of individual GO terms
head(keys(GO.db),4)

### GO OBJECT COLUMNS
### for each key, information corresponding to a number of columns is stored
columns(GO.db)
#we can use the select command to extract information from the databased object
keylist <- head(keys(GO.db),3)
select(GO.db, keys=keylist, columns="TERM") # the name of regulated pathways
select(GO.db, keys=keylist, columns="DEFINITION") # definition of regulated pathways
select(GO.db, keys=keylist, columns="GOID") # GO ID
select(GO.db, keys=keylist, columns="ONTOLOGY") # ontology information i.g MF, molecular function; BP, biological process; CC, cellular components

### GO TERM COUNTS
#let's check the total number of terms in the ontology version in GO.db.
allterms <- select(GO.db, keys=keys(GO.db), columns="ONTOLOGY")
allterms=select(GO.db,keys=keys(GO.db),columns="ONTOLOGY")
table(allterms[,2])
allterms[allterms[,2]=="universal",]

### GO TERM RELATIONSHIP
# GO.db can be used to find offspring, children (direct offspring), ancestr or parents (direct ancestors)
# To access all ancestors of GO:0000002
as.list(GOBPANCESTOR["GO:0000002"])
# notice that the output is a list such that if we give several ontology terms as the input, each list elements will be a vector containing
#the ancestors of single input term

ll=as.list(GOBPANCESTOR["GO:0000002"])
select(GO.db,keys=ll[[1]],columns="TERM")
# To list all the functions in a package we can use ls
ls("package:GO.db")


### ANNOTATIONS
#the annotations connect genes or gene products to the term of the ontology (and give other information on what was measured)
# a collection of annotations can be donwloaded from the ontology site or accessed outsource online brower or through databased package in R
# each annotated item is typically associated with several terms
# annotations are based on varied sources of information, the reliability of which needs to be considered

### ANNOTATIONS IN R
library(hu6800.db)
length(keys(hu6800.db))
head(keys(hu6800.db))
columns(hu6800.db) # GO contains the references to gene ontology terms and notice the number of references to external databases
# the annotation also contains information on the location of the measured items in the genome (CHR, CHRLOC, MAP)
# We can again use select to extract data as needed
xx = select(hu6800.db, keys=keys(hu6800.db)[20], columns="GO")
xx
yy = select(GO.db, keys=xx[["GO"]], columns = "TERM")
yy


### ENRICHMENT ANALYSIS
#integrating proteomic information form different organisms and assigning functions to protein domains
# verifying models of genetic, metabolic, and protein interaction networks
#developming autmotated ways of deriving information about gene function from the literature
#finding functional similarities in gens that are show up as a list of results in differential expression or clustering study
# for example, let's find all annotations associated with GO:0022617 from the genes of the Golub data matrix
library(multtest)
data("golub")
library(GO.db)
library(hu6800.db)
lo <- as.list(GOBPOFFSPRING["GO:0022617"])[[1]] # object change to list of data
lo <- c("GO:0022617", lo) # add the term itself. We build lo by converting from oject to a list and then extracting the first and only component of the list
al <- select(hu6800.db, keys=lo, columns="PROBEID", keytype="GO")$PROBEID
length(al)

# example of a single GO term
intersect(al, golub.gnames[,3])[1:16] # show the 16 first of rows in data matrix
length(intersect(al, golub.gnames[,3])) #30: here we find 30 annotation to our GO term, ni=30. Notice that we use intersect to find the overlap between annotated probe names and Golub matrix probe names.
# this function assumes that each items on the vectors is unique.

# to check if any of the probes is annotated to more than one term in the vector, we can alternatively use table to get the count
tl <- table(al)
sum(tl[golub.gnames[,3]], na.rm=TRUE) # indexing the table with a probe name which is not one of the names in the table produces a NA (no available) element.
#na.rm parameter makes sume drop all these NA elements and add all the real counts corresponding to Golub matrix gene names. As a result, we see that there is a single probe with two annotations.

# a single GO term
m_ind=which.max(tl[golub.gnames[,3]])
m_ind # This probe appears twice, max==2
mi_anns=select(hu6800.db,keys=names(m_ind),columns="GO")
go_ind= intersect(mi_anns$GO,lo) #Which are the two terms?
go_ind
select(GO.db,keys=go_ind,columns="TERM")
as.list(GOBPOFFSPRING["GO:0010716"]) #???
as.list(GOBPOFFSPRING["GO:0010715"]) #???
mi_anns$EVIDENCE # tell us that 1 is TAS (author statement), 2 ISS (less reliable). The annotations thus represent different reliability levels of information


### how manny of the 31 annotations to GO:0022617 are on the list of top 10 genes withe the smallest t-test and p-value
library(genefilter)
tstats = rowttests(golub, factor(golub.cl))
hh = head(order(tstats[,3]), 10)
hh
golub.gnames[hh,3]
intersect(golub.gnames[hh,3],al)


### TOTAL NUMBER OF ANNOTATIONS
# to quantify the level of enrichment, we also need the total amount of annotations N and the total amount of annotations from within the top 10 genes, M. 
# let's focus here on BP ontology annotations (GO:0022617 is BP: biological process)
ss=(select(hu6800.db,keys=golub.gnames[hh,3],
           keytype="PROBEID",columns="GO")$ONTOLOGY)
length(ss)      
ss = (select(hu6800.db, keys=golub.gnames[,3], keytype = "PROBEID", columns="GO")$ONTOLOGY) # 269
length(grep("BP", ss)) # 37791
# this means that in our sample, we have a total of N = 37791 annotations, out of which M =269 are associated with the genes on our top list. Given this, is ki=3, out of ni=31 a lot or not?


### HEYPERGEOMETRIC DISTRIBUTION
#to continue with our example GO:0022617, let's see what the hypergeometric distribution looks like here. Notice that in R, the distribution is parametrized slightly differently, with N, N-M, ni, and ki
dd = dhyper(0:30, 269, 37791-269, 30)
dd
plot(0:30, dd) # notice that the proportion of "intersting" annotations is so small that the distribution is strictly decreasing. This is not necessarily the case in general

# In the cased of a two tailed test, we would like to calculated the probability mass in both tails of the distributions
# a simple and commonly used way of defining the tails for the asymmetrix hypergeometric distribution is to define "more extreme" by all points at which the probability density is lower than or equal to the observed point
sum(dd[dd<=dd[4]])
sum(dhyper(3:30, 269, 37791-269, 30))
1-phyper(2,269,37791-269,30)


### TESTING USING fisher.test()
# To calculate the same p-value with fisher.test() we need a contigency table
mmm <- matrix(nrow=2, c(37791-269-(30-3), 269-3, 30-3,3))
rownames(mmm) = c("Not significant", "Significant")
colnames(mmm) = c("Other terms", "GO:0022617")
mmm
# In other words, Fisher's test is now testing whether GO:0022617 annotations are more likely to come from the significant genes than the other annotations. or Whether the row and column variables are independent
fisher.test(mmm)
# in this case, it seems like the "oods ration"is not equal to 1. That is, the probability of observing a significant annotation is higher than backgroun levels if we are looking G=:0022617.
#GO:002617 is therefore enriched in our gene list



### ENRICHMENT ANALYSIS IN R
# Fo Go enrichment analysis in R, we use the package topGO which is containing the gene list, annotations, and ontology.
library(multtest)
library(genefilter)
data("golub")
library(topGO)

tstats = rowtests(golub, factor(golub.cl))
golubGenes = tstats[,3]
names(golubGenes) = golub.gnames[,3]

golubGOdata <- new("topGOdata",
                   description = "Example analysis", ontology = "BP",
                   allGenes = golubGenes, geneSel = golubSel,
                   nodeSize = 10,
                   annot = annFUN.db, affyLib = "hu6800.db")

golubSel = function(geneList)
{
  res=vector(mode ="logical",length = length(geneList));
  idx=head(order(geneList),10);
  res[idx]=TRUE;
  return(res);
}

source("golubSel.R"); # What does this do?
source("testEnrichment.R")
res <- runTest(golubGOdata
               statistic = "fisher")
res
pvals=score(res)
pp=pvals[head(order(pvals),4)]
pp
