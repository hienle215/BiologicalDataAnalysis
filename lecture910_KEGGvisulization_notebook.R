#### PLOTTING GRAPHS IN R

### A simple graphNEL object
#letäs load th graph package first
library(graph)

# Read about the graphNEL class
?graphNEL

#Let's take a look at some R constants to learn about LETTERS
?Constants


# Creating a simple example with 4 nodes
V = LETTERS[1:4]
V
E = list(c(LETTERS[2:4]), LETTERS[c(1,3,4)], LETTERS[c(1,2,4)], LETTERS[c(1:3)])
E
names(E) = V
E
G = graphNEL(V, E, edgemode = "undirected")
G
plot(G)


### An alternative
# using simple indices to describe the connections instead of the node names
V = LETTERS[1:4]
E=list()
V
E
conns = rbind(c(2:4), c(1,3,4), c(1,2,4), c(1:3)) # notice that all vectors of the sam length
conns
for(i in 1:4) {
  E[[i]] = list(edges=conns[i,])
}
names(E) = V
E # notice that E is list of lists
G = graphNEL(V, E, edgemode ="undirected")
G
plot(G)


#### Plotting GO terms as graphNEL
# Lets draw a graph showing the GO term GO:0000001 and all of its ancestors, as well as the connections between these terms
# First, lets find the term and their connectivity

library(GO.db)
ll = as.list(GOBPANCESTOR["GO:0000001"]) # all ancestors
ll

V = c("GO:0000001", ll[[1]]) # Include the term itself of our intersted GO
V = head(V, -1) # exclude "all" which is the last term
V

E = as.list(GOBPPARENTS[V])
E

bp_idx = grep("all", E) # Find "BP" wih parent "all
E[[bp_idx]] = character(0) # remove "all" connection

G = graphNEL(V, E, edgemode="directed") # "directed" in here will show the trend of interaction
G
plot(G)

G1 = graphNEL(V, E, edgemode = "indrected")
G1
plot(G1)

# Let's turn the graph such that more specific terms are at a lower level. This can be done by inverting the direction of edges in G
library(topGO)
G_reverse = reverseArch(G)
plot(G_reverse)


### CHANGING NODE LABELS
V2 = select(GO.db, keys = V, columns = "TERM")$TERM
V
E2=list()
E
for(i in 1:length(E)){
  E2[[i]] = list(edges=match(E[[i]], V))
}
names(E2) = V2
G2 = graphNEL(V2, E2, edgemode="directed")
G2
plot(reverseArch(G2))

# An alternative way to above is to use remaining capability built-in the graph class
V2 = select(GO.db, keys = V, columns = "TERM")$TERM
G3=G
nodes(G3) = V2
plot(reverseArch(G3))

library(Rgraphviz)
plot(reverseArch(G2), attrs=list(node=list(fontsize=80)))
plot(reverseArch(G3), attrs=list(node=list(fontsize=80)))

# Individual nodes have attributes too, and these can be set using nodeAttrs
nAttrs = list()
nAttrs$shape = rep("rect",16)
nAttrs
names(nAttrs$shape)=V2
nAttrs$shape

nAttrs$height = rep(5,16)
names(nAttrs$height) = V2
nAttrs$height

nAttrs$width = nchar(V2)*2 # nchar is used to count the number of characters
nAttrs$width
?nchar
V2
names(nAttrs$width)=V2
plot(reverseArch(G2), attrs=list(node=list(fontsize=10)), nodeAttrs = nAttrs)

# lets' give some arbirary color to the Graph nodes
nAttrs$fillcolor = rgb(seq(1, 0.5, length.out=16), 0, 0)
names(nAttrs$fillcolor) = V2
nAttrs$fillcolor
plot(reverseArch(G2), attrs=list(nodes=list(fontsize=10)), nodeAttrs=nAttrs)
plot(reverseArch(G2), attrs=list(node=list(fontsize=10)), nodeAttrs=nAttrs)

nAttrs$fillcolor = rgb(seq(1, 0.5, length.out=16), 0, 0)
nAttrs$fillcolor
names(nAttrs$fillcolor) = V2
plot(reverseArch(G2), attrs=list(node=list(fontsize=10)), nodeAttrs=nAttrs)


# another for drawing pathway enrichment
library(genefilter) # loading the library for rowttests function
library(multtest) # loading the data golub
data(golub) #loading packing golub

tstats = rowttests(golub, factor(golub.cl))
tstats # tstats is 
?rowttests
factor(golub.cl)
?factor # the function factor is used to encode a vector as a factor (the term's category and enumered type are also used for factors. If argument ordered is TRUE, the factor levels are assumed to be ordered)

golubGenes =tstats[,3]
golubGenes
names(golubGenes) = golub.gnames[,3]
head(golubGenes)

pvals=tstats[,3]
geneList=golub.gnames[head(order(pvals),50),3]
geneList
res = vector(mode="logical", length=length(geneList))
res

golubSel = function(geneList) {
  res = vector(mode ="logical", length = length(geneList))
  idx = head(order(geneList), 10)
  res[idx] = TRUE
  return(res)
}

golubSel # in here, known as a new function

golubGOdata <- new("topGOdata",
                   description = "Example analysis", 
                   ontology = "BP",
                   allGenes = golubGenes, 
                   geneSel = golubSel,
                   nodeSize = 10,
                   annot = annFUN.db, 
                   affyLib = "hu6800.db")

res <- runTest(golubGOdata, algorithm="classic", statistic="fisher")
?runTest
Gres = showSigOfNodes(golubGOdata, score(res), firstSigNodes=7, useInfo="all") # This plots top 7 terms; no need to use plot()

#let's modify gress e.g to change the font size. Before doing that, let's learn about an operation that we need to use. Pay special attention to the @ operator
?slot # these functions return or set information about the individual slots in an object

Gr =Gres$complete.dag #this is an Ragraph object
for(i in 1:length(AgNode(Gr))){
  AgNode(Gr)[[i]]@txtLabel@labelFontsize =22
  AgNode(Gr)[[i]]@rWidth = as.integer(30)
  AgNode(Gr)[[i]]@lWidth = as.integer(40)
  AgNode(Gr)[[i]]@height = as.integer(40)
}
plot(Gr)

### PATHWAY ENRICHMENT
# connecting example data with KEGG
# Skim through the following help pages as we are going to use them to extract more information about Golub data set in the next code chunk and the upcoming code chunks
?hu6800ENTREZID # is an R object that provides mappings between manufacturer identifiers and Entrez Gene identifiers
#?hu6800PATH
#?hu6800SYMBOL

# To connect our example Golub dataset with KEGG information, we need to provide gene identifers in a format known to KEGG
# at the same time, we can directly find KEGG pathway identifiers from the annotation package for each probe id (what is golubPath like?)

library(multtest)
data(golub)
library(hu6800.db)
golubEntrez = as.list(hu6800ENTREZID[golub.gnames[,3]])
golubPath = as.list(hu6800PATH[golub.gnames[,3]])
mean(is.na(golubEntrez)) # Probes with no ENtrez ID
#as we can see, there are "not available"-values for some probe id:s, meaning that the data from these will be lost in the example analysis

#let's build a new data matrix containing only the probes for which we have the Entrex id
unlist(golubEntrez[1:3]) 

# the same id occurs multiple times. For this example, let's just take the first of each
#let's build a new data matrix containing only the probes fir which we have the Entrez id
selectedIdx = !is.na(golubEntrez) & !duplicated(golubEntrez) # in here, we also remove NA data as well as the duplicated gene id
selectedIdx # including the interested names ID
golub2 = golub[selectedIdx,]
golub2 # including the expression level
rownames(golub2) = golubEntrez[selectedIdx] # setting the names f golub2 based on idx of gene names in golubEntrez
dim(golub2) # here, dim shows the size of data matrix
golubPath = golubPath[selectedIdx]

# let's use the R aggregate function to calculate a new, summarized data matrix from golub data matrix
# in genral, if we have several measurements from the same gene, we may want to take information from all of them and calculate e.g sum, means, or medians for each samples
# such a process is called summarization. In R, we can e.g use aggregate to calculate a new, summarized data matrix from our familar golub
agg_variable = list(unlist(golubEntrez))
golub3 = aggregate(by = agg_variable, x = golub, FUN = median)
dim(golub3) # 1st col is Entrez id
rownames(golub3) = golub3[,1] # name the rows
golub3 = golub3[,-1] # remove the column
golub3["6772", 1:5]

# we can check the result by doing this for the same Entrez id in a simpler way
indexvec = grep("6772", golubEntrez)
indexvec
res = apply(golub[indexvec,],2,median)
res[1:5]

# The above matches aggregate results, but doing this we only got on id at a time and would have to do a for loop to perform the whole task
# for finidng pathways of interest based on our data, we can perform a similar enrichment analysis as we did for GO term. Let's redo a simple t-test again to use a gene list of top 50 genes from there
library(genefilter)
tstats = rowttests(golub, factor(golub.cl))
pvals = tstats[,3]
genelist = golub.gnames[head(order(pvals), 50), 3]                   
head(genelist)

# NExt, we will calculate the numbers of annotations to different KEGG termsn in total and from our gene list
ugolubPath = unlist(golubPath) # palin vector
N = sum(!is.na(ugolubPath)) # total annotation
anns = table(ugolubPath) # quick search n:s for each pathway # get out the number of n of each pathway, n is the vector of annotations
anns2 = rep(0, length(anns)) # for the k:s # get out the k values in here
res = table(unlist(golubPath[genelist])) # int.genes
M = sum(res) #total annotation from gene list
names(anns2) = names(anns)
anns2
anns2[names(res)] = res # note here, res is vector including all the number of gene here
# we can now look at the most highly enriched KEGG pathways and study selected ones in more detail
# we will do it in one-tailed test
head(sort(1-phyper(anns2-1, M, N-M, anns)),4)

# let's have a closer look at KEGG pathway 04062: chemokine signaling pathway. Chemokines are signalling molecules related to immuno_response_and_cell_migration


#### VISUALIZING GOLUB DATA ON A KEGG PATHWAY
# let's calculate the difference in median expression levels for each gene of the pathway we have from the Golub data
golDiff = apply(golub2[, 1:27], 1, median) - apply(golub2[, 28:38], 1, median)
hist(golDiff)
# let's use the Bioconductor package pathview to draw the pathway with colors indicating changes between AML and ALL. After running the command, we have the original 
#...pathway figures as file hsa04062.png and  version with Golub data as hsa04062.golDiff.png
library(pathview)
pv.out = pathview(gene.data = golDiff,
                  pathway.id = "04062", # KEGG ID for chemokine signaling pathway
                  species = "hsa", #id for human
                  limit = list(gene=2, cpd=1), # scale cpd for compounds
                  out.suffix = "golDiff", # output file prefix
                  kegg.native = TRUE) # we want a native KEGG-type graph

# Note that the above code chunk save the .png results on the storage devide. Here, we have manually added the result
# A number of differentially expressed genes can be seen, clustered in two locations on the pathway. Red means higher in ALL, green higher in AML
# for comparison, let's list the pathway genes with the smallest p-values
chkgenes = grep("04062", golubPath)
chkgenes
chkgenes = names(golubPath)[chkgenes]
chkgenes
chkgene_idx = match(chkgenes, golub.gnames[,3])
chkgene_idx
chkpvals = pvals[chkgene_idx]
chkpvals
names(chkpvals)=as.list(hu6800SYMBOL[chkgenes])
chkpvals
topchk=chkpvals[head(order(chkpvals), 10)]
topchk
# many are chemokines such as IL8, CCL23, PPBP, CXCL2, LYN is at Src (tyrosine-protein kinase), CXCR4 is a chemocine receptor, PRKCD at PKC. In case of several genes per mode, pathview adds the values together
# let's add information on the change of median levels for these genes on the list
golDiff2 = apply(golub[chkgene_idx, 1:27],1,median) - apply(golub[chkgene_idx, 28:38],1,median)
golDiff2
order_ind = head(order(chkpvals), 10)
order_ind
topchk = rbind(topchk, golDiff2[order_ind])
topchk
rownames(topchk) = c("pval", "diff")
topchk
