### How to do PCA in R
### How to use the prcomp() function to do PCA
### How to draw a PCA plot using based graphics and ggplot2
### How to determine how much variation each principal component accounts for 
### How to examine the loading scores to determine that variables have the largest effect on the graph

### dataset
datamatrix <- matrix(rnow=100, ncol=10)
colnames(datamatrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep="")
)# This is where we name the samples

rownames(data.matrix) <- paste("gene", 1:100, sep="") #this is where we name the genes (usually we have specific names such as Sox9, Cdk, but since this is a fake dataset, we have gene1, gene2,.... gene100)

for(i in 1:100){
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}

### We call prcomp() to do PCA on our data. The goal is to draw a graph that shows how the samples are related or not related to each other
pca <- prcomp(t(data.matrix), scale=TRUE) 

###NOTE: by default, prcomp() expects the samples to be rows and the gens to be columns
# prcomp() returns three things: x, sdev, and rotation
#x
plot(pca$x[,1], pca$x[,2]) # x contains the principal components Pcs for drawing a graph. here, we are using the first two columns in x to draw a 2-D plot that uses the first two PCs

# since there are 10 samples, there are 10 PCs
# the frist PC accounts for the most variation in the original data(the gene expression across all 10 samples), in the 2nd PCa accounts for the second most variation and so on. To plot a 2-D PCA graph, we usually use the frist 2 PCs. However, sometimes we use PC2 and PC3.
plot(pca$x[,1],pca$x[,2])

pca.var <- pca$sdev^2 # we use the square of sdev which stands for standard deviation, to calculate how much variation in the original data each principal component accounts for.
pca.var.per <- round(pca.var/sum(pca.var)*100,1) #since the percentage of variation that each PC accounts for is way more intersting than the actual value, we calculate the percentages.
barplot(pca.var.per, main="Scree plot", xlab="principal compoments", ylab="percent variation") #plotting the percentage is easy with barplot()

###gplot2
library(ggplot) #we can use the ggplot2 to make a fancy PCA plot that looks nice and also provides us with tons of information
pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2]) #make data frame
pca.data #to see the dataframe looks like. We have one row per sample. each row has a sample ID and X/Y coordinates for that sample.

ggplot(data=pca.data, eas(x=X, y=Y, label=Sample)) +
  geon_test() +
  xlab(paste("Pc1 -", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2-", pca.var.per[2],"%", sep=""))+
  theme_bw()+
  ggtitle("My PCA graph") # here the call is to draw ggplot.
loading_scores <- pca$rotation[,1]
gene_score <- abs(loading_scores)
gene_score_ranked <- names(gene_score_ranked[1:10])
top_10_genes
pca$rotation[top_10_genes,1] #show the scores and +/- sign

#then we see which genes have negative loading scores.
