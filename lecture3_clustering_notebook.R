### practicing the basic concepts with a toy data set

### Let's assume that we have 3 samples s2, s2,s3 and 2 genes g1, g2. Now we want to see how close these samples are from each other

data <- matrix(c(1,4.3,4,1,5.5,5), nrow=2, byrow=TRUE)
colnames(data) <- c("s1","s2","s3")
rownames(data) <- c("g1", "g2")
data

#let's plot the data
par(pty="s")
plot(data[1,], data[2,], xlim=c(0,6), ylim=c(0,6))
text(data[1,] + 0.2, data[2,]-0.2, labels=colnames(data))

#based on the visual inspection, we can see that s2 and s3 are the closest of the samples. However, we will use Euclindean distance to 
#to measure the distance between each pair of samples, we can use the function that calculates the Euclidean distance

euclidean_dist = function(v1, v2){
  return(sqrt(sum((v1-v2)^2)))
}

euclidean_dist(data[, 1], data[, 2]) # Euclidean distance between s1 and s2
euclidean_dist(data[, 1], data[, 3]) # Euclidean distance between s1 and s3
euclidean_dist(data[, 2], data[, 3]) # Euclidean distance between s2 and s3

# The output of euclidean_list() function supports our observation
### We can also use R dist() function to calculate the distance for us.
?dist
#dist() function finds the distance between rows(in here rows are our gens) as it can be seen from below:
dist(data, method = "euclidean", diag=TRUE)

#but we want to find the distance between samples, thus we need to transpose our dta using t() function
dist(t(data), method = "euclidean", diag=TRUE)


### After finishing the calculation of distance, we move on to learn about hierarchical clustering. In this case, hclust() function
?hclust()
hc = hclust(dist(t(data), method = "euclidean"), method = "complete")
plot(hc)
# abline(h=0.58, col="red", lty=2, lwd=0.4)
# abline(h=5.58, col="red", lty=2, lwd=0.4)
#Let's also check the type, class, and structure of the object that is retured by hclust() fucntion
typeof(hc)
class(hc)
str(hc)

### Hierarchical clustering of the Golub data set
install.packages("BioGenerics")
install.packages("parallel")
library(multtest) 
data(golub)
#use hclust() to perform hierarchical clustering overGolub data
hc <- hclust(dist(t(golub), method = "euclidean"), method = "complete") #clustering hc
hc
plot(hc) # hc is an object ready to plot

### HEATMAP
#ploting heat map for our toy data set
heatmap(data) #Do the clustering and draw the heatmap
# we can also set the locor paletter such that the difference are more visible
#heatmap(data, col=fc(20))

# heatmap of genes with high variance in the Golub data set
row_vars = apply(golub, 1, var) # Calculate row variances
Golub_sorted = golub[order(row_vars),] # Find the order
Golub_sorted_subset = Golub_sorted[2800:3051,] # Filtering
heatmap(Golub_sorted_subset) # Do the clustering and draw the heatmap
