### Multi-dimensional scaling (MDS) in R

library(multtest)
data("golub")

#now we have Golub data loaded, let's perform MDS analysis using cmdscale() function
result = cmdscale(dist(t(golub)))
dim(result)

# now let's plot the results of MDS analysis
par(pty="s") # This is used to make the plot square 
plot(result)

#inspecting the component contributions. Note that GOF and results$GOF is short for goodness of fit
result = cmdscale(dist(t(golub)), eig=TRUE)
result$GOF
par(pty="s") # This is used to make the plot square
plot(result$eig)

# adding the class lables from the Golub data set to the resutls of MDS analysis
par(pty="s") # This is used to make the plot square 
plot(result$points, col="red")
points(result$points[28:38,], col="green")

### PRINCIPAL COMPOENT ANALYSIS -PCA in R
#In R, prcomp() function can be used to perform PCA analysis
PCA_result = prcomp(t(golub))

dim(PCA_result$rotation) # 38 basis vectors
#let's take a peek at waht is inside PCA_results$rotation
head(PCA_result$rotation, 2)
head(PCA_result$rotation,3)

#Extracting the frist two principal components through a projection
points = t(t(PCA_result$rotation[, 1:2]) %*% golub) # Projection
dim(points)
head(points, 2)

par(pty="s") # This is used to make the plot square
plot(points, col="red")
points(points[28:38,], col="green")
