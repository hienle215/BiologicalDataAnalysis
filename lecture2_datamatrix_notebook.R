### Lecutre 2 Notebook
#This notebook contains the commands that I saw during the lecture for basic commands in R 

###LOADING GOLUB DATA MATRIX
data(golub)
head(golub, n=1)

### INSPECTING GOLUB DIMENSIONS ANS CLASS LABELS
#now, that we have loaded the GOlub data matrix, we can inspect its dimensions using nrow(), ncol(), and dim() functions
nrow(golub)
ncol(golub)
dim(golub)

#Golub's class labels are stored in golub.cl variable
golub.cl

# inspecting the gene information at index 123 of Golub data
#golub.ganmes is a variable storin a matrix that contains the gene names and microarry probeset identifiers
golub.gnames[123,]


### PLOTTING data for CD33 gene
### fist, we need to find the gene's index using grep() and then use plot() function to plot the extracted data
?grep

gene_index <- grep("CD33", golub.gnames[,2])
plot(golub[gene_index,], xlab="Sample number", ylab="Expression level", main="Golub data for CD33")

#Plotting more than one variable in a plot using matplot() function
ind <- c(182,194,195) #selected genes
datatoplot = t(golub[ind,]) # Notice the transpose
matplot(datatoplot, type="l") # Type specifies a line plot
matplot(datatoplot, type="1") # Type specifies a line plot

#Plotting expression level histograms for CD33 gene using hist() fucntion
par(mfrow=c(1,3)) # 3 histograms side by side
hist(golub[808,], breaks=seq(-1.5,1.5,length=10), ylim=c(0,10)) # All
hist(golub[808,1:27], breaks=seq(-1.5,1.5,length=10), ylim=c(0,10)) # ALL
hist(golub[808,28:38], breaks=seq(-1.5,1.5,length=10), ylim=c(0,10)) # AML


###CONVERTING Golub class lables into factor using factor() function
gol.fac <- factor(x=golub.cl, levels = 0:1, labels=c("ALL", "AML"))
print(gol.fac)

### CALCULATING the mean of the CD33 gene expression values for the AML group using means() function
meanAML <- mean(golub[808, gol.fac=="AML"])
meanAML

meanALL <- mean(golub[808, gol.fac=="ALL"])
meanALL

### VISUALIZATON OF THE DISTRIBUTION of the CD33 gene accross two different groups using boxplot() function

boxplot(golub[808,] ~ gol.fac, boxwex=0.2) ## note the use of tilde ~ operator

###Plotting scatter plots of a set of selected genes using plot() function
# note that we use as.data.frame() function to transform the matrix to data frame
ind <- c(182, 194, 195) #selected genes
df <- as.data.frame(t(golub[ind,])) #select the data to plot
names(df) <- golub.gnames[ind,2] #attach the gene names
plot(df) #the scatter plot itself

### plotting a density function for a GAUssian distribution with ??=1 and ??2=1.5
#before we proceed, let's skim through the dnorm() function's manual page to learn wht it dose
?dnorm

mu = 1
sigma2 = 1.5
x = seq(mu-5, mu+5, len=1000)
plot(x, dnorm(x, mean=mu, sd=sqrt(sigma2)), type="l")

###plotting the quantile-quantile plot for two Gaussians with different population parameters
p <- seq(0,1,len=1000)
p=seq(0, 1, len=1000)
plot(qnorm(p,
           mean=0,
           sd=1),
     qnorm(p,
           mean=2,
           sd=2))

### GENERATING random data from Poisson model ??=5
#note that we are using some new functions var() and table(). If they sound unfimilar to you, skim throught their corresponding namual page
set.seed(111) # setting the seed so that every run results in the same set of random values
data=rpois(1000, lambda=5) # vector of 1000 values
mean(data) # mean
var(data) # Variance
table(data) # what does this function do?

#comparing the random data with a probability model
x <- 0:13
model <- dpois(x, lambda=5)
plot(table(data)/length(data))
points(x, model, col="red")

model <- dpois(x, lambda=4)#when changing the lambda value, we will get the differnt plot. Explain why?

#Exmining the central limit theorem
x <- runif(1e7, -0.5, 0.5)
y <- matrix(x, nrow=1e4)

#for sumes
sums <- colSums(y)
hist(sums)

# for means
means <- colMeans(y)
hist(means)
