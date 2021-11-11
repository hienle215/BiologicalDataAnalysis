### t-test in R
library(multtest)
data(golub)
dim(golub)

#let's perform t-test
t.test(golub[808,1:27], golub[808, 28:38])

### QQ-plot
#let's first read about qqnorm() and qqline
?qqnorm # produces a QQ plot of two datasets
?qqline
?qqplot

#QQ-plot for normal data
#Let's create some random normal/Gaussian data with different means (μ) and standard deviations (σ) and see how the QQ-plot looks like
#Different mean (μ) values, same standard deviation (σ)
par(mfrow=c(1,3), pty="s")
for(mu in c(-3,0,3)){
  set.seed(111)
  sd=1
  random_normal = rnorm(40, mean =mu, sd=sd)
  title = sprintf("mu=%s, sd=%s", mu, sd)
  qqnorm(random_normal, xlim=c(-5,5), ylim=c(-5,5), main=title)
  qqline(random_normal)
  abline(h=0, col="red")
}
#Different standard deviation (σ) values , same mean (μ)
par(mfrow=c(2,3), pty="s")
for (sd in c(0, 0.1, 0.7, 1, 2, 7)) {
  set.seed(111)
  mu = 0
  random_normal = rnorm(40, mean = mu, sd = sd)
  title = sprintf("mu=%s, sd=%s", mu, sd)
  qqnorm(random_normal, xlim=c(-5, 5), ylim=c(-5, 5), main=title) 
  qqline(random_normal)
  abline(v=0, col="red")
}

#QQ-plot for uniform data
set.seed(111)
random_uniform = runif(n=40, min=-3, max=3)
par(pty="s")
qqnorm(random_uniform, xlim=c(-5,5), ylim=c(-5,5))
qqline(random_uniform)
abline(v=0, h=0, col="red")

#A group of data that is from two different uniform distributions
set.seed(111)
random_uniform_1 <- runif(n=20, min=-3, max=-2)
random_uniform_2 <- runif(n= 20, min= 2, max= 3)
random_uniform <- c(random_uniform_1, random_uniform_2)
par(pty="s")
qqnorm(random_uniform, xlim=c(-5,5), ylim=c(-5,5))
qqline(random_uniform)
abline(v=0, h=0, col="red")

#QQ-plot for poisson data()
#Note that Possion data is descrete
sed.seed(111)
random_pois <- rpois(n=40, lambda=1)
par(pty="s")
qqnorm(random_pois, xlim=c(-5,5), ylim=c(-5,5))
qqline(random_pois)
abline(v=0, h=0, col="red")

#QQ-plot for CD33
#for both AML and ALL case
cd33data = golub[808,]
par(pty="s")
qqnorm(cd33data, xlim=c(-3,3), ylim=c(-3,3), main="Both ALL and AML cases")
qqline(cd33data)
abline(v=0, h=0, col="red")

#for each AML and ALL separately cases
par(mfrow = c(1,2))
qqnorm(cd33data[1:27], main ="ALL", xlim=c(-3,3), ylim=c(-3,3))
qqline(cd33data[1:27])
qqnorm(cd33data[28:38], main ="AML", xlim=c(-3,3), ylim=c(-3,3))
qqline(cd33data[28:38])

###PERMUTATION testing in R
set.seed(547) # If this is commented, every run gives you different p-value. why?
# the set.seed() functions sets the starting number used to generated a sequence of random numbers-it ensures that we get the same results if we start wit the same seed each time we run the same process. 
# seed is a base function that it is able to generate together other functions (rnorm, runif, sample) the same random value. We need to set the seed everytime we do some random stuff if we want consistency.
# it is a random number (seed) and 547 has no special meaning there. Using the same seed helps in reproducing the results of a code. 
gene <- 800 # Notice that we pick another gene

NUMPERMS <- 1e4

permstats <- c() # Need to initialize

for (i in 1: NUMPERMS){
  persample <- sample(1:38, size = 27, replace = FALSE)
  data1 <- golub[gene, persample]
  data2 <- golub[gene, -persample]
  permstats[i] <- t.test(data1, data2)$statistic
}

## Actual t-statistic for gene at position 800
ts = t.test(golub[800, 1:27], golub[800, 28:38])$statistic
ts

hist(permstats)
abline(v=ts, col='red')

## Calculating the p-value
p_value = mean(ts < abs(permstats))
p_value

### GENOME-WIDE testing in R
tstats <- mt.teststat(golub, golub.cl, test="t")

par(mfrow=c(1,2))
hist(tstats)
qqnorm(tstats, main="t-statistics")
qqline(tstats)
#Note that, for example, the 1st items in the tstats variable in the above code chunk is the t-statistic results from performing a t-test over the 1st gene in Golub data set. We can see this in the following code chunk.

#Note that mt.testat expects matrix or data frame, thus we use as.matrix() below
t_stat_gene_1 = mt.teststat(t(as.matrix(golub[1,])), golub.cl, test="t")
t_stat_gene_1
t_stat_gene_1 == tstats[1]


### GENOME-WIDE p-values
# approximation using the standard Gaussian distribution
p0 = 2 * (1 - pnorm(abs(tstats)))
p0 <- 2 * (1 - pnorm(abs(tstats)))
hist(p0)

#NOTE that, for example, the p-value for the 1st item in the tstats can be calculated as follows
p_value_gene_1 <- 2 * (1 - pnorm(abs(tstats[1])))
p_value_gene_1
p_value_gene_1 == p0[1]

#using the t-distribution
p1 <- 2 * (1 - pt(abs(tstats), df=38-2))
hist(p1)

### P-values usign rowttest()
#rowttests() function comes from the genefilter package available in Bioconductor
library(genefilter)
gol.fac = factor(golub.cl, levels=0:1, labels=c("ALL","AML"))
tstats = rowttests(golub, gol.fac)
tstats[1,]

#inspecting the results
idx <- order(tstats[,3]) # Smallest first
res <- data.frame(Name=golub.gnames[idx,2], pval=tstats[idx,3])
res[,1] <- substr(res[,1],1,15) #Shorter names
head(res,4) #notice how we get top 4 genes

#creating a gene list
res2 <- res[res[,2] < 1e-4,]
dim(res2)
res2 <- res[res[,2] < 1e-3,]
dim(res2)

### VOLCANO plot
#let's use the effect size for each gene and p-values to plot the volcano plot
par(pty="s")
plot(tstats[,2], -log10(tstats[,3])) #notice the logarithm

#let's set some criteria and select genes such that their p-value is less than 0.001 (above 4 on the Y-axis) and the effect size more than 1 (beyon[-1,1] on X-axis)
statcrit <- which(tstats[,3] < 1e-4)
effcrit <- which(abs(tstats[,2]) > 1)
totcrit <- which(tstats[,3] < 1e-4 & (abs(tstats[,2]) > 1))
length(statcrit)
length(effcrit)
length(totcrit)

#let's do a new volcano plot, coloring the points under different criteria with different colors
par(pty="s")
plot(tstats[,2], -log10(tstats[,3]))
points(tstats[statcrit, 2], -log10(tstats[statcrit, 3]), col="green")
points(tstats[effcrit, 2], -log10(tstats[effcrit, 3]), col="blue")
points(tstats[totcrit, 2], -log10(tstats[totcrit, 3]), col="red")


#### MULTIPLE testing hypothesis
### MULTIPLE testing correction
pcor <- p.adjust(res[,2], method="BH")
res$pcor <- pcor
head(res,4)

#visualizing the correction
par(pty="s")
plot(res$pval, res$pcor, type="l")

# counting the number of significant genes affter multiple testing correction
res2 <- res[res$pcor < 1e-4,]
dim(res2)
res2 <- res[res$pcor < 1e-3,]
dim(res2)

#p-value distribution when there is no true differences
set.seed(152)
library(genefilter) #for rowttests
m = 2e4
p = 30
data = matrix(rnorm(m*p, 0, 1), nrow=m)
labels = factor(c(rep(0,15), rep(1,15)))
res = rowttests(data, labels)
hist(res[,3])


res_ordered = sort(res[,3])
k_per_m = 1:m/m
par(pty="s")
plot(k_per_m, res_ordered)


#p-value distribution when there are true differences
set.seed(152)
m_1 = 1e3
nonnull_ind = sample(m,m_1)
data[nonnull_ind, 1:15] = matrix(rnorm(m_1*15, 0.9, 1), nrow=m_1)
res = rowttests(data, labels)
res_ordered = sort(res[,3])
hist(res[, 3], 40)

#h histogram of generated non-null genes
hist(res[nonnull_ind,3])

k_per_m = 1:m/m
par(pty="s")
plot(res_ordered / k_per_m)


##FDR estimation
#p-value threshold of α = 0.005
# these values may differ from the different time of running as there are some randomnes involved
alpha = 0.005
positives = which(res[,3] < alpha)
lambda = 0.5
pi_0 = length(which(res[,3] > lambda)) / m / (1-lambda)
falsepos = pi_0 * m * alpha
fdr = falsepos / length(positives)
pi_0
falsepos
fdr

# How good is the estimate?
called = vector("logical", m)
null = called
called[positives] = TRUE
null[-nonnull_ind] = TRUE
res_table = table(called, null)
addmargins(res_table)
sum(res_table[2,]) # positives
fdr = res_table[2,2] / sum(res_table[2,])
fdr


# Comparing FDR estimation and B-H correction
length(which(res_ordered) / k_per_m < fdr)
