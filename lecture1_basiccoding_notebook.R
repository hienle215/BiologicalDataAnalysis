###Week1_lecture 1

### NATIONAL CONVENTIONS
print(1:10)

###?Reading R help
?seq #sequence generation

### Vector
x <- c(0.4, -0.2, 1)
x

y <- c(1, -1,2)
y
x*y

y[2]
y[1:2]
y > 0.5
y[y> 0.5]

### List
b <- list(x, 2000)
b

### MATRIX
A <- matrix(c(x,y), nrow=2, ncol=3, byrow=TRUE)
A

#Index matrix
A[2,]
A[1, 2:3]
A[, c(1,3)]

#OPERATION on ROWS and COLUMNS
sum <- c()
for(i in 1: nrow(A))
sum[i] <- mean (A[i,])  
print(sum)

### apply() function
apply(A, 1, mean)


ORDER rows or columns
colmeans <- apply(A,2, mean)
colmeans

ordered_index <- order(colmeans, decreasing = FALSE)
ordered_index

A[, ordered_index]


### DATA FAME
dd = data.frame(xvar=A[1,], yvar=A[2,])
dd[2,]

###WRITING FUNCTIONS
my_function = function(parameter1=0, parameter2=0) {
  return(parameter1+parameter2^2)
}
my_function(1,2)
my_function(parameter2=2, parameter1=1)
