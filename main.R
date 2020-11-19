#loading packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(pacman,OpenImageR,reshape,tidyr,readr,splitstackshape,ggplot2)

#loading functions defined in another file
#setting working dir
setwd("~/Desktop/fisher-discriminant-analysis")
source('functions.r')

#loading data
data = loading_img(paste(getwd(),"Training", sep = "/"))
A = data$data
B = data$labels
data.labeled = cbind(A,B)
dim(data.labeled)

#within class covariance
data <- data.labeled
a <- table(data[,ncol(data)])
class_mean = matrix(data = c(0), nrow = dim(a), ncol = ncol(data) - 1)
global_mean = matrix(data = c(0), nrow = 1, ncol = ncol(data) - 1)
for(i in 1:nrow(data))
{
  class_mean[data[i,ncol(data)],] = class_mean[data[i,ncol(data)],] + data[i,1:ncol(data)-1]
}
for(i in 1:dim(a))
{
  global_mean = global_mean + class_mean[i,]
  class_mean[i,] = class_mean[i,]/as.vector(a)[i]
}
global_mean = global_mean / dim(data)[1]
#normalizin the data, by substracting the mean of each class
for(i in 1:nrow(data))
{
  data[i,1:ncol(data)-1] = data[i,1:ncol(data)-1] - class_mean[data[i,ncol(data)],]
}
#calculating within-class variance Sw

Sw = 0
index = 1
for(i in 1:dim(a))
{
  for(j in index:(index+as.vector(a)[i]-1))
  {
    print(paste("i =",i,", j=",j,"index=",index, sep=" "))
    Sw = Sw + matrix(data = data[index,], nrow = 1) %*% t(matrix(data = data[index,], nrow = 1))
  }
  index = index + as.vector(a)[i]
  
}





