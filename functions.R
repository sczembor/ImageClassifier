loading_img <- function(directory)
{
  original_directory <- getwd()
  setwd(directory)
  
  #loading file names
  filenames <- list.files()
  length(filenames)
  #loading each image to one row in matrix A
  img <- readImage(filenames[1])
  dim(img)
  red <- as.vector(t(img[,,1]))
  green <- as.vector(t(img[,,2]))
  blue <- as.vector(t(img[,,3]))
  A <- matrix(c(red,green,blue),nrow = 1)
  B <- matrix(c(parse_number(filenames[1])),nrow = 1)
  for(i in 2:length(filenames))
  {
    img <- readImage(filenames[i])
    red <- as.vector(t(img[,,1]))
    green <- as.vector(t(img[,,2]))
    blue <- as.vector(t(img[,,3]))
    A <- rbind(A, c(red,green,blue))
    B <- rbind(B, c(parse_number(filenames[i])))
  }
  results <- list("data" = A, "labels" = B)
  setwd(original_directory)
  return(results)
}

FDA <- function(data)
{
  #within class covariance
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
  #TODO calculatin between-class variance Sb
  # Sw = 0
  # index = 1
  # for(i in 1:dim(a))
  # {
  #   Sw = Sw + data[index:index+as.vector(a)[i],] %*% t(data[index:index+as.vector(a)[i],])
  #   index = index + as.vector(a)[i]
  # }
  
  
  results <- list('mean' = mean)
  
}