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
pca <- function(data)
{
  #step 1. calculationg the mean of each column
  vector_ones <- matrix(1.0, nrow = nrow(data), ncol = 1)
  col_means <- matrix(colMeans(data),nrow = 1)
  data.mean <- vector_ones %*% col_means
  #substracting the mean
  data.scaled = data - data.mean
  #computing the cov matrix
  Sigma = data.scaled%*%t(data.scaled)/(ncol(data.scaled)-1)
  Eigen = eigen(Sigma)
  D = Eigen$values
  P = Eigen$vectors
  pca_results <- list('data.mean' = data.mean, 'P' = P, 'D' = D)
  return(pca_results)
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
  Sw = matrix(c(0),nrow= (ncol(data)-1),ncol = (ncol(data)-1))
  index = 1
  for(i in 1:dim(a))
  {
    for(j in index:(index+as.vector(a)[i]-1))
    {
      Sw = Sw + (t(matrix(data = data[index,-dim(data)[2]], nrow = 1)) %*% matrix(data = data[index,-dim(data)[2]], nrow = 1))
    }
    index = index + as.vector(a)[i]
    
  }
  Sb = matrix(c(0),nrow=24,ncol = 24)
  for(i in 1:dim(a))
  {
    Sb = Sb + as.vector(a)[i] *( t((class_mean[i,] - global_mean)) %*% ((class_mean[i,] - global_mean)) )
  }
  eigen = eigen(solve(Sw)%*%Sb)
  P = eigen$vectors
  D = eigen$values

  fda_results <- list('data.mean' = global_mean, 'P' = P, 'D' = D)
  return(fda_results)
}
threshold <- function(data, metric)
{
  data = data[,-ncol(data)]
  distances <- matrix(0,nrow(data),nrow(data))
  if(metric == "euclidan")
  {
    for(i in 1:nrow(data))
    {
      for(j in 1:nrow(data))
      {
        distances[i,j] = euclidan(data[i,],data[j,])
      }
    }
  }
  else if(metric == "manhattan")
  {
    for(i in 1:nrow(data))
    {
      for(j in 1:nrow(data))
      {
        distances[i,j] = manhattan(data[i,],data[j,])
      }
    }
  }
  else if(metric == "cosine")
  {
    for(i in 1:nrow(data))
    {
      for(j in 1:nrow(data))
      {
        distances[i,j] = cosine(data[i,],data[j,])
      }
    }
  }
  else
  {
    print('Wrong metric argument! Choose between: "euclidan", "manhattan", "cosine"')
  }
  return(distances)
}
euclidan <- function(x1,x2)
{
  return(sqrt(sum((x1 - x2) ^ 2)))
}
manhattan <- function(x1,x2)
{
  return(sum(abs(x1-x2)))
}
cosine <- function(x1,x2)
{
  return(sum(x1*x2)/(sqrt(sum(x1^2))*sqrt(sum(x2^2))))
}
my_knn <- function(img, training_data, k, threshold)
{
  training_data = matrix(unlist(training_data), nrow = nrow(training_data))
  #step1. transforming an img to a 1 row matrix
  vector_ones <- matrix(c(1), nrow = nrow(training_data), ncol = 1)
  img_matrix <- vector_ones%*%matrix(unlist(img), nrow = 1)
  result1 <- img_matrix-training_data[,1:ncol(training_data)-1]
  result2 <- result1 %*% t(result1)
  distances <- matrix(nrow = nrow(training_data), ncol = 1)
  for(i in 1:nrow(training_data))
  {
    distances[i,] = sqrt(result2[i,i])
  }
  distances = cbind(distances, training_data[,ncol(training_data)])
  distances.sorted = distances[order(distances[,1],decreasing = F),]
  k_nn = distances.sorted[1:k, 2]
  img.label = names(which.max(table(k_nn)))
  if(distances[1] > threshold) img.label = 0#threshold
  results <- list('img.label' = img.label)
  return(results)
  
}
KNN <- function(img, training_data, k = 4, threshold = 1000 , metric = "euclidan")
{
  
  neighbours <- matrix(Inf, nrow = nrow(training_data), ncol = 1)
  if(metric == "euclidan")
  {
    for(i in 1:nrow(training_data))
    {
      neighbours[i,] = euclidan(img, training_data[i,-ncol(training_data)])
    }
  }
  else if(metric == "manhattan")
  {
    for(i in 1:nrow(data))
    {
      neighbours[i,] = manhattan(img, training_data[i,-ncol(training_data)])
    }
  }
  else if(metric == "cosine")
  {
    for(i in 1:nrow(data))
    {
      neighbours[i,] = cosine(img, training_data[i,-ncol(training_data)])
    }
  }
  else
  {
    print('Wrong metric argument! Choose between: "euclidan", "manhattan", "cosine"')
  }
  neighbours = cbind(neighbours, training_data[,ncol(training_data)])
  neighbours_sorted = neighbours[order(neighbours[,1],decreasing = F),]
  k_nn = neighbours_sorted[1:k,2]
  img_label = names(which.max(table(k_nn)))
  print(neighbours_sorted[1,1])
  if(neighbours_sorted[1,1] > threshold) img_label = 0#threshold
  
  return(img_label)
}

