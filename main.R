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
data_labeled = data.frame(cbind(A,B))
acc = matrix(0, 2, 10)
#obtaining the optimal K and similarity metric
for(i in 1:4)
{
  #sampling the data
  #it is very important to run the FDA and PCA only on training data so we
  #do not affect the results of the test
  str = stratified(data_labeled, paste('X',108001, sep=''),4,bothSets = T)
  train_data = str$SAMP1[,-108001]
  train_labels = str$SAMP1[,108001]
  test_data = str$SAMP2[,-108001]
  test_labels = str$SAMP2[,108001]
  
  #casting to matrix data type
  train_data = matrix(unlist(train_data), nrow = nrow(train_data))
  train_labels = matrix(unlist(train_labels), nrow = nrow(train_labels))
  test_data = matrix(unlist(test_data), nrow = nrow(test_data))
  test_labels = matrix(unlist(test_labels), nrow = nrow(test_labels))
  
  #running pca on only the training set
  pca_results = pca(train_data)
  pca_eigenvalues = pca_results$D
  pca_eigenvectors = pca_results$P
  pca_data_mean = pca_results$data.mean
  
  ######centering and projecting train data
  train_data_scaled = train_data - pca_data_mean
  pca_cum_var = cumsum(pca_eigenvalues)/sum(pca_eigenvalues)
  pca_cum_var = data.frame(
    'pca' = 1:length(pca_cum_var),
    'cum.variance' = pca_cum_var
  )
  pca_cum_var
  pca_eigenvectors = t(train_data_scaled)%*%pca_eigenvectors
  train_data_projected = train_data_scaled %*% pca_eigenvectors[,1:24]#nubmer of calss - 1
  
  #plotting the variance
  ggplot(data=pca_cum_var, aes(x=pca,y=cum.variance))+geom_line()+geom_point()+scale_x_continuous(breaks = seq(0, 150, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  
  #fda
  data <- cbind(train_data_projected,train_labels)
  fda_results <- FDA(data)
  fda_eigenvalues = fda_results$D
  fda_eigenvectors = fda_results$P
  fda_data_mean = fda_results$data.mean
  
  ######centering and projecting the data
  #fda_data_mean_broadcasted = matrix(1, nrow = dim(data)[1], ncol = 1) %*% fda_data_mean
  #data_scaled = data[,-ncol(data)] - fda_data_mean_broadcasted
  fda_cum_var  = cumsum(fda_eigenvalues)/sum(fda_eigenvalues)
  train_data_projected = data[,-ncol(data)] %*% fda_eigenvectors[,1:8]
  train_data_projected = cbind(train_data_projected, data[,ncol(data)])
  
  fda_cum_var  = cumsum(fda_eigenvalues)/sum(fda_eigenvalues)
  fda_cum_var = data.frame(
    'fda' = 1:length(fda_cum_var),
    'cum.variance' = fda_cum_var
  )
  fda_cum_var
  #plotting the variance
  ggplot(data=fda_cum_var, aes(x=fda,y=cum.variance))+geom_line()+geom_point()+scale_x_continuous(breaks = seq(0, 150, by = 4)) +
    scale_y_continuous()
  plot(fda_cum_var)
  
  #transforming the test data with the same parameters as train data
  test_data_scaled = test_data - pca_data_mean[1:nrow(test_data),]
  test_data_projected = test_data_scaled %*% pca_eigenvectors[,1:24]
  test_data_projected = test_data_projected %*% fda_eigenvectors[,1:8]
  
  #threshold
  metrics <- c('euclidan', 'manhattan', 'cosine')
  for(p in 1:2)
  {
    distance_matrix = threshold(train_data_projected,metrics[p])
    threshold1 = mean(distance_matrix)/3
    print(threshold1)
    
    hits = matrix(0, nrow = 1, ncol = 10)
    for(j in 1:nrow(test_data))
    {
      for(k in 1:10)
      {
        result = KNN(test_data_projected[j,], train_data_projected, k, threshold1, metrics[p])
        if(strtoi(result) == test_labels[j,])
        {
          hits[,k] = hits[,k] + 1
        }
      }
    }
    acc[p,] = acc[p,] + hits
  }
}
acc
acc = acc / 200

#plotting with ggplot
eucildan_acc = data.frame(
  'K' = 1:dim(acc)[2],
  'acc' = acc[1,]
)

p <- ggplot(data=eucildan_acc , aes(x=K,y=acc))+geom_line()+geom_point()+scale_x_continuous(breaks = seq(0, 10, by = 1))
p <- p + labs(title = "Accuracy for given K nearest neighbours",
              subtitle = "Similarity metric: Euclidan Distance")
p
manhattan_acc = data.frame(
  'K' = 1:dim(acc)[2],
  'acc' = acc[2,]
)

b <- ggplot(data=manhattan_acc , aes(x=K,y=acc))+geom_line()+geom_point()+scale_x_continuous(breaks = seq(0, 10, by = 1))
b <- b + labs(title = "Accuracy for given K nearest neighbours",
              subtitle = "Similarity metric: Manhattan Distance")
b

#ploting with standard function
par(mfrow=c(2,1))
plot(
  pch = 19,         # PEŁNE KÓŁKO
  cex = 1.0,        # SKALOWANIE
  col = "#cc0000",  # CZEROWNY
  acc[1,], 
  type = 'b',
  xlim = c(0, 10), 
  ylim = c(0, 1),
  xlab = 'K',
  ylab = 'accuracy',
  main = "euclidan distance"
)
plot(
  pch = 19,         # PEŁNE KÓŁKO
  cex = 1.0,        # SKALOWANIE
  col = "#0000cc",  # NIEBIESKI
  acc[2,], 
  type = 'b',
  xlim = c(0, 10), 
  ylim = c(0, 1),
  xlab = 'K',
  ylab = 'accuracy',
  main = "manhattan distance"
)
par(mfrow=c(1,1))

#testing the threshold
  #sampling the data
  
  threshold.data = cbind(A,B)
  train_data = threshold.train.data = threshold.data[7:150,-108001]
  train_labels = threshold.test.label = threshold.data[1:6, 108001]
  test_data = threshold.data[1:6, -108001]
  test_labels = threshold.data[1:6, 108001]
  
  #running pca on only the training set
  pca_results = pca(train_data)
  pca_eigenvalues = pca_results$D
  pca_eigenvectors = pca_results$P
  pca_data_mean = pca_results$data.mean
  
  ######centering and projecting train data
  train_data_scaled = train_data - pca_data_mean
  pca_eigenvectors = t(train_data_scaled)%*%pca_eigenvectors
  train_data_projected = train_data_scaled %*% pca_eigenvectors[,1:24]#nubmer of calss - 1
  
  
  #fda
  data <- cbind(train_data_projected,train_labels)
  fda_results <- FDA(data)
  fda_eigenvalues = fda_results$D
  fda_eigenvectors = fda_results$P
  fda_data_mean = fda_results$data.mean
  
  ######centering and projecting the data
  #fda_data_mean_broadcasted = matrix(1, nrow = dim(data)[1], ncol = 1) %*% fda_data_mean
  #data_scaled = data[,-ncol(data)] - fda_data_mean_broadcasted
  fda_cum_var  = cumsum(fda_eigenvalues)/sum(fda_eigenvalues)
  train_data_projected = data[,-ncol(data)] %*% fda_eigenvectors[,1:8]
  train_data_projected = cbind(train_data_projected, data[,ncol(data)])
  
  #transforming the test data with the same parameters as train data
  test_data_scaled = test_data - pca_data_mean[1:nrow(test_data),]
  test_data_projected = test_data_scaled %*% pca_eigenvectors[,1:24]
  test_data_projected = test_data_projected %*% fda_eigenvectors[,1:8]
  
  #threshold
  denominator <- seq(2,5,0.5)
  acc = c(rep(0,length(denominator)))
  thresshold_value = c(rep(0,length(denominator)))
  for(d in 1:length(denominator))#different denominators
  {
    distance_matrix = threshold(train_data_projected,"euclidan")
    threshold1 = mean(distance_matrix)/denominator[d]
    thresshold_value[d] = threshold1
    
    hits = 0;
    for(j in 1:nrow(test_data))
    {
        result = KNN(test_data_projected[j,], train_data_projected, 4, threshold1, "euclidan")
        if(strtoi(result) == 0)
        {
          hits = hits + 1
        }
    }
    acc[d] = acc[d] + hits
  }
  acc = acc *6
acc = acc / 6
acc
thresshold_acc = data.frame(
  'thresshold' =  1:length(denominator),
  'acc' = acc
)

c <- ggplot(data=thresshold_acc , aes(x=thresshold,y=acc))+geom_line()+geom_point()+scale_x_continuous(breaks = seq(0, length(denominator), by = 1))
c <- c + labs(title = "Accuracy for different thresshold values",
              subtitle = "Similarity metric: Euclidan Distance; K = 4")
c




