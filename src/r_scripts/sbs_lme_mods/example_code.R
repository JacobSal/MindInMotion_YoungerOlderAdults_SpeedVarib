library(R.matlab)
library(R.matlab)

# Read the .mat file
data <- readMat("ark_dat_flasso.mat")

# View the variable names
names(data)


mat <- abs(data$tmp)

mat <- mat ##Took a subset for fast exploration

y <- array(mat)

p <- prod(dim(mat)[1:2])
X <- rbind(diag(p),diag(p),diag(p),diag(p)) 

library(genlasso)
#-- images are stacked and the X determines the uniqueness of these images from each other
# in the example he simply uses diagnols of 1's, but these values could theoretically be 
# replaced wiht others to perform weighting of the images to each other.
out <- fusedlasso2d(y, X, dim(mat)[1],dim(mat)[2])
#-- 
  
beta_hat <- coef(out,lambda = out$lambda[10])$beta
img_est <- matrix(beta_hat,dim(mat)[1],dim(mat)[2])

image(img_est)
