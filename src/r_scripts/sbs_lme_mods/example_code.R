library(R.matlab)
library(R.matlab)

# Read the .mat file
data <- readMat("C:\\Users\\ark007\\Downloads\\ark_dat_flasso.mat")

# View the variable names
names(data)


mat <- abs(data$tmp)

mat <- mat ##Took a subset for fast exploration

y <- array(mat)

p <- prod(dim(mat)[1:2])
X <- rbind(diag(p),diag(p),diag(p),diag(p))

library(genlasso)

out <- fusedlasso2d(y, X, dim(mat)[1],dim(mat)[2])

beta_hat <- coef(out,lambda = out$lambda[10])$beta
img_est <- matrix(beta_hat,dim(mat)[1],dim(mat)[2])

image(img_est)
