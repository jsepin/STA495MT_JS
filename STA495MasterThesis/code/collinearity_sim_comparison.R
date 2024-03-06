 # Small simulation to compare the scaling factor method with the rmvnorm method
library(tidyverse);library(mvtnorm);library(Collinearity)

# Load parameters defined earlier
boston_parameters <- readRDS("../data/boston_parameters.rds")
range_x1 <- boston_parameters$range_x1
range_x2 <- boston_parameters$range_x2
n_explore <- boston_parameters$n_explore
no_runs <- boston_parameters$no_runs
no_coll_magnitude <- boston_parameters$no_coll_magnitude
gamma_true <- boston_parameters$gamma_true

# Beta all the first condition
beta_true <- c(boston_parameters$beta_0, boston_parameters$beta_1[1],
               boston_parameters$beta_2[1])
s_y <- boston_parameters$s_y[1]
mean_x1 <- boston_parameters$mean_x1
mean_x2 <- boston_parameters$mean_x2
sd_x1 <- boston_parameters$sd_x1
sd_x2 <- boston_parameters$sd_x2

##################################################################################
# function to equilibrate
##################################################################################
my_equilibrate <- function(object, dataframe = TRUE){
  # Input is a design matrix X
  object <- as.matrix(object)
  d_ele <- 1 / sqrt( diag(t(object) %*% object) )
  d_mat <- diag(nrow = length(d_ele), ncol = length(d_ele))
  diag(d_mat) <- d_ele
  eX <- object %*% d_mat
  if(dataframe){
    eX <- as.data.frame(eX)
    colnames(eX) <- colnames(object)
  }
  return(eX)
}


##################################################################################
# function to calculate the properties needed
##################################################################################
prop_X <- function(X){
  m <- lm(data = X, y ~ x1 + x2)
  coef <- coef(m);names(coef) <- c("const", "x1", "x2")
  se <- sqrt(diag(vcov(m)));names(se) <- c("const", "x1", "x2")
  X <- X[,-4]
  sd_X <- apply(X,2, sd)
  cond_nu <- max(Collinearity::Var_decom_mat.matrix(as.matrix(X) )[,"cond_ind"])
  E <- my_equilibrate(X, dataframe = FALSE)
  trouble <- diag(solve(t(E)%*%E))
  rho <- cor(X[,-1])[1,2]
  names(trouble) <- c("const", "x1", "x2")
  df <- data.frame(t(c("trouble"=trouble, "sd"=sd_X,"cond_nu"=cond_nu ,"rho"=rho,
                       "beta"= coef, "se" = se)))
  #df$x <- deparse(substitute(X))
  return(df)
}

##################################################################################
# Generating data sets over rmvnorm
##################################################################################
set.seed(1234324)
correlation <- seq(from=-0.99,to = 0, length.out = no_coll_magnitude)# if over rmvnorm
df_rmv <- list()
for(run in 1:no_runs){
  for(i in 1:length(correlation)){
    rho <- correlation[i]
    eps_y <- rnorm(n = n_explore)#maybe has to be added in the inner loop!!
    # over rmvnorm
    X <- mvtnorm::rmvnorm(n =  n_explore, mean = c(0,0),
                sigma = matrix(c(1,rho,rho,1),ncol = 2),method="eigen" )
    X[,1] <- diff(range_x1) *pnorm(X[,1]) + range_x1[1]
    X[,2] <- diff(range_x2) *pnorm(X[,2]) + range_x2[1]
    X <- cbind("const"=1, "x1"=X[,1], "x2"=X[,2])
    X <- data.frame(X)
    X$y <- as.matrix(X)%*%beta_true + s_y*eps_y
    attr(X, "rho_fix") <- rho
    attr(X, "run") <- run
    df_rmv[[length(df_rmv)+1]] <- X
    }
}
# # save it
# saveRDS(df_rmv, file = "../data/df_rmv.rds")

# Generating properties
est_rmv <- cbind(do.call("rbind",lapply(X= df_rmv, FUN = function(x) {
  prop_X(x)
  } )
),
do.call("rbind",lapply(X= df_rmv, FUN = function(x) {
  c("fix"=attributes(x)$rho_fix, "run"=attributes(x)$run)
} )
)
)

# save it
saveRDS(est_rmv, file = "../data/est_rmv.rds")

##################################################################################
# Generating data sets over rmvnorm
# Stick with Normal distribution
##################################################################################
set.seed(1234324)
correlation <- seq(from=-0.99,to = 0, length.out = no_coll_magnitude)# if over rmvnorm
df_rmv_normal <- list()
for(run in 1:no_runs){
  for(i in 1:length(correlation)){
    rho <- correlation[i]
    eps_y <- rnorm(n = n_explore)#maybe has to be added in the inner loop!!
    # over rmvnorm
    X <- mvtnorm::rmvnorm(n =  n_explore, mean = c(0,0),
                          sigma = matrix(c(1,rho,rho,1),ncol = 2),method="eigen" )
    X[,1] <- X[,1]*sd_x1 + mean_x1
    X[,2] <- X[,2]*sd_x2 + mean_x2
    X <- cbind("const"=1, "x1"=X[,1], "x2"=X[,2])
    X <- data.frame(X)
    X$y <- as.matrix(X)%*%beta_true + s_y*eps_y
    attr(X, "rho_fix") <- rho
    attr(X, "run") <- run
    df_rmv_normal[[length(df_rmv_normal)+1]] <- X
  }
}
# save it
# saveRDS(df_rmv_normal, file = "../data/df_rmv_normal.rds")

# Generating properties
est_rmv_normal <- cbind(do.call("rbind",lapply(X= df_rmv_normal, FUN = function(x) {
  prop_X(x)
} )
),
do.call("rbind",lapply(X= df_rmv_normal, FUN = function(x) {
  c("fix"=attributes(x)$rho_fix, "run"=attributes(x)$run)
} )
)
)

# save it
saveRDS(est_rmv_normal, file = "../data/est_rmv_normal.rds")


##################################################################################
#scale finder
##################################################################################
scale_finder <- function(cond_numb_fix, x1, eps_x, par){
  f <- function(scale_f, cond_numb_fix, x1, eps_x, par){
    x2 <- par[1] + par[2]*x1 + eps_x*scale_f
    df <- cbind("const"=1, "x1"=x1, "x2"=x2)
    cond_number_current <- max(Collinearity::Var_decom_mat(df, equilibration = T)[,"cond_ind"])
    return(cond_numb_fix - cond_number_current)
  }
  scale_f <- uniroot(f, interval = c(0,1000), cond_numb_fix = cond_numb_fix, x1 = x1, eps_x=eps_x, par = par)$root
  return(scale_f)
}



##################################################################################
# Scale factor method
##################################################################################

set.seed(1234324)
df_scalefac <- list()
for(run in 1:no_runs){
  x1 <- runif(n = n_explore, min = range_x1[1], max = range_x1[2])
  eps_x <- rnorm(n = length(x1),mean = 0, sd = 1)
  eps_y <- rnorm(n = length(x1),mean = 0, sd = 1)
  
  s_x_grid <- seq(from=scale_finder(cond_numb_fix =max(est_rmv$cond_nu),
                                      x1=x1, eps_x= eps_x, par=gamma_true),
                    to = scale_finder(cond_numb_fix =min(est_rmv$cond_nu),
                                      x1=x1, eps_x= eps_x, par=gamma_true),
                    length.out = no_coll_magnitude)
  
  
  for(i in 1:length(s_x_grid)){
    s_x <- s_x_grid[i]
    x2 <- gamma_true[1] + x1*gamma_true[2] +  eps_x*s_x
    X <- cbind("const"=1, "x1"=x1, "x2"=x2)
    X <- data.frame(X)
    X$y <- as.matrix(X)%*%beta_true + s_y*eps_y# all have the same error!
    attr(X, "s_x") <- s_x
    attr(X, "run") <- run
    df_scalefac[[length(df_scalefac)+1]] <- X
  }
}

# save it
# saveRDS(df_scalefac, file = "../data/df_scalefac.rds")

# Generating properties
est_scalefac<- cbind(do.call("rbind",lapply(X= df_scalefac, FUN = function(x) {
  prop_X(x)
} )
),
do.call("rbind",lapply(X= df_scalefac, FUN = function(x) {
  c("fix"=attributes(x)$s_x, "run"=attributes(x)$run)
} )
)
)

# save it
saveRDS(est_scalefac, file = "../data/est_scalefac.rds")



################################################################################
# Plotting demonstration data frames ----
################################################################################

wanted_cond <- c(20,100)
df_show <- list()

#===============================================================================
# Scale factor method:
#===============================================================================
set.seed(12343249)
x1 <- runif(n = n_explore, min = range_x1[1], max = range_x1[2])
eps_x <- rnorm(n = length(x1),mean = 0, sd = 1)
for(cond in wanted_cond ){
  s_x <- scale_finder(cond_numb_fix =cond,x1=x1, eps_x= eps_x, par=gamma_true)
  x2 <- gamma_true[1] + x1*gamma_true[2] +  eps_x*s_x
  X <- cbind("const"=1, "x1"=x1, "x2"=x2)
  X <- data.frame(X)
  cond_nu <- max(Collinearity::Var_decom_mat(as.matrix(X))[,"cond_ind"])
  rho <- cor(X[,-1])[1,2]
  attr(X, "cond_nu") <- cond_nu
  attr(X, "rho") <- rho
  X$method <- paste0("scale_",cond)
  df_show[[length(df_show)+1]] <- X
}

#===============================================================================
# rmvnorm Normal distribution:
#===============================================================================

set.seed(12343249)
scale_finder_norm <- function(cond_numb_fix, sd_x1,sd_x2, mean_x1, mean_x2){
  f <- function(rho, cond_numb_fix, sd_x1,sd_x2, mean_x1, mean_x2){
    X <- mvtnorm::rmvnorm(n =  n_explore, mean = c(0,0),
                          sigma = matrix(c(1,-abs(rho),-abs(rho),1),ncol = 2),
                          method="eigen" )
    X[,1] <- X[,1]*sd_x1 + mean_x1
    X[,2] <- X[,2]*sd_x2 + mean_x2
    X <- cbind("const"=1, "x1"=X[,1], "x2"=X[,2])
    cond_number_current <- max(Collinearity::Var_decom_mat(X,
                                                equilibration = T)[,"cond_ind"])
    return(cond_numb_fix - cond_number_current)
  }
  scale_f <- uniroot(f, interval = c(0,1), cond_numb_fix = cond_numb_fix,
              sd_x1 = sd_x1, sd_x2=sd_x2,mean_x1=mean_x1 ,mean_x2=mean_x2)$root
  return(scale_f)
}

for(cond in wanted_cond ){
  rho <- scale_finder_norm(cond_numb_fix=cond,  sd_x1,sd_x2, mean_x1, mean_x2)
  X <- mvtnorm::rmvnorm(n =  n_explore, mean = c(0,0),
                        sigma = matrix(c(1,-abs(rho),-abs(rho),1),ncol = 2),method="eigen" )
  X[,1] <- X[,1]*sd_x1 + mean_x1
  X[,2] <- X[,2]*sd_x2 + mean_x2
  X <- cbind("const"=1, "x1"=X[,1], "x2"=X[,2])
  X <- data.frame(X)
  cond_nu <- max(Collinearity::Var_decom_mat(as.matrix(X))[,"cond_ind"])
  rho <- cor(X[,-1])[1,2]
  attr(X, "cond_nu") <- cond_nu
  attr(X, "rho") <- rho
  X$method <- paste0("rmvnormal_",cond)
  df_show[[length(df_show)+1]] <- X
}

#===============================================================================
# rmvnorm Uniform distribution:
#===============================================================================

set.seed(12343249)
scale_finder_norm <- function(cond_numb_fix, sd_x1,sd_x2, mean_x1, mean_x2){
  f <- function(rho, cond_numb_fix, sd_x1,sd_x2, mean_x1, mean_x2){
    X <- mvtnorm::rmvnorm(n =  n_explore, mean = c(0,0),
                          sigma = matrix(c(1,-abs(rho),
                                           -abs(rho),1),ncol = 2),method="eigen" )
    X[,1] <- diff(range_x1) *pnorm(X[,1]) + range_x1[1]
    X[,2] <- diff(range_x2) *pnorm(X[,2]) + range_x2[1]
    X <- cbind("const"=1, "x1"=X[,1], "x2"=X[,2])
    cond_number_current <- max(Collinearity::Var_decom_mat(X, 
                                                equilibration = T)[,"cond_ind"])
    return(cond_numb_fix - cond_number_current)
  }
  scale_f <- uniroot(f, interval = c(0,1), cond_numb_fix = cond_numb_fix,
                sd_x1 = sd_x1, sd_x2=sd_x2,mean_x1=mean_x1 ,mean_x2=mean_x2)$root
  return(scale_f)
}

for(cond in wanted_cond ){
  rho <- scale_finder_norm(cond_numb_fix=cond,  sd_x1,sd_x2, mean_x1, mean_x2)
  X <- mvtnorm::rmvnorm(n =  n_explore, mean = c(0,0),
                        sigma = matrix(c(1,-abs(rho),-abs(rho),1),ncol = 2),
                        method="eigen" )
  X[,1] <- diff(range_x1) *pnorm(X[,1]) + range_x1[1]
  X[,2] <- diff(range_x2) *pnorm(X[,2]) + range_x2[1]
  X <- cbind("const"=1, "x1"=X[,1], "x2"=X[,2])
  X <- data.frame(X)
  cond_nu <- max(Collinearity::Var_decom_mat(as.matrix(X))[,"cond_ind"])
  rho <- cor(X[,-1])[1,2]
  attr(X, "cond_nu") <- cond_nu
  attr(X, "rho") <- rho
  X$method <- paste0("rmv_",cond)
  df_show[[length(df_show)+1]] <- X
}

################################################################################
# Combining and Saving
df_to_plot <- do.call("rbind", df_show)
df_to_plot_prop <- do.call("rbind",lapply(df_show, function(x){
c("cond_nu"=attributes(x)$cond_nu,
"rho"=attributes(x)$rho,
"method"=unique(x$method)
)
}
))

# save it
saveRDS(list("df" = df_to_plot, "prop"= df_to_plot_prop ), 
        file = "../data/cross_section_df.rds")



################################################################################
# What is the lowest trouble value? ----
################################################################################

library(stringr); library(Collinearity)
set.seed(324324)
df_min_trouble <- matrix(NA,ncol = 5, nrow = 200)
colnames(df_min_trouble) <- c("cond_nu", "trouble_const", "trouble_x1",
                              "trouble_x2", "rho")
for(i in 1:nrow(df_min_trouble)){
  X <- mvtnorm::rmvnorm(n =  n_explore, mean = c(mean_x1,mean_x2),
                        sigma = matrix(c(sd_x1^2,0,0,sd_x2^2),ncol = 2)
  )
  X <- cbind("const"=1, "x1"=X[,1], "x2"=X[,2])
  X <- data.frame(X)
  cond_nu <- max(Collinearity::Var_decom_mat(as.matrix(X))[,"cond_ind"])
  rho <- cor(X[,-1])[1,2]
  E <- Collinearity::equilibrate_matrix(as.matrix(X))
  trouble <- diag(solve(t(E)%*%E))
  df_min_trouble[i,] <- c("cond_nu"=cond_nu, "trouble"=trouble, "rho" = rho)
}
df_min_trouble <- data.frame(df_min_trouble)
df_min_trouble <- df_min_trouble %>% gather(-c(cond_nu,rho),key = "variable",
                                            value="trouble" )
df_min_trouble$variable <- str_replace(df_min_trouble$variable, "trouble_","")

# save it
saveRDS(df_min_trouble, file = "../data/df_min_trouble.rds")

################################################################################
################################################################################

# clearing stuff
rm(df_min_trouble, boston_parameters,df_rmv, df_rmv_normal, df_scalefac, df_show,
   df_to_plot, df_to_plot_prop, est_rmv, est_rmv_normal,est_scalefac, X)




