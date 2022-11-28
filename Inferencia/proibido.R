#Código proibido, cuidado meu jovem.
seed = 123
set.seed(seed)
# Generate a random exponential variable, using the quantilic function
r_exp <- function(n, beta){
  vec<-numeric(n)
  for(i in 1:n){
    q <- runif(1, 0, 1)
    vec[i] <- -log(1-q)/beta
  }
  return(vec)
}

#Example 
#n = 10, n = 30 e n = 50,
x <- r_exp(n = 1000, beta = 0.5) #I remove the "seed" that had in the function
hist(x)

est_emv_exp<- 1/mean(x)
mean_emv_exp <- n*beta/(n-1)*beta
var_emv_exp <- ((n*beta)^2)/(((n-1)^2)*(n-2))
sd_em_exp <- sqrt(var_emv_exp)


# Making the first script about bootstrap.
# After that i will start to make a general function.
B = 50
n = 10
beta = 0.5
seed = 123
vec_est_B <- numeric(B)
x <- r_exp(n = n, beta = beta, seed)
for (i in 1:B){
  resample_x <- sample(x, replace = T)
  vec_est_B[i] <- 1/mean(resample_x)
}
est_B <- mean(vec_est_B)
sd_B <- sd(vec_est_B)


#General Function
bootstrap_exp <- function(nobs = c(10, 30, 50), B = c(50, 100, 500), beta = c(1, 0.5), seed = 123, 
                          type = 'non-parametric', view_results = T, export = T){
  
  if(!is.numeric(B) | !is.numeric(nobs) | !is.numeric(beta) |
     !is.numeric(seed) |!is.logical(view_results) |
     !is.logical(export)){
    stop('Some arguments entered in the function are wrong. Follow the instructions: \n
     * "nobs" -> These are the sample sizes that will be evaluated in the simulations, it must be a numerical vector.\n
     * "B" -> Represents the number of bootstrap sample size, must be a numeric argument.\n
     * "beta" -> These are the parameters of Exponential distribution that you want to estimate in you simulation. It must be a numerical vector.\n
     * "view_results" -> If you prefer to print the simulation results for each sample size and parameterization.\n
     * "export" -> Whether to save the results in .txt file.\n')
  }
  
  set.seed(seed)
  
  res_est_b_list <- list()
  res_est_b_list_full<-list()
  for(i in 1:length(beta)){
    for(j in 1:length(B)){
      res_est_b_list[[j]] <- `colnames<-`(matrix(c(0), nrow = B[j], ncol = length(nobs)), c(paste0('Nobs =  ', nobs)))
      
    }
    names(res_est_b_list) <- c(paste0('B = ', B))
    res_est_b_list_full[[i]] <- res_est_b_list
  }
  names(res_est_b_list_full) <- c(paste0('Beta = ', beta))
  
  # Setting the environment to receive the results
  res_est_array <- array(c(0), dim = c(length(beta), length(B)+1, length(nobs)), 
                         dimnames = list(paste0('Beta = ', beta),
                                         c('EMV', paste0('Est.B ', B)),
                                         paste0('nobs = ', nobs)))
  
  res_bias_array <- array(c(0), dim = c(length(beta), length(B)+1, length(nobs)), 
                          dimnames = list(paste0('Beta = ', beta),
                                          c('Bias EMV', paste0('Bias Est.B ', B)),
                                          paste0('nobs = ', nobs)))
  
  res_sd_array <- array(c(0), dim = c(length(beta), length(B)+1, length(nobs)), 
                        dimnames = list(paste0('Beta = ', beta),
                                        c('Sd EMV', paste0('Sd Est.B ', B)),
                                        paste0('nobs = ', nobs)))
  #Starts the simulation
  for(i in 1:length(nobs)){
    for(j in 1:length(beta)){
      x <- rexp(n = nobs[i], rate = beta[j])
      n <-length(x)
      est_emv_exp <- 1/mean(x)
      sd_emv_exp <- sqrt(((n*beta[j])^2)/(((n-1)^2)*(n-2)))
      res_est_array[j, 1, i] <- est_emv_exp
      res_bias_array[j, 1, i] <- beta[j] - est_emv_exp
      res_sd_array[j, 1, i] <- sd_emv_exp
      for(k in 1:length(B)){
        vec_est_B <- numeric(B[k])
        for(l in 1:B[k]){
          resample_x <- sample(x, replace = T)
          vec_est_B[l] <- 1/mean(resample_x)
        }
        res_est_b_list_full[[j]][[k]][,i] <- vec_est_B
        res_est_array[j, k+1, i] <- mean(vec_est_B)
        res_bias_array[j, k+1, i] <- beta[j] - mean(vec_est_B)
        res_sd_array[j, k+1, i] <-  sd(vec_est_B)
        
      }    
    }
  }
  
  # if you want to print all the results in the terminal
  if(view_results){
    cat('\n')
    cat("Exponential Distribution \n")
    cat('\n')
    for(i in 1:length(nobs)){
      cat("####################################################\n")
      cat("Sample size:", nobs[i], '\n')
      cat("####################################################\n")
      for(j in 1:length(beta)){
        cat('\n')
        cat("Results for parameter", rownames(res_est_array)[j],"\n")
        cat("Beta:", beta[j],'\n')
        cat("Estimate EMV:", res_est_array[j, 1, i],' \n')
        cat("BIAS EMV:", res_bias_array[j, 1, i],' \n')
        cat("Standard deviation EMV:", res_sd_array[j, 1, i],' \n')
        cat('\n')
        for(k in 1:length(B)){
          cat("Bootstrap sample size:", B[k], '\n')
          cat("Est.B", B[k], ':', res_est_array[j, k, i],' \n')
          cat("BIAS Est.B", B[k], ':', res_bias_array[j, k, i],' \n')
          cat("Sd Est.B", B[k],':', res_sd_array[j, k, i],' \n')
          cat('\n')
        }
        cat("####################################################\n")
      }
      cat('\n')
    }
  }
  
  # if you want to export the results saving as .txt file
  if(export){
    name_file_export<- paste0('Simulation results.txt')
    file_sink <- file(name_file_export)
    sink(name_file_export, append = T)
    cat('\n')
    cat("Exponential Distribution \n")
    cat('\n')
    for(i in 1:length(nobs)){
      cat("####################################################\n")
      cat("Sample size:", nobs[i], '\n')
      cat("####################################################\n")
      for(j in 1:length(beta)){
        cat('\n')
        cat("Results for parameter", rownames(res_est_array)[j],"\n")
        cat("Beta:", beta[j],'\n')
        cat("Estimate EMV:", res_est_array[j, 1, i],' \n')
        cat("BIAS EMV:", res_bias_array[j, 1, i],' \n')
        cat("Standard deviation EMV:", res_sd_array[j, 1, i],' \n')
        cat('\n')
        for(k in 1:length(B)){
          cat("Bootstrap sample size:", B[k], '\n')
          cat("Est.B", B[k], ':', res_est_array[j, k, i],' \n')
          cat("BIAS Est.B", B[k], ':', res_bias_array[j, k, i],' \n')
          cat("Sd Est.B", B[k],':', res_sd_array[j, k, i],' \n')
        }
        cat("####################################################\n")
      }
      cat('\n')
    }
    sink()
    close(file_sink)
  }
 
  #Storing within a list.
  list(res_est_b_list_full = res_est_b_list_full,
       res_est_array = res_est_array,
       res_bias_array = res_bias_array,
       res_sd_array = res_sd_array,
       nobs = nobs, B = B, beta = beta)
  
}

#Testa aqui, mizera.
res <- bootstrap_exp(nobs = c(10, 30, 50), B = c(50, 100, 500), beta = c(1, 0.5), seed = 1987, 
                 type = 'non-parametric', view_results = F, export = F)

res$res_est_array

##########################################
####### SECTION FOR THE GRAPHICS #########
##########################################

library(tidyverse)
beta<-res$beta
B<-res$B
nobs<-res$nobs

# Graphics of distributions of the estimatives for beta = 1
par(mfrow = c(3,3))
for(j in 1:length(B)){
  for(k in 1:length(nobs)){
    hist(res$res_est_b_list_full[[1]][[j]][,k], freq = F,
         main = paste0('Beta = ', beta[1],' | B = ', B[j], ' | N = ', nobs[k]), 
         xlab = NULL)
  }
}

# Graphics of distributions of the estimatives for beta = 0.5
par(mfrow = c(3,3))
for(j in 1:length(B)){
  for(k in 1:length(nobs)){
    hist(res$res_est_b_list_full[[2]][[j]][,k], freq = F,
         main = paste0('Beta = ', beta[2],' | B = ', B[j], ' | N = ', nobs[k]), 
         xlab = NULL)
  }
}

par(mfrow = c(2,3))
seed = 1987
set.seed(seed)
# Nobs = 10 and Beta = 1
x <- rexp(n = 10, rate = 1)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 10 | Beta = ', beta[1]), col = "gray", ylim = c(0, max(fun)), xlab = '')
lines(x2, dexp(x2, rate = res$res_est_array[1,1,1]), col = 1, lty = 1, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,2,1]), col = 2, lty = 2, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,3,1]), col = 3, lty = 3, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,4,1]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 0.6, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)


set.seed(seed)
# Nobs = 30 and Beta = 1
x <- rexp(n = 30, rate = 1)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 30 | Beta = ', beta[1]), col = "gray", ylim = c(0, max(fun)), xlab = '')
lines(x2, dexp(x2, rate = res$res_est_array[1,1,2]), col = 1, lty = 1, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,2,2]), col = 2, lty = 2, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,3,2]), col = 3, lty = 3, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,4,2]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 0.6, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

set.seed(seed)
# Nobs = 50 and Beta = 1
x <- rexp(n = 50, rate = 1)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 50 | Beta = ', beta[1]), col = "gray", ylim = c(0, max(fun)), xlab = '')
lines(x2, dexp(x2, rate = res$res_est_array[1,1,3]), col = 1, lty = 1, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,2,3]), col = 2, lty = 2, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,3,3]), col = 3, lty = 3, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[1,4,3]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 0.6, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

####
set.seed(seed)
# Nobs = 10 and Beta = 0.5
x <- rexp(n = 10, rate = 0.5)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 10 | Beta = ', beta[2]), col = "gray", ylim = c(0, max(fun)), xlab = '')
lines(x2, dexp(x2, rate = res$res_est_array[2,1,1]), col = 1, lty = 1, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,2,1]), col = 2, lty = 2, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,3,1]), col = 3, lty = 3, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,4,1]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 0.6, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

set.seed(seed)
# Nobs = 30 and Beta = 0.5
x <- rexp(n = 30, rate = 0.5)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 30 | Beta = ', beta[2]), col = "gray", ylim = c(0, max(fun)), xlab = '')
lines(x2, dexp(x2, rate = res$res_est_array[2,1,2]), col = 1, lty = 1, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,2,2]), col = 2, lty = 2, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,3,2]), col = 3, lty = 3, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,4,2]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 0.6, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

set.seed(seed)
# Nobs = 50 and Beta = 0.5
x <- rexp(n = 50, rate = 0.5)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 50 | Beta = ', beta[2]), col = "gray", ylim = c(0, max(fun)), xlab = '')
lines(x2, dexp(x2, rate = res$res_est_array[2,1,3]), col = 1, lty = 1, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,2,3]), col = 2, lty = 2, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,3,3]), col = 3, lty = 3, lwd = 1)
lines(x2, dexp(x2, rate = res$res_est_array[2,4,3]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 0.6, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)







