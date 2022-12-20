bootstrap_uni <- function(nobs = c(10, 30, 50), B = c(50, 100, 500), theta = c(1, 0.5), seed = 1987, 
                          type = 'non-parametric', view_results = T, export = T){
  
  if(!is.numeric(B) | !is.numeric(nobs) | !is.numeric(theta) |
     !is.numeric(seed) |!is.logical(view_results) |
     !is.logical(export) | (length(which(type == c('non-parametric', 'parametric')))==0 | length(type) > 1)){
    stop('Some arguments entered in the function are wrong. Follow the instructions: \n
     * "nobs" -> These are the sample sizes that will be evaluated in the simulations, it must be a numerical vector.\n
     * "B" -> Represents the number of bootstrap sample size, must be a numeric argument.\n
     * "theta" -> These are the parameters of Uniform distribution that you want to estimate in you simulation. It must be a numerical vector.\n
     * "view_results" -> If you prefer to print the simulation results for each sample size and parameterization.\n
     * "export" -> Whether to save the results in .txt file.\n')
  }
  
  if(!require(MASS)){install.packages('MASS')}
  
  res_est_b_list <- list()
  res_est_b_list_full<-list()
  for(i in 1:length(theta)){
    for(j in 1:length(B)){
      res_est_b_list[[j]] <- `colnames<-`(matrix(c(0), nrow = B[j], ncol = length(nobs)), c(paste0('Nobs =  ', nobs)))
      
    }
    names(res_est_b_list) <- c(paste0('B = ', B))
    res_est_b_list_full[[i]] <- res_est_b_list
  }
  names(res_est_b_list_full) <- c(paste0('theta = ', theta))
  
  # Setting the environment to receive the results
  res_est_array <- array(c(0), dim = c(length(theta), length(B)+1, length(nobs)), 
                         dimnames = list(paste0('theta = ', theta),
                                         c('EMV', paste0('Est.B ', B)),
                                         paste0('nobs = ', nobs)))
  
  res_bias_array <- array(c(0), dim = c(length(theta), length(B)+1, length(nobs)), 
                          dimnames = list(paste0('theta = ', theta),
                                          c('Bias EMV', paste0('Bias Est.B ', B)),
                                          paste0('nobs = ', nobs)))
  
  res_sd_array <- array(c(0), dim = c(length(theta), length(B)+1, length(nobs)), 
                        dimnames = list(paste0('theta = ', theta),
                                        c('Sd EMV', paste0('Sd Est.B ', B)),
                                        paste0('nobs = ', nobs)))
  #Starts the simulation
  for(i in 1:length(nobs)){
    for(j in 1:length(theta)){
      set.seed(seed)
      x <- runif(n = nobs[i], min = 0, max = theta[j])
      n <-length(x)
      est_emv_uni <- max(x)
      sd_emv_uni <- sqrt(((n*theta[j])^2)/(((n+1)^2)*(n+2)))
      res_est_array[j, 1, i] <- est_emv_uni
      res_bias_array[j, 1, i] <- theta[j] - est_emv_uni
      res_sd_array[j, 1, i] <- sd_emv_uni
      for(k in 1:length(B)){
        vec_est_B <- numeric(B[k])
        for(l in 1:B[k]){
          if(type == 'non-parametric'){
            resample_x <- sample(x, replace = T)
          }else{
            est_x <- fitdistr(x, dunif, list(min=0, max = 2), lower = c(-Inf, -Inf))
            resample_x <- runif(length(x), min=0, max = est_x$estimate[2])
          }
          vec_est_B[l] <- max(resample_x)
        }
        res_est_b_list_full[[j]][[k]][,i] <- vec_est_B
        res_est_array[j, k+1, i] <- mean(vec_est_B)
        res_bias_array[j, k+1, i] <- theta[j] - mean(vec_est_B)
        res_sd_array[j, k+1, i] <-  sd(vec_est_B)
        
      }    
    }
  }
  
  # if you want to print all the results in the terminal
  if(view_results){
    cat('\n')
    cat("Uniform Distribution \n")
    cat('\n')
    cat(type)
    cat('\n')
    for(i in 1:length(nobs)){
      cat("####################################################\n")
      cat("Sample size:", nobs[i], '\n')
      cat("####################################################\n")
      for(j in 1:length(theta)){
        cat('\n')
        cat("Results for parameter", rownames(res_est_array)[j],"\n")
        cat("theta:", theta[j],'\n')
        cat("Estimate EMV:", res_est_array[j, 1, i],' \n')
        cat("BIAS EMV:", res_bias_array[j, 1, i],' \n')
        cat("Standard deviation EMV:", res_sd_array[j, 1, i],' \n')
        cat('\n')
        for(k in 1:length(B)){
          cat("Bootstrap sample size:", B[k], '\n')
          cat("Est.B", B[k], ':', res_est_array[j, k+1, i],' \n')
          cat("BIAS Est.B", B[k], ':', res_bias_array[j, k+1, i],' \n')
          cat("Sd Est.B", B[k],':', res_sd_array[j, k+1, i],' \n')
          cat('\n')
        }
        cat("####################################################\n")
      }
      cat('\n')
    }
  }
  
  # if you want to export the results saving as .txt file
  if(export){
    name_file_export <- paste0('Simulation results uniform', ' - ', type, '.txt')
    file_sink <- file(name_file_export)
    sink(name_file_export, append = T)
    cat('\n')
    cat("Uniform Distribution \n")
    cat('\n')
    cat(type)
    cat('\n')
    for(i in 1:length(nobs)){
      cat("####################################################\n")
      cat("Sample size:", nobs[i], '\n')
      cat("####################################################\n")
      for(j in 1:length(theta)){
        cat('\n')
        cat("Results for parameter", rownames(res_est_array)[j],"\n")
        cat("theta:", theta[j],'\n')
        cat("Estimate EMV:", res_est_array[j, 1, i],' \n')
        cat("BIAS EMV:", res_bias_array[j, 1, i],' \n')
        cat("Standard deviation EMV:", res_sd_array[j, 1, i],' \n')
        cat('\n')
        for(k in 1:length(B)){
          cat("Bootstrap sample size:", B[k], '\n')
          cat("Est.B", B[k], ':', res_est_array[j, k+1, i],' \n')
          cat("BIAS Est.B", B[k], ':', res_bias_array[j, k+1, i],' \n')
          cat("Sd Est.B", B[k],':', res_sd_array[j, k+1, i],' \n')
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
       nobs = nobs, B = B, theta = theta)
  
}

res<-bootstrap_uni(nobs = c(10, 30, 50), B = c(50, 100, 500), theta = c(1, 2), seed = 1987, 
                   type = 'non-parametric', view_results = F, export = F)

res$res_est_array

##########################################
####### SECTION FOR THE GRAPHICS #########
##########################################

library(tidyverse)
theta<-res$theta
B<-res$B
nobs<-res$nobs

# Graphics of distributions of the estimatives for theta = 1
par(mfrow = c(3,3))
for(j in 1:length(B)){
  for(k in 1:length(nobs)){
    hist(res$res_est_b_list_full[[1]][[j]][,k], freq = F,
         main = paste0('Theta = ', theta[1],' | B = ', B[j], ' | N = ', nobs[k]), 
         xlab = NULL)
  }
}

# Graphics of distributions of the estimatives for theta = 0.5
par(mfrow = c(3,3))
for(j in 1:length(B)){
  for(k in 1:length(nobs)){
    hist(res$res_est_b_list_full[[2]][[j]][,k], freq = F,
         main = paste0('Theta = ', theta[2],' | B = ', B[j], ' | N = ', nobs[k]), 
         xlab = NULL)
  }
}

par(mfrow = c(2,3))
seed = 1987
set.seed(seed)
# Nobs = 10 and Theta = 1
x <- runif(n = 10, min = 0, max = 1)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 10 | Theta = ', theta[1]), col = "gray", xlab = '')
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,1,1]), col = 1, lty = 1, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,2,1]), col = 2, lty = 2, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,3,1]), col = 3, lty = 3, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,4,1]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 1, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)


set.seed(seed)
# Nobs = 30 and Theta = 1
x <- runif(n = 30, min = 0, max = 1)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 30 | Theta = ', theta[1]), col = "gray", xlab = '')
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,1,2]), col = 1, lty = 1, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,2,2]), col = 2, lty = 2, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,3,2]), col = 3, lty = 3, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,4,2]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 1, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

set.seed(seed)
# Nobs = 50 and Theta = 1
x <- runif(n = 50, min = 0, max = 1)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 50 | Theta = ', theta[1]), col = "gray", xlab = '')
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,1,3]), col = 1, lty = 1, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,2,3]), col = 2, lty = 2, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,3,3]), col = 3, lty = 3, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[1,4,3]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 1, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

####
set.seed(seed)
# Nobs = 10 and Theta = 2
x <- runif(n = 10, min = 0, max = 2)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 10 | Theta = ', theta[2]), col = "gray", xlab = '')
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,1,1]), col = 1, lty = 1, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,2,1]), col = 2, lty = 2, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,3,1]), col = 3, lty = 3, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,4,1]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 1, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

set.seed(seed)
# Nobs = 30 and Theta = 2
x <- runif(n = 30, min = 0, max = 2)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 30 | Theta = ', theta[2]), col = "gray", xlab = '')
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,1,2]), col = 1, lty = 1, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,2,2]), col = 2, lty = 2, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,3,2]), col = 3, lty = 3, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,4,2]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 1, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

set.seed(seed)
# Nobs = 50 and Theta = 2
x <- runif(n = 50, min = 0, max = 2)
x2 <- seq(min(x), max(x), length = length(x))
hist(x, prob = T,  main = paste0('N = 50 | Theta = ', theta[2]), col = "gray", xlab = '')
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,1,3]), col = 1, lty = 1, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,2,3]), col = 2, lty = 2, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,3,3]), col = 3, lty = 3, lwd = 1)
lines(x2, dunif(x2, min = 0, max = res$res_est_array[2,4,3]), col = 4, lty = 4, lwd = 1)
legend('topright', cex= 1, lty = c(1:4), col= c(1:4),
       legend=c("EMV", "EMV B = 100", "EMV B = 300", "EMV B = 500"),
       box.lty=0)

