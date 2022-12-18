bootstrap_uni <- function(nobs = c(10, 30, 50), B = c(50, 100, 500), theta = c(1, 2), seed = 123, 
                          type = 'non-parametric', view_results = T, export = T){
  
  if(!is.numeric(B) | !is.numeric(nobs) | !is.numeric(theta) |
     !is.numeric(seed) |!is.logical(view_results) |
     !is.logical(export)){
    stop('Some arguments entered in the function are wrong. Follow the instructions: \n
     * "nobs" -> These are the sample sizes that will be evaluated in the simulations, it must be a numerical vector.\n
     * "B" -> Represents the number of bootstrap sample size, must be a numeric argument.\n
     * "theta" -> These are the parameters of Uniform distribution that you want to estimate in you simulation. It must be a numerical vector.\n
     * "view_results" -> If you prefer to print the simulation results for each sample size and parameterization.\n
     * "export" -> Whether to save the results in .txt file.\n')
  }
  
  set.seed(seed)
  
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
      x <- rexp(n = nobs[i], rate = theta[j])
      n <-length(x)
      est_emv_exp <- max(x)
      sd_emv_exp <- sqrt(((n*theta[j])^2)/(((n-1)^2)*(n-2)))
      res_est_array[j, 1, i] <- est_emv_exp
      res_bias_array[j, 1, i] <- theta[j] - est_emv_exp
      res_sd_array[j, 1, i] <- sd_emv_exp
      for(k in 1:length(B)){
        vec_est_B <- numeric(B[k])
        for(l in 1:B[k]){
          resample_x <- sample(x, replace = T)
          vec_est_B[l] <- 1/mean(resample_x)
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
    name_file_export<- paste0('Simulation results uniform.txt')
    file_sink <- file(name_file_export)
    sink(name_file_export, append = T)
    cat('\n')
    cat("Uniform Distribution \n")
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

res <- bootstrap_uni(nobs = c(10, 30, 50), B = c(50, 100, 500), theta = c(1, 2), seed = 1987, 
                     type = 'non-parametric', view_results = F, export = T)

res$res_est_array
