### @@@ (half) SurvPHBinSim 

library(GJRM)
library(VineCopula)
library(gamlss.dist)
library(MASS)
library(dplyr)


SurvPHBinSim <- function(n,cs,pb){
  
  sigma <- 1  # standard error
  
  Sigma <- matrix(0.5, 10, 10)
  diag(Sigma) <- 1
  cov <- MASS::mvrnorm(n, rep(0, 10), Sigma)
  cov <- pnorm(cov)  # 
  
  
  x1 <- abs(cov[, 1])  # 
  x2 <- abs(cov[, 2])  #
  x3 <- abs(cov[, 3])  # 
  x4 <- round(cov[, 4]) #
  x5 <- cut(cov[, 5],   
            breaks = 3, 
            labels = c("1", "2", "3"), 
            include.lowest = TRUE)
  x5 <- as.numeric(x5, levels = c("1", "2", "3"))
 
  x6 <- round(cov[, 6])
  x7 <- abs(cov[, 7])  # 
  x8 <- abs(cov[, 8])  # 
  x9 <- cov[, 9]      # 
  x10 <- cov[, 10]  #
  
  # define function f1 for variable x2
  f1 <- function(x) cos(pi * 2 * x) + sin(pi * x)
  
  
  #######################################################
  # simulated error term: u
  # log survival time t_i has a normal distribution
  # BiCopSim(n,1,0.5)
  # 1 is normal; 0.5 is correlation between error terms
  ########################################################
  

  u <- BiCopSim(n, 1, 0.5) 
  
  ####################################
  # simulatied binary variable: y1 (LVI)
  ######################################
  
  
  
  y1 <- ifelse(0.9 * x1 + f1(x2) + 
                 1.4 * x3 + ## 1.2
                 0.7*x4 + #clinical
                 #-1.08 * x5_2 + #####-1.08
                 #-0.7 * x5_3 +#### -0.7
                 -1.08 * x5+
                 0.6 * x6 +  
                 + qnorm(u[, 1]) > 0, 1, 0 )
  #y1 <- ifelse(-1.55  +0 * x1 + f1(x2) + 
  #             0 *x4  +
  #           -1.08 * x5 + #####-1.08
  #         0 * x6 +  
  #       + qnorm(u[, 1]) > 0, 1, 0)
  
  
  ####################################
  # simulated binary variable: y1 (LVI)
  ######################################
  
  beta <- c(0, -1.39, -1.31, -1.4, 1.12, 1.6,-1.23)
  #beta <- c(0, -1.39, -1.31, 0, 1.12, 0,0)
  
  ##aft model
  #time
  #error2=log(-log(1 - u[, 2])) #cloglog
  log_t <- beta[1] + beta[2] * y1  + beta[3]*x3 + beta[4]*x7 + beta[5] * x8 + beta[6]*x9 + beta[7]*x10 + qnorm(u[, 2])
  
  ########
  ### LVI: the increase of y1 will decrese the probability of survival 
  ########
  t<-exp(log_t)##### check this line
  
  
  # control censoring rate
  # 0.1 = 0.7
  # 0.3 = 0.5 
  # 0.5 = 0.3
  
  
  ## censoring independent from t
  c <- runif(n, quantile(t, cs), quantile(t, cs + .2)) ########censoring rate
  
  
  time = apply(cbind(c,t),1,min)#use diff name
  
  delta1 = as.numeric(t <= c)###delta1 ==1 not censored 
  
  
  df <- data.frame(
    id = 1:n,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x5 = x5,
    x6 = x6,
    x7 = x7,
    x8 = x8,
    x9 = x9,
    x10 = x10,
    y1 = y1,
    logt=log_t,
    time = time,  
    status = delta1
  )
  
  
  t<- df$time
  cutoff <- as.numeric(quantile(t, probs = pb))
  
  
  # Initialize t1c with default values (assuming similar to right-censoring)
  y2 <- rep(0, length(t))
  
  
  #####################
  ###Generate y2
  ####################
  
  
  # Alive (not an event) and time at or beyond cutoff
  y2[df$status == 0 & df$time >= cutoff] <- 1  # survived past cutoff
  
  # Alive (not an event) and time before cutoff
  y2[df$status == 0 & df$time < cutoff] <- NA  # left censoring
  
  # Deceased (event occurred) and time beyond cutoff
  y2[df$status == 1 & df$time > cutoff] <- 1  # Censored, survived past cutoff
  
  # Deceased (event occurred) and time at or before cutoff
  y2[df$status == 1 & df$time <= cutoff] <- 0  # Event observed within study period
  
  
  
  
  df <- data.frame(
    time = time,
    status = delta1,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x5 = x5,
    x6 = x6,
    x7 = x7,
    x8 = x8,
    x9 = x9,
    x10 = x10,
    y1 = y1,
    y2 = y2,
    logt=log_t,
    e1 = qnorm(u[,1]),
    e2 = sigma*qnorm(u[, 2])
  )

  return(df)
}






## calculate effect


calculate_effect <- function(y1,x3,x7,x8,x9,x10) {
  #Define beta coefficient
  beta <- c(0, -1.39, -1.31, -1.4, 1.12, 1.6,-1.23)
  #beta <- c(0, -1.39, -1.31, 0, 1.12, 0,0)
  
  beta[1] + beta[2] * y1  + beta[3]*x3 + beta[4]*x7 + beta[5] * x8 + beta[6]*x9 + beta[7]*x10
  
}




Surv_aft <- function(t, meanlog) {
  plnorm(t, meanlog = meanlog, sdlog = 1, lower.tail = FALSE)
}

check_convergence <- function(out) {
  conv_output <- capture.output(conv.check(out))
  return(any(grepl("Observed information matrix is positive definite", conv_output)))
}


### @@!!@@

Sumfunc <- function(j, f.list, marg, n, cs, pb) {
  
  success <- FALSE
  out <- NULL
  results <- NULL
  
  while (!success) {
    
    
    rm(list = c("df_samp", "df_clean", "df_copula", "df"), envir = .GlobalEnv)
    
    
    df_samp <- SurvPHBinSim(120000, cs, pb)
    assign("df_samp", df_samp, envir = .GlobalEnv)
    
    df_clean <- df_samp %>% filter(!is.na(y2))
    assign("df_clean", df_clean, envir = .GlobalEnv)
    
    df_copula <- df_clean %>% sample_n(n)
    assign("df_copula", df_copula, envir = .GlobalEnv)
    
    df <- df_clean %>% sample_n(30000) # oracle用的大样本
    assign("df", df, envir = .GlobalEnv)
    
    # 2) Copula 
    Cop <- c("N", "C0", "C90", "C180", "C270", 
             "GAL0", "GAL90", "GAL180", "GAL270", 
             "J0", "J90", "J180", "J270", 
             "G0", "G90", "G180", "G270", 
             "F", "AMH", "FGM", "T", "PL", "HO")[j]
    
    # 3) 
    out_try <- try(
      gjrm(f.list, data = df_copula, margins = marg, model = "B", copula = Cop, uni.fit = TRUE),
      silent = TRUE
    )
    
    if (!inherits(out_try, "try-error") && check_convergence(out_try)) {
      out <- out_try
      
      
      model_summary <- summary(out)
      
      y1_est <- model_summary$tableP2[2, 1]
      y1_se  <- model_summary$tableP2[2, 2]
      y1_value <- round(y1_est, 3)
      y1_ci <- c(round(y1_est - 1.96 * y1_se, 4), round(y1_est + 1.96 * y1_se, 4))
      
      kkk <- k.tau(out)
      tau_value <- as.numeric(kkk[1])
      low_tau   <- as.numeric(kkk[2])
      up_tau    <- as.numeric(kkk[3])
      
      at <- ATE(out, 'y1', n.sim = 1000)
      ate_value <- as.numeric(at$res["ATE"])
      low_ate   <- as.numeric(at$res["2.5%"])
      up_ate    <- as.numeric(at$res["97.5%"])
      
      # --- Oracle 部分：ATT + CI ---
      t0 <- as.numeric(quantile(df$time, pb))
      
      effect_1 <- calculate_effect(1, df$x3, df$x7, df$x8, df$x9, df$x10)
      effect_0 <- calculate_effect(0, df$x3, df$x7, df$x8, df$x9, df$x10)
      
      S_t0_y1_1 <- Surv_aft(t0, meanlog = effect_1) 
      S_t0_y1_0 <- Surv_aft(t0, meanlog = effect_0)
      
      att_est <- mean(S_t0_y1_1 - S_t0_y1_0)
      se_att_est <- sd(S_t0_y1_1 - S_t0_y1_0) / sqrt(nrow(df))
      z_critical <- qnorm(0.975)
      ci_lower <- att_est - z_critical * se_att_est
      ci_upper <- att_est + z_critical * se_att_est
      
     
      results <- c(
        ate_value, low_ate, up_ate,        
        att_est,  ci_lower, ci_upper,      
        AIC(out), BIC(out),                
        tau_value, low_tau, up_tau,        
        y1_value, y1_ci,                   
        sum(df$y2 == 0, na.rm = TRUE),     
        sum(df$y2 == 1, na.rm = TRUE),     
        sum(is.na(df$y2))                  
      )
      
      success <- TRUE  
    } else {
      warning("The model did not converge; regenerate the data and retry...")
    }
  }
  
  return(round(results, 4))
}





# Define the function to run Sumfunc multiple times and calculate average results
run_and_average_sumfunc <- function(m, j, f.list, marg, n, cs, pb) {
  # Determine the expected length of the result
  expected_length <- 17
  
  # Initialize a list to store the valid results
  results_list <- vector("list", m)
  valid_count <- 0
  iterations <- 0
  
  # Run Sumfunc to collect m valid results
  while (valid_count < m && iterations < m * 300) {  # Allow up to 2*m iterations to find m valid results
    result <- tryCatch({
      # Capture warnings and errors, return NULL for invalid results
      withCallingHandlers({
        res <- Sumfunc(j, f.list, marg, n, cs, pb)
        if (length(res) == expected_length) {
          res
        } else {
          NULL  # Skip this iteration if length mismatch
        }
      }, warning = function(w) {
        invokeRestart("muffleWarning")
        return(NULL)  # Skip this iteration if warning occurs
      })
    }, error = function(e) {
      return(NULL)  # Skip this iteration if error occurs
    })
    
    if (!is.null(result)) {
      valid_count <- valid_count + 1
      results_list[[valid_count]] <- result
      
      # Print progress every 10 valid results
      if (valid_count %% 10 == 0) {
        message(sprintf("Progress: %d valid results collected out of %d", valid_count, m))
      }
      

    }
    iterations <- iterations + 1
  }
  
  # Check if we have exactly m valid results, if not, extend the list with NA rows
  if (valid_count < m) {
    warning("Not enough valid results available. Padding with NAs to reach the required number of results.")
    for (i in (valid_count + 1):m) {
      results_list[[i]] <- rep(NA, expected_length)
    }
  }
  
  valid_results <- do.call(rbind, results_list)
  
  # Add true_tau variable based on the interval condition
  true_tau <- ifelse(valid_results[, 10] <= 0.33 & valid_results[, 11] >= 0.33, 1, 0)
  valid_results <- cbind(valid_results, true_tau)  # Append true_tau as a new column
  
  # Calculate the average values for each column in the valid_results matrix
  average_results <- colMeans(valid_results, na.rm = TRUE)
  
  # Define column names, including the new true_tau column and the iteration column
  col_names <- c("AT", "low_AT_ci", "up_AT_ci", "AAT", "low_AAT_ci", "up_AAT_ci",
                 "AIC", "BIC", "tau", "low_tau_ci", "up_tau_ci", "y1", "low_y1_ci", "up_y1_ci",
                 "0_y2", "1_y2", "na_y2", "true_tau")
  
  # Return the average results as a named vector with the column names
  named_average_results <- setNames(round(average_results, 4), col_names)  # Exclude the iteration column from averaging
  #named_average_results
  
  return(named_average_results)
  
}



j <- 1 # normal copula


n <- 200 # sample size
m <- 200# number of runs
marg <- c("probit", "probit") # marginals

f.list <- list(
  eq1 <- y1 ~ s(x1) + s(x2) + s(x3) + x4 + x5 + x6,
  eq2 <- y2 ~ y1 + s(x3) + s(x7) + s(x8) + s(x9) +s(x10)
)


cs <- 0.80 ## adjusting for censoring rate

# Initialize an empty data frame to store the results
results_df_200_30_0_all <- data.frame(pb = numeric(), average_results = numeric())

# Function to run and store results for different pb values
for (pb in seq(0.1, 0.9, by = 0.1)) {
  average_results <- run_and_average_sumfunc(m, j, f.list, marg, n, cs, pb)
  results_df_200_30_0_all <- rbind(results_df_200_30_0_all, data.frame(pb = pb, average_results =t(average_results) ))
}

# Print the results data frame
print(results_df_200_0_all)
