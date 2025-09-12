# Load required libraries
library(dplyr)
library(mgcv)      # For s() smooth terms
library(survival)  # For coxph() survival model
library(GJRM)
library(VineCopula)
library(gamlss.dist)
library(MASS)

# ------------------------------------------------------------------------------
# Function: SurvPHBinSim
# Simulate covariates, LVI (binary treatment), survival time, and survival status
# ------------------------------------------------------------------------------

SurvPHBinSim <- function(n, cs, pb, e) {
  sigma <- 1  # Standard deviation for error
  
  # Simulate multivariate normal covariates
  Sigma <- matrix(0.5, 10, 10)
  diag(Sigma) <- 1
  cov <- MASS::mvrnorm(n, rep(0, 10), Sigma)
  cov <- pnorm(cov)  # Transform to uniform via normal CDF
  
  # Assign covariates
  x1 <- abs(cov[, 1])
  x2 <- abs(cov[, 2])
  x3 <- abs(cov[, 3])
  x4 <- round(cov[, 4])
  x5 <- cut(cov[, 5], breaks = 3, labels = c("1", "2", "3"), include.lowest = TRUE)
  x5 <- as.numeric(x5)
  x6 <- round(cov[, 6])
  x7 <- abs(cov[, 7])
  x8 <- abs(cov[, 8])
  x9 <- cov[, 9]
  x10 <- cov[, 10]
  
  # Nonlinear function for x2
  f1 <- function(x) cos(pi * 2 * x) + sin(pi * x)
  
  # Simulate errors with specified copula dependence
  u <- BiCopSim(n, 1, e)
  
  # Generate binary treatment y1 (LVI)
  y1 <- ifelse(0.9 * x1 + f1(x2) + 1.4 * x3 + 0.7 * x4 - 1.4 * x5 + 0.6 * x6 + qnorm(u[, 1]) > 0, 1, 0)
  
  # True coefficients for log survival time
  beta <- c(-1.39, -1.31, -1.4, 1.12, 1.6, -1.23)
  log_t <- beta[1] * y1 + beta[2] * x3 + beta[3] * x7 + beta[4] * x8 + beta[5] * x9 + beta[6] * x10 + qnorm(u[, 2])
  t <- exp(log_t)  # Actual survival time
  
  # Generate right-censoring variable to control censoring rate
  c <- runif(n, quantile(t, cs), quantile(t, cs + 0.2))
  time <- pmin(c, t)
  status <- as.numeric(t <= c)
  
  # Generate binary outcome y2 based on survival beyond a fixed cutoff
  cutoff <- as.numeric(quantile(time, probs = pb))
  y2 <- rep(0, n)
  y2[status == 0 & time >= cutoff] <- 1
  y2[status == 0 & time < cutoff] <- NA
  y2[status == 1 & time > cutoff] <- 1
  y2[status == 1 & time <= cutoff] <- 0
  
  # Return data frame
  df <- data.frame(
    time = time, status = status,
    x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6,
    x7 = x7, x8 = x8, x9 = x9, x10 = x10,
    y1 = y1, y2 = y2, logt = log_t,
    e1 = qnorm(u[, 1]), e2 = sigma * qnorm(u[, 2])
  )
  return(df)
}

# ------------------------------------------------------------------------------
# Check convergence of GJRM model
# ------------------------------------------------------------------------------

check_convergence <- function(out) {
  conv_output <- capture.output(conv.check(out))
  return(any(grepl("Observed information matrix is positive definite", conv_output)))
}

# ------------------------------------------------------------------------------
# Simulation parameters
# ------------------------------------------------------------------------------

n_sim <- 500
sample_sizes <- c(200, 500, 1000)
pb_values <- c(0.25, 0.5, 0.75)
cs_values <- c(0.8, 0.75, 0.67, 0.58, 0.48, 0.39, 0.29)
true_value <- -1.39
results_df <- data.frame()

# ------------------------------------------------------------------------------
# Main simulation loop
# ------------------------------------------------------------------------------

for (n in sample_sizes) {
  for (pb in pb_values) {
    for (cs in cs_values) {
      
      valid_sim_two_stage <- valid_sim_copula <- valid_sim_survival <- 0
      y1_coefficients_two_stage <- y1_coefficients_copula <- y1_coefficients_survival <- numeric(n_sim)
      current_sim <- 1
      
      while (current_sim <= n_sim) {
        simulation_result <- tryCatch({
          df_sample <- SurvPHBinSim(100000, cs, pb, e)  ########################### e = 0，0.5 ，0.7 for adjusting error terms correlation 
          df <- df_sample %>% filter(!is.na(y2)) %>% sample_n(n)
          
          # --- Copula-based method --- 
          f.list <- list(
            eq1 = y1 ~ s(x1) + s(x2) + s(x3) + x4 + x5 + x6,
            eq2 = y2 ~ y1 + s(x3) + s(x7) + s(x8) + s(x9) + s(x10)
          )
          out <- gjrm(f.list, data = df, margins = c("probit", "probit"),
                      model = "B", copula = "N", uni.fit = TRUE)
          
          if (!check_convergence(out)) stop("Model not converged")
          coef_y1_copula <- summary(out)$tableP2[2, 1]
          if (coef_y1_copula < -5 || coef_y1_copula > 5) stop("y1 estimate out of bounds")
          
          # --- Two-stage method ---
          model_eq1 <- glm(y1 ~ x1 + x2 + x3 + x4 + x5 + x6, family = binomial(), data = df)
          df$y1_hat <- predict(model_eq1, type = "response")
          model_eq2 <- glm(y2 ~ y1_hat + x3 + x7 + x8 + x9 + x10, family = binomial(), data = df)
          coef_y1_two_stage <- coef(model_eq2)["y1_hat"]
          
          # --- One-stage survival (Cox) ---
          surv_model <- coxph(Surv(time, status) ~ y1 + x3 + x7 + x8 + x9 + x10, data = df)
          coef_y1_survival <- coef(surv_model)["y1"]
          
          list(
            coef_two_stage = coef_y1_two_stage,
            coef_copula = coef_y1_copula,
            coef_survival = coef_y1_survival
          )
        }, error = function(e) {
          cat("Simulation", current_sim, "skipped:", e$message, "\n")
          return(NULL)
        })
        
        if (!is.null(simulation_result)) {
          y1_coefficients_two_stage[current_sim] <- simulation_result$coef_two_stage
          y1_coefficients_copula[current_sim] <- simulation_result$coef_copula
          y1_coefficients_survival[current_sim] <- simulation_result$coef_survival
          
          valid_sim_two_stage <- valid_sim_two_stage + 1
          valid_sim_copula <- valid_sim_copula + 1
          valid_sim_survival <- valid_sim_survival + 1
          
          cat("Simulation", current_sim, "done | n =", n, "pb =", pb, "cs =", cs, "\n")
          current_sim <- current_sim + 1
        } else {
          cat("Simulation", current_sim, "failed, retrying...\n")
        }
      }
      
      # Compute bias
      bias_two_stage <- mean(y1_coefficients_two_stage, na.rm = TRUE) - true_value
      bias_copula <- mean(y1_coefficients_copula, na.rm = TRUE) - true_value
      bias_survival <- -mean(y1_coefficients_survival, na.rm = TRUE) - true_value
      
      # Append results
      results_df <- rbind(results_df, data.frame(
        Sample_Size = n,
        PB = pb,
        cs_value = cs,
        Mean_Coef_Two_Stage = mean(y1_coefficients_two_stage, na.rm = TRUE),
        Bias_Two_Stage = bias_two_stage,
        Abs_Bias_Two_Stage = abs(bias_two_stage),
        Mean_Coef_Copula = mean(y1_coefficients_copula, na.rm = TRUE),
        Bias_Copula = bias_copula,
        Abs_Bias_Copula = abs(bias_copula),
        Mean_Coef_Survival = mean(y1_coefficients_survival, na.rm = TRUE),
        Bias_Survival = bias_survival,
        Abs_Bias_Survival = abs(bias_survival)
      ))
    }
  }
}

# Print and save results
print(results_df)

