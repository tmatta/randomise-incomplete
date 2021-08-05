#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~ 2_mcar1_sim.R
#~~~
#~~~ Simulation design:  
#~~~ Dropout is covariate-dependent missing based on time and group
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(mvtnorm)

rm(list = ls())
set.seed(675635)

reps <- 100
miss_output <- array(NA, dim = c(reps, 2, 2), 
                    dimnames = list(1:reps, c("x0", "x1"), c("t1", "t2")))

mod_output <- matrix(NA, nrow = reps, ncol = 3)

#--- Sample -------------------------------------------------------------------#
# i = 1, ..., N subjects
# j = 1, ..., n_i observations per subject i
#------------------------------------------------------------------------------#
n <- 80*2
j <- 3    

#--- Fixed effects ------------------------------------------------------------#
b0 <- 20      # intercept (same for both groups) 
b1 <-  2      # linear slope group 1 
b2 <-  1      # treatment difference 
b_vec <- c(b0, b1, b2)   # vector of fixed effects

#--- Residual Variances -------------------------------------------------------#
#--- e ~ N (0, s2)
sd_e1 <- 1
s2_mat <- diag(sd_e1^2, j)  

#--- Fixed effect design matrix -----------------------------------------------#
#                  b0         b1           b2              
x_mat_1 <- cbind(rep(1, j), seq(1:j)-1, seq(1:j)-1)
x_mat_0 <- cbind(rep(1, j), seq(1:j)-1, rep(0, j))

#--- Z Matrix -----------------------------------------------------------------#
z_mat <- cbind(rep(1, j), seq(1:j)-1)

#--- tau Matrix ---------------------------------------------------------------#
sd_i <-  1.25
sd_s <-  0.75
r_is <- -0.3  # cor(y_int, y_slp) 

tau <- diag(c(sd_i^2, sd_s^2))      # variances
lower_tau <- c(r_is * sd_i * sd_s)  # covariance
tau[lower.tri(tau)] <- lower_tau
tau <- tau + t(tau) - diag(diag(tau))

#--- Means and variances
y_var <- (z_mat %*% tau %*% t(z_mat)) + s2_mat
y_mu_1 <- x_mat_1 %*% b_vec
y_mu_0 <- x_mat_0 %*% b_vec

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~ BEGIN SIMULATION                                                       ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for(rr in 1:reps){

  #=== Generate Longitudinal Data ===========================================#
  y_1 <- rmvnorm(n / 2, mean = y_mu_1, sigma = y_var)
  y_0 <- rmvnorm(n / 2, mean = y_mu_0, sigma = y_var)
          
  #--- Name data 
  y_1_w <- data.frame(y_1)
  y_0_w <- data.frame(y_0)
      
  colnames(y_1_w) <- ynames <- paste0("y", seq(1:ncol(y_1_w))-1)
  colnames(y_0_w) <- ynames <- paste0("y", seq(1:ncol(y_0_w))-1)
      
  #--- ID variable, works for balanced groups
  id_var <- 1:n
  n_g = 2
  id_group <- split(id_var, sort(id_var %% n_g))
          
  y_1_w$id <- id_group[[1]]
  y_0_w$id <- id_group[[2]]
      
  #--- Covariates
  y_1_w$x <- 1     
  y_0_w$x <- 0     
      
  y_w <- rbind(y_1_w, y_0_w)
  rownames(y_w) <- NULL
      
  #--- Reshape to long & reorder 
  y_l <- reshape(y_w, varying = ynames, v.names = "y", timevar = "time", 
                 times = (1 : j) - 1, direction = "long") 
  
  y_l <- y_l[order(y_l$id, y_l$time), ]
  rownames(y_l) <- NULL
  
  y_l$x_t <- y_l$x * y_l$time 
  y_l$index <- rownames(y_l)
  
  #--- Final full data data.frame
  y_l <- y_l[, c("index", "id", "time", "x", "x_t", "y")]
  
  #=== Generate Dropout =====================================================#
  a0 <- log(0.10 / (1 - 0.10)) # Log-odds of event at time 1
  a1 <- log(0.15 / (1 - 0.15)) # Log-odds of event at time 2
  a2 <- log(0.99 / (1 - 0.99)) # probability of event at time 3 (1)
  a3 <- log(.65 / (1 - .65))   # Log-odds of dropping out for x = 1
  
  #--- Create a factor wave variable
  d <- y_l
  d$t <- as.factor(d$time)
  d <- cbind.data.frame(d[ , -3],  model.matrix(~ -1  + time + t, data = d))
  d$t <- NULL
  
  #--- Hazard is conditional on the event not occurring at a prior time.
  #--- haz(t_j) = Pr(T_i = j | T_i >= j)
  d$haz <- 1 / (1 + exp(-(a0*d$t0 + a1*d$t1 + a2*d$t2 + a3*d$x)))
  
  d$r <- 0
  d$d <- NA
  #--- For each subject, conditional on not dropping out previously,
  #--- We do it like this to limit the runif() statement to only those left in the risk set.
  #--- For subject 1
  if (d$time[1] == 0) {
    d$r[1] <- runif(1); d$d[1] <- ifelse(d$haz[1] > d$r[1], 1, 0)
  } 
  if (d$time[2] == 1 & d$d[1] == 0) {
    d$r[2] <- runif(1); d$d[2] <- ifelse(d$haz[2] > d$r[2], 1, 0)
  } 
  if (d$time[2] == 1 & d$d[1] != 0) {
    d$r[2] <- NA; d$d[2] <- 0
  } 
  if (d$time[3] == 2 & sum(d$d[2:1]) == 0) {
    d$r[3] <- NA; d$d[3] <- 1
  }
  if (d$time[3] == 2 & sum(d$d[2:1]) != 0) {
    d$r[3] <- NA; d$d[3] <- 0
  }
  
  #--- For subject 2:N
  for(i in 4:nrow(d)) {
    if (d$time[i] == 0) {
      d$r[i] <- runif(1); d$d[i] <- ifelse(d$haz[i] > d$r[i], 1, 0)
    }   
    if (d$time[i] == 1 & d$d[i - 1] == 0) {
      d$r[i] <- runif(1); d$d[i] <- ifelse(d$haz[i] > d$r[i], 1, 0)
    } 
    if (d$time[i] == 1 & d$d[i - 1] != 0) {
      d$r[i] <- NA; d$d[i] <- 0
    } 
    if (d$time[i] == 2 & sum(d$d[(i-2):(i-1)]) == 0) {
      d$r[i] <- NA; d$d[i] <- 1
    }
    if (d$time[i] == 2 & sum(d$d[(i-2):(i-1)]) != 0) {
      d$r[i] <- NA; d$d[i] <- 0
    }  
  }  
  
  #--- Create event occurrence variable
  d <- d[order(d$id, d$time), ]
  rownames(d) <- NULL 
  d$event <- ifelse(d$time == 0 & d$d == 1, 1,
               ifelse(d$time == 1 & d$d == 1, 1, 0)) 
  
  #--- Keep 1 observation per id  
  d1 <- d[which(d$d==1), c("id", "x", "time", "haz", "d", "event")]
  d1$time <- ifelse(d1$time == 2, 1, d1$time) 
  d1$time <- d1$time + 1
  
  #--- Expand data frame based on duration
  d2 <- d1[rep(d1$id, d1$time), ] 
  colnames(d2) <- c("id", "x", "dur", "haz", "d", "event")
      
  #--- Create unstructured data
  d2$t1 <- as.numeric(row.names(d2))           # Create variable from row name      
  rownames(d2) <- NULL
  d2$t2 <- floor(d2$t1)                        # Take the integer from the decimal 
  d2$t <-  round((d2$t1 - d2$t2) * 10 + 1, 0)  # Subtract the int from the decimal
  d2$t_num <- d2$t
  d2$t <- as.factor(d2$t)
  d2$r <- ifelse(d2$event == 1 & d2$dur == d2$t_num, 1, 0)  
  d2$time <- d2$t_num
  
  d3 <- d2[, c("id", "x", "haz", "time", "r")]
  d3 <- d3[order(d3$id, d3$time), ]
  d3$t_fac <- as.factor(d3$time)
  
  d_df <- cbind.data.frame(d3[ , -4],  model.matrix(~ -1  + time + t_fac, data = d3))
  d_df$t_fac <- NULL
  
  m1 <- glm(r ~ t_fac2 + x, data = d_df, family = "binomial")
  

  miss_output[rr, 1, 1] <- length(which(d_df$time==1 & d_df$x==0 & d_df$r == 1)) / length(which(d_df$time==1 & d_df$x==0)) 
  miss_output[rr, 1, 2] <- length(which(d_df$time==2 & d_df$x==0 & d_df$r == 1)) / length(which(d_df$time==2 & d_df$x==0)) 

  miss_output[rr, 2, 1] <- length(which(d_df$time==1 & d_df$x==1 & d_df$r == 1)) / length(which(d_df$time==1 & d_df$x==1)) 
  miss_output[rr, 2, 2] <- length(which(d_df$time==2 & d_df$x==1 & d_df$r == 1)) / length(which(d_df$time==2 & d_df$x==1)) 
  
  mod_output[rr, ] <- c(coef(m1)[1], sum(coef(m1)[1:2]), coef(m1)[3])

}

dropout_mis_output <- miss_output
dropout_mod_output <- mod_output

save(dropout_mis_output, file = "1_dropout_mis_output.rda")
save(dropout_mod_output, file = "1_dropout_mod_output.rda")



round(colMeans(mod_output[,1:3]), 2)
round(c(a0, a1, a3), 2)

round(colMeans(dropout_mis_output[,1, 1:2]), 2)
round(1 / (1 + exp(-(a0))), 2)      # probability of event at time 1
round(1 / (1 + exp(-(a1))), 2)      # probability of event at time 2

round(colMeans(dropout_mis_output[,2, 1:2]), 2)
round(1 / (1 + exp(-(a0 + a3))), 2)      # probability of event at time 1
round(1 / (1 + exp(-(a1 + a3))), 2)      # probability of event at time 2





