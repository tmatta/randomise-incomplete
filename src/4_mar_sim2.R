#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~ 4_mar_sim.R
#~~~
#~~~ Simulation design:  
#~~~   3 waves, balanced groups
#~~~   Dropout is dependent on time and previous y (as a deviation from the y-bar) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(mvtnorm)
library(lme4)

rm(list = ls())
set.seed(675635)

reps <- 5000

#--- Dropout Probabilities ----------------------------------------------------#

t_labs <- c("t1", "t2", "t3")
par_labs <- c("a2", "int", "slp", "grp_slp", "int_se", "slp_se", "grp_slp_se", "sd_int", "sd_slp", "int_slp_cor", "resid")

#=== Simulation Arrays ========================================================#

mod_output <- array(NA, dim = c(reps, 2, length(par_labs)), 
                    dimnames = list(1:reps, c("av", "cc"), par_labs))

#--- Sample -------------------------------------------------------------------#
# i = 1, ..., N subjects
# j = 1, ..., n_i observations per subject i
#------------------------------------------------------------------------------#
n <- 100*2
j <- 3    

#--- Fixed effects ------------------------------------------------------------#
b0 <-  5.50      # intercept (same for both groups) 
b1 <- -1.75      # linear slope group 1 
b2 <-  0.75      # treatment difference 
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
sd_s <-  0.50
r_is <- -0.30  # cor(y_int, y_slp) 

tau <- diag(c(sd_i^2, sd_s^2))      # variances
lower_tau <- c(r_is * sd_i * sd_s)  # covariance
tau[lower.tri(tau)] <- lower_tau
tau <- tau + t(tau) - diag(diag(tau))

#--- Means and variances
y_var <- (z_mat %*% tau %*% t(z_mat)) + s2_mat
y_mu_1 <- x_mat_1 %*% b_vec
y_mu_0 <- x_mat_0 %*% b_vec

y_mu <- (y_mu_1 + y_mu_0) / 2

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

  # create z-score based on marginal mean at each time
  y_l$zy <- ifelse(y_l$time == 0, y_l$y - y_mu[1], 
             ifelse(y_l$time == 1, y_l$y - y_mu[2], 
              ifelse(y_l$time == 2, y_l$y - y_mu[3], NA)))
                         
  y_l$index <- rownames(y_l)

  #--- Final full data data.frame
  y_l <- y_l[, c("index", "id", "time", "x", "x_t", "y", "zy")]

  #=== Generate Dropout =====================================================#
  #--- Log-odds, as intercepts 
  i0 <- log(0.10 / (1 - 0.10))
  i1 <- log(0.15 / (1 - 0.15))
  # i2 <- log(or_drop[[jj]])
  i2 <- log(runif(1, min = 1, max = 3.5))

  #--- Logistic Regression Parameters
  a0 <- i0 
  a1 <- i1 - i0 
  a2 <- i2

  d <- y_l

  #--- Create a factor wave variable
  d$t <- as.factor(d$time)
  d_mat <- cbind.data.frame(d,  model.matrix(~ -1  + time + t, data = d))

  #--- Hazard is conditional on the event not occurring at a prior time.
  #--- haz(t_j) = Pr(T_i = j | T_i >= j)
  d$haz <- 1 / (1 + exp(-(a0 + a1*d_mat$t1 + a2*d_mat$zy)))

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

  #--- Collect prportion of missing data at each time per z=0 and z=1
  #miss_output[rr, ,jj] <- prop.table(table(d$time, d$d), 1)[, 2]

  #--- order d and reset row names 
  d <- d[order(d$id, d$time), ]
  rownames(d) <- NULL 

  #--- Keep 1 observation per subject, the wave of dropout
  d1 <- d[, c("id", "time", "d")]
  d2 <- d1[which(d1$d==1), c("id", "time")]
  d2$time <- d2$time + 1 
    
  #--- Expand data frame based on duration
  d3 <- d2[rep(d2$id, d2$time), ] 
  colnames(d3) <- c("id", "dur")
    
  #--- Create time unstructured data
  d3$t1 <- as.numeric(row.names(d3))           # Create variable from row name      
  rownames(d3) <- NULL
  d3$t2 <- floor(d3$t1)                        # Take the integer from the decimal 
  d3$t <-  round((d3$t1 - d3$t2) * 10 + 1, 0)  # Subtract the int from the decimal
  d3$t_num <- d3$t
  d3$t <- as.factor(d3$t)
  d3$r <- ifelse(d3$dur == d3$t_num, 1, 0)  
  d3$time <- d3$t_num - 1 
  d3 <- d3[, c("id", "time", "r")]
  d3 <- d3[order(d3$id, d3$time), ]

  #--- baseline hazard
  y_miss <- merge(y_l, d3, by = c("id", "time"))
  y_miss <- y_miss[order(y_miss$id, y_miss$time), ]

  y_cc <- y_miss$id[which(y_miss$time == 2)]

  #=== Estimation ===========================================================#
  #--- Model 1: Available Data
  m_av <- lmer(y ~ 1 + time + x_t + (1 + time | id), data = y_miss)
  m_av_fe <- as.numeric(fixef(m_av))      
  m_av_se <- round(sqrt(diag(vcov(m_av))), 2)
  m_av_vc <- as.numeric(as.data.frame(VarCorr(m_av))[, 5])  
    
  #--- Model 2: Complete-Case Data
  m_cc <- lmer(y ~ 1 + time + x_t + (1 + time | id), data = y_miss[which(y_miss$id %in% y_cc), ])
  m_cc_fe <- as.numeric(fixef(m_cc))  
  m_cc_se <- round(sqrt(diag(vcov(m_cc))), 2)   
  m_cc_vc <- as.numeric(as.data.frame(VarCorr(m_cc))[, 5])  

  mod_output[rr, 1, ] <- c(a2, m_av_fe, m_av_se, m_av_vc) 
  mod_output[rr, 2, ] <- c(a2, m_cc_fe, m_cc_se, m_cc_vc) 
}


mar_mod_output <- mod_output

save(mar_mod_output, file = "4_mar_mod_output2.rda")



round(colMeans(miss_output[, 1:3,  1]), 2)

round(colMeans(mod_output[, 1, 1]), 2)
round(colMeans(mod_output[, 2, 11:5]), 2)

