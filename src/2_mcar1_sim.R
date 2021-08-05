#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~ 2_mcar1_sim.R
#~~~
#~~~ Simulation design:  
#~~~   3 waves, 2 balanced groups
#~~~   Dropout is dependent on time and x 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(mvtnorm)
library(lme4)

rm(list = ls())
set.seed(675635)
reps <- 1000

#--- Sample -------------------------------------------------------------------#
# i = 1, ..., N subjects
# j = 1, ..., n_i observations per subject i
#------------------------------------------------------------------------------#
n <- 100*2
j <- 3    

#--- Dropout Probabilities ----------------------------------------------------#
pr_drop <- c(.55, .65, .75)
pr_drop_labs <- paste0("pr_", pr_drop*100)

t_labs <- c("t1", "t2")
par_labs <- c("int", "slp", "grp_slp", "int_se", "slp_se", "grp_slp_se", "sd_int", "sd_slp", "int_slp_cor", "resid")
mod_labs <- c("f1", "f2", "cc", "av1", "av2")

#=== Simulation Arrays ========================================================#
miss_output <- array(NA, dim = c(reps, 2, length(t_labs), length(pr_drop_labs)), 
                    dimnames = list(1:reps, c("x0", "x1"), t_labs, pr_drop_labs))

mod_output <- array(NA, dim = c(reps, length(mod_labs), length(par_labs), length(pr_drop_labs)), 
                    dimnames = list(1:reps, mod_labs, par_labs, pr_drop_labs))

df_list <- list()

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~ BEGIN SIMULATION                                                       ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for(jj in 1:length(pr_drop)){

  for(rr in 1:reps){

    #=== Generate Longitudinal Data ===========================================#
    y_1 <- mvtnorm::rmvnorm(n / 2, mean = y_mu_1, sigma = y_var)
    y_0 <- mvtnorm::rmvnorm(n / 2, mean = y_mu_0, sigma = y_var)
        
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
    #==========================================================================#
  
    #=== Generate Dropout =====================================================#
    a0 <- log(0.10 / (1 - 0.10)) # Log-odds of event at time 1
    a1 <- log(0.15 / (1 - 0.15)) # Log-odds of event at time 2
    a2 <- log(0.99 / (1 - 0.99)) # probability of event at time 3 (1)
    a3 <- log(pr_drop[jj] / (1 - pr_drop[jj]))   # Log-odds of dropping out for x = 1
    
    #--- Create a factor wave variable
    d <- y_l
    d$t <- as.factor(d$time)
    d <- cbind.data.frame(d[ , c("index", "id", "x", "x_t")],  model.matrix(~ -1  + time + t, data = d))
    d$t <- NULL
    # > head(d)
    #   index id x x_t time t0 t1 t2
    # 1     1  1 1   0    0  1  0  0
    # 2     2  1 1   1    1  0  1  0
    # 3     3  1 1   2    2  0  0  1

    #--- Hazard is conditional on the event not occurring at a prior time.
    #--- haz(t_j) = Pr(T_i = j | T_i >= j)
    d$haz <- 1 / (1 + exp(-(a0*d$t0 + a1*d$t1 + a2*d$t2 + a3*d$x)))
    d$r <- 0
    d$d <- NA
    # > head(d)
    #   index id x x_t time t0 t1 t2       haz r  d
    # 1     1  1 1   0    0  1  0  0 0.1710526 0 NA
    # 2     2  1 1   1    1  0  1  0 0.2468354 0 NA
    # 3     3  1 1   2    2  0  0  1 0.9945904 0 NA

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

    #--- order d and reset row names 
    d <- d[order(d$id, d$time), ]
    rownames(d) <- NULL 
    d$event <- ifelse(d$time == 0 & d$d == 1, 1,
                 ifelse(d$time == 1 & d$d == 1, 1, 0)) 
    #==========================================================================#


    #==== Calculate proportion dropout at each time ===========================#
    #--- Keep 1 observation per id  
    d1 <- d[which(d$d==1), c("id", "x", "time", "d", "event")]
    d1$time <- ifelse(d1$time == 2, 1, d1$time) 
    d1$time <- d1$time + 1

    #--- Expand data frame based on duration
    d2 <- d1[rep(d1$id, d1$time), ] 
    colnames(d2) <- c("id", "x", "dur", "d", "event")

    #--- Create unstructured data
    d2$t1 <- as.numeric(row.names(d2))           # Create variable from row name      
    rownames(d2) <- NULL
    d2$t2 <- floor(d2$t1)                        # Take the integer from the decimal 
    d2$t <-  round((d2$t1 - d2$t2) * 10 + 1, 0)  # Subtract the int from the decimal
    d2$t_num <- d2$t
    d2$t <- as.factor(d2$t)
    d2$r <- ifelse(d2$event == 1 & d2$dur == d2$t_num, 1, 0)  
    d2$time <- d2$t_num   

    #--- Final event occurrence dataset
    d3 <- d2[, c("id", "x", "time", "r")]
    d3 <- d3[order(d3$id, d3$time), ]
    d3$t_fac <- as.factor(d3$time)

    miss_output[rr, 1, 1, jj] <- length(which(d3$time==1 & d3$x==0 & d3$r == 1)) / length(which(d3$time==1 & d3$x==0)) 
    miss_output[rr, 1, 2, jj] <- length(which(d3$time==2 & d3$x==0 & d3$r == 1)) / length(which(d3$time==2 & d3$x==0)) 
    miss_output[rr, 2, 1, jj] <- length(which(d3$time==1 & d3$x==1 & d3$r == 1)) / length(which(d3$time==1 & d3$x==1)) 
    miss_output[rr, 2, 2, jj] <- length(which(d3$time==2 & d3$x==1 & d3$r == 1)) / length(which(d3$time==2 & d3$x==1)) 
    #==========================================================================#

    #==== Delete data in y_l ==================================================#
    #--- Keep 1 observation per subject, the wave of dropout
    m1 <- d[, c("id", "time", "d")]
    m1 <- m1[which(m1$d==1), c("id", "time")]
    m1$time <- m1$time + 1 
    
    #--- Expand data frame based on duration
    m2 <- m1[rep(m1$id, m1$time), ] 
    colnames(m2) <- c("id", "dur")
    
    #--- Create time unstructured data
    m2$t1 <- as.numeric(row.names(m2))           # Create variable from row name      
    rownames(m2) <- NULL
    m2$t2 <- floor(m2$t1)                        # Take the integer from the decimal 
    m2$t <-  round((m2$t1 - m2$t2) * 10 + 1, 0)  # Subtract the int from the decimal
    m2$t_num <- m2$t
    m2$t <- as.factor(m2$t)
    m2$r <- ifelse(m2$dur == m2$t_num, 1, 0)  
    m2$time <- m2$t_num - 1 
    m2 <- m2[, c("id", "time", "r")]
    m2 <- m2[order(m2$id, m2$time), ]

    #--- Delete data according to m2
    y_miss <- merge(y_l, m2, by = c("id", "time"))
    y_miss <- y_miss[order(y_miss$id, y_miss$time), ]
    y_miss <- y_miss[, c("index", "id", "time", "x", "x_t", "y", "r")]
    row.names(y_miss) <- NULL

    y_cc <- y_miss$id[y_miss$time == 2]

    #=== Estimation ===========================================================#

    #--- Model 1: Full Data properly specified
    m_f1 <- lme4::lmer(y ~ 1 + time + x_t + (1 + time | id), data = y_l)
    m_f1_fe <- as.numeric(fixef(m_f1))  
    m_f1_se <- sqrt(diag(vcov(m_f1)))
    m_f1_vc <- as.numeric(as.data.frame(VarCorr(m_f1))[, 5])  
    mod_output[rr, 1, , jj] <- c(m_f1_fe, m_f1_se, m_f1_vc) 

    #--- Model 2: Full Data misspecified
    m_f2 <- lme4::lmer(y ~ 1 + time + (1 + time | id), data = y_l)
    m_f2_fe <- as.numeric(fixef(m_f2))  
    m_f2_se <- round(sqrt(diag(vcov(m_f2))), 2)
    m_f2_vc <- as.numeric(as.data.frame(VarCorr(m_f2))[, 5])  
    mod_output[rr, 2, 1:2 , jj] <- m_f2_fe 
    mod_output[rr, 2, 3 , jj] <- NA 
    mod_output[rr, 2, 4:5 , jj] <- m_f2_se 
    mod_output[rr, 2, 6 , jj] <- NA 
    mod_output[rr, 2, 7:10 , jj] <- m_f2_vc 

    #--- Model 3: Complete-Case Data
    m_cc <- lme4::lmer(y ~ 1 + time + x_t + (1 + time | id), data = y_miss[which(y_miss$id %in% y_cc), ])
    m_cc_fe <- as.numeric(fixef(m_cc))  
    m_cc_se <- sqrt(diag(vcov(m_cc)))
    m_cc_vc <- as.numeric(as.data.frame(VarCorr(m_cc))[, 5]) 
    mod_output[rr, 3, , jj] <- c(m_cc_fe, m_cc_se, m_cc_vc) 

    #--- Model 4: Available Data
    m_av1 <- lme4::lmer(y ~ 1 + time + x_t + (1 + time | id), data = y_miss)
    m_av1_fe <- as.numeric(fixef(m_av1))  
    m_av1_se <- sqrt(diag(vcov(m_av1)))
    m_av1_vc <- as.numeric(as.data.frame(VarCorr(m_av1))[, 5])  
    mod_output[rr, 4, , jj] <- c(m_av1_fe, m_av1_se, m_av1_vc) 

    #--- Model 5: Available Data misspecified
    m_av2 <- lme4::lmer(y ~ 1 + time + (1 + time | id), data = y_miss)
    m_av2_fe <- as.numeric(fixef(m_av2))  
    m_av2_se <- round(sqrt(diag(vcov(m_av2))), 2)
    m_av2_vc <- as.numeric(as.data.frame(VarCorr(m_av2))[, 5])  
    mod_output[rr, 5, 1:2 , jj] <- m_av2_fe 
    mod_output[rr, 5, 3 , jj] <- NA 
    mod_output[rr, 5, 4:5 , jj] <- m_av2_se 
    mod_output[rr, 5, 6 , jj] <- NA 
    mod_output[rr, 5, 7:10 , jj] <- m_av2_vc 

  }
}

mcar1_mis_output <- miss_output
mcar1_mod_output <- mod_output

save(mcar1_mis_output, file = "2_mcar1_mis_output.rda")
save(mcar1_mod_output, file = "2_mcar1_mod_output.rda")


