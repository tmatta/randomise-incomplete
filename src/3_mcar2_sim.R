#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~ 3_mcar2_sim.R
#~~~
#~~~ Simulation design:  
#~~~   3 waves, 2 balanced groups
#~~~   Dropout is dependent on time and z 
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

#=== Simulation Arrays ========================================================#
miss_output <- array(NA, dim = c(reps, 2, length(t_labs), length(pr_drop_labs)), 
                    dimnames = list(1:reps, c("z0", "z1"), t_labs, pr_drop_labs))

mod_output <- array(NA, dim = c(reps, 2, length(par_labs), length(pr_drop_labs)), 
                    dimnames = list(1:reps, c("av", "cc"), par_labs, pr_drop_labs))

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
    y_11 <- rmvnorm(n / 4, mean = y_mu_1, sigma = y_var)
    y_10 <- rmvnorm(n / 4, mean = y_mu_1, sigma = y_var)
    y_01 <- rmvnorm(n / 4, mean = y_mu_0, sigma = y_var)
    y_00 <- rmvnorm(n / 4, mean = y_mu_0, sigma = y_var)
        
    #--- Name data 
    y_11_w <- data.frame(y_11)
    y_10_w <- data.frame(y_10)
    y_01_w <- data.frame(y_01)
    y_00_w <- data.frame(y_00)
    
    colnames(y_11_w) <- ynames <- paste0("y", seq(1:ncol(y_11_w))-1)
    colnames(y_10_w) <- ynames <- paste0("y", seq(1:ncol(y_10_w))-1)
    colnames(y_01_w) <- ynames <- paste0("y", seq(1:ncol(y_01_w))-1)
    colnames(y_00_w) <- ynames <- paste0("y", seq(1:ncol(y_00_w))-1)
    
    #--- ID variable, works for balanced groups
    id_var <- 1:n
    n_g = 4
    id_group <- split(id_var, sort(id_var %% n_g))
        
    y_11_w$id <- id_group[[1]]
    y_10_w$id <- id_group[[2]]
    y_01_w$id <- id_group[[3]]
    y_00_w$id <- id_group[[4]]
    
    #--- Covariates
    y_11_w$x <- 1;  y_11_w$z <- 1    
    y_10_w$x <- 1;  y_10_w$z <- 0    
    y_01_w$x <- 0;  y_01_w$z <- 1    
    y_00_w$x <- 0;  y_00_w$z <- 0
    
    y_w <- rbind(y_11_w, y_10_w, y_01_w, y_00_w)
    rownames(y_w) <- NULL
    
    #--- Reshape to long & reorder 
    y_l <- reshape(y_w, varying = ynames, v.names = "y", timevar = "time", 
                   times = (1 : j) - 1, direction = "long") 

    y_l <- y_l[order(y_l$id, y_l$time), ]
    rownames(y_l) <- NULL

    y_l$x_t <- y_l$x * y_l$time 
    y_l$index <- rownames(y_l)

    #--- Final full data data.frame
    y_l <- y_l[, c("index", "id", "time", "x", "x_t", "z", "y")]



    #=== Generate Dropout =====================================================#
    a0 <- log(0.10 / (1 - 0.10)) # Log-odds of event at time 1
    a1 <- log(0.15 / (1 - 0.15)) # Log-odds of event at time 2
    a2 <- log(0.99 / (1 - 0.99)) # probability of event at time 3 (1)
    a3 <- log(pr_drop[jj] / (1 - pr_drop[jj]))   # Log-odds of dropping out for x = 1
    
    #--- Create a factor wave variable
    d <- y_l
    d$t <- as.factor(d$time)
    d <- cbind.data.frame(d[ , c("index", "id", "x", "x_t", "z")],  model.matrix(~ -1  + time + t, data = d))
    d$t <- NULL
    # > head(d)
    #   index id x x_t time t0 t1 t2
    # 1     1  1 1   0    0  1  0  0
    # 2     2  1 1   1    1  0  1  0
    # 3     3  1 1   2    2  0  0  1

    #--- Hazard is conditional on the event not occurring at a prior time.
    #--- haz(t_j) = Pr(T_i = j | T_i >= j)
    d$haz <- 1 / (1 + exp(-(a0*d$t0 + a1*d$t1 + a2*d$t2 + a3*d$z)))
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
    d1 <- d[which(d$d==1), c("id", "z", "time", "d", "event")]
    d1$time <- ifelse(d1$time == 2, 1, d1$time) 
    d1$time <- d1$time + 1
    #--- Expand data frame based on duration
    d2 <- d1[rep(d1$id, d1$time), ] 
    colnames(d2) <- c("id", "z", "dur", "d", "event")
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
    d3 <- d2[, c("id", "z", "time", "r")]
    d3 <- d3[order(d3$id, d3$time), ]
    d3$t_fac <- as.factor(d3$time)

    miss_output[rr, 1, 1, jj] <- length(which(d3$time==1 & d3$z==0 & d3$r == 1)) / length(which(d3$time==1 & d3$z==0)) 
    miss_output[rr, 1, 2, jj] <- length(which(d3$time==2 & d3$z==0 & d3$r == 1)) / length(which(d3$time==2 & d3$z==0)) 
    miss_output[rr, 2, 1, jj] <- length(which(d3$time==1 & d3$z==1 & d3$r == 1)) / length(which(d3$time==1 & d3$z==1)) 
    miss_output[rr, 2, 2, jj] <- length(which(d3$time==2 & d3$z==1 & d3$r == 1)) / length(which(d3$time==2 & d3$z==1)) 
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
    y_miss <- y_miss[, c("index", "id", "time", "x", "x_t", "y", "z", "r")]
    row.names(y_miss) <- NULL

    y_cc <- y_miss$id[which(y_miss$time == 2)]
    #==========================================================================#

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

    mod_output[rr, 1, , jj] <- c(m_av_fe, m_av_se, m_av_vc) 
    mod_output[rr, 2, , jj] <- c(m_cc_fe, m_cc_se, m_cc_vc) 


    y_c <- y_l
    y_c$r <- ifelse(y_l$index %in% y_miss$index, 1, 0)
    df_list[[rr]] <- y_c
  }
}



mcar2_mis_output <- miss_output
mcar2_mod_output <- mod_output
mcar2_data <- df_list

save(mcar2_mis_output, file = "3_mcar2_mis_output.rda")
save(mcar2_mod_output, file = "3_mcar2_mod_output.rda")
save(mcar2_data, file = "3_mcar2_data.rda")









test <- 
density(df_list[[1]]$y[which(df_list[[1]]$r == 0 & 
                                     df_list[[1]]$x == 1 & 
                                     df_list[[1]]$z == 1 &
                                     df_list[[1]]$time == 1)]) / 
density(df_list[[1]]$y[which(df_list[[1]]$x == 1 & 
                                     df_list[[1]]$z == 1 &
                                     df_list[[1]]$time == 1)])

plot(NULL, xlim = c(15, 30), ylim = c(0, .5))
for(i in 1:reps){
  lines(density(df_list[[i]]$y[which(df_list[[i]]$x == 1 & 
                                     df_list[[i]]$z == 1 &
                                     df_list[[i]]$time == 1)]), lwd = .5, lty = 1) 
  lines(density(df_list[[i]]$y[which(df_list[[i]]$r == 0 & 
                                     df_list[[i]]$x == 1 & 
                                     df_list[[i]]$z == 1 &
                                     df_list[[i]]$time == 1)]), lwd = .5, lty = 1, col = "red") 

}

hist(df_list[[1]]$y[which(df_list[[1]]$r == 1 & 
                          df_list[[1]]$z == 1 &
                          df_list[[1]]$time == 1)], breaks=20)

hist(df_list[[1]]$y[which(df_list[[1]]$r == 0 & 
                          df_list[[1]]$z == 1 &
                          df_list[[1]]$time == 1)], breaks=20)

hist(df_list[[1]]$y[which(df_list[[1]]$z == 1 &
                          df_list[[1]]$time == 1)], breaks=20)


length(df_list[[i]]$x[which(df_list[[i]]$r == 0 & 
                            df_list[[i]]$x == 0 &
                            df_list[[i]]$z == 1 &
                            df_list[[i]]$time == 1)])

length(df_list[[i]]$x[which(df_list[[i]]$x == 0 &
                            df_list[[i]]$z == 1 &
                            df_list[[i]]$time == 1)])

length(df_list[[i]]$x[which(df_list[[i]]$r == 0 & 
                            df_list[[i]]$x == 1 &
                            df_list[[i]]$z == 1 &
                            df_list[[i]]$time == 1)])

length(df_list[[i]]$x[which(df_list[[i]]$x == 1 &
                            df_list[[i]]$z == 1 &
                            df_list[[i]]$time == 1)])


lines(density(y_complete$y[which(y_complete$r == 0 & 
                                 y_complete$x == 1 & 
                                 y_complete$z == 1 &
                                 y_complete$time == 1)])) 

lines(density(y_complete$y[which(y_complete$r == 1 & 
                                 y_complete$x == 0 & 
                                 y_complete$time == 1)]), lty =2) 
lines(density(y_complete$y[which(y_complete$r == 0 & 
                                 y_complete$x == 0 & 
                                 y_complete$time == 1)])) 

points(x = y_complete$time[which(y_complete$r == 0 & y_complete$x == 1)] + .05, 
       y = y_complete$y[which(y_complete$r == 0 & y_complete$x == 1)],
       col = "blue")


