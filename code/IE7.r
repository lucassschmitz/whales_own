library(readr)
library(MASS)
library(dplyr)
library(fixest)
library(lmtest)
library(ggplot2)
library(np)
library(broom)
library(skimr)
library(sandwich)
library(lmtest)

#library(plm)
rm(list = ls())


script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname(script_dir))
print(getwd())


fig_folder <- "Figures/meeting_250603"

df <- read_csv("Data/masters_voyages_merged.csv")

# clean data and show some summary stats
    df <- df %>% select(product, yearout, duration, vesselid, mastercode, tonn_avg, yearin, age, lay, re_prod, n_voyage, l_re_prod, change_boat, stay_boat)
    

    df <- df %>% mutate(vesselcode = as.numeric(factor(vesselid)))

    df <- df[order(df$mastercode), ]
         
    df <- df %>%
      arrange(mastercode, yearout) %>%
      group_by(mastercode) %>%
      mutate( N = n(), 
       master_last = as.integer(yearout == max(yearout, na.rm = TRUE)), # dummy for master's last voyage. 
       first_year = as.integer(min(yearout, na.rm = TRUE)), # first year of the master. 
       n_voyage2 = yearout -first_year, # time since first voyage. 
       exp_years = lag(cumsum(duration), default = 0) #time at sea before the current voyage 
       ) %>% ungroup()

    df <- df %>% filter(!is.na(product))

#################################################
## Part 1: different measure of experience does not matter much 
##################################################
    m1 <- feols(re_prod ~ 1 + n_voyage , data = df)
    
    # try different experience measures 
    m2 <- feols(re_prod ~ 1 + n_voyage  + master_last   |  mastercode + vesselcode + yearout , data = df)
    m2b <- feols(re_prod ~ 1 + n_voyage  + master_last + l_re_prod  |  mastercode + vesselcode + yearout , data = df)

    m3 <- feols(re_prod ~ 1 + exp_years  + master_last  |  mastercode + vesselcode + yearout , data = df)


  # try includign the lagg with both of the experience measures. 
    m4 <- feols(re_prod ~ 1 + n_voyage2  + master_last + l_re_prod  |  mastercode + vesselcode + yearout , data = df)
    m5 <- feols(re_prod ~ 1 + exp_years  + master_last + l_re_prod |  mastercode + vesselcode + yearout , data = df)


 esttex(   m1,  m2, m2b, m3, m5,  
    dict = c(re_prod = "Production/Duration", n_voyage  = "Experience", exp_years = "Exp. in years"),
    file   = "Writeup/Tables/IE7_experience_on_production.tex",
    replace = TRUE,          
    label  = "tab:n_voyage", 
    title  = "Effect of Experience on Production",
    digits = 3,
    notes = "Exp. in years refers to the cumulative number of years of previous voyages." )


#################################################
## Part 2a:  correlation for switchers and non-switchers [not finished}]
##################################################

View(df)
 
    df <- df %>%
      arrange(mastercode, n_voyage) %>%
      group_by(mastercode) %>%
      mutate( change_boat2 = dplyr::lead(change_boat))  %>% ungroup()

df_changers <- df %>% filter(change_boat == 1 | change_boat2 == 1) # keep only the ones that changed boat at least once.


m4 <- feols(re_prod ~ 1 + n_voyage  |  mastercode + vesselcode + yearout , data = df_changers)
m4b <- feols(re_prod ~ 1 + n_voyage  |  mastercode + vesselcode + yearout , data = df)

# pull out the lists of fixed effects
fe4 <- fixef(m4)
fe4b <- fixef(m4b)

# add them to your df under distinct names
df <- df %>% 
  mutate(
    master_fe_m4 = fe4$mastercode[ match(mastercode, names(fe4$mastercode)) ],
    vessel_fe_m4 = fe4$vesselcode[ match(vesselcode, names(fe4$vesselcode)) ],
    master_fe_m4b = fe4b$mastercode[ match(mastercode, names(fe4b$mastercode)) ],
    vessel_fe_m4b= fe4b$vesselcode[ match(vesselcode, names(fe4b$vesselcode)) ]
  )

# now you get two (potentially) different correlations
cor(df$master_fe_m4, df$vessel_fe_m4, use = "complete.obs")
cor(df$master_fe_m4b, df$vessel_fe_m4b, use = "complete.obs")
 

#################################################
## Part 2b:  do the correction 
##################################################
m4 <- feols(re_prod ~ 1 | mastercode + vesselcode, data = df)
sigma2e <- m4$sigma2


df_aux <- df %>% select(product, yearout, duration, vesselcode, mastercode, tonn_avg, yearin, re_prod, n_voyage, l_re_prod, master_last)


# Get unique mastercodes/vesselcodes and sort them
unique_mastercodes <- sort(unique(df_aux$mastercode))
vesselcodes <- sort(unique(df_aux$vesselcode))


# Create the mastercode fixed effects matrix
D <- model.matrix(~ 0 + factor(mastercode, levels = unique_mastercodes), data = df_aux)
F <- model.matrix(~ 0 + factor(vesselcode, levels = vesselcodes), data = df_aux)
Z <- matrix(1, nrow = nrow(D), ncol = 1) ## we later could include covariates. 


I <- diag(nrow(D))  # Identity matrix 
ones <- matrix(1, nrow = nrow(D), ncol = 1)  # Vector of ones
A <- I - ones %*% solve(t(ones) %*% ones) %*% t(ones)  # Projection matrix for the constant term



M_d <- I - D %*% solve(t(D) %*% D) %*% t(D)  # Projection matrix for mastercode fixed effects
M_f <- I - F %*% solve(t(F) %*% F) %*% t(F)  # Projection matrix for vesselcode fixed effects
M_z <- I - Z %*% solve(t(Z) %*% Z) %*% t(Z)  # Projection matrix for the constant term


F_tilde <- M_z %*% F   
M_f_tilde <- F_tilde %*% solve(t(F_tilde) %*% F_tilde) %*% t(F_tilde) 
D_tilde <- M_z %*% D
M_d_tilde <- D_tilde %*% solve(t(D_tilde) %*% D_tilde) %*% t(D_tilde)


Q_f_tilde <- M_z %*% M_f_tilde  %*% M_z
Q_d_tilde <- M_z %*% M_d_tilde  %*% M_z

theta_hat <- solve(t(D) %*% Q_f_tilde %*% D) %*% t(D) %*% Q_f_tilde %*% df_aux$re_prod # equation 14
phi_hat <- solve(t(F)%*% Q_d_tilde %*% F) %*% t(F) %*% Q_d_tilde %*% df_aux$re_prod # equation 15


# constructing estimator 
part1 <- t(theta_hat) %*% t(D) %*% A %*% F %*% phi_hat
part2a <- t(D) %*% A %*% F %*% solve(t(F) %*% Q_d_tilde %*% F)  
part2b <- t(F) %*% M_z %*% D %*% solve(t(D) %*% M_z %*% D)

consistent_estimator <- (part1 + sigma2e * sum(diag(part2a %*% part2b)) )/nrow(D) 
 




#################################################
## Part 3:  why is it that there is a negative auto-correlation in production?
##################################################
df1 <- df %>% filter(N == 1)
df2 <- df %>% filter(N == 2)
df3 <- df %>% filter(N == 3)
df4 <- df %>% filter(N == 4)
df5 <- df %>% filter(N >= 5)

df_last <- df %>% filter(master_last == 1) 
df_previous <- df %>% filter(master_last == 0)

    m10 <- feols(re_prod ~ 1 + l_re_prod , data = df)
    m11 <- feols(re_prod ~ 1 + l_re_prod  |  mastercode  , data = df)
    m12 <- feols(re_prod ~ 1 + l_re_prod  |   vesselcode   , data = df)
    m13 <- feols(re_prod ~ 1 + l_re_prod   |  mastercode + vesselcode + yearout , data = df)

esttex( m10, m11,  m12, m13, 
    dict = c(re_prod = "Whole sample ", n_voyage  = "Experience"),
    file   = "Writeup/Tables/IE7_aux2.tex",
    replace = TRUE,          
    label  = "tab:aux2", 
    title  = "Whole sample",
    digits = 3,
    notes = "" )



  m10a <- feols(re_prod ~ 1 + l_re_prod , data = df2)
    m10a2 <- feols(re_prod ~ 1 + l_re_prod , data = df3)
    m10a3 <- feols(re_prod ~ 1 + l_re_prod , data = df4)
    m10a4 <- feols(re_prod ~ 1 + l_re_prod , data = df5)
    
esttex( m10a,m10a2, m10a3,  m10a4, 
    dict = c(re_prod = "Whole sample ", n_voyage  = "Experience"),
    file   = "Writeup/Tables/IE7_aux3.tex",
    replace = TRUE,          
    label  = "tab:aux3", 
    title  = "Autocorrelation for different samples. ",
    digits = 3,
    notes = "Column i has the captains that did i+1 voyages" )



    m11a <- feols(re_prod ~ 1 + l_re_prod  |  mastercode  , data = df3)
    m12a <- feols(re_prod ~ 1 + l_re_prod  |   vesselcode   , data = df3)
    m13a <- feols(re_prod ~ 1 + l_re_prod   |  mastercode + vesselcode + yearout , data = df3)

    m14a <- feols(re_prod ~ 1 + l_re_prod  |  mastercode  , data = df4)
    m15a <- feols(re_prod ~ 1 + l_re_prod  |   vesselcode   , data = df4)
    m16a <- feols(re_prod ~ 1 + l_re_prod   |  mastercode + vesselcode + yearout , data = df4)

    m17a <- feols(re_prod ~ 1 + l_re_prod  |  mastercode  , data = df5)
    m18a <- feols(re_prod ~ 1 + l_re_prod  |   vesselcode   , data = df5)
    m19a <- feols(re_prod ~ 1 + l_re_prod   |  mastercode + vesselcode + yearout , data = df5)
esttex(   m11a,  m12a, m13a, m14a, m15a, m16a, m17a, m18a, m19a,
    dict = c(re_prod = "Whole sample ", n_voyage  = "Experience"),
    file   = "Writeup/Tables/IE7_aux4.tex",
    replace = TRUE,          
    label  = "tab:n_voyage", 
    title  = "Effect of Experience on Production",
    digits = 3,
    notes = "The first three columns use the sample with 3 voyages, the next 3 cols use the sample with 4 voyages and so on. " )

   m10_last <- feols(re_prod ~ 1 + l_re_prod , data = df_last)
    m11_last <- feols(re_prod ~ 1 + l_re_prod  |  mastercode  , data = df_last)
    m12_last <- feols(re_prod ~ 1 + l_re_prod  |   vesselcode   , data = df_last)
    m13_last <- feols(re_prod ~ 1 + l_re_prod   |  mastercode + vesselcode + yearout , data = df_last)

m10_previous <- feols(re_prod ~ 1 + l_re_prod , data = df_previous)
    m11_previous <- feols(re_prod ~ 1 + l_re_prod  |  mastercode  , data = df_previous)
    m12_previous <- feols(re_prod ~ 1 + l_re_prod  |   vesselcode   , data = df_previous)
    m13_previous <- feols(re_prod ~ 1 + l_re_prod   |  mastercode + vesselcode + yearout , data = df_previous)
models <- list(m10_last, m11_last, m12_last, m13_last,
               m10_previous, m11_previous, m12_previous, m13_previous)

esttex(models, 
    dict = c(re_prod = "Whole sample ", n_voyage  = "Experience"),
    file   = "Writeup/Tables/IE7_aux5.tex",
    replace = TRUE,          
    label  = "tab:n_voyage", 
    title  = "Last voyage or not last voyage",
    digits = 3,
    notes = "" )

