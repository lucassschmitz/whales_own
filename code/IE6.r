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

rm(list = ls())


script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname(script_dir))
print(getwd())


fig_folder <- "Figures/meeting_250527"

df <- read_csv("Data/masters_voyages_merged.csv")

# clean data and show some summary stats
    df <- df %>% select(product, yearout, duration, vesselid, mastercode, tonn_avg, yearin, age, lay)
    df$duration[df$duration == 0] <- 1 # if yearin = yearout 
    df$re_prod <- df$product/df$duration # real production 
    df <- df %>% mutate(vesselcode = as.numeric(factor(vesselid)))

    df <- df[order(df$mastercode), ]
    df$decade <- (df$yearout %/% 10) * 10
        
    df <- df %>%
      arrange(mastercode, yearout) %>%
      group_by(mastercode) %>%
      mutate(n_voyage = row_number(),
       N = n(), 
       master_last = as.integer(yearout == max(yearout, na.rm = TRUE)),
       first_year = as.integer(min(yearout, na.rm = TRUE)), 
       n_voyage2 = yearout -first_year ) %>%
      ungroup()

    df <- df %>% filter(!is.na(product))


#################################################
## Part 1 
##################################################
  
# lag for production. 
df <- df %>%   arrange(mastercode, n_voyage) %>%     
  group_by(mastercode) %>%  mutate(l_re_prod = lag(re_prod)) %>%  ungroup()

### regressions 
    m1 <- feols(re_prod ~ 1 + n_voyage , data = df)
    m2 <- feols(re_prod ~ 1 + n_voyage  |  mastercode  , data = df)
    m3 <- feols(re_prod ~ 1 + n_voyage  |  vesselcode  , data = df)
    m4 <- feols(re_prod ~ 1 + n_voyage  |  mastercode + vesselcode + yearout , data = df)

   
    m5 <- feols(re_prod ~ 1 + n_voyage + l_re_prod , data = df)
    m6 <- feols(re_prod ~ 1 + n_voyage + l_re_prod |  mastercode  , data = df)
    m7 <- feols(re_prod ~ 1 + n_voyage + l_re_prod |  vesselcode  , data = df)
    m8 <- feols(re_prod ~ 1 + n_voyage + l_re_prod |  mastercode + vesselcode + yearout , data = df)
  

# pull out the lists of fixed effects
fe4 <- fixef(m4)
fe8 <- fixef(m8)

# add them to your df under distinct names
df <- df %>% 
  mutate(
    master_fe_m4 = fe4$mastercode[ match(mastercode, names(fe4$mastercode)) ],
    vessel_fe_m4 = fe4$vesselcode[ match(vesselcode, names(fe4$vesselcode)) ],
    master_fe_m8 = fe8$mastercode[ match(mastercode, names(fe8$mastercode)) ],
    vessel_fe_m8 = fe8$vesselcode[ match(vesselcode, names(fe8$vesselcode)) ]
  )

# now you get two (potentially) different correlations
cor(df$master_fe_m4, df$vessel_fe_m4, use = "complete.obs")
cor(df$master_fe_m8, df$vessel_fe_m8, use = "complete.obs")

    #export 
    esttex(  m4,  m5, m6, m7, m8, 
    keep   = c("%n_voyage", "%l_re_prod" ),
    dict = c(re_prod = "Production/Duration", n_voyage  = "Experience"),
    file   = "Writeup/Tables/IE6_tab_n_voyage2.tex",
    replace = TRUE,          
    label  = "tab:n_voyage", 
    title  = "Effect of Experience on Production",
    digits = 3,
    notes = "" )



 
#################################################
############# ANOVA test rejects null 
##################################################
# restricted: drop vesselcode
m9_res <- lm(re_prod ~ n_voyage + tonn_avg + l_re_prod + factor(mastercode) + factor(yearout),  data = df)
m9_unres <- lm(re_prod ~ n_voyage + tonn_avg + l_re_prod + factor(mastercode) + factor(yearout) + factor(vesselcode), data = df )

adj_r2_res <- summary(m9_res)$adj.r.squared
adj_r2_unres <- summary(m9_unres)$adj.r.squared
anova(m9_res, m9_unres) ## 1.7863 (pvalue 0.000) without l_re_prod,  1.779 (pvalue 0.000) with l_re_prod,

#m10_res <- lm(re_prod ~ n_voyage + tonn_avg  + factor(yearout) + factor(vesselcode), data = df )
#anova(m9_unres, m10_res) ## 1.26 (pvalue )

m11 <- lm(scale(vessel_fe_m8) ~ scale(tonn_avg), data = df)
summary(m11) 

ggplot(df, aes(x = vessel_fe_m8, y = tonn_avg)) +
  geom_bin2d(bins = 50) +
  labs(x = "Vessel FE", y = "Weight (tons) ", title = "") +
  scale_x_continuous(limits = c(-50000, 50000)) + 
  theme_minimal(base_size = 20)
ggsave(filename = file.path(fig_folder,  "IE6_weight_shipFE.png"), dpi = 560, bg = "white") 
 


 
#################################################
############# Ship evolution 
##################################################
 df <- df %>%  group_by(vesselcode) %>%
  mutate(built = min(yearout, na.rm = TRUE)) %>%  ungroup()



ggplot(df, aes(x = built, y = vessel_fe_m8 )) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    x     = "Year Built",
    y     = "Vessel FE",
    title = "Evolution of Vessel Fixed Effects by Build Year", 
    caption = "Year built is the defined as the year of departure of the first voyage."
  ) +
  theme_minimal(base_size = 25)
  
ggsave(filename = file.path(fig_folder,  "IE6_shipFE_built.png"), dpi = 560, bg = "white") 


ggplot(df, aes(x = built, y = tonn_avg)) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    x     = "Year Built",
    y     = "Weight (tons)",
    title = "Evolution of Weight by Build Year"
  ) +
  theme_minimal(base_size = 25)
ggsave(filename = file.path(fig_folder,  "IE6_weight_built.png"), dpi = 560, bg = "white") 

############ Correlation between lay and tonn_avg

 

 # Fit the logit
logit_mod <- glm(master_last ~ re_prod,
                 data   = df,
                 family = binomial(link = "logit"))

# View coefficients, standard errors, and fit statistics
summary(logit_mod)


 

