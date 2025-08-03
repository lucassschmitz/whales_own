library(readr)
library(MASS)
library(dplyr)
library(fixest)
library(ggplot2)
library(np)
library(broom)
library(skimr)
library(sandwich)
library(lmtest)

rm(list = ls())

setwd("C:/Users/lucas/Dropbox/whales_own")

fig_folder <- "Figures/temp2"

df <- read_csv("Data/masters_voyages_merged.csv")

# clean data and show some summary stats
    df <- df %>% select(product, yearout, duration, vesselid, mastercode, tonn_avg, yearin, age)
    df$duration[df$duration == 0] <- 1 # if yearin = yearout 
    df$re_prod <- df$product/df$duration # real production 
    df <- df %>% mutate(vesselcode = as.numeric(factor(vesselid)))

    df <- df[order(df$mastercode), ]
    df$decade <- (df$yearout %/% 10) * 10
        
    df <- df %>%
      arrange(mastercode, yearout) %>%
      group_by(mastercode) %>%
      mutate(n_voyage = row_number(), N = n() ) %>%
      ungroup()

    skim(df$re_prod)

    voyages_per_captain <- df %>%  count(mastercode, name = "voyage_count")
    summary(voyages_per_captain$voyage_count)
    ggplot(voyages_per_captain , aes(x = voyage_count)) +
      geom_histogram(bins = 30, color = "black", fill = "steelblue") +
      labs(title = "Voyages per Captain", x = "Number of Voyages") +
      theme_minimal()


    df <- df %>% filter(!is.na(product))
    #df <- df %>% filter(yearout >= 1790 & yearout <= 1900)

    cat("Number of different mastercode values:", length(unique(df$mastercode))) 
    length(df$mastercode)
    #2247  captains and 5424 obs

    #voyages per captain after filtering. 
    voyages_per_captain <- df %>%  count(mastercode, name = "voyage_count")
    summary(voyages_per_captain$voyage_count)
    ggplot(voyages_per_captain, aes(x = voyage_count)) +
      geom_histogram(bins = 30, color = "black", fill = "steelblue") +
      labs(title = "Voyages per Captain", x = "Number of Voyages") +
      theme_minimal()


##################################################
##### Evolution of real production over time 
##################################################


    df %>%  filter(yearout > 1790 & yearout < 1901) %>%
    ggplot(aes(x = yearout, y = re_prod)) + 
     geom_point() + 
      geom_smooth(method = "loess", se = TRUE, color = "red") +
      labs(x = "Year Out", y = "Production/Duration", 
          title = "Productivity vs Year (1791-1901)")
      ggsave(filename = file.path(fig_folder,  "reprod_year.png"), dpi = 360) 

  df %>%  filter(yearout > 1790 & yearout < 1901 ) %>%
    ggplot(aes(x = yearout, y = product)) +
          geom_point() +
      geom_smooth(method = "loess", se = TRUE, color = "red") +
      labs(x = "Year Out", y = "Production", title = "Productivity vs Year Out (1791- 1901)")
    ggsave(filename = file.path(fig_folder,  "prod_year.png"), dpi = 360)


  # plot 5-year fixed effects. 
 

df <- df %>%   mutate(
    lustrum = ((yearout - 1800) %/% 5) * 5 + 1800,
    lustrum = factor(lustrum, levels = seq(1800, 1895, by = 5))
  )


# Run regression of y on x with robust standard errors
m_yx <- lm(re_prod ~ lustrum, data = df)
 coeftest(m_yx, vcov = vcovHC(m_yx, type = "HC1"))
ct <- coeftest(m_yx, vcov = vcovHC(m_yx, type = "HC1"))


coef_rows <- grep("^lustrum", rownames(ct))
est       <- ct[coef_rows, 1]
se        <- ct[coef_rows, 2]
 
base_level        <- levels(df$lustrum)[1]  
fitted_levels     <- sub("^lustrum", "", rownames(ct)[coef_rows])
lustrum_plot_labs <- c(base_level, fitted_levels)

# 4) make your plotting data-frame
plot_df <- data.frame(
  lustrum  = lustrum_plot_labs,
  estimate = c(0, est),  ci_low = c(0, est - 1.96 * se), ci_high  = c(0, est + 1.96 * se)
)

plot_df$lustrum <- factor(plot_df$lustrum, levels = lustrum_plot_labs)

ggplot(plot_df, aes(x = lustrum, y = estimate, group = 1)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "5-year period", y = "Production FE",
    title = "Lustrum Effects on Production\n(normalized to first bin = 0)") +
  theme_minimal()

ggsave(filename = file.path(fig_folder,  "reprod_yearFE.png"), dpi = 360, bg = "white")




# Run regression of y on x with robust standard errors
m_yx <- lm(product ~ lustrum, data = df)
ct <- coeftest(m_yx, vcov = vcovHC(m_yx, type = "HC1"))


coef_rows <- grep("^lustrum", rownames(ct))
est       <- ct[coef_rows, 1]
se        <- ct[coef_rows, 2]
 
base_level        <- levels(df$lustrum)[1]  
fitted_levels     <- sub("^lustrum", "", rownames(ct)[coef_rows])
lustrum_plot_labs <- c(base_level, fitted_levels)

# 4) make your plotting data-frame
plot_df <- data.frame(
  lustrum  = lustrum_plot_labs,
  estimate = c(0, est),  ci_low = c(0, est - 1.96 * se), ci_high  = c(0, est + 1.96 * se)
)

plot_df$lustrum <- factor(plot_df$lustrum, levels = lustrum_plot_labs)

ggplot(plot_df, aes(x = lustrum, y = estimate, group = 1)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "5-year period", y = "Production FE",
    title = "Lustrum Effects on Production\n(normalized to first bin = 0)") +
  theme_minimal()

ggsave(filename = file.path(fig_folder,  "prod_yearFE.png"), dpi = 360, bg = "white")


##################################################
## Figure 3 
##################################################
  
  # plot original figure 3 weight
  df %>%   filter(n_voyage -1 < 5) %>%                 
      ggplot(aes(x = n_voyage - 1, y = tonn_avg)) +
      geom_smooth(method = "loess",  linewidth = 1.1) +
      stat_summary(fun  = mean, geom = "line", linewidth = 1.1) +
      labs(x = "Experience", y = "Weight(Tons) ")
  ggsave(filename = file.path(fig_folder,  "fig3_weight.png"), dpi = 360)

  # original figure 3 but real production 
  df %>%   filter(n_voyage -1 < 5) %>%                
    ggplot(aes(x = n_voyage - 1, y = re_prod)) +
    geom_smooth(method = "loess",  linewidth = 1.1) +
    stat_summary(fun  = mean, geom = "line", linewidth = 1.1) +
    labs(x = "Experience", y = "Production/Duration")
    ggsave(filename = file.path(fig_folder,  "fig3_reprod.png"), dpi = 360)


### regressions 
    m1 <- feols(re_prod ~ 1 + n_voyage, data = df)
    m2 <- feols(re_prod ~ 1 + n_voyage  |  mastercode  , data = df)
    m3 <- feols(re_prod ~ 1 + n_voyage  |  vesselcode  , data = df)
    m4 <- feols(re_prod ~ 1 + n_voyage  |  mastercode + vesselcode + yearout , data = df)

    # reduce the number of fe by time trend or proxying ship fe with weight 
    m5 <- feols(re_prod ~ 1 + n_voyage + yearout |  mastercode + vesselcode , data = df)
    m6 <- feols(re_prod ~ 1 + n_voyage + tonn_avg |  mastercode + yearout , data = df)
 
    # change the sample
    df2 <- df %>%  add_count(mastercode, name = "N2")
    df2 <- df2 %>% filter(N2 > 2)
    m7 <- feols(re_prod ~ 1 + n_voyage |  mastercode + vesselcode + yearout , data = df2)
    m8 <- feols(re_prod ~ 1 + n_voyage + yearout  |  mastercode + vesselcode , data = df2)
    m9 <- feols(re_prod ~ 1 + n_voyage + tonn_avg  |  mastercode + yearout , data = df2)


    #export 
    esttex( m1, m2, m3,  m4,  m7,
    keep   = "%n_voyage",
    dict = c(re_prod = "Production/Duration", n_voyage  = "Experience"),
    file   = "Writeup/Tables/tab_n_voyage.tex",
    replace = TRUE,          
    label  = "tab:n_voyage", 
    title  = "Effect of Experience on Production",
    digits = 3,
    notes = "The last column uses only voyages with captains with more than 2 voyages, the other columns use the full sample. " )


    # heatmap of captain and ship fixed effects 

    fe <- fixef(m4)
    corr <- cor(fe$mastercode[as.character(df$mastercode)],
              fe$vesselcode[as.character(df$vesselcode)],
              use = "complete.obs")
    print(corr)
    uv <- unique(df[c("mastercode","vesselcode")])
    uv$alpha <- fe$mastercode[as.character(uv$mastercode)]
    uv$gamma <- fe$vesselcode[as.character(uv$vesselcode)]
    qs_alpha <- quantile(uv$alpha, probs = c(0.01, 0.99), na.rm = TRUE)
    qs_gamma <- quantile(uv$gamma, probs = c(0.01, 0.99), na.rm = TRUE)
    uv_trim <- uv %>%
      filter(
        alpha > qs_alpha[1], alpha < qs_alpha[2],
        gamma > qs_gamma[1], gamma < qs_gamma[2]
      )

    ggplot(uv_trim, aes(x = alpha, y = gamma)) +
      stat_bin2d(bins = 50) +                        # counts in 50×50 grid
      scale_fill_viridis_c(option = "D") +           # continuous colour scale
      labs(x  = "Captain FE", y = "Vessel FE", fill  = "Pairs",
        title = "Heatmap of Captain vs. Vessel Fixed Effects") +
      theme_minimal()
    ggsave(filename = file.path(fig_folder,  "heat_vessel&master_fe.png"), dpi = 360, bg = "white")

    # Plottign experience FE 
   df <- df %>% mutate(
      n_voy_f = factor( 
        ifelse(n_voyage >= 5, "5+", as.character(n_voyage)),
        levels = c("1","2","3","4","5+")
      ))

    m2 <- feols(re_prod ~ i(n_voy_f, ref = "1") | mastercode + vesselcode + yearout,data = df)

    tib <- tidy(m2, conf.int = TRUE) %>% 
      filter(grepl("^n_voy_f::", term)) %>%
      mutate(n_voyage = sub("n_voy_f::", "", term),
        estimate = estimate, ci_low   = conf.low, ci_high  = conf.high)

    tib <- bind_rows(
      tibble(n_voyage = "1", estimate = 0, ci_low   = 0, ci_high  = 0),
      tib
    ) %>%
      mutate(n_voyage = factor(n_voyage, levels = c("1","2","3","4","5+")))

    ggplot(tib, aes(x = n_voyage, y = estimate, group = 1)) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
      geom_line() +
      geom_point(size = 2) +
      labs( x = "Experience", y = "Experience FE", title = "Experience Effects") +
      theme_minimal()
        ggsave(filename = file.path(fig_folder,  "experienceFE.png"), dpi = 360, bg = "white")


#################################################
############## Figure 5 ###################### 
##################################################
# Figures 
  ggplot(df, aes(x = age, y = tonn_avg)) +
    geom_smooth(method = "loess",  linewidth = 1.1) +
    labs(x = "Age", y = "Weight(Tons)")
  ggsave(filename = file.path(fig_folder,  "fig5_weight.png"), dpi = 360)
  
    # original figure 5 but real production 
    ggplot(df, aes(x = age, y = re_prod)) +
      geom_smooth(method = "loess",  linewidth = 1.1) +
      labs(x = "Age", y = "Production/Duration")
      ggsave(filename = file.path(fig_folder,  "fig5_reprod.png"), dpi = 360)
  
# regressions 
m10 <- feols(re_prod ~ 1 + age + n_voyage |  mastercode + vesselcode + yearout , data = df)
m11 <-feols(re_prod ~ 1 + age |  mastercode + vesselcode + yearout , data = df)
m12 <- feols(re_prod ~ 1 + age |  mastercode + vesselcode + yearout , data = df2)
m13 <- feols(re_prod ~ 1 + age + yearout  |  mastercode + vesselcode , data = df2)
m14 <- feols(re_prod ~ 1 + age + tonn_avg  |  mastercode + yearout , data = df2)

 #export 
    esttex( m10, m11, m12,  m13,  m14,
    keep   = "%age",
    dict = c(re_prod = "Production/Duration", n_voyage = "Experience", age  = "Age"),
    file   = "Writeup/Tables/tab_age.tex",
    replace = TRUE,          
    label  = "tab:age", 
    title  = "Effect of Experience on Production",
    digits = 3,
    notes = "The last column uses only voyages with captains with more than 2 voyages, the other columns use the full sample. " )


 corr <- cor(df$age, df$n_voyage, use = "complete.obs")

##################### Trip value added 

df <- df %>% select(-c(yearin, age, vesselid, decade, lustrum, n_voy_f))

  m21 <- feols(re_prod ~ 1 |  mastercode + vesselcode + yearout , data = df)
  df <- df %>% 
  mutate(
    master_fe  = fe$mastercode[ match(mastercode, names(fe$mastercode)) ],
    vessel_fe  = fe$vesselcode[ match(vesselcode, names(fe$vesselcode)) ],
    year_fe    = fe$yearout[   match(yearout,   as.integer(names(fe$yearout))) ]
  )
  
  df$VA <- df$re_prod - df$vessel_fe - df$year_fe 

df <- df %>%   arrange(mastercode, n_voyage) %>%     
  group_by(mastercode) %>%              
  mutate(l_VA = lag(VA), 
          del_ship = vessel_fe - lag(vessel_fe),
          avg_l_VA = lag(cummean(VA)),           
          prev_vessel = lag(vesselcode), 
          change = as.integer(!is.na(prev_vessel) & vesselcode != prev_vessel), 
          not_change = as.integer(!is.na(prev_vessel) & vesselcode == prev_vessel)
          ) %>%
  ungroup()
View(df) 

skim(df$change) 
skim(df$not_change)


 m22 <- feols(scale(VA) ~ 1 + scale(l_VA) , data = df)

m23 <-feols(scale(vessel_fe) ~ 1 + scale(l_VA) , data = df)
m24 <-feols(scale(del_ship) ~ 1 + scale(l_VA) , data = df)

m26 <-feols(scale(del_ship) ~ 1 + scale(l_VA) , data = df, subset = ~ change == 1)



esttex( m22 ,  m23 , m24 , m26 ,
    file   = "Writeup/Tables/tab_va.tex",
    keep = "l_VA", 
    dict = c(VA = "VA_{t}"),
    replace = TRUE,          
    label  = "tab:value_added", 
    title  = "",
    digits = 3,
    notes = "Column iv) only invludes observations where the captain changed changed ship. " )

#keep   = "%age",
#    dict = c(re_prod = "Production/Duration", l_n_voyage = "Experience", age  = "Age"),




  ##################### Figure 6 ########################### 

  # replicate figure 6 of last week. 

