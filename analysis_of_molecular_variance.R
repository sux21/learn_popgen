# general linear model: Y = Xb + e
# Y: actual values
# Xb: predicted values
# goal: find b such that || Y-Xb || is minimum. 


# Analysis of molecular variance 
# p_jig = p + a_g + b_ig + c_jig
# g: group, i: population, j: individual
# p_ijg: observed values
# p: mean
# a_g: group effect
# b_ig: population effect
# c_jig: individuals effect


# How p_jig = p + a_g + b_ig + c_jig relates to Y = Xb + e?
# p_jig forms the vector Y
# coefficients of p, a_g, b_ig form matrix X
# p, all a_g, all b_ig form marix b
# c_jig forms the vector e



library(dplyr)




# simulaion

# parameters
num_ind <- 4 # number of individuals in each population (assume each population has same number of individuals)
num_pop <- 3 # number of populations in each group (assume each group has same number of populations)
num_grp <- 2 # number of groups
lambda <- 10 # mean of Poisson distribution

set.seed(1234)

sim <- tibble(group = factor(rep(1:num_grp, each = num_pop*num_ind)) 
              , pop = factor(rep(1:num_pop, each = num_ind, times = num_grp))
              , p = rpois(num_grp*num_pop*num_ind, lambda)
)
print(sim)

means <- summarise(sim, 
          mean = mean(p),
          .by = c(group, pop)
)

tot_ind <- nrow(sim) # total number of individuals
grand_mean <- mean(sim$p) # mean of all individuals


# Calculate SSD(total)
SS_tot <- (sim$p - grand_mean)^2 |> sum() # sum of squares (total) 

# Calculate SSD(WP)
for (g in 1:num_grp) {
  for (i in 1:num_pop) {
    sim_sub <- subset(sim, group == g & pop == i)
    means_sub <- subset(means, group == g & pop == i)
    
    val <- sim_sub$p
    m <- means_sub$mean
    
    ss <- (val-m)^2 |> sum()
    print(ss)
   # SS_wp <- SS_wp + ss
  }
}

# calculate SSD(AP/WG)
grp <- levels(sim$group)
pop <- levels(sim$pop)

for (g in seq_along(grp)) {
  for (i in seq_along(pop)) {
    for (ii in seq_along(pop)[-i]) {
      sim_sub <- subset(sim, group == g & pop == i)
      means_sub <- subset(means, group == g & pop == ii)
      
      val <- sim_sub$p
      m <- means_sub$mean
      
      ss <- (val-m)^2 |> sum()
      print(ss)
    }
  }
}

# 
# 
# print(sim
#       |> summarize(m = mean(d2), .by=c(group, pop))
# )
# 
# 
# 
# fit.aov <- aov(p ~ group + pop, data = sim)
# summary(fit.aov)
# 
# 
# model <- (lm(p ~ group+pop, data = sim))
# print(model)
# summary(model)
# 
# ## predict(m)
# dat <- (dat
#         |> mutate(
#           fit = predict(m)
#           , resid=val-fit
#         )
# )
# 
# dat <- summarize(dat,
#                  sP=var(val)
#                  , sG = var(fit)
#                  , sI = var(resid)
# )
# print(dat)
# 
# quit()
# summary(m)
# anova(m)
# car::Anova(m)
# 
# ## Explore residuals and variance
# 
# # this may be interesting to repeat the analysis in the article
# # Table 2 data
# region <- c("Asia", "West_Africa", "America", "Europe", "Middle_East")
# num_region <- length(region)
# 
# pop <- c("Tharu", "Oriental", "Wolof", "Peul", "Pima", "Maya", "Finnish",
#          "Sicilian", "Israeli_Jews", "Israeli_Arabs")
# num_pop <- length(pop)
# 
# geo_dat <- data.frame(region = rep(region, each = 2)
#                       , pop = pop)
# 
# hap <- c(1, 2, 6, 7, 8, 9, 
#          10, 11, 12, 13, 17, 18, 
#          21, 22, 23, 27, 28, 29, 
#          31, 34, 36, 37, 38, 39, 
#          40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 
#          50, 51, 52, 53, 54, 56, 57, 
#          64, 65, 66, 67, 68, 69, 
#          71, 72, 73, 75, 76, 77, 
#          82, 83, 
#          95)
# num_hap <- length(hap)
# 
# as.vector(sapply(c(1,2,3,4,5), rep, 10))
# 
# hap_dat <- data.frame(hap = c(rep(1, 10), 
# )
# fre = c(),)