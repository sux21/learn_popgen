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



# parameters
num_ind <- 4 # number of individuals in each population (assume each population has same number of individuals)
num_pop <- 3 # number of populations in each group (assume each group has same number of populations)
num_grp <- 2 # number of groups
lambda <- 10 # mean of Poisson distribution

set.seed(1545)

sim <- tibble(group = factor(rep(1:num_grp, each = num_pop*num_ind)) 
              , pop = factor(rep(1:num_pop, each = num_ind, times = num_grp))
              , p = rpois(num_grp*num_pop*num_ind, lambda)
)
#contrasts(sim$group) <- "contr.sum"

print(sim)

print(sim
      |> summarize(m = mean(p), .by=c(group, pop))
)


model <- (lm(p ~ group+pop, data = sim))
print(model)
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