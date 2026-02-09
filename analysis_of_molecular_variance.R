# Excoffier L., Smouse PE., Quattro JM. 1992. 
# Analysis of molecular variance inferred from metric distances 
# among dna haplotypes: Application to human mitochondrial dna 
# restriction data. Genetics 131:479–491. 
# Available at: http://www.genetics.org/content/131/2/479.abstract


# The main goal is to illustrate how to manually perform analysis of molecular variance

library(dplyr)

# Simulation 

# parameters
num_ind <- 4 # number of individuals in each population (assume each population has same number of individuals)
num_pop <- 3 # number of populations in each group (assume each group has same number of populations)
num_grp <- 2 # number of groups
num_sit <- 10 # number of restriction sites

tot_ind <- num_ind*num_pop*num_grp # total number of individuals

# simulate a haplotype for each individual
hap_dat <- vector("list", tot_ind)

set.seed(1234)
rand_seeds <- sample(1000:2000, tot_ind, replace = FALSE)

for (i in 1:tot_ind) {
  set.seed(rand_seeds[i])
  hap_dat[[i]] <- sample(c(0,1), num_sit, replace = TRUE, prob = c(0.5, 0.5))
}

# group information for each individual
grp_info <- data.frame(ind = seq(tot_ind)
                       , group = factor(rep(1:num_grp
                                            , each = num_pop*num_ind)) 
                       , pop = factor(rep(1:num_pop
                                          , each = num_ind
                                          , times = num_grp)
                       )
)


# Calculate Euclidean distance between each pair of haplotypes
# Set weight matrix W to identity matrix

hap_pairs <- expand.grid(grp_info$ind, grp_info$ind)

d2 <- vector("list", nrow(hap_pairs))

for (i in seq_along(hap_pairs$Var1)) {
  j <- hap_pairs[i,"Var1"]
  k <- hap_pairs[i,"Var2"]
  
  pj <- hap_dat[[j]]
  pk <- hap_dat[[k]]
  
  d2[[i]] <- (pj - pk) %*% (pj - pk) # equation 3a
}

d2_dat <- ( hap_pairs 
            |> mutate(d2 = unlist(d2)) 
            |> rename(ind1 = Var1, ind2 = Var2)
            |> left_join(grp_info, by = join_by(ind1 == ind))
            |> rename(group1 = group, pop1 = pop)
            |> left_join(grp_info, by = join_by(ind2 == ind))
            |> rename(group2 = group, pop2 = pop)
            )

# Calculate SSD(total)
SSD_tot <- sum(d2_dat$d2)/(2*tot_ind) # equation 5b
print(SSD_tot)

# Calculate SSD(within populations), equation 8a
SSD_wp <- 0

for (g in 1:num_grp) {
  d2_grp <- subset(d2_dat, group1 == g & group2 == g) # both haplotypes in same group
  for (i in 1:num_pop) {
    d2_pop <- subset(d2_grp, pop1 == i & pop2 == i) # both haplotypes in same population
    N_pop <-  length(unique(d2_pop$ind1)) # number of haplotypes in a pop
    
    current_sum <- sum(d2_pop$d2)/(2*N_pop)
    SSD_wp <- SSD_wp + current_sum
  }
}
print(SSD_wp)

# Calculate SSD(among populations/within groups), equation 8b
SSD_apwg <- 0

for (g in 1:num_grp) {
  d2_grp <- subset(d2_dat, group1 == g) # both haplotypes in same group
  N_grp <- length(unique(d2_grp$ind1))
  
  pop_tot <- 0
  for (i in 1:num_pop) {
    d2_pop <- subset(d2_grp, pop1 == i & pop2 == i)
    N_pop <-  length(unique(d2_pop$ind1))
    
    pop_sum <- sum(d2_pop$d2)/(2*N_pop)
    pop_tot <- pop_tot + pop_sum
  }
  grp_sum <- sum(d2_grp$d2)/(2*N_grp) # 2 haplotypes can in any populations
  
  SSD_apwg <- SSD_apwg + (grp_sum - pop_tot)
}
print(SSD_apwg)

# Calculate SSD(among groups), equation 8c
for (g in 1:num_grp) {
  d2_grp <- subset(d2_dat, group1 == g) # both haplotypes in same group
  N_grp <- length(unique(d2_grp$ind1))
  
  grp_sum <- sum(d2_grp$d2)/(2*N_grp)
}
SSD_ag <- SSD_tot - grp_sum
print(SSD_ag)

# Show sum of squares add up to total sum of squares
print(SSD_tot == SSD_wp + SSD_apwg + SSD_ag) # equation 7



# How to do it using general linear model?
# How to do it with different different number of populations in each groups? 
# Different number of individuals in each population?
