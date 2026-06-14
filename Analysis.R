# Sebastian Nöth
# 09 May 2026

######################################################################################################################################################################################
#### LOAD AND CLEAN DATA ####
# Set working directory
setwd("D:/Uni/Master/2. Semester/Project/Data")

# Read in the foraminifera morphology data 

foram_morph <- read.csv(file = "Just_spinose_foram_data.csv", sep = ";")

# Load in the Triton data (Available at https://doi.org/10.1038/s41597-021-00942-7)
load("triton_no_dup.RData")

# make a working copy of the data that we can modify
triton_working <- triton.no.dup.pres

# Adding a column where the genus and species names are separated by an underscore instead of a space
# This can sometimes help avoid computation issues or the computer "getting confused"

triton_working$gen_sp <- gsub(" ", "_", triton_working$species)

# Sometimes, microfossils can become "reworked", i.e. the sediment can get mixed up, so fossils appear higher or lower in the sediment
# than where they should really be. This can lead falsely-assigned age estimates

# First, subset to include only fossils that occur before the species went extinct
triton_working <- triton_working[triton_working$age >= triton_working$Extinction,]

# Then, subset to include only fossils that occur after the species originated
triton_working <- triton_working[triton_working$age <= triton_working$Speciation,]


######################################################################################################################################################################################
#### ASSIGN FOSSIL DATA TO TEMPORAL BINS ####

# For now, use bin size of one million years (you can experiment with this, for example 0.1 for 100,000 years etc...)

binsize <- 1

triton_working$age_bins <- cut(triton_working$age, breaks = seq(0, max(triton_working$age)+binsize, by=binsize),
                               labels = FALSE, include.lowest = TRUE)


######################################################################################################################################################################################
#### SUBSET INTO SPINOSE AND NON-SPINOSE ####
# Make smaller lists of just spinose and just non-spinose species:

spines <- foram_morph$Morphospecies[foram_morph$Spinose == 1]

no_spines <- foram_morph$Morphospecies[foram_morph$Spinose == 0]

# Using those lists of species, subset the Triton occurence data into two subsets

triton_spinose <- triton_working[triton_working$species %in% spines,]

triton_non <- triton_working[triton_working$species %in% no_spines,]


######################################################################################################################################################################################
#### DIVERSITY ANALYSIS WITH DIVDYN ####

# install.packages("divDyn")
library(divDyn)


divdyn_spinose <- divDyn::divDyn(triton_spinose, 
                                 bin="age_bins", 
                                 tax="gen_sp",
                                 revtime = T)

divdyn_non <- divDyn::divDyn(triton_non, 
                             bin="age_bins", 
                             tax="gen_sp",
                             revtime = T)


######################################################################################################################################################################################
#### MORPHOLOGY ####
# extinction rate per capita
mean(divdyn_non$extPC, na.rm = T)

mean(divdyn_spinose$extPC, na.rm = T)


# standard deviation and standard error
sd_non <- sd(divdyn_non$extPC, na.rm = TRUE)
sd_spinose <- sd(divdyn_spinose$extPC, na.rm = TRUE)

se_non <- sd_non / sqrt(length(no_spines))
se_spinose <- sd_spinose / sqrt(length(spines))


# speciation rate
mean(divdyn_non$oriPC, na.rm = T)
mean(divdyn_spinose$oriPC, na.rm = T)


# turnover rate (speciation rate + extinction rate)
turnover_non <- divdyn_non$extPC + divdyn_non$oriPC
turnover_spinose <- divdyn_spinose$extPC + divdyn_spinose$oriPC
mean(turnover_non, na.rm = T)
mean(turnover_spinose, na.rm = T)


# three-timer sampling completeness
mean(divdyn_non$samp3t, na.rm = T)
mean(divdyn_spinose$samp3t, na.rm = T)


# three-timer turnover
tt_turnover_non <- divdyn_non$ext3t + divdyn_non$ori3t
tt_turnover_spinose <- divdyn_spinose$ext3t + divdyn_spinose$ori3t
mean(tt_turnover_non, na.rm = T)
mean(tt_turnover_spinose, na.rm = T)


######################################################################################################################################################################################
#### ECOLOGY ####
# Extracting the different eco-groups
open_ocean_mixed_symb <- foram_morph$Morphospecies[foram_morph$EcogroupList == 1]
open_ocean_mixed_nosymb <- foram_morph$Morphospecies[foram_morph$EcogroupList == 2]
open_ocean_thermocline <- foram_morph$Morphospecies[foram_morph$EcogroupList == 3]
open_ocean_subthermocline <- foram_morph$Morphospecies[foram_morph$EcogroupList == 4]
Highlatitude <- foram_morph$Morphospecies[foram_morph$EcogroupList == 5]
Upwelling <- foram_morph$Morphospecies[foram_morph$EcogroupList == 6]


# Subsetting the Triton Data into the different eco-groups
triton_open_ocean_mixed_symb <- triton_working[triton_working$species %in% open_ocean_mixed_symb,]
triton_open_ocean_mixed_nosymb <- triton_working[triton_working$species %in% open_ocean_mixed_nosymb,]
triton_open_ocean_thermocline <- triton_working[triton_working$species %in% open_ocean_thermocline,]
triton_open_ocean_subthermocline <- triton_working[triton_working$species %in% open_ocean_subthermocline,]
triton_Highlatitude <- triton_working[triton_working$species %in% Highlatitude,]
triton_Upwelling <- triton_working[triton_working$species %in% Upwelling,]


# Extracting diversity parameters through divdyn
divdyn_open_ocean_mixed_symb <- divDyn::divDyn(triton_open_ocean_mixed_symb, 
                                               bin="age_bins", 
                                               tax="gen_sp",
                                               revtime = T)
divdyn_open_ocean_mixed_nosymb <- divDyn::divDyn(triton_open_ocean_mixed_nosymb, 
                                                 bin="age_bins", 
                                                 tax="gen_sp",
                                                 revtime = T)
divdyn_open_ocean_thermocline <- divDyn::divDyn(triton_open_ocean_thermocline, 
                                                bin="age_bins", 
                                                tax="gen_sp",
                                                revtime = T)
divdyn_open_ocean_subthermocline <- divDyn::divDyn(triton_open_ocean_subthermocline, 
                                                   bin="age_bins", 
                                                   tax="gen_sp",
                                                   revtime = T)
divdyn_Highlatitude <- divDyn::divDyn(triton_Highlatitude, 
                                      bin="age_bins", 
                                      tax="gen_sp",
                                      revtime = T)
divdyn_Upwelling <- divDyn::divDyn(triton_Upwelling, 
                                   bin="age_bins", 
                                   tax="gen_sp",
                                   revtime = T)

# Extinction Rate by eco-group
mean(divdyn_open_ocean_mixed_symb$extPC, na.rm = T)
mean(divdyn_open_ocean_mixed_nosymb$extPC, na.rm = T)
mean(divdyn_open_ocean_thermocline$extPC, na.rm = T)
mean(divdyn_open_ocean_subthermocline$extPC, na.rm = T)
mean(divdyn_Highlatitude$extPC, na.rm = T)
mean(divdyn_Upwelling$extPC, na.rm = T)


# Extinction Rate by eco-group for spinose and non-spinose respectively
# Ecogroup 1: Open ocean mixed layer tropical/subtropical,with symbionts
open_ocean_mixed_symb_spinose <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 1 &
    foram_morph$Spinose == 1
]

open_ocean_mixed_symb_non <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 1 &
    foram_morph$Spinose == 0
]

triton_open_ocean_mixed_symb_spinose <- triton_working[triton_working$species %in% open_ocean_mixed_symb_spinose,]
triton_open_ocean_mixed_symb_non <- triton_working[triton_working$species %in% open_ocean_mixed_symb_non,]

divdyn_open_ocean_mixed_symb_spinose <- divDyn::divDyn(triton_open_ocean_mixed_symb_spinose, 
                                                  bin="age_bins", 
                                                  tax="gen_sp",
                                                  revtime = T)
divdyn_open_ocean_mixed_symb_non<- divDyn::divDyn(triton_open_ocean_mixed_symb_non, 
                                             bin="age_bins", 
                                             tax="gen_sp",
                                             revtime = T)

mean(divdyn_open_ocean_mixed_symb_spinose$extPC, na.rm = T)
mean(divdyn_open_ocean_mixed_symb_non$extPC, na.rm = T)


#Ecogroup 2: Open ocean mixed layer tropical/subtropical,without symbionts
open_ocean_mixed_nosymb_spinose <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 2 &
    foram_morph$Spinose == 1
]

open_ocean_mixed_nosymb_non <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 2 &
    foram_morph$Spinose == 0
]


triton_open_ocean_mixed_nosymb_spinose <- triton_working[triton_working$species %in% open_ocean_mixed_nosymb_spinose,]
triton_open_ocean_mixed_nosymb_non <- triton_working[triton_working$species %in% open_ocean_mixed_nosymb_non,]

divdyn_open_ocean_mixed_nosymb_spinose <- divDyn::divDyn(triton_open_ocean_mixed_nosymb_spinose, 
                                                  bin="age_bins", 
                                                  tax="gen_sp",
                                                  revtime = T)
divdyn_open_ocean_mixed_nosymb_non<- divDyn::divDyn(triton_open_ocean_mixed_nosymb_non, 
                                             bin="age_bins", 
                                             tax="gen_sp",
                                             revtime = T)

mean(divdyn_open_ocean_mixed_nosymb_spinose$extPC, na.rm = T)
mean(divdyn_open_ocean_mixed_nosymb_non$extPC, na.rm = T)


#Ecogroup 3: Open ocean thermocline
open_ocean_thermocline_spinose <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 3 &
    foram_morph$Spinose == 1
]

open_ocean_thermocline_non <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 3 &
    foram_morph$Spinose == 0
]

triton_open_ocean_thermocline_spinose <- triton_working[triton_working$species %in% open_ocean_thermocline_spinose,]
triton_open_ocean_thermocline_non <- triton_working[triton_working$species %in% open_ocean_thermocline_non,]

divdyn_open_ocean_thermocline_spinose <- divDyn::divDyn(triton_open_ocean_thermocline_spinose, 
                                                  bin="age_bins", 
                                                  tax="gen_sp",
                                                  revtime = T)
divdyn_open_ocean_thermocline_non<- divDyn::divDyn(triton_open_ocean_thermocline_non, 
                                             bin="age_bins", 
                                             tax="gen_sp",
                                             revtime = T)

mean(divdyn_open_ocean_thermocline_spinose$extPC, na.rm = T)
mean(divdyn_open_ocean_thermocline_non$extPC, na.rm = T)
#We can still see a difference but the values are much closer together at 0.098 and 0.111


#Ecogroup 4: Open ocean sub-thermocline
open_ocean_subthermocline_spinose <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 4 &
    foram_morph$Spinose == 1
]

open_ocean_subthermocline_non <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 4 &
    foram_morph$Spinose == 0
]

triton_open_ocean_subthermocline_spinose <- triton_working[triton_working$species %in% open_ocean_subthermocline_spinose,]
triton_open_ocean_subthermocline_non <- triton_working[triton_working$species %in% open_ocean_subthermocline_non,]

divdyn_open_ocean_subthermocline_spinose <- divDyn::divDyn(triton_open_ocean_subthermocline_spinose, 
                                                  bin="age_bins", 
                                                  tax="gen_sp",
                                                  revtime = T)
divdyn_open_ocean_subthermocline_non<- divDyn::divDyn(triton_open_ocean_subthermocline_non, 
                                             bin="age_bins", 
                                             tax="gen_sp",
                                             revtime = T)

mean(divdyn_open_ocean_subthermocline_spinose$extPC, na.rm = T)
mean(divdyn_open_ocean_subthermocline_non$extPC, na.rm = T)
#There is a distinct difference but the values are much lower than in other 


#Ecogroup 5: High latitude
#Not enough data

#Ecogroup 6: Upwelling/highproductivity
#Very low amount of data, likely skewing the result
Upwelling_spinose <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 6 &
    foram_morph$Spinose == 1
]

Upwelling_non <- foram_morph$Morphospecies[
  foram_morph$EcogroupList == 6 &
    foram_morph$Spinose == 0
]


triton_Upwelling_spinose <- triton_working[triton_working$species %in% Upwelling_spinose,]
triton_Upwelling_non <- triton_working[triton_working$species %in% Upwelling_non,]

divdyn_Upwelling_spinose <- divDyn::divDyn(triton_Upwelling_spinose, 
                                                  bin="age_bins", 
                                                  tax="gen_sp",
                                                  revtime = T)
divdyn_Upwelling_non<- divDyn::divDyn(triton_Upwelling_non, 
                                             bin="age_bins", 
                                             tax="gen_sp",
                                             revtime = T)

mean(divdyn_Upwelling_spinose$extPC, na.rm = T)
mean(divdyn_Upwelling_non$extPC, na.rm = T) #not enough data to calculate


# Turnover by eco-group in per capita
turnover_open_ocean_mixed_symb <- divdyn_open_ocean_mixed_symb$extPC + divdyn_open_ocean_mixed_symb$oriPC
turnover_open_ocean_mixed_nosymb <- divdyn_open_ocean_mixed_nosymb$extPC + divdyn_open_ocean_mixed_nosymb$oriPC
turnover_open_ocean_thermocline <- divdyn_open_ocean_thermocline$extPC + divdyn_open_ocean_thermocline$oriPC
turnover_open_ocean_subthermocline <- divdyn_open_ocean_subthermocline$extPC + divdyn_open_ocean_subthermocline$oriPC
turnover_Highlatitude <- divdyn_Highlatitude$extPC + divdyn_Highlatitude$oriPC
turnover_Upwelling <- divdyn_Upwelling$extPC + divdyn_Upwelling$oriPC

mean(turnover_open_ocean_mixed_symb, na.rm = T)
mean(turnover_open_ocean_mixed_nosymb, na.rm = T)
mean(turnover_open_ocean_thermocline, na.rm = T)
mean(turnover_open_ocean_subthermocline, na.rm = T)
mean(turnover_Highlatitude, na.rm = T)
mean(turnover_Upwelling, na.rm = T)


# three-timer turnover by eco-group
tt_turnover_open_ocean_mixed_symb <- divdyn_open_ocean_mixed_symb$ext3t + divdyn_open_ocean_mixed_symb$ori3t
tt_turnover_open_ocean_mixed_nosymb <- divdyn_open_ocean_mixed_nosymb$ext3t + divdyn_open_ocean_mixed_nosymb$ori3t
tt_turnover_open_ocean_thermocline <- divdyn_open_ocean_thermocline$ext3t + divdyn_open_ocean_thermocline$ori3t
tt_turnover_open_ocean_subthermocline <- divdyn_open_ocean_subthermocline$ext3t + divdyn_open_ocean_subthermocline$ori3t
tt_turnover_Highlatitude <- divdyn_Highlatitude$ext3t + divdyn_Highlatitude$ori3t
tt_turnover_Upwelling <- divdyn_Upwelling$ext3t + divdyn_Upwelling$ori3t

mean(tt_turnover_open_ocean_mixed_symb, na.rm = T)
mean(tt_turnover_open_ocean_mixed_nosymb, na.rm = T)
mean(tt_turnover_open_ocean_thermocline, na.rm = T)
mean(tt_turnover_open_ocean_subthermocline, na.rm = T)
mean(tt_turnover_Highlatitude, na.rm = T)
mean(tt_turnover_Upwelling, na.rm = T)


# three-timer sampling completeness by eco-group
mean(divdyn_open_ocean_mixed_symb$samp3t, na.rm = T)
mean(divdyn_open_ocean_mixed_nosymb$samp3t, na.rm = T)
mean(divdyn_open_ocean_thermocline$samp3t, na.rm = T)
mean(divdyn_open_ocean_subthermocline$samp3t, na.rm = T)
mean(divdyn_Highlatitude$samp3t, na.rm = T)
mean(divdyn_Upwelling$samp3t, na.rm = T)


######################################################################################################################################################################################
#### ECOLOGY/MORPHOLOGY PROPORTIONS ####
# Count spinose and non-spinose species within each eco-group
sum_foraminifer_df <- as.data.frame.matrix(
  table(
    foram_morph$EcogroupList,
    foram_morph$Spinose
  )
)

colnames(sum_foraminifer_df) <- c(
  "Non-spinose",
  "Spinose"
)


# Assign/Add eco-groups
sum_foraminifer_df$Ecogroup <- c("EG1","EG2","EG3","EG4","EG5","EG6")


# Add necessary columns
sum_foraminifer_df <- sum_foraminifer_df[
  , c("Ecogroup", "Spinose", "Non-spinose")
]


# Calculate and add percentages
sum_foraminifer_df$Total <- 
  sum_foraminifer_df$Spinose +
  sum_foraminifer_df$`Non-spinose`

sum_foraminifer_df$Spinose_Percent <- 
  sum_foraminifer_df$Spinose /
  sum_foraminifer_df$Total * 100

sum_foraminifer_df$Nonspinose_Percent <- 
  sum_foraminifer_df$`Non-spinose` /
  sum_foraminifer_df$Total * 100

