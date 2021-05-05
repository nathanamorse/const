################################################################################
#--------------- DATA: CHARACTERISTICS of CONSTITUTIONAL DESIGN ---------------#

# 1) Load and clean EGM's replication data (constitutions)
# 2) Load and clean full longitudinal data (constitution-years)
# 3) Save data

################################################################################

# Packages
library(haven) # loading stata data
library(tidyverse) # data manipulation and plots
library(countrycode) # converting country names and codes
library(fastDummies) # converting categorical variables into dummies


#----- (1) DESIGN VARIABLES from EGM (2009) -----#

# Load EGM replication data on constitutional design (CCP) [tidyverse]
egm = read.delim("../data/ccp/egm/design_variables/design_variables.txt") %>%
  select(const=system_num, code=cowcode, everything())

# Load chronology data and calculate constitution ages (CCP) [tidyverse]
chron = read.delim("../data/ccp/ccpcce_v1_3/ccpcce_v1_3.txt") %>%
  select(const=systid, year, code=cowcode, evntid, evnttype) %>%
  arrange(const, year) %>% # put rows in order
  group_by(const) %>% # group by constitution
  mutate(prom = min(year), # year promulgated
         age = 1:n(), # count age of constitution
         failure = (age==max(age) & year<2018)) # indicates if const was replaced

# Combine data on constitutional design and age [tidyverse, countrycode]
const = chron %>%
  select(-year) %>% # remove 'year' to avoid confusion with 'prom' (year promulgated)
  filter(age==max(age)) %>% # retain just the latest row for each constitution
  merge(egm, by=c("const", "code")) %>% # merge with design data
  
  # Code region and era variables
  mutate(region = countrycode(code, "cown", "region"),
         region_lam = (region=="Latin America & Caribbean"), 
         region_eur = (region=="Europe & Central Asia"),
         region_mid = (region=="Middle East & North Africa"),
         region_afr = (region=="Sub-Saharan Africa"),
         region_sas = (region=="South Asia"),
         region_eas = (region=="East Asia & Pacific"),
         era_mid = (prom>=1914 & prom<1946),
         era_late = (prom>1946)) %>%
  
  # Organize dataset
  select(code, year=prom, const, age:failure, legacy, interim:ppi, 
         region_lam:era_late) %>%
  arrange(code, year)


#----- (2) FULL CONSTITUTIONAL CHARACTERISTICS DATA -----#

# Numeric variables
j_num = c("length", "docs", "preambw", "rghtwrds", "demonum")

# Numeric variables to be removed
j_num_remove = c("hosterm", "hogterm", "agterm", "lhseats", "uhseats", 
                 "uhterm", "quorumw","chfterm", "ordterm", "adterm", 
                 "conterm", "ecterm")

# Character variables to be removed
j_char = c("regions", "hosname", "depname", "lhname", "uhname", "supname",
           "ordname", "exinst1", "langoffw", "langnatw")

# Constitution data (CCP) [haven, tidyverse]
cnc_full = read_dta("../data/ccp/ccpcnc_v3.dta") %>%
  
  # Select, rename, and remove columns
  select(code=cowcode, year, const=systid, prom=systyear, 
         all_of(j_num), -all_of(j_num_remove), -all_of(j_char), 
         ends_with("age") & !marriage, model:achighed & !all_of(j_num), 
         -ends_with("comments"), -ends_with("article")) %>%
  
  # Recode variables
  mutate(
    # Recode age limits: 1 means age of adulthood. Extract and replace with 18
    across(ends_with("age") & !marriage, as.numeric), # numericize age variables
    across(ends_with("age") & !marriage, ~.==1, .names="{.col}_maj"), # extract values
    across(ends_with("age") & !marriage, ~ifelse(.==1,18,.)), # replace with 18
    
    # Replace values with NA
    across(all_of(j_num), ~na_if(., "")), # replace blank with NA
    across(model:achighed, ~na_if(., 97)), # replace 97 with NA
    across(model:achighed, ~na_if(., 98)), # replace 98 with NA
    across(length:achighed, ~na_if(., 99)), # replace 99 with NA
    
    # Recode numerical variables
    across(all_of(j_num), ~as.numeric(as.character(.)))
  ) %>%
  
  # Calculate duration and survival variables
  group_by(const) %>%
  mutate(
    age = 1:n(), # age of constitution
    failure = ifelse(age==max(age) & year<2019,1,0) # indicates ended
  ) %>%
  
  # Organize dataset
  select(code:prom, age, failure, everything()) %>%
  arrange(code, year)

# Categorical variable names
j_cat = names(cnc_full)[which(names(cnc_full)=="model"):which(names(cnc_full)=="achighed")]

# Convert categorical variables into dummies [fastDummies]
cnc_wide1 = dummy_cols(cnc_full, j_cat, ignore_na=TRUE, remove_selected_columns=TRUE)

# Make sure all columns are numeric [tidyverse]
cnc_wide = mutate(cnc_wide1, across(everything(), as.numeric))


#----- (3) SAVE DATA-----#
save(egm, chron, const, cnc_wide, file="../data/const.Rdata")



