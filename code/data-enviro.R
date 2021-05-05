################################################################################
#--------------------- DATA: ENVIRONMENTS of CONSTITUTIONS --------------------#

# 1) Load and clean external data on environmental variables (country-years)
# 3) Combine environmental data into single dataset (country-years)
# 3) Save data

################################################################################

# Packages
library(haven) # loading stata data
library(tidyverse) # data manipulation and plots
library(zoo) # data manipulation
library(countrycode) # converting country names and codes
load("../data/const.Rdata") # load constitution data


#----- (1) LOAD SOURCE DATASETS -----#

# Calculate number of new constitutions in the world each year [tidyverse]
global = chron %>%
  group_by(year) %>% # collapse by year
  summarise(global_const = sum(age==1)) # count new constitutions per year

# Regime types and transitions (Polity) [tidyverse]
polity = read.csv("../data/polity/p5v2018.csv", na.strings="") %>%
  mutate_all(funs(replace(., .<(-10), NA))) %>% # recode all NAs
  group_by(ccode) %>% # collapse by country
  mutate(change = polity - lag(polity), # change in score from previous year
         prev3 = rollsumr(change, 3, fill=NA), # sum change over last 3 years
         dem_trans = (prev3>2 & change!=0), # dem trans if score increased 3+
         aut_trans = (prev3<0 & change!=0), # aut trans if score decreased 3+
         dem = (polity>5)) %>% # democratic regime if score is 6+
  select(code=ccode, year, dem_trans, aut_trans, dem, polity) %>%
  arrange(code, year)

# Territorial changes (COW) [tidyverse]
tc = read.csv("../data/cow/terr-changes-v6/tc2018.csv") %>%
  filter(gainer>1, loser>1) # remove NAs
tgain = tc %>% 
  select(year, code=gainer) %>% # focus on countries that gained territory
  mutate(terr_gain=TRUE) %>% # indicator of gain
  arrange(year, code) 
tloss = tc %>% 
  select(year, code=loser) %>% # focus on countries that lost territory
  mutate(terr_loss=TRUE) %>%  # indicator of loss
  arrange(year, code)

# Interstate conflict (COW) [tidyverse]
mida = read.csv("../data/cow/mids/MIDA 5.0.csv") %>% # dispute-level data
  filter(hostlev>3, outcome==2) # filter to wars won by side B
midb = read.csv("../data/cow/mids/MIDB 5.0.csv") %>% # country-level data
  mutate(defeat = (dispnum %in% mida$dispnum) & (sidea==1)) %>%
  # determine if country was side A in a war won by side B
  filter(defeat) %>%
  select(code=ccode, year=endyear, defeat) %>%
  arrange(code, year)

# Intrastate conflict (COW) [tidyverse]
intra1 = read.csv("../data/cow/intra/INTRA-STATE WARS v5.1 CSV.csv") %>%
  select(code=CcodeA, contains("Yr")) %>%
  mutate_all(funs(replace(., .<0, NA))) %>% # recode negatives as NA
  rowwise() %>%
  mutate(start = min(c_across(2:9), na.rm=TRUE), # first year of dispute
         end = max(c_across(2:9), na.rm=TRUE), # last year of dispute
         years = list(seq(start, end))) %>% # list of all years
  select(code, years) %>% 
  arrange(code)
intra2 = polity %>% select(code, year) %>% # initialize country-year dataframe
  mutate(dom_crisis=FALSE) # initialize indicator of domestic crisis
for (i in 1:nrow(intra2)) {
  if (intra2$code[i] %in% intra1$code) { # if country had any disputes
    x = which(intra1$code==intra2$code[i]) # indexes of country's disputes
    intra2$dom_crisis[i] = (intra2$year[i] %in% unlist(intra1$years[x]))
    # determine if a dispute occurred in that country-year
  }
}

# Economic development and crisis (COW) [tidyverse]
nmc = read.csv("../data/cow/nmc/NMC_5_0.csv") %>%
  group_by(ccode) %>% # collapse by country
  mutate(ecpc = pec/tpop, # energy consumption per capita (ECPC)
         ecpc_prev = lag(ecpc), # previous year
         diff = ecpc/ecpc_prev, # ECPC growth ratio
         econ_crisis = (diff<.9)) %>% # crisis if ECPC decreases by >10%) %>% 
  select(code=ccode, year, ecpc, econ_crisis)

# Ethnic fractionalization (Drazanova 2020) [tidyverse]
hief = read.csv("../data/drazanova/HIEF_data.csv") %>%
  mutate(Country = ifelse(Country=="Democratic Republic of Vietnam", 
                          "Vietnam", Country), # rewrite name to prevent error
         code = countrycode(Country, "country.name", "cown", nomatch=340)) %>%
  # recode countries from names to COW codes
  select(code, year=Year, ethnic=EFindex) %>%
  arrange(code, year)

# Leadership changes (Archigos) [haven, tidyverse]
arch = read_dta("../data/archigos/Archigos_4.1_stata14.dta") %>%
  mutate(endyear = format(as.Date(enddate), "%Y"), # extract exit year
         startyear = format(as.Date(startdate), "%Y")) # extract entry year
arch_entries = arch %>% 
  select(code=ccode, year=startyear, entry) # store new leaders
arch_exits = arch %>% 
  select(code=ccode, year=endyear, exit) # store exiting leaders
arch_year = merge(arch_entries, arch_exits, by=c("code", "year"), all=TRUE) %>%
  # merge entries and exits into condensed data frame
  group_by(code, year) %>% # collapse by country-year
  summarise(leader_intra = sum((entry=="Regular") | (exit=="Regular"))>1,
            leader_extra = sum((entry=="Irregular") | (entry=="Foreign Imposition") |
                                 (exit=="Irregular") | (exit=="Foreign"))>1) %>%
  # determine if at least 1 regular and at least 1 irregular transition occurred
  arrange(code, year)


#----- (2) MERGING ENVIRONMENTAL DATA -----#

# Merge into single dataset [tidyverse, zoo]
enviro = select(polity, code, year) %>%
  
  # Global constitutional events
  merge(global, all.x=TRUE) %>%
  mutate(global_const = replace_na(global_const, 0)) %>%
  
  # Territory changes
  merge(tgain, all.x=TRUE) %>% 
  merge(tloss, all.x=TRUE) %>%
  mutate(terr_gain = replace_na(terr_gain, FALSE),
         terr_loss = replace_na(terr_loss, FALSE)) %>%
  
  # Defeat in war
  merge(midb, all.x=TRUE) %>% 
  mutate(defeat = replace_na(defeat, FALSE)) %>%
  
  # Domestic crisis
  merge(intra2, all.x=TRUE) %>% 
  
  # Economic crisis
  merge(select(nmc, !ecpc), all.x=TRUE) %>%
  mutate(econ_crisis = replace_na(econ_crisis, FALSE)) %>%
  
  # Regime transition
  merge(select(polity, code:aut_trans), all.x=TRUE) %>%
  
  # Leadership changes
  merge(arch_year, all.x=TRUE) %>% 
  mutate(leader_extra = replace_na(leader_extra, FALSE),
         leader_intra = replace_na(leader_intra, FALSE)) %>%
  
  # Democracy
  merge(select(polity, code, year, dem), all.x=TRUE) %>%
  
  # Ethnic heterogeneity
  merge(hief, all.x=TRUE) %>% 
  
  # Economic development
  merge(select(nmc, !econ_crisis), all.x=TRUE) %>%
  
  # Calculate 5-year rolling averages for each year [zoo]
  group_by(code) %>% # group by country
  mutate(across(global_const:dem, rollmaxr, k=5, fill=NA), # sum last 5 years
         across(ethnic:ecpc, rollmeanr, k=5, fill=NA)) %>% # average last 5 years
  
  # Organize dataset
  arrange(code, year) %>%
  select(code, year, global_const:ecpc) %>%
  distinct()


#----- (3) SAVE DATA-----#
save(enviro, polity, file="../data/enviro.Rdata")




#----- (4) PREPARE CONSTITUTIONAL SURVIVAL DATASET -----#

# Combine data [tidyverse]
surv = merge(enviro, const, by=c("code", "year")) %>% # merge with constitution data
  select(const, code:year, age:failure, global_const:ecpc, legacy:era_late) %>%
  arrange(const)

# Keep track of predictors and formulas for imputation
surv_vars = paste(names(select(surv, global_const:era_late)), collapse=" + ")
surv_form = as.formula(paste("Surv(age, failure)", surv_vars, sep=" ~ "))

# Impute missing data [randomForestSRC]
set.seed(2000)
surv_imp = impute(surv_form, data=surv, # formula and data
                  nodesize = 1, nsplit = 10) # hyperparameters


#----- (5) PREPARE DEMOCRATIC BACKSLIDING DATASET -----#

# Democratic backsliding data [tidyverse]
bs = polity %>%
  
  # Merge polity and constitution chronology data
  merge(chron, by=c("code", "year")) %>%
  select(const, code, age, prom, year, aut_trans) %>%
  arrange(const, year) %>%
  
  # Determine first occurrence of democratic backsliding
  group_by(const) %>% # group by constitution
  mutate(count = cumsum(ifelse(is.na(aut_trans),0,aut_trans))) %>%
  # count events for each constitution
  filter((aut_trans & count==1) | count==0) %>% # remove repeated events
  filter(age==max(age)) %>% # retain just the latest row for each constitution
  
  # Organize dataset and merge with survival dataset
  select(const:age, backslid=aut_trans) %>%
  merge(select(surv, !year:failure), by=c("const", "code")) %>%
  arrange(const)

# Keep track of predictors and formulas for imputation
bs_vars = paste(names(select(bs, global_const:era_late)), collapse=" + ")
bs_form = as.formula(paste("Surv(age, backslid)", bs_vars, sep=" ~ "))

# Impute missing data [randomForestSRC]
set.seed(2000)
bs_imp = impute(bs_form, data=bs, # formula and data
                nodesize = 1, nsplit = 10) # hyperparameters

