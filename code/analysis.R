################################################################################
#------------- REPLICATION of ELKINS, GINSBURG, and MELTON (2009) -------------#

# 1) Load data, run Cox model, store residuals
# 2) Run Cox model
# 3) Run preliminary random forests to determine the 9 most important
#    environmental variables. The data inclues 9 design variables, so only the 
#    top 9 environmental variables are retained to keep comparisons balanced.
# 4) Define formulas for models
# 5) Split data and run random forests and Cox models
# 6) Calculate model performance metrics
# 7) Calculate variable importance metrics

################################################################################


#----- (1) PREPARE CONSTITUTIONAL SURVIVAL DATASET -----#

# Combine data [tidyverse]
surv = merge(enviro, const, by=c("code", "year")) %>% # merge data
  select(const, code:year, age:failure, global_const:ecpc, legacy:era_late) %>%
  arrange(const)

# Keep track of predictors and formula for imputation
surv_vars = paste(names(select(surv, global_const:era_late)), collapse=" + ")
surv_form = as.formula(paste("Surv(age, failure)", surv_vars, sep=" ~ "))

# Impute missing data [randomForestSRC]
set.seed(2000)
surv_imp = impute(surv_form, data=surv, # formula and data
                  nodesize = 1, nsplit = 10) # hyperparameters


#----- (2) COX MODEL -----#

# Prepare data [tidyverse]
surv_cox = surv_imp %>% 
  select(!const:year) %>% # remove ID cols
  mutate(amend_rate2 = amend_rate^2) # square amendment rate

# Cox model [survival]
cox_full = coxph(Surv(age, failure) ~., surv_cox)


#----- (3) DETERMINE 9 BEST ENVIRONMENTAL VARIABLES -----#

# Load variable information
vars = read.csv("data/variables.csv") %>% filter(!grepl("region", variable))
vars_enviro = vars$variable[vars$type=="enviro"]

# Set formulas for forests
y_s0 = "Surv(age, failure)"
y_b0 = "Surv(age, backslid)"
x_s0 = paste(vars$variable, collapse=" + ")
x_b0 = paste(filter(vars, bs)$variable, collapse=" + ")
form_s0 = as.formula(paste(y_s0, x_s0, sep=" ~ "))
form_b0 = as.formula(paste(y_b0, x_b0, sep=" ~ "))

# Train models [randomForestSRC]
set.seed(2000)
rf_s0 = rfsrc(form_s0, data=surv_imp)
rf_b0 = rfsrc(form_b0, data=bs_imp)

# Extract variable importance [randomForestSRC]
vimp_s0 = vimp(rf_s0)$importance
vimp_b0 = vimp(rf_b0)$importance

# Determine top environmental variables for survival [tidyverse]
vs0 = data.frame(variable=names(vimp_s0), vimp=vimp_s0*100) %>%
  filter(variable %in% vars_enviro) %>%
  arrange(-vimp)

# Determine top environmental variables for backsliding [tidyverse]
vb0 = data.frame(variable=names(vimp_b0), vimp=vimp_b0*100) %>%
  filter(variable %in% vars_enviro) %>%
  arrange(-vimp)


#----- (4) DEFINE MODEL SPECIFICATIONS -----#

# Save number of variables of each type
k_d = length(which(vars$type=="design"))
k_e = length(which(vars$type=="enviro"))
k = min(k_d, k_e)

# Define variables
vars_all = vars$variable
vars_design = vars$variable[vars$type=="design"]
vars_s_env = vs0$variable[1:k]
vars_b_env = vb0$variable[1:k]
vars_control = vars$variable[vars$type=="control"]

# Paste variables together for formulas
x_all = paste(vars_all, collapse=" + ")
x_design = paste(c(vars_design, vars_control), collapse=" + ")
x_s_env = paste(c(vars_s_env, vars_control), collapse=" + ")
x_b_env = paste(c(vars_b_env, vars_control), collapse=" + ")
y_surv = "Surv(age, failure)"
y_bs = "Surv(age, backslid)"

# Set formulas for models
form_surv = as.formula(paste(y_surv, x_all, sep=" ~ "))
form_surv_design = as.formula(paste(y_surv, x_design, sep=" ~ "))
form_surv_enviro = as.formula(paste(y_surv, x_s_env, sep=" ~ "))
form_bs = as.formula(paste(y_bs, x_all, sep=" ~ "))
form_bs_design = as.formula(paste(y_bs, x_design, sep=" ~ "))
form_bs_enviro = as.formula(paste(y_bs, x_b_env, sep=" ~ "))


#----- (5) RANDOM SURVIVAL FORESTS and COX MODELS -----#

# Sample indices to split data
n = nrow(surv_imp) # sample size
p = 0.7 # proportion to put in training data
set.seed(1500) # set seed for random sample
train = sample(n, p*n) # randomly select cases for training data

# Cox models [survival]
cox_surv = coxph(form_surv, surv_imp[train,], x=TRUE)
cox_bs = coxph(form_bs, bs_imp[train,], x=TRUE)

# Train random forests [randomForestSRC]
set.seed(2000)
rf_surv = rfsrc(form_surv, data=surv_imp[train,])
rf_surv_design = rfsrc(form_surv_design, data=surv_imp[train,])
rf_surv_enviro = rfsrc(form_surv_enviro, data=surv_imp[train,])
rf_bs = rfsrc(form_bs, data=bs_imp[train,])
rf_bs_design = rfsrc(form_bs_design, data=bs_imp[train,])
rf_bs_enviro = rfsrc(form_bs_enviro, data=bs_imp[train,])


#----- (6) MODEL PERFORMANCE METRICS with CROSS VALIDATION -----#

# Bootstrap cross-validation on constitutional survival models [pec]
c_surv = cindex(list(`Cox (Full)`=cox_surv, # list models
                     `RF (Full)`=rf_surv, 
                     `RF (Design)`=rf_surv_design, 
                     `RF (Enviro)`=rf_surv_enviro),
                formula = Surv(age, failure)~1, # model form
                data = surv_imp[-train,], # run on test data
                splitMethod="BootCv", B=100) # bootstrapping specs

# Bootstrap cross-validation on democratic backsliding models [pec]
c_bs = cindex(list(`Cox (Full)`=cox_bs, # list models
                   `RF (Full)`=rf_bs, 
                   `RF (Design)`=rf_bs_design, 
                   `RF (Enviro)`=rf_bs_enviro),
              formula = Surv(age, backslid)~1, # model form
              data = bs_imp[-train,], # run on test data
              splitMethod="BootCv", B=100) # bootstrapping specs

# Store c-indexes in dataframe [tidyverse]
perf = data.frame(model = names(unlist(c_surv$BootCvCindex)),
                  surv = unlist(c_surv$BootCvCindex),
                  bs = unlist(c_bs$BootCvCindex)) %>%
  
  # Organize dataframe
  pivot_longer(surv:bs, names_to="outcome", values_to="c") %>% # reshape
  arrange(desc(outcome), c) %>% # order rows
  mutate(outcome = factor(outcome, levels=unique(outcome)), # turn into factor
         model = factor(model, levels=unique(model)), # turn into factor
         order = row_number()) # keep track of order


#----- (7) COMPARE VARIABLE IMPORTANCE -----#

# Variable importance [randomForestSRC, tidyverse]
rf_surv_vimp = vimp(rf_surv)
surv_vimp = data.frame(variable=names(rf_surv_vimp$importance), # variable names
                       vimp=as.numeric(rf_surv_vimp$importance)*100) %>% # vimp scores
  merge(vars, by="variable") %>% # merge with variable information
  filter(type != "control") %>% # remove control variables
  mutate(index = sapply(variable, grep, names(data))) %>% 
  # get the column number from data corresponding to variable
  arrange(-vimp)

# Summary statistics for variable importance [tidyverse]
surv_vimp_sum = surv_vimp %>%
  group_by(type) %>% # group by type (design/enviro)
  summarize(sum = round(sum(vimp),2), # sum scores
            avg = round(mean(vimp),2)) %>% # average scores
  pivot_longer(sum:avg, names_to="variable", values_to="vimp") %>% # reshape data
  mutate(name = c("Design (sum)", # add labels
                  "Design (average)",
                  "Environment (sum)", 
                  "Environment (average)")) %>%
  select(variable, vimp, name, type)

# Variable importance for graph [tidyverse]
surv_vimp2 = surv_vimp %>%
  add_row(surv_vimp_sum) %>% # append summary statistics to importance data
  filter(variable != "sum") %>% # remove sums
  arrange(vimp) %>% # put rows in order
  mutate(avg = (variable=="avg"), # indicator of summary statistic
         name = factor(name, levels=name)) # make factor to retain order



################################################################################
#----------------------------- REVISION INTERVALS -----------------------------#

# 1) Calculate revision intervals
# 2) Merge with constitution data
# 3) Principal component analysis
# 4) Train random survival forests

################################################################################


#----- (1) REVISION INTERVALS -----#

# Revisions
rev = chron %>% ungroup() %>%
  
  # Recode events
  mutate(type = ifelse(evnttype %in% c("new", "reinstated"), "new", evnttype),
         type = ifelse(type=="samendment", "amendment", type),
         type = ifelse(type=="non-event", "non", type)) %>%
  
  # Remove non-event rows except last one to indicate right-censoring
  group_by(code) %>% arrange(code, year) %>%
  filter(type %in% c("new", "amendment") | row_number()==n()) %>%
  
  # Country-level statistics
  mutate(events = length(which(type!="non")),
         state_age = year - min(year),
         final_age = max(year) - min(year),
         avg_dur = final_age/events,
         avg_const = final_age/length(which(type=="new")),
         num_const = cumsum(type=="new")) %>%
  
  # Constitution-level statistics
  ungroup() %>% group_by(code, const) %>%
  mutate(num_amends = 1:n()) %>%
  
  # Significant amendments
  mutate(type = ifelse(!(type %in% c("new", "non")) & (avg_const>10) &
                         (year-lag(year) > avg_dur), 
                       "sig", type),
         type = factor(type, levels=c("new", "amendment", "sig", "non"))) %>%
  select(code, year, const, const_age=age, type:num_amends)


# Significant revisions
sigrev = rev %>% ungroup() %>% group_by(code) %>%
  filter(type != "amendment") %>%
  mutate(dur = year - lag(year),
         type = ifelse(type=="non", 0, ifelse(type=="sig", 1, 2))) %>%
  filter(!is.na(dur)) %>%
  select(code:const, dur, type, avg_dur, const_age, state_age, 
         num_const, num_amends)


#----- (2) CONSTITUTION DATA -----#

# Merge data on constitutional events and characteristics
sr = sigrev %>%
  merge(cnc_wide, by=c("code", "year", "const")) %>% # merge data
  
  # Remove unneeded rows and columns
  select(where(~mean(is.na(.))<.3)) %>% # remove variables with >40% NAs
  select(which(apply(., 2, var, na.rm=TRUE) != 0), type) %>% # remove constants
  select(code:dur, type, everything(), -age, -failure)

# Merge data with EGM replication data
sr_egm = sigrev %>%
  merge(select(surv_imp, !year:failure), by=c("code", "const"))

# Impute missing data [randomForestSRC]
sr_imp = cbind(select(sr, code:const),
               impute.rfsrc(Surv(dur ~ type) ~ ., # formula
                            data=select(sr, !code:const), # data
                            nsplit=3))  # hyperparameters


#----- (3) DIMENSION REDUCTION -----#

# Principal component analysis
pca = select(sr_imp, !code:type) %>%
  prcomp(scale=TRUE, rank=25)

# Store component scores and add environmental data
sr_pc = cbind(select(sr, code:type), predict(pca)) %>%
  merge(enviro, by=c("code", "year"), all.x=TRUE)


#----- (4) RANDOM FORESTS -----#

# Sample indices to split data
n = nrow(sr_pc) # sample size
p = 0.7 # proportion to put in training data
set.seed(1000) # set seed for random sample
train = sample(n, p*n) # randomly select cases for training data

# Remove ID cols
sr_rf_data = select(sr_pc, !code:const)
sr_egm_data = select(sr_egm, !code:year)
sr_cox_data = mutate(sr_egm_data, type = ifelse(type>0, 1, 0))

# Train forests
cox_sr = coxph(Surv(dur, type) ~., x=TRUE,  # formula
               data=sr_cox_data[train,]) # data
rf_sr_egm = rfsrc(Surv(dur, type) ~., # formula
                  data = sr_egm_data[train,]) # data
rf_sr = rfsrc(Surv(dur, type) ~., # formula
              data = sr_rf_data[train,], # data
              na.action = "na.impute") # impute missing data

# Bootstrap cross-validation on constitutional survival models [pec]
c_sr1 = cindex(list(`Cox`=cox_sr), formula = Surv(dur, type)~1, 
               data = sr_cox_data[-train,], splitMethod="BootCv", B=10)
c_sr2 = cindex(list(`RF (Original)`=rf_sr_egm), formula = Surv(dur, failure)~1, 
               data = sr_egm_data[-train,], splitMethod="BootCv", B=10)
c_sr3 = cindex(list(`RF (PCA)`=rf_sr), formula = Surv(dur, failure)~1, 
               data = sr_rf_data[-train,], splitMethod="BootCv", B=10)

