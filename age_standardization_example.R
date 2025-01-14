library(tidyverse)



these_agecat <- paste0('Age',c('0-4','5-9','10-14','15-19','20-24','25-29',
                               '30-34','35-44','45-54','55-64','65-74','75-84','85+'))

# read SEER standard population data
seer_age19 <- read_fwf("stdpop.19ages.txt",
                       fwf_widths(c(3,3,8),
                                  c("standard","age","std_raw"))) %>%
  # In this example, we will use standard "201", which corresponds to
  # 2000 U.S. Std Million (19 age groups) according to the documentation
  filter(standard=="201") |>
  # In this example, we need to collapse some of the categories to match what we have available in the mortality data
  mutate(agecat=recode(age,
                       '000'="Age0-4",                   
                       '001'="Age0-4",
                       '002'="Age5-9",
                       '003'="Age10-14",
                       '004'="Age15-19",
                       '005'="Age20-24",
                       '006'="Age25-29",
                       '007'="Age30-34",
                       '008'="Age35-44",
                       '009'="Age35-44",
                       '010'="Age45-54",
                       '011'="Age45-54",
                       '012'="Age55-64",
                       '013'="Age55-64",
                       '014'="Age65-74",
                       '015'="Age65-74",
                       '016'="Age75-84",
                       '017'="Age75-84",
                       '018'="Age85+"),
         std.pop=as.numeric(std_raw)) |>
  group_by(agecat) |>
  summarise(std=sum(std.pop)) |>
  mutate(agecat = factor(agecat, levels = these_agecat))



# Create the example data frame -- this is simulated data
set.seed(123)  # for reproducibility

df_mortality_example <- data.frame(
  group  = rep(1:3, each = length(these_agecat)),
  agecat = rep(these_agecat, 3)
) |>
  dplyr::mutate(
    # Larger Poisson mean for population
    population = rpois(n = dplyr::n(), lambda = 10000),
    # Smaller Poisson mean for deaths
    deaths     = rpois(n = dplyr::n(), lambda = 10)
  )



# Now we do age standardization.
# As a reminder, look at the formulas for age-standardized rates,
# rate differences, and rate ratios from the Public Health Geocoding Project Monograph
# https://hsph.harvard.edu/wp-content/uploads/2024/10/2004_The-Public-Health-Disparities-Geocoding-Project-Monograph_krieger-et-al_v2_final_06-21-24.pdf

# We have to merge the dataset with the age standard, by age category
df_age_standardized <- df_mortality_example |>
  # merge with the SEER age standard
  left_join(seer_age19, by = "agecat") |>
  # calculate stratum-specific rates
  # and the variance of the stratum-specific rates
  mutate(rate = deaths / population,
         var_rate = deaths / (population^2)) |>
  # to aggregate over age strata, we need to group
  group_by(group) |>
  # summarize by taking weighted sums, where the weights are the 
  # weights from the SEER standard
  # Don't forget to divide by the sum of the weights
  dplyr::summarize(age_std_rate = sum(std * rate)/sum(std),
                   # Note that the variance of the age-standardized rate is itself 
                   # a weighted sum, where the weights  are 
                   var_age_std_rate = sum(std^2 * var_rate)/(sum(std)^2)) |>
  ungroup()


# Suppose that we want to compute age-standardized rate differences
# or rate ratios comparing groups 2 and 3 to group 1

# make a dataset with just the reference group
df_age_standardized_reference_rates <- df_age_standardized |>
    filter(group == 1) |>
    rename(ref_rate = age_std_rate,
           ref_var_rate = var_age_std_rate) |>
    # create a dummy variable to facilitate merging
    mutate(dummy = 1) |>
    select(-group)

df_age_standardized_comparison_groups <- df_age_standardized |>
  filter(!group == 1) |>
  mutate(dummy = 1) |>
  left_join(df_age_standardized_reference_rates, by = "dummy") |>
  select(-dummy) |>
  # compute age-standardized rate differences and age-standardized rate ratios
  mutate(rd_std = age_std_rate - ref_rate,
         var_rd_std = var_age_std_rate + ref_var_rate,
         rr_std = age_std_rate / ref_rate,
         var_log_irr_std = (var_age_std_rate / age_std_rate^2) + (ref_var_rate / ref_rate^2)) |>
  # compute 95% CI and express rate differences as rate differences per 100,000
  mutate(rd_std_per_100k = 1e+5 * rd_std,
         rd_std_per_100k_lo95 = 1e+5 * (rd_std - 1.96*var_rd_std),
         rd_std_per_100k_up95 = 1e+5 * (rd_std + 1.96*var_rd_std),
         rr_std_lo95 = exp(log(rr_std) - 1.96 * var_log_irr_std),
         rr_std_up95 = exp(log(rr_std) + 1.96 * var_log_irr_std))
         
         


