library(tidyverse)
library(readxl)
library(refineR)
library(clipr)
 

 

m <- function() {
  

# ptz_data ----------------------------------------------------------------
# sample_id sex   age_year result_num age_tam_sayi
  
  ptz_ref_data <- read_excel("data/coa_results/PT-2 2025.10.17.xlsx")  |> 
    janitor::clean_names() |> 
    mutate(age_tam_sayi  = ordered(floor(age_year))) |> 
    select(-age_group) |> 
    filter(age_tam_sayi >= 1 & age_tam_sayi <18)
  
  
  
  ptz_ref_data <- read_excel("data/coa_results/tce 2026 1-18 pt.xlsx")  |> 
    janitor::clean_names() |> 
    select(sample_id, sex, age_year, result_num, ) |> 
    mutate(age_tam_sayi  = ordered(floor(age_year))) |> 
    select(-age_year) #|> 
    # filter(age_tam_sayi >= 1 & age_tam_sayi <18)
  
  
  ptz_ref_data |> count(age_tam_sayi)

 
  aptt_ref_data <- read_excel("data/coa_results/aPTT 2025.10.17.xlsx")  |> 
    janitor::clean_names() |> 
    select(sample_id:result_num) |> 
    mutate(age_tam_sayi  = ordered(floor(age_year))) |> 
    select(-age_group) |> 
    filter(age_tam_sayi >= 1 & age_tam_sayi <18)
  
  
  
  aptt_ref_data <- read_excel("data/coa_results/tce 2026 1-18 aptt.xlsx")  |> 
    janitor::clean_names() |> 
    select(sample_id, sex, age_year, result_num, ) |> 
    mutate(age_tam_sayi  = ordered(floor(age_year))) |> 
    select(-age_year) #|> 
  # filter(age_tam_sayi >= 1 & age_tam_sayi <18)
  
  aptt_ref_data |> count(age_tam_sayi)
  
  
  fib_ref_data <- read_excel("data/coa_results/Fib-2 2025.10.17.xlsx")  |> 
    janitor::clean_names() |> 
    select(sample_id:result_num) |> 
    mutate(age_tam_sayi  = ordered(floor(age_year))) |> 
    select(-age_group) |> 
    filter(age_tam_sayi >= 1)
  
  
  fib_ref_data <- read_excel("data/coa_results/tce 2026 1-18 fib.xlsx")  |> 
    janitor::clean_names() |> 
    select(sample_id, sex, age_year, result_num, ) |> 
    mutate(age_tam_sayi  = ordered(floor(age_year))) |> 
    select(-age_year) 
  
  
  
  ptz_ref_data |> 
    count(age_tam_sayi )
  
  aptt_ref_data |> 
    count(age_tam_sayi )
  
  fib_ref_data |> 
    count(age_tam_sayi )
  
  
  ptz_ref_data |> 
    filter(result_num < 30) |> 
    group_by(age_tam_sayi, sex   ) |> 
    mutate(med = median(result_num )) |> 
    ungroup() |> 
    ggplot(aes(x =age_tam_sayi ,  y = result_num , color = sex)) +
    geom_point(position = position_jitterdodge(), alpha=0.5) +
    geom_boxplot(aes(color = sex))  +
    geom_line(aes(y=med, color = sex, group = sex), linewidth = 2)    + 
    theme_bw() 
  
 
  
  # 
  # openxlsx::write.xlsx(ptz_ref_data, "ptz_ref_data.xlsx")
  # 
  
  # ptz_ref_1_18 <- refineR::findRI(ptz_ref_data$result_num, NBootstrap = 200)
  # 
  # saveRDS(ptz_ref_1_18, "ptz_ref_1_18.RDS")
  # 
  
  ptz_ref_1_18 <- readRDS("data/processed/ptz_ref_1_18.RDS")
  ptz_ref_1_18 #8.56 - 10.9  ///Üretci 8.4 - 10.6
  

  refineR::getRI(ptz_ref_1_18, RIperc = c(0.025,0.50, 0.975) )
  
  ptz_ref_data_1_12 <- ptz_ref_data |>
    filter(age_tam_sayi >= 1, age_tam_sayi <=12 )
  ptz_ref_1_12 <- refineR::findRI(ptz_ref_data_1_12$result_num, NBootstrap = 200)
  ptz_ref_1_12 # 8.55 - 10.8
  saveRDS(ptz_ref_1_12, "data/processed/ptz_ref_1_12.RDS")
  
  
  ptz_ref_data_1_12 <- readRDS("data/processed/ptz_ref_1_12.RDS")
  
  
  ptz_ref_by_age_raw <- ptz_ref_data |>  
    group_by(age_tam_sayi) |> 
    nest() |> 
    ungroup() |> 
    arrange(age_tam_sayi) |> 
    # slice(1:1) |> 
    mutate(n = map_dbl(data, nrow)) |> 
    filter(n>10) |> 
    ungroup() |> 
    mutate(ri = map(data, ~findRI(.x$result_num, NBootstrap = 200))  ) 
  
  saveRDS(ptz_ref_by_age_raw, "data/processed/ptz_ref_by_age_raw.RDS")
  
  ptz_ref_by_age_raw <- readRDS("data/processed/ptz_ref_by_age_raw.RDS")
  
  
  ptz_ref_by_age_raw 
  
  

  ptz_ref_by_age <- ptz_ref_by_age_raw |> 
    mutate(ric = map(ri,  function(x) getRI(x, RIperc = c(0.025,0.50, 0.975)) )) |> 
    unnest(ric) |> 
    select(age_tam_sayi, n, Percentile:CIHigh) |> 
    mutate(RI_type = "refiner") |> 
    select(RI_type, everything()) |> 
    pivot_wider(names_from = Percentile, values_from = c(PointEst, CILow, CIHigh),names_prefix = "l_" ) |> 
    arrange(age_tam_sayi) |> 
    select(RI_type, age_tam_sayi, n, contains("l_0.025"), contains("l_0.5"), contains("l_0.975") )
  
  ptz_ref_by_age |> write_clip()
  
  
  ptz_ref_by_age |> 
    ggplot(aes(x = age_tam_sayi)) +
    geom_errorbar(aes(ymin = CILow_l_0.025   , ymax = CIHigh_l_0.025     ), linewidth =1) +
    geom_errorbar(aes(ymin = CILow_l_0.975    , ymax = CIHigh_l_0.975    ), linewidth =1) +
    geom_errorbar(aes(ymin = PointEst_l_0.025     , ymax = PointEst_l_0.975    ), linewidth =1, color = "red") +
    geom_point(aes(y=PointEst_l_0.5   ),  color = "blue", size = 2) +
    # geom_line(aes(y=PointEst_l_0.975   ),  color = "orange", size = 2, group = 1) +
    ggtitle("ptz")
  
  
  ptz_ref_by_age |> write_clip()
  
 
  
  # aptt ----------------------------------------------------------------
  
  
  
-
  
  
  aptt_ref_1_18 <- refineR::findRI(aptt_ref_data$result_num, NBootstrap = 100)
  
  aptt_ref_1_18 #  23.8 - 35.2  ///Üretici   UL:30.6    // Mevcut 32.4
  
  saveRDS(aptt_ref_1_18, "data/processed/aptt_ref_1_18.RDS")
  
    aptt_ref_1_18 <- readRDS("data/processed/aptt_ref_1_18.RDS")
  
  
  aptt_ref_data_1_11 <- aptt_ref_data |> 
    filter(age_tam_sayi >= 1, age_tam_sayi < 12 )
  
  
  aptt_ref_data_1_11 |> count(age_tam_sayi)
  
  aptt_ref_1_11 <- refineR::findRI(aptt_ref_data_1_12$result_num, NBootstrap = 200)
  saveRDS(aptt_ref_1_11, "data/processed/aptt_ref_1_11.RDS")
  
  aptt_ref_1_11 <- readRDS("data/processed/aptt_ref_1_11.RDS")
  
  aptt_ref_1_11 #  24.3-35.3
  refineR::getRI(aptt_ref_1_11, RIperc = c(0.025,0.50, 0.975) )
  
  aptt_ref_data_12_18 <- aptt_ref_data |> 
    filter(age_tam_sayi >= 12)
  
  
  aptt_ref_data_12_18 |> count(age_tam_sayi)
  
  aptt_ref_data_12_18 <- refineR::findRI(aptt_ref_data_12_18$result_num, NBootstrap = 200)
  saveRDS(aptt_ref_data_12_18, "data/processed/aptt_ref_data_12_18.RDS")
  
  aptt_ref_data_12_18 <- readRDS("data/processed/aptt_ref_data_12_18.RDS")
  
  aptt_ref_data_12_18 #  24-35.5
  

 
  
  
  
  aptt_ref_by_age_raw <- aptt_ref_data |>  
    group_by(age_tam_sayi) |> 
    nest() |> 
    ungroup() |> 
    arrange(age_tam_sayi) |> 
    # slice(1:1) |> 
    mutate(n = map_dbl(data, nrow)) |> 
    filter(n>10) |> 
    ungroup() |> 
    mutate(ri = map(data, ~findRI(.x$result_num, NBootstrap = 200))  ) 
  
  saveRDS(aptt_ref_by_age_raw, "data/processed/aptt_ref_by_age_raw.RDS")
  
  aptt_ref_by_age_raw <- readRDS("data/processed/aptt_ref_by_age_raw.RDS")
  
 
  
  aptt_ref_by_age <- aptt_ref_by_age_raw |> 
    mutate(ric = map(ri,  function(x) getRI(x, RIperc = c(0.025,0.50, 0.975)) )) |> 
    unnest(ric) |> 
    select(age_tam_sayi, n, Percentile:CIHigh) |> 
    mutate(RI_type = "refiner") |> 
    select(RI_type, everything()) |> 
    pivot_wider(names_from = Percentile, values_from = c(PointEst, CILow, CIHigh),names_prefix = "l_" ) |> 
    arrange(age_tam_sayi) |> 
    select(RI_type, age_tam_sayi, n, contains("l_0.025"), contains("l_0.5"), contains("l_0.975") )
  
  
  
  aptt_ref_by_age
  
  aptt_ref_by_age |> clipr::write_clip()
  
  aptt_ref_by_age |> 
    ggplot(aes(x = age_tam_sayi)) +
    geom_errorbar(aes(ymin = CILow_l_0.025   , ymax = CIHigh_l_0.025     ), linewidth =1) +
    geom_errorbar(aes(ymin = CILow_l_0.975    , ymax = CIHigh_l_0.975    ), linewidth =1) +
    geom_errorbar(aes(ymin = PointEst_l_0.025     , ymax = PointEst_l_0.975    ), linewidth =1, color = "red") +
    geom_point(aes(y=PointEst_l_0.5   ),  color = "blue", size = 2) +
    # geom_line(aes(y=PointEst_l_0.975   ),  color = "orange", size = 2, group = 1) +
    ggtitle("aptt")
  
  aptt_ref_by_age |> 
    ggplot(aes(x = age_tam_sayi)) +
    geom_errorbar(aes(ymin = CILow_l_0.975    , ymax = CIHigh_l_0.975    ), linewidth =1) 
  

  


#   # Sex -----------------------------------------------------------------

  
  
   ptz_ref_by_age_sex_raw <- ptz_ref_data |>  
     group_by(sex, age_tam_sayi) |> 
     nest() |> 
     ungroup() |> 
     arrange(age_tam_sayi, sex) |> 
     # slice(1:1) |> 
     mutate(n = map_dbl(data, nrow)) |> 
     filter(n>10) |> 
     ungroup() |> 
     mutate(ri = map(data, ~findRI(.x$result_num, NBootstrap = 200))  ) 
  
  saveRDS(ptz_ref_by_age_sex_raw, "data/processed/ptz_ref_by_age_sex_raw.RDS")
  
  ptz_ref_by_age_sex_raw <- readRDS("data/processed/ptz_ref_by_age_sex_raw.RDS")
  

  
 
  
  ptz_ref_by_age_sex <- ptz_ref_by_age_sex_raw  |> 
    mutate(ric = map(ri,  function(x) getRI(x, RIperc = c(0.025,0.50, 0.975), CIprop = 0.90) )) |> 
    unnest(ric) |> 
    select(age_tam_sayi,sex, n, Percentile:CIHigh) |> 
    mutate(RI_type = "refiner") |> 
    select(RI_type, everything()) |> 
    pivot_wider(names_from = Percentile, values_from = c(PointEst, CILow, CIHigh),names_prefix = "l_" ) |> 
    arrange(age_tam_sayi) |> 
    select(RI_type, age_tam_sayi, sex, n, contains("l_0.025"), contains("l_0.5"), contains("l_0.975") )
  
  
  ptz_ref_by_age_sex |> write_clip()
  
  
  ptz_ref_by_age_sex |> 
    ggplot(aes(x = age_tam_sayi, color = sex)) +
    geom_errorbar(aes(ymin = CILow_l_0.025, ymax = CIHigh_l_0.025 ), linewidth =1) +
    geom_errorbar(aes(ymin = CILow_l_0.975 , ymax = CIHigh_l_0.975  ), linewidth =1)
  
  
  ptz_ref_by_age_sex |> 
    ggplot(aes(x = age_tam_sayi, color = sex)) +
    geom_errorbar(aes(ymin = CILow_l_0.975 , ymax = CIHigh_l_0.975  ), linewidth =1)
  
  
  
  
  
  
  plot(ptz_ref_by_age_sex_raw$ri[[1]])
  
  
  
  aptt_ref_by_age_sex_raw <- aptt_ref_data |>  
    group_by(sex, age_tam_sayi) |> 
    nest() |> 
    ungroup() |> 
    arrange(age_tam_sayi, sex) |> 
    # slice(1:1) |> 
    mutate(n = map_dbl(data, nrow)) |> 
    filter(n>10) |> 
    ungroup() |> 
    mutate(ri = map(data, ~findRI(.x$result_num, NBootstrap = 200))  ) 
  
  saveRDS(aptt_ref_by_age_sex_raw, "data/processed/aptt_ref_by_age_sex_raw.RDS")
  
   aptt_ref_by_age_sex_raw <- readRDS("data/processed/aptt_ref_by_age_sex_raw.RDS")
  
  aptt_ref_by_age_sex_raw 
  
  
  aptt_ref_by_age_sex <- aptt_ref_by_age_sex_raw  |> 
    mutate(ric = map(ri,  function(x) getRI(x, RIperc = c(0.025,0.50, 0.975), CIprop = 0.90) )) |> 
    unnest(ric) |> 
    select(age_tam_sayi,sex, n, Percentile:CIHigh) |> 
    mutate(RI_type = "refiner") |> 
    select(RI_type, everything()) |> 
    pivot_wider(names_from = Percentile, values_from = c(PointEst, CILow, CIHigh),names_prefix = "l_" ) |> 
    arrange(age_tam_sayi) |> 
    select(RI_type, age_tam_sayi, sex, n, contains("l_0.025"), contains("l_0.5"), contains("l_0.975") )
  
  
  aptt_ref_by_age_sex |> write_clip()
  
  
  
  aptt_ref_by_age_sex |> 
    ggplot(aes(x = age_tam_sayi, color = sex)) +
    geom_errorbar(aes(ymin = CILow_l_0.025   , ymax = CIHigh_l_0.025     ), linewidth =1)
  
  
  aptt_ref_by_age_sex |> 
    ggplot(aes(x = age_tam_sayi, color = sex)) +
    geom_errorbar(aes(ymin = CILow_l_0.975    , ymax = CIHigh_l_0.975    ), linewidth =1)
  
  
  aptt_ref_by_age_sex |> 
  ggplot(aes(x = age_tam_sayi, color = sex)) +
    geom_errorbar(aes(ymin = CILow_l_0.025, ymax = CIHigh_l_0.025 ), linewidth =1) +
    geom_errorbar(aes(ymin = CILow_l_0.975 , ymax = CIHigh_l_0.975  ), linewidth =1)
  
  
  # fib  ----------------------------------------------------------------
  
  
  
  fib_ref_1_18 <- refineR::findRI(fib_ref_data$result_num, NBootstrap = 200)
  
  fib_ref_1_18 #  1.88 - 3.81  ///Üretici   1.96-4.28   // Direct 2.00 -  4.13 
  
  saveRDS(fib_ref_1_18, "data/processed/fib_ref_1_18.RDS")
  
 
  fib_ref_1_18 <- readRDS("data/processed/fib_ref_1_18.RDS")
  
  
  fib_ref_by_age_raw <- fib_ref_data |>  
    group_by(age_tam_sayi) |> 
    nest() |> 
    ungroup() |> 
    arrange(age_tam_sayi) |> 
    # slice(1:1) |> 
    mutate(n = map_dbl(data, nrow)) |> 
    filter(n>10) |> 
    ungroup() |> 
    mutate(ri = map(data, ~findRI(.x$result_num, NBootstrap = 200))  ) 
  
   saveRDS(fib_ref_by_age_raw, "data/processed/fib_ref_by_age_raw.RDS")
  
   fib_ref_by_age_raw <-  readRDS("data/processed/fib_ref_by_age_raw.RDS")
  
  fib_ref_by_age <- fib_ref_by_age_raw |> 
    mutate(ric = map(ri,  function(x) getRI(x, RIperc = c(0.025,0.50, 0.975)) )) |> 
    unnest(ric) |> 
    select(age_tam_sayi, n, Percentile:CIHigh) |> 
    mutate(RI_type = "refiner") |> 
    select(RI_type, everything()) |> 
    pivot_wider(names_from = Percentile, values_from = c(PointEst, CILow, CIHigh),names_prefix = "l_" ) |> 
    arrange(age_tam_sayi) |> 
    select(RI_type, age_tam_sayi, n, contains("l_0.025"), contains("l_0.5"), contains("l_0.975") )
  
  
  
  fib_ref_by_age |> clipr::write_clip()
  
  
  
  
  install.packages(c("gamlss","gamlss.add"))
  
  
  
}
