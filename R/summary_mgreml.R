library(tidyverse)
library(data.table)

read_vc <- function(file) {
  vc <- fread(file)

  vg <- vc %>%
    filter(V1 == "genetic covariance") %>%
    select(V2, V3, estimate) %>%
    spread(V2, estimate) %>%
    column_to_rownames("V3")

  vg[upper.tri(vg)] <- 0
  tmp_vg <- vg
  diag(tmp_vg) <- 0
  final_vg <- vg + t(tmp_vg)

  # vg_se <- vc %>%
  #   filter(V1 == "genetic covariance") %>% 
  #   select(V2, V3, `standard error`) %>%
  #   spread(V2, `standard error`) %>%
  #   column_to_rownames("V3")

  # vg_se[upper.tri(vg_se)] <- 0
  # tmp_vg_se <- vg_se
  # diag(tmp_vg_se) <- 0
  # final_vg_se <- vg_se + t(tmp_vg_se)

  ve <- vc %>%
    filter(V1 == "environment covariance") %>%
    select(V2, V3, estimate) %>%
    spread(V2, estimate) %>%
    column_to_rownames("V3")

  ve[upper.tri(ve)] <- 0
  tmp_ve <- ve
  diag(tmp_ve) <- 0
  final_ve <- ve + t(tmp_ve)

  # ve_se <- vc %>%
  #   filter(V1 == "environment covariance") %>%
  #   select(V2, V3, `standard error`) %>%
  #   spread(V2, `standard error`) %>%
  #   column_to_rownames("V3")

  # ve_se[upper.tri(ve_se)] <- 0
  # tmp_ve_se <- ve_se
  # diag(tmp_ve_se) <- 0
  # final_ve_se <- vg_se + t(tmp_ve_se)
  
  return(list(vg = as.matrix(final_vg), ve = as.matrix(final_ve)))
}