rm(list = ls())
library(tidyverse)
library(data.table)
library(jjutil)

read_vc <- function(file) {
  vc <- fread(file)

  vg <- vc %>%
    filter(V1 == "genetic covariance") %>%
    select(V2, V3, estimate) %>%
    spread(V2, estimate) %>%
    column_to_rownames("V3")
  vg[upper.tri(vg)] <- vg[lower.tri(vg)]

  vg_se <- vc %>%
    filter(V1 == "genetic covariance") %>% 
    select(V2, V3, `standard error`) %>%
    spread(V2, `standard error`) %>%
    column_to_rownames("V3")
  vg_se[upper.tri(vg_se)] <- vg_se[lower.tri(vg_se)]


  ve <- vc %>%
    filter(V1 == "environment covariance") %>%
    select(V2, V3, estimate) %>%
    spread(V2, estimate) %>%
    column_to_rownames("V3")
  ve[upper.tri(ve)] <- ve[lower.tri(ve)]


  ve_se <- vc %>%
    filter(V1 == "environment covariance") %>%
    select(V2, V3, `standard error`) %>%
    spread(V2, `standard error`) %>%
    column_to_rownames("V3")
  ve_se[upper.tri(ve_se)] <- ve_se[lower.tri(ve_se)]

  # browser()

  list(vg = vg, vg_se = vg_se, ve = ve, ve_se = ve_se)
}

read_vc("mgreml/tutorial/components.VCs.out")

read_vc("mgreml/tutorial/components.VCs.out")






fread("mgreml/tutorial/pheno.txt")
