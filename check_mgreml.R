rm(list = ls())
library(tidyverse)
library(data.table)
library(gemma2) # https://github.com/fboehm/gemma2

library(jjutil) # devtools::install_github("junghyunJJ/jjuitl")
source("R/summary_mgreml.R")


# Run mgreml (https://github.com/devlaming/mgreml?tab=readme-ov-file#variance-components)
# python ./mgreml.py --grm ./tutorial/data --pheno ./tutorial/pheno.txt \
#                    --covar ./tutorial/covar.txt \
#                    --variance-components \
#                    --out ./tutorial/components


vc <- read_vc("data/components.VCs.out")

h2 <- fread("data/components.HSq.out") %>%
  column_to_rownames("V1")

gcor <- fread("data/components.RhoG.out") %>%
  column_to_rownames("V1")

