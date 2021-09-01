# Parsing  Data
# obtained from: https://github.com/Riselya/Prevalence-Thresholds-Metaanalysis

# Load RDS object from data folder
## contains human, meerkat, deer, spiny rat, carollia, mouse lemur, flamingo, stint
## we want only the human
phylo <- readRDS("data/phylo_list_unrarefied_frontiers.RDS")

# Rarefy as done in the original paper (10,000 across all data sets)
