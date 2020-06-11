# Data_preparation
library(ape) 

# Data_loading
dat_wf <- read.table("data/dat_wf.csv", sep = ",", header = TRUE)
dat_sf <- read.table("data/dat_sf.csv", sep = ",", header = TRUE)
dat_he <- read.table("data/dat_he.csv", sep = ",", header = TRUE)
spec_codes <- read.table(file = "data/spec_codes.csv", sep = ",", header = TRUE)
tree <- read.tree("data/selected_tree.tre")

# Add half-minimum to zero values (for logarithmic transformation)
dat_wf$nleaves[which(dat_wf$nleaves == 0)] <- 
  min(dat_wf$nleaves[dat_wf$nleaves > 0], na.rm = T) / 2
dat_wf$bion[which(dat_wf$bion == 0)] <- 
  min(dat_wf$bion[dat_wf$bion > 0], na.rm = T) / 2

## Auxiliary function to prep. data for Stan
prep_dat <- function(dataset, variable, tree, spec_codes, startval = FALSE, 
  t_log = FALSE){
  if (startval) {
    sel_var <- cbind(dataset[, variable], 
      dataset[, paste0(variable, "_start")])
  }else{
    sel_var <- dataset[, variable, drop = FALSE]
  }
  sel <- !is.na(rowSums(sel_var))
  if (t_log) sel_var <- log(sel_var)
  spec_codes$Species <- gsub(" ", "_", spec_codes$Species)
  sel_tree <- keep.tip(tree, 
    spec_codes$Species[spec_codes$Shortcut %in% dataset$spec[sel]])
  tree_vcv <- vcv(sel_tree)
  res <- list(
    Nspec = length(unique(dataset$spec[sel])),
    Ntot = length(dataset$spec[sel]),
    FinSize = sel_var[sel, 1],
    Treat = dataset$treat[sel],
    Woody = spec_codes$Woody[match(unique(dataset$spec[sel]), spec_codes$Shortcut)],
    phy = tree_vcv[order(rownames(tree_vcv)), order(colnames(tree_vcv))],
    DruhP = table(dataset$spec[sel])
  )
  if (startval) res <- c(res, list(StaSize = sel_var[sel, 2]))
  res
}

#_
