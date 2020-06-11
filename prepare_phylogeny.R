# Prepare_phylogeny
# Requires phylogenetic tree from Zanne et al. 2014 (supplements) to be in the folder "data"
library(ape)

z_tree <- read.tree("data/Vascular_Plants_rooted.dated.tre")
spec_codes <- read.table(file = "data/spec_codes.csv", sep = ",", header = TRUE)

## Adjustments (3 species - not present in z_tree)
# Parietaria officinalis - Parietarias represent monophyletic clade in our 
#   dataset. Thus we take another species from the same genera
# Ribes uva-crispa - different spelling
# Epilobium angustifolium - in z_tree is synonym Chamerion angustifolium 
#   (based on The Plant List)
z_tree$tip.label[z_tree$tip.label == "Parietaria_debilis"] <-
  "Parietaria_officinalis"  
z_tree$tip.label[z_tree$tip.label == "Ribes_uva"] <-
  "Ribes_uva-crispa"  
z_tree$tip.label[z_tree$tip.label == "Chamerion_angustifolium"] <-
  "Epilobium_angustifolium"  

# Tree selection
tree <- keep.tip(z_tree, gsub(" ", "_", spec_codes$Species))
tree$node.label <- NULL
write.tree(tree, file = "data/selected_tree.tre")

#_