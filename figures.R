# Figures
# Follows script "data_preparation.R"

cols <- c("lightblue", "lightpink")
cols2 <- palette()[c(4,2)]
calc_mort <- function(dat, spec_codes) {
  aux_prop <- function(x) sum(x) / length(x)
  list(mort = tapply(!dat$alive, 
    list(spec_codes$Woody[match(dat$spec, spec_codes$Shortcut)], 
    dat$treat), aux_prop),
    se = tapply(!dat$alive, 
    list(spec_codes$Woody[match(dat$spec, spec_codes$Shortcut)], 
    dat$treat), function(x) sqrt((aux_prop(x) * (1 - aux_prop(x))) / length(x)))
  )
}
aux_bar <- function(dat, coln, space){
  barplot(dat$mort[, coln], add = TRUE, col = cols2, space = space, 
    names.arg = "", axes = FALSE)
  arrows(x0 = c(0.5, 1.5) + space[1], y0 = dat$mort[, coln], 
    y1 = dat$mort[, coln] + dat$se[, coln], length = 0.1, angle = 90)
}
calc_means <- function(dat, vec){
  res <- tapply(log(dat[, vec]), list(dat$spec, dat$treat), mean, na.rm = TRUE)
  cbind(res, W = spec_codes$Woody[match(rownames(res), spec_codes$Shortcut)])
}
sp <- c(2,4,6)
aux_lines <- function(dat, vec, sp, shift = 0) {
  aux <- calc_means(dat, vec)
  aux <- aux[!is.na(rowSums(aux)), ]
  apply(aux, 1, function(x, sp, shift) {
    lines(1:2 + sp, x[1:2] + shift, col = cols[x[3]+1])
  }, sp, shift)
  aux_m <- apply(aux, 2, tapply, aux[, 3], mean)
  apply(aux_m, 1, function(x, sp, shift) {
    lines(1:2 + sp, x[1:2] + shift, col = cols2[x[3]+1], lwd = 2)
  }, sp, shift)
  invisible()
}
aux_axis <- function(x_shifts, y_shift, labs = c(1,10,100), 
  unts = c("cm", "cm", "")){
  unts <- unts[1:length(x_shifts)]
  Map(function(x, unts) {
    lines(c(x+0.8,x+2.2), rep(log(0.2),2) + y_shift)
    lines(rep(x+0.8, 2), log(c(0.2,300)) + y_shift)
    sapply(log(c(1,10,100)), function(y) 
      lines(c(x+0.8,x+0.7), rep(y, 2) + y_shift))
    text(rep(x+0.65,3), log(c(1,10,100)) + y_shift, paste(labs, unts), 
      cex = 0.7, adj = c(1,0.5))
  }, x_shifts, unts)
  invisible()
}
s_y1 <- -8
s_y2 <- -16

#___________________________________________________________________________
# Figure 1 - phylogenetic tree
library(phytools)
tree_spec <- function(dat, codes, tree){
  sp <- table(gsub(" ", "_", 
    codes$Species[match(dat$spec, codes$Shortcut)]))
  sp_num <- sp[match(tree$tip.label, names(sp))]
  sp_num[is.na(sp_num)] <- 0
  sp_num
}
spec_mat <- cbind(
  tree_spec(dat_wf, spec_codes, tree), 
  tree_spec(dat_sf, spec_codes, tree),
  tree_spec(dat_he, spec_codes, tree)
)
rownames(spec_mat) <- tree$tip.label
tree_h <- as.numeric(!spec_codes$Woody[match(tree$tip.label, gsub(" ", "_", spec_codes$Species))])
spec_mat_b <- (spec_mat > 0) + 0
spec_mat[spec_mat == 0] <- ""

# png("figures/Fig_1_tree.png", height=480*10, width=480*10, res=72*10)
# jpeg("figures/Fig_1_tree.jpg", height=480*10, width=480*10, res=72*10, quality = 100)
par(mai=c(0,0,0.3,0))
plot(tree, label.offset = 32)
invisible(mapply(plotrix::draw.circle, x = 125, y = 1:27, radius = spec_mat_b[,1]*4, col = tree_h*2+2, lty = 0))
invisible(mapply(plotrix::draw.circle, x = 135, y = 1:27, radius = spec_mat_b[,2]*4, col = tree_h*2+2, lty = 0))
invisible(mapply(plotrix::draw.circle, x = 145, y = 1:27, radius = spec_mat_b[,3]*4, col = tree_h*2+2, lty = 0))
text(125, 1:27, spec_mat[, 1], cex = 0.75, adj = c(0.5,0.4), col = "white")
text(135, 1:27, spec_mat[, 2], cex = 0.75, adj = c(0.5,0.4), col = "white")
text(145, 1:27, spec_mat[, 3], cex = 0.75, adj = c(0.5,0.4), col = "white")
par(new = T, mai=c(0,0,0,0))
plot(c(0,250), c(0,28), type = "n", axes = F, ann = F)
text(c(125,135,145), 27.8, c("Wf.", "Sf.", "He."))

#___________________________________________________________________________
# Figure 2 - mortality
# png("figures/Fig_2_mortality.png", height=480*10, width=480*10, res=72*10)
# jpeg("figures/Fig_2_mortality.jpg", height=480*10, width=480*10, res=72*10, quality = 100)
par(mai=c(0.4,0.82,0,0))
mort <- calc_mort(dat_wf, spec_codes)
plot(c(0, 15.5),c(0, 0.8), type = "n", axes = FALSE, xlab = "", 
  ylab = "Mortality")
aux_bar(mort, coln = 1, space = c(0, 0))
aux_bar(mort, coln = 2, space = c(2.5, 0))

mort <- calc_mort(dat_sf, spec_codes)
aux_bar(mort, coln = 1, space = c(5.5, 0))
aux_bar(mort, coln = 2, space = c(8, 0))

mort <- calc_mort(dat_he, spec_codes)
aux_bar(mort, coln = 1, space = c(11, 0))
aux_bar(mort, coln = 2, space = c(13.5, 0))

abline(h = 0)
axis(2, labels = paste0(c(0:4*20), "%"), at = c(0:4*20)/100, las = 2)
axis(1, labels = rep(c("Contr", "Treat"), 3), at = c(1,3.5,6.5,9,12,14.5), 
  tick = FALSE, cex.axis = 0.8, line = -2)
axis(1, labels = c("Winter freezing","Spring freezing","Herbivory"), 
  at = c(2.25,7.75,13.25), tick = FALSE, line = -0.5)
legend(0, 0.8, fill = cols2, legend = c("Herbaceous","Woody"), bty = "n")

#___________________________________________________________________________
# Figure 3 - size
# png("figures/Fig_3_size.png", height=480*10, width=480*10, res=72*10)
# jpeg("figures/Fig_3_size.jpg", height=480*10, width=480*10, res=72*10, quality = 100)
par(mai = c(0.6,1.0,0,0))
plot(c(0.3, 8.2), c(-16.8, 5.5), type = "n", axes = FALSE, xlab = "", ylab = "")
aux_lines(dat_wf, "plainw", 0)
aux_lines(dat_wf, "length", sp[1])
aux_lines(dat_wf, "nleaves", sp[2])
aux_lines(dat_wf, "bio", sp[3], shift = log(10))

aux_lines(dat_sf, "plainw", 0, shift = s_y1)
aux_lines(dat_sf, "length", sp[1], shift = s_y1)
aux_lines(dat_sf, "nleaves", sp[2], shift = s_y1)

aux_lines(dat_he, "plainw", 0, shift = s_y2)
aux_lines(dat_he, "length", sp[1], shift = s_y2)

axis(1, labels = rep(c("Contr", "Treat"), 4), line = -0.8,
  at = sort(c(1.1,1.9,sp+1.1,sp+1.9)), tick = FALSE, cex.axis = 0.8)
axis(1, labels = c("Plain view", "Height", "Leaves", "Biomass"), 
  at = c(0,sp)+1.5, tick = FALSE, line = 0.7)
axis(2, labels = c("Herbivory","Spring \nfreezing","Winter \nfreezing"), 
  at = c(s_y2,s_y1,0)+log(10), line = -0.5, tick = FALSE, las = 2, 
  cex.axis = 1.2)
legend(sp[2]+2, s_y2+2, fill = cols2, legend = c("Herbaceous","Woody"), 
  bty = "n")

aux_axis(c(0,sp[1:2]), 0)
aux_axis(sp[3], 0, labs = c(0.1,1,10), unts = "g")
aux_axis(c(0,sp[1:2]), s_y1)
aux_axis(c(0,sp[1]), s_y2)

#_