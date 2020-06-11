# Figures
# Follows script "data_preparation.R"

cols <- c("lightblue", "lightpink")
cols2 <- c("blue", "red")
calc_mort <- function(dat, spec_codes) {
  tapply(!dat$alive, 
    list(spec_codes$Woody[match(dat$spec, spec_codes$Shortcut)], 
    dat$treat), function(x) sum(x) / length(x))
}
aux_bar <- function(dat, space){
  barplot(dat, add = TRUE, col = cols2, space = space, 
    names.arg = "", axes = FALSE)
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
aux_axis <- function(x_shifts, y_shift, labs = c(1,10,100)){
  sapply(x_shifts, function(x) axis(2, labels = labs, 
    at = log(c(1,10,100)) + y_shift, line = x, cex.axis = 0.8))
  invisible()
}
s_y1 <- -8
s_y2 <- -16

#___________________________________________________________________________
# Figure 1 - mortality
# png("figures/Fig_1_mortality.png", height=480*10, width=480*10, res=72*10)
par(mai=c(0.92,0.82,0,0))
mort <- calc_mort(dat_wf, spec_codes)
plot(c(0, 15.5),c(0, 0.8), type = "n", axes = FALSE, xlab = "", 
  ylab = "Mortality [%]")
aux_bar(mort[, 1], space = c(0, 0))
aux_bar(mort[, 2], space = c(2.5, 0))

mort <- calc_mort(dat_sf, spec_codes)
aux_bar(mort[, 1], space = c(5.5, 0))
aux_bar(mort[, 2], space = c(8, 0))

mort <- calc_mort(dat_he, spec_codes)
aux_bar(mort[, 1], space = c(11, 0))
aux_bar(mort[, 2], space = c(13.5, 0))

axis(2, labels = c(0:4*20), at = c(0:4*20)/100)
axis(1, labels = rep(c("Contr", "Treat"), 3), at = c(1,3.5,6.5,9,12,14.5), 
  tick = FALSE, cex.axis = 0.8, line = -2)
axis(1, labels = c("Winter\nfreezing","Spring\nfreezing","Herbivory\n"), 
  at = c(2.25,7.75,13.25), tick = FALSE, line = 0.5)
legend(0, 0.8, fill = cols2, legend = c("Herbaceous","Woody"), bty = "n")

par(xpd = NA)
arrows(x0 = 0.5,y0 = -0.13, x1 = 15,y1 = -0.13, code = 3, length = 0.15)
text(2.5, -0.145, "Predictable", cex = 0.8)
text(13, -0.145, "Unpredictable", cex = 0.8)

#___________________________________________________________________________
# Figure 2 - size
# png("figures/Fig_2_size.png", height=480*10, width=480*10, res=72*8)
par(mai = c(1.02,1.72,0,0))
plot(c(1, 8.2), c(-16, 5.5), type = "n", axes = FALSE, xlab = "", ylab = "")
aux_lines(dat_wf, "plainw", 0)
aux_lines(dat_wf, "length", sp[1])
aux_lines(dat_wf, "nleaves", sp[2])
aux_lines(dat_wf, "bio", sp[3], shift = log(10))

aux_lines(dat_sf, "plainw", 0, shift = s_y1)
aux_lines(dat_sf, "length", sp[1], shift = s_y1)
aux_lines(dat_sf, "nleaves", sp[2], shift = s_y1)

aux_lines(dat_he, "plainw", 0, shift = s_y2)
aux_lines(dat_he, "length", sp[1], shift = s_y2)

axis(1, labels = rep(c("Contr", "Treat"), 4), 
  at = sort(c(1.1,1.9,sp+1.1,sp+1.9)), tick = FALSE, cex.axis = 0.8)
axis(1, labels = c("Plain view\n[cm]", "Height\n[cm]", "Leaves\n[number]", "Biomass\n[g]"), 
  at = c(0,sp)+1.5, tick = FALSE, line = 2.5)
axis(2, labels = c("Herbivory","Spring \nfreezing","Winter \nfreezing"), 
  at = c(s_y2,s_y1,0)+log(10), line = 1.5, tick = FALSE, las = 2, 
  cex.axis = 1.2)
legend(sp[2]+2.2, s_y2+2, fill = cols2, legend = c("Herbaceous","Woody"), 
  bty = "n")

cf <- 4.27
# cf <- 3.15  # shift for console printing
aux_axis(c(0,-sp[1:2] * cf), 0)
aux_axis(-sp[3] * cf, 0, labs = c(0.1,1,10))
aux_axis(c(0,-sp[1:2] * cf), s_y1)
aux_axis(c(0,-sp[1] * cf), s_y2)

par(xpd = NA)
ssh <- 0.5
# ssh <- -0.2  # shift for console printing
arrows(x0 = -1.5+ssh, y0 = -17, x1 = -1.5+ssh, y1 = 5, code = 3, length = 0.2)
text(-1.65+ssh, 2, "Predictable", srt = 90)
text(-1.65+ssh, -14, "Unpredictable", srt = 90)
#_