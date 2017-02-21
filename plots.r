## Generate Figure 3a.
data <- read.csv("2cat_random/results.txt", header=TRUE, sep="\t")
z <- seq(0,0.999,0.001)
par(mar=c(5,5,1,1), mgp=c(3,.7,0))
qqplot(z, data$muni.p, xlim=c(0,1), ylim=c(0,1),  cex=0.75, pch=1,
       xlab="", ylab="p-value", cex.lab=2,  cex.axis=1.3, col="#ff7f00")
points(sort(z), sort(data$mantel.eucl.p),  col="#e41a1c", pch=16)
points(sort(z), sort(data$HSIC.p), col="#4daf4a", pch=15)
points(sort(z), sort(data$muni.median.p), col="#377eb8", pch=17)
points(sort(z), sort(data$HSIC.opt.p), col="#a65628", pch=3)
points(sort(z), sort(data$Delaunay.p), col="#984ea3", pch=18)
abline(0,1)
legend("topleft", 
       legend = c("Join counts",
                  "HSIC: median",
                  "HSIC: sweep (H)",
                  "Moran's I: median",
                  "Moran's I: sweep (MI)",
                  "Mantel"
       ),
       col = c("#984ea3",
               "#4daf4a",
               "#a65628",
               "#377eb8",
               "#ff7f00",
               "#e41a1c"),
       pch = c(18,15,3,17,1, 16),
       bty = "n",
       cex = 2)

text(x=c(0.8,0.61), y=c(0.4,0.45), pos=4, labels=c('MI','H'), cex=2)

## Generate Figure 3b
data <- read.csv("multi_random/results.txt", header=TRUE, sep="\t")
z <- seq(0,0.999,0.001)
par(mar=c(5,5,1,1), mgp=c(3,.7,0))
qqplot(z, data$mantel.eucl.p, xlim=c(0,1), ylim=c(0,1),  cex=0.75, pch=16,
       xlab="", ylab="p-value", cex.lab=2,  cex.axis=1.3, col="#e41a1c")
points(sort(z), sort(data$HSIC.p), col="#4daf4a", pch=15)
points(sort(z), sort(data$HSIC.opt.p), col="#a65628", pch=3)
points(sort(z), sort(data$Delaunay.p), col="#984ea3", pch=18)
abline(0,1)
legend("topleft", 
       legend = c("Join counts",
                  "HSIC: median",
                  "HSIC: sweep (H)",
                  "Mantel"
       ),
       col = c("#984ea3",
               "#4daf4a",
               "#a65628",
               "#e41a1c"),
       pch = c(18,15,3,16),
       bty = "n",
       cex = 2)

text(x=c(0.55), y=c(0.4), pos=4, labels=c('H'), cex=2)


## Generate Figure 3c
data <- read.csv("freq_random/results.txt", header=TRUE, sep="\t")
z <- seq(0,0.999,0.001)
par(mar=c(5,5,1,1), mgp=c(3,.7,0))
qqplot(z, data$muni.p, xlim=c(0,1), ylim=c(0,1),  cex=0.75, pch=1,
       xlab="", ylab="p-value", cex.lab=2,  cex.axis=1.3, col="#ff7f00")
points(sort(z), sort(data$mantel.eucl.freq.p),  col="#e41a1c", pch=16)
points(sort(z), sort(data$HSIC.freq.p), col="#4daf4a", pch=15)
points(sort(z), sort(data$muni.median.p), col="#377eb8", pch=17)
points(sort(z), sort(data$HSIC.opt.freq.p), col="#a65628", pch=3)
abline(0,1)
legend("topleft", 
       legend = c("HSIC: median",
                  "HSIC: sweep (H)",
                  "Moran's I: median",
                  "Moran's I: sweep (MI)",
                  "Mantel"
       ),
       col = c("#4daf4a",
               "#a65628",
               "#377eb8",
               "#ff7f00",
               "#e41a1c"),
       pch = c(15,3,17,1, 16),
       bty = "n",
       cex = 2)
text(x=c(0.70,0.55), y=c(0.3,0.40), pos=4, labels=c('MI','H'), cex=1.3)

## Figure 4a
data <- read.csv("2cat_angles_bandwidth_obs/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(HSIC.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(HSIC_bandwidth_median_abs_0.25.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.25.p <- tmp$HSIC_bandwidth_median_abs_0.25.p
tmp <- aggregate(HSIC_bandwidth_median_abs_0.5.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.5.p  <- tmp$HSIC_bandwidth_median_abs_0.5.p
tmp <- aggregate(HSIC_bandwidth_median_abs_0.75.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.75.p <- tmp$HSIC_bandwidth_median_abs_0.75.p
tmp <- aggregate(HSIC_bandwidth_median_abs_0.8.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.8.p  <- tmp$HSIC_bandwidth_median_abs_0.8.p 
tmp <- aggregate(HSIC_bandwidth_median_abs_1.25.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_1.25.p  <- tmp$HSIC_bandwidth_median_abs_1.25.p 
tmp <- aggregate(HSIC_bandwidth_median_abs_1.5.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_1.5.p  <- tmp$HSIC_bandwidth_median_abs_1.5.p 
tmp <- aggregate(HSIC_bandwidth_median_abs_2.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_2.p  <- tmp$HSIC_bandwidth_median_abs_2.p 



options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(log10(new_data$Obs.mu), new_data$HSIC.p, col="#377eb8", type="o", pch=17, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_0.5.p, type="o", pch=1,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_0.8.p, type="o", pch=18,  col="#4daf4a", lwd=4, lty=3, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_1.25.p, type="o", pch=16, col="#984ea3", lwd=4, lty=4, cex=2)
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_2.p, type="o", pch=19, col="#000000", lwd=4, lty=5, cex=2)


legend("bottomright", 
       legend = c(expression(paste("HSIC: 0.5 median")),
                  expression(paste("HSIC: 0.8 median")),
                  "HSIC: median",
                  expression(paste("HSIC: 1.25 median")),
                  expression(paste("HSIC: 2 median"))
       ), 
       lty = c(2,3,1,6,7),
       pch = c(1,18,17,16,19),
       col= c("#e41a1c", "#4daf4a",  "#377eb8", "#984ea3", "#000000" ),
       bty = "n", pt.cex = 2, cex=2, lwd=4, seg.len=4)

## Figure 4b
data <- read.csv("2cat_centers_bandwidth_obs/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(HSIC.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(HSIC_bandwidth_median_abs_0.25.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.25.p <- tmp$HSIC_bandwidth_median_abs_0.25.p
tmp <- aggregate(HSIC_bandwidth_median_abs_0.5.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.5.p  <- tmp$HSIC_bandwidth_median_abs_0.5.p
tmp <- aggregate(HSIC_bandwidth_median_abs_0.75.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.75.p <- tmp$HSIC_bandwidth_median_abs_0.75.p
tmp <- aggregate(HSIC_bandwidth_median_abs_0.8.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_0.8.p  <- tmp$HSIC_bandwidth_median_abs_0.8.p 
tmp <- aggregate(HSIC_bandwidth_median_abs_1.25.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_1.25.p  <- tmp$HSIC_bandwidth_median_abs_1.25.p 
tmp <- aggregate(HSIC_bandwidth_median_abs_1.5.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_1.5.p  <- tmp$HSIC_bandwidth_median_abs_1.5.p 
tmp <- aggregate(HSIC_bandwidth_median_abs_2.p ~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC_bandwidth_median_abs_2.p  <- tmp$HSIC_bandwidth_median_abs_2.p 



options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(log10(new_data$Obs.mu), new_data$HSIC.p, col="#377eb8", type="o", pch=17, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_0.5.p, type="o", pch=1,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_0.8.p, type="o", pch=18,  col="#4daf4a", lwd=4, lty=3, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_1.25.p, type="o", pch=16, col="#984ea3", lwd=4, lty=4, cex=2)
lines(log10(new_data$Obs.mu),new_data$HSIC_bandwidth_median_abs_2.p, type="o", pch=19, col="#000000", lwd=4, lty=5, cex=2)


legend("topleft", 
       legend = c(expression(paste("HSIC: 0.5 median")),
                  expression(paste("HSIC: 0.8 median")),
                  "HSIC: median",
                  expression(paste("HSIC: 1.25 median")),
                  expression(paste("HSIC: 2 median"))
       ), 
       lty = c(2,3,1,6,7),
       pch = c(1,18,17,16,19),
       col= c("#e41a1c", "#4daf4a",  "#377eb8", "#984ea3", "#000000" ),
       bty = "n", pt.cex = 2, cex=2, lwd=4, seg.len=4)

## Figure 4c
data <- read.csv("2cat_angles_bandwidth_obs/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(muni.median.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(mI_cutoff_median_abs0.25.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.25.p <- tmp$mI_cutoff_median_abs0.25.p
tmp <- aggregate(mI_cutoff_median_abs0.5.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.5.p <- tmp$mI_cutoff_median_abs0.5.p
tmp <- aggregate(mI_cutoff_median_abs0.75.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.75.p <- tmp$mI_cutoff_median_abs0.75.p
tmp <- aggregate(mI_cutoff_median_abs0.8.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.8.p <- tmp$mI_cutoff_median_abs0.8.p
tmp <- aggregate(mI_cutoff_median_abs1.25.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs1.25.p <- tmp$mI_cutoff_median_abs1.25.p
tmp <- aggregate(mI_cutoff_median_abs1.5.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs1.5.p <- tmp$mI_cutoff_median_abs1.5.p
tmp <- aggregate(mI_cutoff_median_abs2.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs2.p <- tmp$mI_cutoff_median_abs2.p

options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))

plot(log10(new_data$Obs.mu), new_data$muni.median.p, col="#377eb8", type="o", pch=17, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs0.5.p, type="o", pch=1,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs0.8.p, type="o", pch=18,  col="#4daf4a", lwd=4, lty=3, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs1.25.p, type="o", pch=16, col="#984ea3", lwd=4, lty=4, cex=2)
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs2.p, type="o", pch=19, col="#000000", lwd=4, lty=5, cex=2)



legend("bottomright", 
       legend = c(expression(paste("Moran's I: 0.5 median")),
                  expression(paste("Moran's I: 0.8 median")),
                  "Moran's I: median",
                  expression(paste("Moran's I: 1.25 median")),
                  expression(paste("Moran's I: 2 median"))
       ), 
       lty = c(2,3,1,6,7),
       pch = c(1,18,17,16,19),
       col= c("#e41a1c", "#4daf4a",  "#377eb8", "#984ea3", "#000000" ),
       bty = "n", pt.cex = 2, cex=2, lwd=4, seg.len=4)

## Figure 4d
data <- read.csv("2cat_centers_bandwidth_obs/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(muni.median.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(mI_cutoff_median_abs0.25.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.25.p <- tmp$mI_cutoff_median_abs0.25.p
tmp <- aggregate(mI_cutoff_median_abs0.5.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.5.p <- tmp$mI_cutoff_median_abs0.5.p
tmp <- aggregate(mI_cutoff_median_abs0.75.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.75.p <- tmp$mI_cutoff_median_abs0.75.p
tmp <- aggregate(mI_cutoff_median_abs0.8.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs0.8.p <- tmp$mI_cutoff_median_abs0.8.p
tmp <- aggregate(mI_cutoff_median_abs1.25.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs1.25.p <- tmp$mI_cutoff_median_abs1.25.p
tmp <- aggregate(mI_cutoff_median_abs1.5.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs1.5.p <- tmp$mI_cutoff_median_abs1.5.p
tmp <- aggregate(mI_cutoff_median_abs2.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mI_cutoff_median_abs2.p <- tmp$mI_cutoff_median_abs2.p

options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))

plot(log10(new_data$Obs.mu), new_data$muni.median.p, col="#377eb8", type="o", pch=17, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs0.5.p, type="o", pch=1,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs0.8.p, type="o", pch=18,  col="#4daf4a", lwd=4, lty=3, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs1.25.p, type="o", pch=16, col="#984ea3", lwd=4, lty=4, cex=2)
lines(log10(new_data$Obs.mu),new_data$mI_cutoff_median_abs2.p, type="o", pch=19, col="#000000", lwd=4, lty=5, cex=2)



legend("topleft", 
       legend = c(expression(paste("Moran's I: 0.5 median")),
                  expression(paste("Moran's I: 0.8 median")),
                  "Moran's I: median",
                  expression(paste("Moran's I: 1.25 median")),
                  expression(paste("Moran's I: 2 median"))
       ), 
       lty = c(2,3,1,6,7),
       pch = c(1,18,17,16,19),
       col= c("#e41a1c", "#4daf4a",  "#377eb8", "#984ea3", "#000000" ),
       bty = "n", pt.cex = 2, cex=2, lwd=4, seg.len=4)



## Figure 5
data <- read.csv("2cat_angle_example/results.txt", head=TRUE, sep="\t")
data$angle_new <- floor(data$angle/20) * 20 + 10

new_data <- aggregate(muni.median.p~angle_new, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(HSIC.p~angle_new, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC.p <- tmp$HSIC.p
tmp <- aggregate(Delaunay.p~angle_new, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$Delaunay.p <- tmp$Delaunay.p
tmp <- aggregate(mantel.eucl.p~angle_new, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.p <- tmp$mantel.eucl.p

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(new_data$angle_new, new_data$HSIC.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), cex=0.75, 
     xlab="Angle", ylab="power", cex.lab=2,  cex.axis=1.3, lwd=5) 
lines(new_data$angle_new,new_data$mantel.eucl.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(new_data$angle_new,new_data$muni.median.p, type="o", pch=17,  col="#377eb8", lwd=4, lty=3, cex=1.7)
lines(new_data$angle_new,new_data$Delaunay.p, type="o", pch=18, col="#984ea3", lwd=4, lty=4, cex=2)



legend("bottomleft", 
       legend = c("Mantel",
                  "HSIC",
                  "Moran's I",
                  "Join counts"
       ),
       col = c("#e41a1c",
               "#4daf4a",
               "#377eb8",
               "#984ea3"),
       pch= c(16, 15, 17, 18),
       bty = "n",
       cex = 2)


## Figure 6a
data <- read.csv("freq_angles_noise/results.txt", header=TRUE, sep="\t")

new_data <- aggregate(HSIC.freq.p~Noise, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(muni.median.p~Noise, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$muni.median.p <- tmp$muni.median.p
tmp <- aggregate(mantel.eucl.freq.p~Noise, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.freq.p <- tmp$mantel.eucl.freq.p

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(new_data$Noise, new_data$HSIC.freq.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab="Outlier probability", ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(new_data$Noise,new_data$mantel.eucl.freq.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(new_data$Noise,new_data$muni.median.p, type="o", pch=17,  col="#377eb8", lwd=4, lty=3, cex=1.7)

legend("bottomleft", 
       legend = c(
         "HSIC",
         "Mantel", 
         "Moran's I"
       ), 
       lty = c(1,2,3),
       col = c("#4daf4a", "#e41a1c", "#377eb8"),
       pch = c(15,16,17),
       bty = "n", pt.cex = 2, cex=2, lwd=4)

## Figure 6b
data <- read.csv("freq_centers_noise/results.txt", header=TRUE, sep="\t")

new_data <- aggregate(HSIC.freq.p~Noise, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(muni.median.p~Noise, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$muni.median.p <- tmp$muni.median.p
tmp <- aggregate(mantel.eucl.freq.p~Noise, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.freq.p <- tmp$mantel.eucl.freq.p

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(new_data$Noise, new_data$HSIC.freq.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab="Outlier probability", ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(new_data$Noise,new_data$mantel.eucl.freq.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(new_data$Noise,new_data$muni.median.p, type="o", pch=17,  col="#377eb8", lwd=4, lty=3, cex=1.7)

legend("topright", 
       legend = c(
         "HSIC",
         "Mantel", 
         "Moran's I"
       ), 
       lty = c(1,2,3),
       col = c("#4daf4a", "#e41a1c", "#377eb8"),
       pch = c(15,16,17),
       bty = "n", pt.cex = 2, cex=2, lwd=4)

# Figure 7a
data <- read.csv("2cat_angles_obs/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(muni.median.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(HSIC.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC.p <- tmp$HSIC.p
tmp <- aggregate(mantel.eucl.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.p <- tmp$mantel.eucl.p
tmp <- aggregate(Delaunay.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$Delaunay.p <- tmp$Delaunay.p

options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(log10(new_data$Obs.mu), new_data$HSIC.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$mantel.eucl.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$muni.median.p, type="o", pch=17,  col="#377eb8", lwd=4, lty=3, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$Delaunay.p, type="o", pch=18, col="#984ea3", lwd=4, lty=4, cex=2)

legend("bottomright", 
       legend = c("Join counts",
                  "HSIC",
                  "Mantel", 
                  "Moran's I"
       ), 
       lty = c(4,1,2,3),
       col = c("#984ea3", "#4daf4a", "#e41a1c", "#377eb8"),
       pch = c(18,15,16,17),
       bty = "n", pt.cex = 2, cex=2, lwd=4)


## Figure 7b
data <- read.csv("multi_angles_obs/results.txt", header=TRUE, sep="\t")
new_data <-aggregate(HSIC.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC.p <- tmp$HSIC.p
tmp <- aggregate(mantel.eucl.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.p <- tmp$mantel.eucl.p
tmp <- aggregate(Delaunay.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$Delaunay.p <- tmp$Delaunay.p

options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(log10(new_data$Obs.mu), new_data$HSIC.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$mantel.eucl.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$Delaunay.p, type="o", pch=18, col="#984ea3", lwd=4, lty=4, cex=2)

legend("bottomright", 
       legend = c("Join counts",
                  "HSIC",
                  "Mantel"
       ), 
       lty = c(4,1,2),
       col = c("#984ea3", "#4daf4a", "#e41a1c"),
       pch = c(18,15,16),
       bty = "n", pt.cex = 2, cex=2, lwd=4)

## Figure 7c
data <- read.csv("freq_angles_sigma/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(HSIC.freq.p~Sigma, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(muni.median.p~Sigma, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$muni.median.p <- tmp$muni.median.p
tmp <- aggregate(mantel.eucl.freq.p~Sigma, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.freq.p <- tmp$mantel.eucl.freq.p

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(new_data$Sigma, new_data$HSIC.freq.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab=expression(sigma), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(new_data$Sigma,new_data$mantel.eucl.freq.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(new_data$Sigma,new_data$muni.median.p, type="o", pch=17,  col="#377eb8", lwd=4, lty=3, cex=1.7)

legend("bottomleft", 
       legend = c(
         "HSIC",
         "Mantel", 
         "Moran's I"
       ), 
       lty = c(1,2,3),
       col = c("#4daf4a", "#e41a1c", "#377eb8"),
       pch = c(15,16,17),
       bty = "n", pt.cex = 2, cex=2, lwd=4)

## Figure 8a
data <- read.csv("2cat_centers_obs/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(muni.median.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(HSIC.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC.p <- tmp$HSIC.p
tmp <- aggregate(mantel.eucl.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.p <- tmp$mantel.eucl.p
tmp <- aggregate(Delaunay.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$Delaunay.p <- tmp$Delaunay.p

options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(log10(new_data$Obs.mu), new_data$HSIC.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$mantel.eucl.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$muni.median.p, type="o", pch=17,  col="#377eb8", lwd=4, lty=3, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$Delaunay.p, type="o", pch=18, col="#984ea3", lwd=4, lty=4, cex=2)

legend("topleft", 
       legend = c("Join counts",
                  "HSIC",
                  "Mantel", 
                  "Moran's I"
       ), 
       lty = c(4,1,2,3),
       col = c("#984ea3", "#4daf4a", "#e41a1c", "#377eb8"),
       pch = c(18,15,16,17),
       bty = "n", pt.cex = 2, cex=2, lwd=4)

## Figure 8b
data <- read.csv("multi_centers_obs/results.txt", header=TRUE, sep="\t")
new_data <-aggregate(HSIC.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$HSIC.p <- tmp$HSIC.p
tmp <- aggregate(mantel.eucl.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.p <- tmp$mantel.eucl.p
tmp <- aggregate(Delaunay.p~Obs.mu, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$Delaunay.p <- tmp$Delaunay.p

options(scipen=999)

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(log10(new_data$Obs.mu), new_data$HSIC.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab=expression(log(mu[obs])), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(log10(new_data$Obs.mu),new_data$mantel.eucl.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(log10(new_data$Obs.mu),new_data$Delaunay.p, type="o", pch=18, col="#984ea3", lwd=4, lty=4, cex=2)

legend("topleft", 
       legend = c("Join counts",
                  "HSIC",
                  "Mantel"
       ), 
       lty = c(4,1,2),
       col = c("#984ea3", "#4daf4a", "#e41a1c"),
       pch = c(18,15,16),
       bty = "n", pt.cex = 2, cex=2, lwd=4)

## Figure 8c
data <- read.csv("freq_centers_sigma/results.txt", header=TRUE, sep="\t")
new_data <- aggregate(HSIC.freq.p~Sigma, data=data, FUN= function(x) sum(x<0.05)/length(x))
tmp <- aggregate(muni.median.p~Sigma, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$muni.median.p <- tmp$muni.median.p
tmp <- aggregate(mantel.eucl.freq.p~Sigma, data=data, FUN= function(x) sum(x<0.05)/length(x))
new_data$mantel.eucl.freq.p <- tmp$mantel.eucl.freq.p

par(mar=c(5,5,1,1), mgp=c(3,.7,0))
plot(new_data$Sigma, new_data$HSIC.freq.p, col="#4daf4a", type="o", pch=15, ylim=c(0,1), 
     xlab=expression(sigma), ylab="power", cex.lab=3,  cex.axis=1.3, lwd=4, lty=1, cex=1.7) 
lines(new_data$Sigma,new_data$mantel.eucl.freq.p, type="o", pch=16,   col="#e41a1c", lwd=4, lty=2, cex=1.7)
lines(new_data$Sigma,new_data$muni.median.p, type="o", pch=17,  col="#377eb8", lwd=4, lty=3, cex=1.7)

legend("topright", 
       legend = c(
         "HSIC",
         "Mantel", 
         "Moran's I"
       ), 
       lty = c(1,2,3),
       col = c("#4daf4a", "#e41a1c", "#377eb8"),
       pch = c(15,16,17),
       bty = "n", pt.cex = 2, cex=2, lwd=4)