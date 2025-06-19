# create lattice dotplots to illustrate where missmatches or false matches are along chromosome
# command line:
# R CMD BATCH --no-save --no-restore '--args filename' create_lattice_dotplots.R

# NOTE: This script is currently run as part of a container, not from the bin folder
# Keep here for reference, and eventually for extraction from the container

args <- commandArgs(TRUE)

pod_table <- read.table(args[1])

plot_dir = args[2]

name = args[3]

pdf(paste(plot_dir,'/', name,'_POD_karyotype.pdf', sep=""), width=20, heigh=10)

library(karyoploteR)

pod_table_F <- pod_table[which(pod_table$V3 == 'F-AAB' | pod_table$V3 == 'F-ABB'),]
pod_table_M <- pod_table[which(pod_table$V3 == 'M-AAB' | pod_table$V3 == 'M-ABB'),]

pod_table_F1 <- pod_table_F[which(pod_table_F$V1 == '1'),]
pod_table_F2 <- pod_table_F[which(pod_table_F$V1 == '2'),]
pod_table_F3 <- pod_table_F[which(pod_table_F$V1 == '3'),]
pod_table_F4 <- pod_table_F[which(pod_table_F$V1 == '4'),]
pod_table_F5 <- pod_table_F[which(pod_table_F$V1 == '5'),]
pod_table_F6 <- pod_table_F[which(pod_table_F$V1 == '6'),]
pod_table_F7 <- pod_table_F[which(pod_table_F$V1 == '7'),]
pod_table_F8 <- pod_table_F[which(pod_table_F$V1 == '8'),]
pod_table_F9 <- pod_table_F[which(pod_table_F$V1 == '9'),]
pod_table_F10 <- pod_table_F[which(pod_table_F$V1 == '10'),]
pod_table_F11 <- pod_table_F[which(pod_table_F$V1 == '11'),]
pod_table_F12 <- pod_table_F[which(pod_table_F$V1 == '12'),]
pod_table_F13 <- pod_table_F[which(pod_table_F$V1 == '13'),]
pod_table_F14 <- pod_table_F[which(pod_table_F$V1 == '14'),]
pod_table_F15 <- pod_table_F[which(pod_table_F$V1 == '15'),]
pod_table_F16 <- pod_table_F[which(pod_table_F$V1 == '16'),]
pod_table_F17 <- pod_table_F[which(pod_table_F$V1 == '17'),]
pod_table_F18 <- pod_table_F[which(pod_table_F$V1 == '18'),]
pod_table_F19 <- pod_table_F[which(pod_table_F$V1 == '19'),]
pod_table_F20 <- pod_table_F[which(pod_table_F$V1 == '20'),]
pod_table_F21 <- pod_table_F[which(pod_table_F$V1 == '21'),]
pod_table_F22 <- pod_table_F[which(pod_table_F$V1 == '22'),]
pod_table_FX <- pod_table_F[which(pod_table_F$V1 == 'X'),]
pod_table_FY <- pod_table_F[which(pod_table_F$V1 == 'Y'),]

pod_table_M1 <- pod_table_M[which(pod_table_M$V1 == '1'),]
pod_table_M2 <- pod_table_M[which(pod_table_M$V1 == '2'),]
pod_table_M3 <- pod_table_M[which(pod_table_M$V1 == '3'),]
pod_table_M4 <- pod_table_M[which(pod_table_M$V1 == '4'),]
pod_table_M5 <- pod_table_M[which(pod_table_M$V1 == '5'),]
pod_table_M6 <- pod_table_M[which(pod_table_M$V1 == '6'),]
pod_table_M7 <- pod_table_M[which(pod_table_M$V1 == '7'),]
pod_table_M8 <- pod_table_M[which(pod_table_M$V1 == '8'),]
pod_table_M9 <- pod_table_M[which(pod_table_M$V1 == '9'),]
pod_table_M10 <- pod_table_M[which(pod_table_M$V1 == '10'),]
pod_table_M11 <- pod_table_M[which(pod_table_M$V1 == '11'),]
pod_table_M12 <- pod_table_M[which(pod_table_M$V1 == '12'),]
pod_table_M13 <- pod_table_M[which(pod_table_M$V1 == '13'),]
pod_table_M14 <- pod_table_M[which(pod_table_M$V1 == '14'),]
pod_table_M15 <- pod_table_M[which(pod_table_M$V1 == '15'),]
pod_table_M16 <- pod_table_M[which(pod_table_M$V1 == '16'),]
pod_table_M17 <- pod_table_M[which(pod_table_M$V1 == '17'),]
pod_table_M18 <- pod_table_M[which(pod_table_M$V1 == '18'),]
pod_table_M19 <- pod_table_M[which(pod_table_M$V1 == '19'),]
pod_table_M20 <- pod_table_M[which(pod_table_M$V1 == '20'),]
pod_table_M21 <- pod_table_M[which(pod_table_M$V1 == '21'),]
pod_table_M22 <- pod_table_M[which(pod_table_M$V1 == '22'),]
pod_table_MX <- pod_table_M[which(pod_table_M$V1 == 'X'),]
pod_table_MY <- pod_table_M[which(pod_table_M$V1 == 'Y'),]


pp <- getDefaultPlotParams(plot.type=2) 
pp$ideogramheight <- 150
pp$data1outmargin <- 50
pp$data2outmargin <- 50
pp$topmargin <- 20
pp$bottommargin <- 20
pp$data1height <- 250
pp$data2height <- 250


kp <- plotKaryotype(genome="hg19", plot.type=2, plot.params=pp)

kpAddCytobandLabels(kp, cex=0.45)

kpDataBackground(kp, data.panel = 1, col="#FFAACB")
kpDataBackground(kp, data.panel = 2, col="#AACBFF")

kpAxis(kp, ymin=-2.32, ymax=2.32, numticks=2, tick.pos = c(-1.5,1.5), labels = c("Mat ABB","Mat AAB"), cex=0.45, data.panel=1)
kpAxis(kp, ymin=2.32, ymax=-2.32, numticks=2, tick.pos = c(-1.5,1.5), labels = c("Pat ABB","Pat AAB"), cex=0.45, data.panel=2)

kpPoints(kp, chr="chr1", ymin=-2.32, ymax=2.32, x=pod_table_M1$V2, y=log2(pod_table_M1$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr1", ymin=2.32, ymax=-2.32, x=pod_table_F1$V2, y=log2(pod_table_F1$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr2", ymin=-2.32, ymax=2.32, x=pod_table_M2$V2, y=log2(pod_table_M2$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr2", ymin=2.32, ymax=-2.32, x=pod_table_F2$V2, y=log2(pod_table_F2$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr3", ymin=-2.32, ymax=2.32, x=pod_table_M3$V2, y=log2(pod_table_M3$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr3", ymin=2.32, ymax=-2.32, x=pod_table_F3$V2, y=log2(pod_table_F3$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr4", ymin=-2.32, ymax=2.32, x=pod_table_M4$V2, y=log2(pod_table_M4$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr4", ymin=2.32, ymax=-2.32, x=pod_table_F4$V2, y=log2(pod_table_F4$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr5", ymin=-2.32, ymax=2.32, x=pod_table_M5$V2, y=log2(pod_table_M5$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr5", ymin=2.32, ymax=-2.32, x=pod_table_F5$V2, y=log2(pod_table_F5$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr6", ymin=-2.32, ymax=2.32, x=pod_table_M6$V2, y=log2(pod_table_M6$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr6", ymin=2.32, ymax=-2.32, x=pod_table_F6$V2, y=log2(pod_table_F6$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr7", ymin=-2.32, ymax=2.32, x=pod_table_M7$V2, y=log2(pod_table_M7$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr7", ymin=2.32, ymax=-2.32, x=pod_table_F7$V2, y=log2(pod_table_F7$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr8", ymin=-2.32, ymax=2.32, x=pod_table_M8$V2, y=log2(pod_table_M8$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr8", ymin=2.32, ymax=-2.32, x=pod_table_F8$V2, y=log2(pod_table_F8$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr9", ymin=-2.32, ymax=2.32, x=pod_table_M9$V2, y=log2(pod_table_M9$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr9", ymin=2.32, ymax=-2.32, x=pod_table_F9$V2, y=log2(pod_table_F9$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr10", ymin=-2.32, ymax=2.32, x=pod_table_M10$V2, y=log2(pod_table_M10$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr10", ymin=2.32, ymax=-2.32, x=pod_table_F10$V2, y=log2(pod_table_F10$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr11", ymin=-2.32, ymax=2.32, x=pod_table_M11$V2, y=log2(pod_table_M11$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr11", ymin=2.32, ymax=-2.32, x=pod_table_F11$V2, y=log2(pod_table_F11$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr12", ymin=-2.32, ymax=2.32, x=pod_table_M12$V2, y=log2(pod_table_M12$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr12", ymin=2.32, ymax=-2.32, x=pod_table_F12$V2, y=log2(pod_table_F12$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr13", ymin=-2.32, ymax=2.32, x=pod_table_M13$V2, y=log2(pod_table_M13$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr13", ymin=2.32, ymax=-2.32, x=pod_table_F13$V2, y=log2(pod_table_F13$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr14", ymin=-2.32, ymax=2.32, x=pod_table_M14$V2, y=log2(pod_table_M14$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr14", ymin=2.32, ymax=-2.32, x=pod_table_F14$V2, y=log2(pod_table_F14$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr15", ymin=-2.32, ymax=2.32, x=pod_table_M15$V2, y=log2(pod_table_M15$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr15", ymin=2.32, ymax=-2.32, x=pod_table_F15$V2, y=log2(pod_table_F15$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr16", ymin=-2.32, ymax=2.32, x=pod_table_M16$V2, y=log2(pod_table_M16$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr16", ymin=2.32, ymax=-2.32, x=pod_table_F16$V2, y=log2(pod_table_F16$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr17", ymin=-2.32, ymax=2.32, x=pod_table_M17$V2, y=log2(pod_table_M17$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr17", ymin=2.32, ymax=-2.32, x=pod_table_F17$V2, y=log2(pod_table_F17$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr18", ymin=-2.32, ymax=2.32, x=pod_table_M18$V2, y=log2(pod_table_M18$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr18", ymin=2.32, ymax=-2.32, x=pod_table_F18$V2, y=log2(pod_table_F18$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr19", ymin=-2.32, ymax=2.32, x=pod_table_M19$V2, y=log2(pod_table_M19$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr19", ymin=2.32, ymax=-2.32, x=pod_table_F19$V2, y=log2(pod_table_F19$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr20", ymin=-2.32, ymax=2.32, x=pod_table_M20$V2, y=log2(pod_table_M20$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr20", ymin=2.32, ymax=-2.32, x=pod_table_F20$V2, y=log2(pod_table_F20$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr21", ymin=-2.32, ymax=2.32, x=pod_table_M21$V2, y=log2(pod_table_M21$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr21", ymin=2.32, ymax=-2.32, x=pod_table_F21$V2, y=log2(pod_table_F21$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chr22", ymin=-2.32, ymax=2.32, x=pod_table_M22$V2, y=log2(pod_table_M22$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chr22", ymin=2.32, ymax=-2.32, x=pod_table_F22$V2, y=log2(pod_table_F22$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chrX", ymin=-2.32, ymax=2.32, x=pod_table_MX$V2, y=log2(pod_table_MX$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chrX", ymin=2.32, ymax=-2.32, x=pod_table_FX$V2, y=log2(pod_table_FX$V4), cex=0.25, data.panel = 2)

kpPoints(kp, chr="chrY", ymin=-2.32, ymax=2.32, x=pod_table_MY$V2, y=log2(pod_table_MY$V4), cex=0.25, data.panel = 1)
kpPoints(kp, chr="chrY", ymin=2.32, ymax=-2.32, x=pod_table_FY$V2, y=log2(pod_table_FY$V4), cex=0.25, data.panel = 2)

dev.off()
