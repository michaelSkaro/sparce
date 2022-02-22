if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")

library(karyoploteR)


# plot the values of the CPG islands for the human genome
library(AnnotationHub)
ahub <- AnnotationHub()
ahub["AH5086"]

# All cpg islands
library(karyoploteR)
kp <- plotKaryotype()
kpPlotRegions(kp, data=cpgs)

# Density of CpG islands
kp <- plotKaryotype()
kpPlotRegions(kp, data=cpgs, r0=0, r1=0.5)
kpPlotDensity(kp, data=cpgs, r0=0.5, r1=1)


# window sizes in chromosome 12

pp <- getDefaultPlotParams(plot.type = 2)
pp$data1height <- 50
kp <- plotKaryotype(chromosomes="chr10", plot.type=2, plot.params = pp)
kpPlotRegions(kp, data=cpgs, col="#AAAAAA", border="#AAAAAA")
kpPlotDensity(kp, data=cpgs, data.panel=2, col="#AA88FF", r0=0.5, r1=1)
kpPlotDensity(kp, data=cpgs, data.panel=2, col="#FF88AA", window.size = 100000, r0=0.5, r1=0)

# variable window sizes for chromosome 10

pp <- getDefaultPlotParams(plot.type = 2)
pp$data1height <- 50
kp <- plotKaryotype(chromosomes="chr10", plot.type=1)
kpPlotRegions(kp, data=cpgs, col="#CCCCCC44", border="#CCCCCC44")
kpPlotDensity(kp, data=cpgs, data.panel=1, col="#8844FF", window.size= 1000000, r0=0, r1=1)
kpPlotDensity(kp, data=cpgs, data.panel=1, col="#AA66FF", window.size = 500000, r0=0.25, r1=1)
kpPlotDensity(kp, data=cpgs, data.panel=1, col="#CC88FF", window.size = 200000, r0=0.5, r1=1)
kpPlotDensity(kp, data=cpgs, data.panel=1, col="#EEAAFF", window.size = 100000, r0=0.75, r1=1)
