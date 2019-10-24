#! /usr/bin/env Rscript

## MI flow cytometry processing
##########################################
### Panel 2 - Regulatory T cell subsets ##
##########################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 2C")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P2.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P02")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P2.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P2.flow.data) <- sample.ids

# subset samples into group C
P2.flow.data <- P2.flow.data[529:792]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P2.flow.data))){
  events <- nrow(P2.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P2.fs_filt <- P2.flow.data[seq(along=P2.flow.data)]
if(length(extract)){P2.fs_filt <- P2.fs_filt[-(extract),]}

# remove channels that aren't interesting
param.desc <-  pData(parameters(P2.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P2.flow.data[[1]])
param.desc[is.na(param.desc)] <- names(param.desc)[is.na(param.desc)]
param.meta <- do.call(cbind.data.frame,
                      list("Param"=names(param.desc),
                           "Marker"=param.desc))

# transform the data
arcTrans <- arcsinhTransform()
logTrans <- logTransform()
biexpTrans <- biexponentialTransform()
scaleTrans <- scaleTransform()
logicleTrans <- logicleTransform()
linearTrans <- linearTransform()

# transform each paramter
transform1 <- transformList("V1-A", biexpTrans)
transform2 <- transformList("V2-A", biexpTrans)
transform3 <- transformList("B1-A", biexpTrans)
transform4 <- transformList("B2-A", biexpTrans)
transform5 <- transformList("B3-A", biexpTrans)
transform6 <- transformList("B4-A", biexpTrans)
transform7 <- transformList("R1-A", biexpTrans)
transform8 <- transformList("R2-A", biexpTrans)
transformForwardH <- transformList("FSC-H", biexpTrans)
transformSideA <- transformList("SSC-A", biexpTrans)
transformForwardA <- transformList("FSC-A", biexpTrans)

P2.fs_filt <- transform(P2.fs_filt, transform1)
P2.fs_filt <- transform(P2.fs_filt, transform2)
P2.fs_filt <- transform(P2.fs_filt, transform3)
P2.fs_filt <- transform(P2.fs_filt, transform4)
P2.fs_filt <- transform(P2.fs_filt, transform5)
P2.fs_filt <- transform(P2.fs_filt, transform6)
P2.fs_filt <- transform(P2.fs_filt, transform7)
P2.fs_filt <- transform(P2.fs_filt, transform8)

# gate on leukocytes from FSC and SSC
# ggcyto(P2.fs_filt[1:6], aes(x=`FSC-A`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 1e5))

# set up a CD3 gate
# setup side scatter and CD3 gates
ssc.points <- c(0, 5e4, 6e4, 7.5e4, 7.5e4, 0)
fsc.points <- c(0, 1e3, 1e5, 1.5e5, 2.75e5, 2.75e5)

leuko.gate <- polygonGate(filterId="Leukocytes", .gate=cbind("SSC-A"=ssc.points, "FSC-A"=fsc.points))
leuko.filter <- filter(P2.fs_filt, leuko.gate)
leuko.gates <- lapply(leuko.filter, 
                      function(res) filterDetails(res, "Leukocytes")$filter)

# subset all CD3+ cells before all other analyses
P2.leuko <- Subset(P2.fs_filt, leuko.gate)

# ggcyto(P2.fs_filt[1:6], aes(x=`FSC-A`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(leuko.gates)

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P2.leuko)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("B1-A", "B4-A")
third.warp <- c("R1-A")
P2.warping <- warpSet(P2.leuko, stains=setdiff(pars, c(second.warp, third.warp)))
P2.warping <- warpSet(P2.warping, stains=second.warp, peakNr=2)
P2.warping <- warpSet(P2.warping, stains=third.warp, bwFac = 1, nbreaks=20, peakNr=2)

densityplot(x=~`R1-A`,
            data=P2.warping,
            channel=c("R1-A"))

# set up gating scheme using a GatingSet object
P2.gs <- GatingSet(P2.warping)

rm(list=c("P2.warping", "P2.flow.data", "P2.fs_filt",
          "P2.leuko"))
gc()

# set up a singlets gate
# setup forward and side scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))

flowWorkspace::add(P2.gs, singlet.gate, parent="root")
recompute(P2.gs)

# forward and side scater gate
ssc.cov <- matrix(c(2e9, 1e6, 1.5e9, 1e10), ncol=2,
                  dimnames=list(c("SSC-A", "FSC-H"), c("SSC-A", "FSC-H")))
ssc.mean <- c("SSC-A"=2e4, "FSC-H"=2.5e4)
ssc.ellipsoid <- ellipsoidGate(filterId="SSC.Ellipse", .gate=ssc.cov, mean=ssc.mean)

flowWorkspace::add(P2.gs, ssc.ellipsoid, parent="Singlet")
recompute(P2.gs)

# set up a viability gate
# set up a singlets gate
# setup forward and side scatter gates
side.points <- c(0, 7.5e4, 7.5e4, 0, 0)
viability.points <- c(0, 0, 9, 9, 0)
viability.gate <- polygonGate(filterId="Live", .gate=cbind("SSC-A"=side.points, "V2-A"=viability.points))

flowWorkspace::add(P2.gs, viability.gate, parent="SSC.Ellipse")
recompute(P2.gs)

# # gate on singlets
# ggcyto(P2.gs[1:6], aes(x=`SSC-A`, y=`V2-A`), subset="SSC.Ellipse") +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 12)) +
#   geom_gate(viability.gate)

# plot CD4 Vs CD8 to get major T cell subsets
#ggcyto(P2.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="Live") +
#  geom_hex(bins=256) +
#  lims(x=c(0, 12), y=c(0, 12))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7.5, 12, 12, 7.5, 7.5)
cd4.points <- c(0, 0, 7, 7, 0)
cd8.gate <- polygonGate(filterId="CD8", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(8, 8, 0, 0, 8)
cd4.points <- c(7.5, 11, 11, 7.5, 7.5)
cd4.gate <- polygonGate(filterId="CD4", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(0, 7.5, 7.5, 0, 0)
cd4.points <- c(0, 0, 7, 7, 0)
DN.gate <- polygonGate(filterId="DN", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(8, 12, 12, 8, 8)
cd4.points <- c(7, 7, 11, 11, 7)
DP.gate <- polygonGate(filterId="DP", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# # # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P2.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="Live") +
#   geom_hex(bins=256) +
#   labs(x="CD4", y="CD8") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cd8.gate) +
#   geom_gate(cd4.gate) +
#   geom_gate(DN.gate) +
#   geom_gate(DP.gate)

# add the T cell gates
flowWorkspace::add(P2.gs, cd4.gate, parent="Live")
recompute(P2.gs)

flowWorkspace::add(P2.gs, cd8.gate, parent="Live")
recompute(P2.gs)

flowWorkspace::add(P2.gs, DN.gate, parent="Live")
recompute(P2.gs)

flowWorkspace::add(P2.gs, DP.gate, parent="Live")
recompute(P2.gs)

# T cell groups defined by this population:
# naive T regs
# memory T regs
# activated T regs

## Select the CD127- CD25+ cells T regs first
cd127.points <- c(4, 4, 4, 5, 8, 7.5, 7.5, 6.75, 6.5, 6, 4)
cd25.points <- c(6, 7.75, 10, 11, 11, 8.75, 8.25, 8.25, 7.75, 6.75, 6)
treg.gate <- polygonGate(filterId="Treg", .gate=cbind("B2-A"=cd25.points, "R1-A"=cd127.points))

flowWorkspace::add(P2.gs, treg.gate, parent="CD4")
recompute(P2.gs)

flowWorkspace::add(P2.gs, treg.gate, parent="CD8")
recompute(P2.gs)

#######################
# Naive HLADR-CD45RA+ #
#######################
# naive T reg cells express CD45RA but lack HLA-DR
cd45ra.points <- c(7.5, 7.5, 12, 12, 7.5)
hladr.points <- c(0, 7, 7, 0, 0)
tregnaive.gate <- polygonGate(filterId="naiveTreg", .gate=cbind("B3-A"=hladr.points, "B1-A"=cd45ra.points))

flowWorkspace::add(P2.gs, tregnaive.gate, parent="CD4/Treg")
recompute(P2.gs)

flowWorkspace::add(P2.gs, tregnaive.gate, parent="CD8/Treg")
recompute(P2.gs)

###########################
# activated HLADR+CD45RA- #
###########################
# central memmory T cells express HLA-DR but lack CD45RA
cd45ra.points <- c(2.5, 2.5, 7.5, 7.5, 2.5)
hladr.points <- c(7, 11, 11, 7, 7)
tregact.gate <- polygonGate(filterId="ActivatedTreg", .gate=cbind("B3-A"=hladr.points, "B1-A"=cd45ra.points))

flowWorkspace::add(P2.gs, tregact.gate, parent="CD4/Treg")
recompute(P2.gs)

flowWorkspace::add(P2.gs, tregact.gate, parent="CD8/Treg")
recompute(P2.gs)

#######################
# memory HLADR-CD45RA- #
#######################
#  memmory T reg cells lack HLADR and CD45RA
cd45ra.points <- c(2.5, 2.5, 7.5, 7.5, 2.5)
hladr.points <- c(0, 7, 7, 0, 0)
tregem.gate <- polygonGate(filterId="memTreg", .gate=cbind("B3-A"=hladr.points, "B1-A"=cd45ra.points))

flowWorkspace::add(P2.gs, tregem.gate, parent="CD4/Treg")
recompute(P2.gs)

flowWorkspace::add(P2.gs, tregem.gate, parent="CD8/Treg")
recompute(P2.gs)

# ggcyto(P2.gs[1:6], aes(x=`B3-A`, y=`B1-A`), subset="CD4/Treg") +
#   geom_hex(bins=128) +
#   labs(x="HLA-DR", y="CD45RA") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(tregnaive.gate) +
#   geom_gate(tregact.gate) +
#   geom_gate(tregem.gate)

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P2.gs)[7:length(getNodes(P2.gs))]
samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P2.gs, y=node)
  samp.names <- sampleNames(node.flow)
  node.list <- list()
  for(x in seq_along(samp.names)){
    samp <- samp.names[x]
    samp.flow <- node.flow[[samp]]
    samp.mat <- exprs(samp.flow)
    if(nrow(samp.mat) > 1){
      samp.mat <- samp.mat[, c("FSC-A", "SSC-A", "FSC-H", pars)]
      cell.prefix <- paste(samp, gsub(node, pattern="/", replacement="_"), sep=".")
      cell.names <- paste(cell.prefix, paste0("Cell", 1:nrow(samp.mat)), sep=".")
      rownames(samp.mat) <- cell.names
      node.list[[samp]] <- samp.mat
    } else if(nrow(samp.mat) < 1 & nrow(samp.mat) > 0){
      samp.mat <- as.matrix(samp.mat[, c("FSC-A", "SSC-A", "FSC-H", pars)], ncol=length(pars))
      cell.prefix <- paste(samp, gsub(node, pattern="/", replacement="_"), sep=".")
      cell.names <- paste(cell.prefix, paste0("Cell", 1:nrow(samp.mat)), sep=".")
      rownames(samp.mat) <- cell.names
      node.list[[samp]] <- samp.mat
    }
  }
  samp.list[[node]] <- do.call(rbind.data.frame,
                               node.list)
}

# turn this into one massive data set?
P2.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P2.df <- P2.df[!duplicated(P2.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "ICOS", "Viability", "CD45RA",
                     "CD25", "HLADR", "CD8B", "CD127", "CD4")
P2.df <- P2.df[, parm.goi]
colnames(P2.df) <- names(parm.goi)

celltype <- unlist(lapply(strsplit(rownames(P2.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P2.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P2.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P2.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P2.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
individuals <- levels(P2.meta$Individual)
celltypes <- levels(P2.meta$CellType)

cell.summary.list <- list()

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P2.df[grepl(rownames(P2.df), pattern=celltype), ]
  indiv.summary.list <- list()
  
  for(j in seq_along(individuals)){
    indiv <- individuals[j]
    indiv.df <- celltype.df[grepl(rownames(celltype.df), pattern=indiv), ]
    
    if(nrow(indiv.df) > 0){
      # truncate negative values at 0
      indiv.df[indiv.df < 0] <- 0
      
      # normalize the expression by the cell volume, then compute summaries
      # normalized values are very small, need to work on a log scale
      #norm.df <- apply(indiv.df, 2, FUN=function(X) log10(X)/log10(indiv.df[, "FSC.A"]^3))
      norm.df <- indiv.df
      if(typeof(norm.df) != "list"){
        norm.df <- as.data.frame(matrix(norm.df, ncol=ncol(indiv.df)))
        colnames(norm.df) <- colnames(indiv.df)
      }
      mean.param <- apply(norm.df, 2, mean)
      var.param <- apply(norm.df, 2, var)
      cv2.param <- apply(norm.df, 2, FUN=function(X) var(X)/(mean(X)^2))
      cv.param <- apply(norm.df, 2, FUN=function(X) sd(X)/mean(X))
      fano.param <- apply(norm.df, 2, FUN=function(X) var(X)/mean(X))
      # calculate cell sizes and correlations on original cell size, not normalized
      meanSize.param <- mean(indiv.df[, "FSC.A"]^3)
      sizeCor.param <- apply(norm.df, 2, FUN=function(X) cor(indiv.df[, "FSC.A"]^3, X))
      
      indiv.summary <- do.call(cbind.data.frame,
                               list("Variance"=var.param,
                                    "Mean"=mean.param,
                                    "CV2"=cv2.param,
                                    "CV"=cv.param,
                                    "Fano"=fano.param,
                                    "MeanSize"=meanSize.param,
                                    "SizeCor"=sizeCor.param))
      indiv.summary$Individual <- indiv
      # need to record numbers of cells, might need to remove very small samples
      indiv.summary$NCells <- nrow(indiv.df)
      indiv.summary$Parameters <- rownames(indiv.summary)
      indiv.summary.list[[indiv]] <- indiv.summary
    }
  }
  
  cell.summ.df <- do.call(rbind.data.frame,
                          indiv.summary.list)
  
  cell.summ.df$CellType <- celltype
  cell.summary.list[[celltype]] <- cell.summ.df
}

all.cell.df <- do.call(rbind.data.frame,
                       cell.summary.list)

# remove observations with < 100 cells
all.cell.df <- all.cell.df[all.cell.df$NCells > 100, ]

write.table(all.cell.df,
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel2C_summary.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

