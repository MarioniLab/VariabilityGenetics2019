#! /usr/bin/env Rscript

## TwinsUK flow cytometry processing
############################
### Panel 9 - Th subsets ##
###########################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 9C")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P9.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P09")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P9.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P9.flow.data) <- sample.ids

P9.flow.data <- P9.flow.data[634:951]
# there is no compensation matrix in these files
# comp.matrix <- keyword(P9.flow.data[[1]])$`$SPILLOVER`
# P9.flow.data <- compensate(P9.flow.data, spillover <- comp.matrix)

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P9.flow.data))){
  events <- nrow(P9.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P9.fs_filt <- P9.flow.data[seq(along=P9.flow.data)]
if(length(extract)){P9.fs_filt <- P9.fs_filt[-(extract),]}

## subset P9.fs_filt for testing
#P9.fs_filt <- P9.fs_filt[1:20]

# remove channels that aren't interesting
param.desc <-  pData(parameters(P9.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P9.flow.data[[1]])
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

P9.fs_filt <- transform(P9.fs_filt, transform1)
P9.fs_filt <- transform(P9.fs_filt, transform2)
P9.fs_filt <- transform(P9.fs_filt, transform3)
P9.fs_filt <- transform(P9.fs_filt, transform4)
P9.fs_filt <- transform(P9.fs_filt, transform5)
P9.fs_filt <- transform(P9.fs_filt, transform6)
P9.fs_filt <- transform(P9.fs_filt, transform7)
P9.fs_filt <- transform(P9.fs_filt, transform8)

extract <- numeric(0)
for(x in seq_len(length(P9.fs_filt))){
  events <- nrow(P9.fs_filt[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P9.fs_filt <- P9.fs_filt[-(extract),]}
#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P9.fs_filt)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("B2-A")
third.warp <- c("B1-A")
P9.warping <- warpSet(P9.fs_filt, stains=setdiff(pars, c(second.warp, third.warp)))
P9.warping <- warpSet(P9.warping, stains=second.warp, peakNr=2, bwFac = 0.75, nbreaks=20)
P9.warping <- warpSet(P9.warping, stains=third.warp, bwFac = 2.5, nbreaks=20, peakNr=1)

densityplot(x=~`B1-A`,
            data=P9.warping,
            channel=c("B1-A"))

# set up gating scheme using a GatingSet object
P9.gs <- GatingSet(P9.warping)
rm(list=c("P9.warping", "P9.flow.data", "P9.fs_filt"))
gc()

# # # gate on TCRgd- cells
# ggcyto(P9.fs_filt[1:6], aes(x=`B2-A`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 2.7e5))

# set up gates to select gamma-delta TCR + and - cells
# setup side scatter and CD3 gates
side.points <- c(0, 0, 1.0e5, 1.0e5, 0)
tcrgd.points <- c(2.5, 6, 6, 2.5, 2.5)

gdtneg.gate <- polygonGate(filterId="gdTneg", .gate=cbind("SSC-A"=side.points, "B2-A"=tcrgd.points))

# set up a CD3 gate
# setup side scatter and CD3 gates
side.points <- c(0, 0, 1.0e5, 1.0e5, 0)
tcrgd.points <- c(6.5, 11, 11, 6.5, 6.5)

gdtpos.gate <- polygonGate(filterId="gdT", .gate=cbind("SSC-A"=side.points, "B2-A"=tcrgd.points))

# # gate on TCRgd- cells
# ggcyto(P9.gs[1:6], aes(x=`B2-A`, y=`SSC-A`), subset="root") +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 2.7e5)) +
#   geom_gate(gdtneg.gate) +
#   geom_gate(gdtpos.gate)

flowWorkspace::add(P9.gs, gdtneg.gate, parent="root")
recompute(P9.gs)

flowWorkspace::add(P9.gs, gdtpos.gate, parent="root")
recompute(P9.gs)

# set up a singlets gate
# setup forward and side scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))

# # # gate on singlets
# ggcyto(P9.gs[1:6], aes(x=`FSC-H`, y=`FSC-W`), subset="gdTneg") +
#   geom_hex(bins=128) +
#   lims(x=c(0, 3e5), y=c(0, 3e5)) +
#   geom_gate(singlet.gate)

flowWorkspace::add(P9.gs, singlet.gate, parent="gdT")
recompute(P9.gs)

flowWorkspace::add(P9.gs, singlet.gate, parent="gdTneg")
recompute(P9.gs)

# # gate on side scatter
# set up a side scatter gate requries a covariance matrix
# setup covariance gate
ssc.cov <- matrix(c(2e9, 1e6, 1.5e9, 1e10), ncol=2,
                  dimnames=list(c("SSC-H", "SSC-A"), c("SSC-H", "SSC-A")))
ssc.mean <- c("SSC-H"=2e4, "SSC-A"=2.5e4)
ssc.ellipsoid <- ellipsoidGate(filterId="SSC.Ellipse", .gate=ssc.cov, mean=ssc.mean)
# 
# ggcyto(P9.gs[1:6], aes(x=`SSC-H`, y=`SSC-A`), subset="/gdT/Singlet") +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.ellipsoid)

flowWorkspace::add(P9.gs, ssc.ellipsoid, parent="/gdT/Singlet")
recompute(P9.gs)

flowWorkspace::add(P9.gs, ssc.ellipsoid, parent="/gdTneg/Singlet")
recompute(P9.gs)

P9.prewarp <- getData(P9.gs, "/gdTneg/Singlet/SSC.Ellipse")
extract <- numeric(0)
for(x in seq_len(length(P9.prewarp))){
  events <- nrow(P9.prewarp[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P9.prewarp <- P9.prewarp[-(extract),]}

#######################################
# Do I need to do a second warping at this point?
# normalisation of these data tends to be quite specific to the parameters
second.warp <- c("V2-A", "R2-A")
#third.warp <- c("B1-A")
P9.warping <- warpSet(P9.prewarp, 
                      stains=setdiff(pars, second.warp))
P9.warping <- warpSet(P9.warping, stains=second.warp, peakNr=3, bwFac=0.55, nbreaks=15)
P9.warping <- warpSet(P9.warping, stains=third.warp, bwFac = 1, nbreaks=20, peakNr=1)
# 
# densityplot(x=~`R2-A`,
#             data=P9.warping,
#             channel=c("R2-A"))

P9.gs <- GatingSet(P9.warping)

############################################
# First the CD3+ alpha-beta T cell subsets #
############################################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7.5, 12, 12, 7.5, 7.5)
cd4.points <- c(0, 0, 7.5, 7.5, 0)
cd8.gate <- polygonGate(filterId="CD8", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7.5, 7.5, 0, 0, 7.5)
cd4.points <- c(7.75, 11, 11, 7.75, 7.75)
cd4.gate <- polygonGate(filterId="CD4", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(0, 7.5, 7.5, 0, 0)
cd4.points <- c(0, 0, 7.5, 7.5, 0)
DN.gate <- polygonGate(filterId="DN", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(8, 12, 12, 8, 8)
cd4.points <- c(7.75, 7.75, 11, 11, 7.75)
DP.gate <- polygonGate(filterId="DP", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# # # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P9.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="root") +
#   geom_hex(bins=256) +
#   labs(x="CD4", y="CD8") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cd8.gate) +
#   geom_gate(cd4.gate) +
#   geom_gate(DN.gate) +
#   geom_gate(DP.gate)

# add the T cell gates
flowWorkspace::add(P9.gs, cd4.gate, parent="root")
recompute(P9.gs)

flowWorkspace::add(P9.gs, cd8.gate, parent="root")
recompute(P9.gs)

flowWorkspace::add(P9.gs, DN.gate, parent="root")
recompute(P9.gs)

flowWorkspace::add(P9.gs, DP.gate, parent="root")
recompute(P9.gs)

# T cell groups defined by this population:
# CD8+CCR6-CXCR3-
# CD8+CCR6+CXCR3+
# CD4+CCR6+CXCR3-
# CD4+CCR6-CXCR3-CXCR5+
# CD4+CCR6-CXCR3-CRTh2+

###################
# CD8+CCR6-CXCR3- #
###################
# What T cell subset are these? CD8+ Th2?
ccr6.points <- c(0, 0, 7, 7, 0)
cxcr3.points <- c(0, 7.5, 7.5, 0, 0)
cd8th1.gate <- polygonGate(filterId="Th2", .gate=cbind("B3-A"=ccr6.points, "R1-A"=cxcr3.points))

flowWorkspace::add(P9.gs, cd8th1.gate, parent="/CD8")
recompute(P9.gs)

###################
# CD8+CCR6+CXCR3- #
###################
# What T cell subset are these? CD8+ Th1?
ccr6.points <- c(7.5, 7.5, 12, 12, 7.5)
cxcr3.points <- c(0, 7.5, 7.5, 0, 0)
cd8ccr6.gate <- polygonGate(filterId="Th1", .gate=cbind("B3-A"=ccr6.points, "R1-A"=cxcr3.points))

flowWorkspace::add(P9.gs, cd8ccr6.gate, parent="/CD8")
recompute(P9.gs)

###################
# CD4+CCR6-CXCR3- #
###################
# What T cell subset are these? CD8+ Th1?
ccr6.points <- c(0, 0, 7, 7, 0)
cxcr3.points <- c(0, 7.5, 7.5, 0, 0)
cd4ccr6neg.gate <- polygonGate(filterId="CCR6neg", .gate=cbind("B3-A"=ccr6.points, "R1-A"=cxcr3.points))

flowWorkspace::add(P9.gs, cd4ccr6neg.gate, parent="/CD4")
recompute(P9.gs)

###################
# CD4+CCR6+CXCR3- #
###################
# What T cell subset are these? CD8+ Th1?
ccr6.points <- c(7.5, 7.5, 12, 12, 7.5)
cxcr3.points <- c(0, 7.5, 7.5, 0, 0)
cd4ccr6.gate <- polygonGate(filterId="Th17", .gate=cbind("B3-A"=ccr6.points, "R1-A"=cxcr3.points))

flowWorkspace::add(P9.gs, cd4ccr6.gate, parent="/CD4")
recompute(P9.gs)

###################
# CD4+CCR6-CXCR3+ #
###################
# What T cell subset are these? CD8+ Th1?
ccr6.points <- c(0, 0, 7, 7, 0)
cxcr3.points <- c(7.5, 11, 11, 7.5, 7.5)
cd4th1.gate <- polygonGate(filterId="Th1", .gate=cbind("B3-A"=ccr6.points, "R1-A"=cxcr3.points))

flowWorkspace::add(P9.gs, cd4th1.gate, parent="/CD4")
recompute(P9.gs)

#########################
# CD4+CCR6-CXCR3-CXCR5+ #
#########################
# I think these are Tfh cells
crth2.points <- c(0, 0, 7.5, 7.5, 0)
cxcr5.points <- c(7.5, 11, 11, 7.5, 7.5)
tfh.gate <- polygonGate(filterId="Tfh", .gate=cbind("V2-A"=crth2.points, "B1-A"=cxcr5.points))

flowWorkspace::add(P9.gs, tfh.gate, parent="/CD4/CCR6neg")
recompute(P9.gs)

#########################
# CD4+CCR6-CXCR3-CRTh2+ #
#########################
# I think these are Tfh cells?
crth2.points <- c(7.5, 7.5, 11, 11, 7.5)
cxcr5.points <- c(0, 7.5, 7.5, 0, 0)
th2.gate <- polygonGate(filterId="Th2", .gate=cbind("V2-A"=crth2.points, "B1-A"=cxcr5.points))

flowWorkspace::add(P9.gs, th2.gate, parent="/CD4/CCR6neg")
recompute(P9.gs)

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P9.gs)[4:length(getNodes(P9.gs))]

# remove the CCR6neg gate  as it contains Tfh and Th2 cells which
# are defined later
all.nodes <- all.nodes[!grepl(all.nodes, pattern="CCR6neg$")]

samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P9.gs, y=node)
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
P9.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P9.df <- P9.df[!duplicated(P9.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CCR4", "CRTh2", "CXCR5", 
                     "TCRgd", "CCR6", "CD8B", "CXCR3", "CD4")
P9.df <- P9.df[, parm.goi]
colnames(P9.df) <- names(parm.goi)

celltype <- unlist(lapply(strsplit(rownames(P9.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P9.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P9.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P9.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P9.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
rm(list=c("P9.gs"))
gc()

individuals <- levels(P9.meta$Individual)
celltypes <- levels(P9.meta$CellType)

print("Performing bootstraps")
cell.summary.list <- list()
n.boots <- 100
boot.prop <- 0.66

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P9.df[grepl(rownames(P9.df), pattern=celltype), ]
  indiv.summary.list <- list()
  
  for(j in seq_along(individuals)){
    indiv <- individuals[j]
    indiv.df <- celltype.df[grepl(rownames(celltype.df), pattern=indiv), ]
    
    if(nrow(indiv.df) > 0){
      # truncate negative values at 0
      indiv.df[indiv.df < 0] <- 0
      
      # normalize the expression by the cell volume, then compute summaries
      # normalized values are very small, need to work on a log scale
      norm.df <- indiv.df
      
      # subsample the number of cells and calculate the summary statistics for each of these
      # use this to calculate the mean and variance of these bootstraps
      # I also need to save each boot strap for latter to bootstrap the 
      # mean-variance dependence fitting
      
      if(typeof(norm.df) != "list"){
        norm.df <- as.data.frame(matrix(norm.df, ncol=ncol(indiv.df)))
        colnames(norm.df) <- colnames(indiv.df)
      }
      boot.list <- list()
      for(k in seq_along(1:n.boots)){
        boot.rows <- sample(seq_len(nrow(norm.df)), size=floor(boot.prop*nrow(norm.df)))
        
        boot.df <- norm.df[boot.rows, ]
        boot.mean <- apply(boot.df, 2, mean)
        boot.var <- apply(boot.df, 2, var)
        boot.cv2 <- apply(boot.df, 2, FUN=function(X) var(X)/(mean(X)^2))
        
        boot.list[[as.character(k)]] <- do.call(cbind.data.frame,
                                                list("Variance"=boot.var,
                                                     "Mean"=boot.mean,
                                                     "CV2"=boot.cv2,
                                                     "Resample"=k))
      }
    }
    
    boot.summ <- do.call(rbind.data.frame,
                         boot.list)
    boot.summ$Parameter <- gsub(rownames(boot.summ), pattern="^[0-9]+.", replacement="")
    boot.summ$Individual <- indiv
    indiv.summary.list[[indiv]] <- boot.summ
  }
  
  cell.summ.df <- do.call(rbind.data.frame,
                          indiv.summary.list)
  cell.summ.df$CellType <- celltype
  cell.summary.list[[celltype]] <- cell.summ.df
}

all.cell.df <- do.call(rbind.data.frame,
                       cell.summary.list)

write.table(all.cell.df,
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel9C_bootstrap.txt",
            quote=FALSE, row.names=FALSE, sep="\t")
