#! /usr/bin/env Rscript

## MI flow cytometry processing
##############################################
### Panel 1 - Naive & Memory T cell subsets ##
##############################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 1B")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
P1.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P01")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P1.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P1.flow.data) <- sample.ids

# subset into group B
P1.flow.data <- P1.flow.data[327:652]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P1.flow.data))){
  events <- nrow(P1.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P1.fs_filt <- P1.flow.data[seq(along=P1.flow.data)]
if(length(extract)){P1.fs_filt <- P1.fs_filt[-(extract),]}

## subset P1.fs_filt for testing
#P1.fs_filt <- P1.fs_filt[1:20]

# remove channels that aren't interesting
param.desc <-  pData(parameters(P1.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P1.flow.data[[1]])
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

P1.fs_filt <- transform(P1.fs_filt, transform1)
P1.fs_filt <- transform(P1.fs_filt, transform2)
P1.fs_filt <- transform(P1.fs_filt, transform3)
P1.fs_filt <- transform(P1.fs_filt, transform4)
P1.fs_filt <- transform(P1.fs_filt, transform5)
P1.fs_filt <- transform(P1.fs_filt, transform6)
P1.fs_filt <- transform(P1.fs_filt, transform7)
P1.fs_filt <- transform(P1.fs_filt, transform8)

# P1.fs_filt <- transform(P1.fs_filt, transformForwardH)
# P1.fs_filt <- transform(P1.fs_filt, transformSideA)
# P1.fs_filt <- transform(P1.fs_filt, transformForwardA)

# # gate on CD3+ cells
# ggcyto(P1.fs_filt[1:6], aes(x=`CD3-A`, y=`SSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 2.7e5))

# set up a CD3 gate
# setup side scatter and CD3 gates
side.points <- c(0, 0, 1.0e5, 1.0e5, 0)
cd3.points <- c(7.9, 12, 12, 7.9, 7.9)

cd3.gate <- polygonGate(filterId="CD3", .gate=cbind("SSC-A"=side.points, "V1-A"=cd3.points))
cd3.filter <- filter(P1.fs_filt, cd3.gate)
cd3.gates <- lapply(cd3.filter, 
                    function(res) filterDetails(res, "CD3")$filter)

# subset all CD3+ cells before all other analyses
P1.cd3 <- Subset(P1.fs_filt, cd3.gate)

# ggcyto(P1.cd3[1:6], aes(x=`CD3-A`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 2.7e5)) +
#   geom_gate(cd3.gates)

# gate on singlets
# ggcyto(P1.fs_filt[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 3e5), y=c(0, 3e5))

# set up a singlets gate
# setup forward and side scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))
singlet.filter <- filter(P1.cd3, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# # gate on singlets
# ggcyto(P1.cd3[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 3e5), y=c(0, 3e5)) +
#   geom_gate(singlet.gates)

# subset all singlet cells before all other analyses
P1.singlet <- Subset(P1.cd3, singlet.gate)

# # gate on side scatter
# ggcyto(P1.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5))

# set up a side scatter gate requries a covariance matrix
# setup covariance gate
ssc.cov <- matrix(c(2e9, 1e6, 1.5e9, 1e10), ncol=2,
                  dimnames=list(c("SSC-H", "SSC-A"), c("SSC-H", "SSC-A")))
ssc.mean <- c("SSC-H"=2e4, "SSC-A"=2.5e4)
ssc.ellipsoid <- ellipsoidGate(filterId="SSC.Ellipse", .gate=ssc.cov, mean=ssc.mean)
ssc.filter <- filter(P1.singlet, ssc.ellipsoid)
ssc.gates <- lapply(ssc.filter, 
                    function(res) filterDetails(res, "SSC.Ellipse")$filter)

# # gate on side scatter
# ggcyto(P1.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.ellipsoid)

# subset all singlet cells before all other analyses
P1.ssc <- Subset(P1.singlet, ssc.gates)
extract <- numeric(0)
for(x in seq_len(length(P1.ssc))){
  events <- nrow(P1.ssc[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P1.ssc <- P1.ssc[-(extract),]}

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P1.ssc)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("V2-A", "R2-A")
third.warp <- c("R1-A")
P1.warping <- warpSet(P1.ssc, stains=setdiff(pars, c(second.warp, third.warp)))
P1.warping <- warpSet(P1.warping, stains=second.warp, peakNr=2)
P1.warping <- warpSet(P1.warping, stains=third.warp, bwFac = 0.75, nbreaks=20)

# densityplot(x=~`R1-A`,
#             data=P1.warping,
#             channel=c("R1-A"))

# set up gating scheme using a GatingSet object
P1.gs <- GatingSet(P1.warping)

rm(list=c("P1.warping", "P1.flow.data", "P1.fs_filt",
          "P1.singlet", "P1.ssc", "P1.cd3"))
gc()

# # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P1.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="root") +
#   geom_hex(bins=256) +
#   lims(x=c(0, 12), y=c(0, 12))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7.5, 11, 11, 7.5, 7.5)
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
cd8.points <- c(8, 11, 11, 8, 8)
cd4.points <- c(7, 7, 11, 11, 7)
DP.gate <- polygonGate(filterId="DP", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# # # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P1.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="root") +
#   geom_hex(bins=256) +
#   labs(x="CD4", y="CD8") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cd8.gate) +
#   geom_gate(cd4.gate) +
#   geom_gate(DN.gate) +
#   geom_gate(DP.gate)

# add the T cell gates
flowWorkspace::add(P1.gs, cd4.gate, parent="root")
recompute(P1.gs)

flowWorkspace::add(P1.gs, cd8.gate, parent="root")
recompute(P1.gs)

flowWorkspace::add(P1.gs, DN.gate, parent="root")
recompute(P1.gs)

flowWorkspace::add(P1.gs, DP.gate, parent="root")
recompute(P1.gs)

# T cell groups defined by this population:
# T naive
# T central memory
# T effector memory
# T effector memory RA

######################
# Naive CD27+CD45RA+ #
######################
# naive T cells co-express CD45RA and CD27 (also CCR7+)
cd45ra.points <- c(7.5, 7.5, 12, 12, 7.5)
cd27.points <- c(8, 11, 11, 8, 8)
tnaive.gate <- polygonGate(filterId="Tnaive", .gate=cbind("B3-A"=cd27.points, "B1-A"=cd45ra.points))

flowWorkspace::add(P1.gs, tnaive.gate, parent="/CD4")
recompute(P1.gs)

flowWorkspace::add(P1.gs, tnaive.gate, parent="/CD8")
recompute(P1.gs)

###############################
# Central memory CD27+CD45RA- #
###############################
# central memmory T cells express CD27 but lack CD45RA and CCR7
cd45ra.points <- c(2.5, 2.5, 7.5, 7.5, 2.5)
cd27.points <- c(7, 11, 11, 7, 7)
tcm.gate <- polygonGate(filterId="TCM", .gate=cbind("B3-A"=cd27.points, "B1-A"=cd45ra.points))

flowWorkspace::add(P1.gs, tcm.gate, parent="/CD4")
recompute(P1.gs)

flowWorkspace::add(P1.gs, tcm.gate, parent="/CD8")
recompute(P1.gs)

################################
# Effector memory CD27-CD45RA- #
################################
# Effector memmory T cells lack CD27, CD45RA and CCR7
cd45ra.points <- c(2.5, 2.5, 7.5, 7.5, 2.5)
cd27.points <- c(2.5, 7, 7, 2.5, 2.5)
tem.gate <- polygonGate(filterId="TEM", .gate=cbind("B3-A"=cd27.points, "B1-A"=cd45ra.points))

flowWorkspace::add(P1.gs, tem.gate, parent="/CD4")
recompute(P1.gs)

flowWorkspace::add(P1.gs, tem.gate, parent="/CD8")
recompute(P1.gs)

#######################
# T EMRA CD27-CD45RA+ #
#######################
# Effector memmory T cells express CD45RA but lack CD27 and CCR7
cd45ra.points <- c(7.5, 7.5, 11, 11, 7.5)
cd27.points <- c(2.5, 7, 7, 2.5, 2.5)
temra.gate <- polygonGate(filterId="TEMRA", .gate=cbind("B3-A"=cd27.points, "B1-A"=cd45ra.points))

flowWorkspace::add(P1.gs, temra.gate, parent="/CD4")
recompute(P1.gs)

flowWorkspace::add(P1.gs, temra.gate, parent="/CD8")
recompute(P1.gs)

# # plot CD27 Vs CD45RA to get major T cell subsets
# ggcyto(P1.gs[1:6], aes(x=`B3-A`, y=`B1-A`), subset="CD8") +
#   geom_hex(bins=128) +
#   labs(x="CD27", y="CD45RA") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(tnaive.gate) +
#   geom_gate(tcm.gate) +
#   geom_gate(tem.gate) +
#   geom_gate(temra.gate)

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P1.gs)[4:length(getNodes(P1.gs))]
samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P1.gs, y=node)
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
P1.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P1.df <- P1.df[!duplicated(P1.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD3", "HLADR", "CD45RA", 
                     "CD8A", "CD27", "CD8B", "CCR7", "CD4")
P1.df <- P1.df[, parm.goi]
colnames(P1.df) <- names(parm.goi)

celltype <- unlist(lapply(strsplit(rownames(P1.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P1.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P1.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P1.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P1.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
rm(list=c("P1.gs"))
gc()

individuals <- levels(P1.meta$Individual)
celltypes <- levels(P1.meta$CellType)

print("Performing bootstraps")
cell.summary.list <- list()
n.boots <- 100
boot.prop <- 0.66

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P1.df[grepl(rownames(P1.df), pattern=celltype), ]
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
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel1B_bootstrap.txt",
            quote=FALSE, row.names=FALSE, sep="\t")