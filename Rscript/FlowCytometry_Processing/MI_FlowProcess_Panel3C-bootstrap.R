#! /usr/bin/env Rscript

## TwinsUK flow cytometry processing
###############################
### Panel 3 - MAIT/NKT cells ##
###############################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 3C")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P3.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P03")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P3.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P3.flow.data) <- sample.ids

# subset samples into group C
P3.flow.data <- P3.flow.data[641:961]

# there is no compensation matrix in these files
# comp.matrix <- keyword(P3.flow.data[[1]])$`$SPILLOVER`
# P3.flow.data <- compensate(P3.flow.data, spillover <- comp.matrix)

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P3.flow.data))){
  events <- nrow(P3.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P3.fs_filt <- P3.flow.data[seq(along=P3.flow.data)]
if(length(extract)){P3.fs_filt <- P3.fs_filt[-(extract),]}

## subset P3.fs_filt for testing
#P3.fs_filt <- P3.fs_filt[1:20]

# remove channels that aren't interesting
param.desc <-  pData(parameters(P3.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P3.flow.data[[1]])
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

P3.fs_filt <- transform(P3.fs_filt, transform1)
P3.fs_filt <- transform(P3.fs_filt, transform2)
P3.fs_filt <- transform(P3.fs_filt, transform3)
P3.fs_filt <- transform(P3.fs_filt, transform4)
P3.fs_filt <- transform(P3.fs_filt, transform5)
P3.fs_filt <- transform(P3.fs_filt, transform6)
P3.fs_filt <- transform(P3.fs_filt, transform7)
P3.fs_filt <- transform(P3.fs_filt, transform8)

# gate on CD3+ cells
# ggcyto(P3.fs_filt[1:6], aes(x=`CD3-A`, y=`SSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 2.7e5))

# set up a CD3 gate
# setup side scatter and CD3 gates
side.points <- c(0, 0, 1.0e5, 1.0e5, 0)
cd3.points <- c(7.9, 12, 12, 7.9, 7.9)

cd3.gate <- polygonGate(filterId="CD3", .gate=cbind("SSC-A"=side.points, "V1-A"=cd3.points))
cd3.filter <- filter(P3.fs_filt, cd3.gate)
cd3.gates <- lapply(cd3.filter, 
                    function(res) filterDetails(res, "CD3")$filter)

# subset all CD3+ cells before all other analyses
P3.cd3 <- Subset(P3.fs_filt, cd3.gate)

# ggcyto(P3.fs_filt[1:6], aes(x=`CD3-A`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 2.7e5)) +
#   geom_gate(cd3.gates)

# # gate on singlets
# ggcyto(P3.fs_filt[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 3e5), y=c(0, 3e5))

# set up a singlets gate
# setup forward and side scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))
singlet.filter <- filter(P3.cd3, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# # gate on singlets
# ggcyto(P3.cd3[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 3e5), y=c(0, 3e5)) +
#   geom_gate(singlet.gates)

# subset all singlet cells before all other analyses
P3.singlet <- Subset(P3.cd3, singlet.gate)

# # gate on side scatter
# ggcyto(P3.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5))

# set up a side scatter gate requries a covariance matrix
# setup covariance gate
ssc.cov <- matrix(c(2e9, 1e6, 1.5e9, 1e10), ncol=2,
                  dimnames=list(c("SSC-H", "SSC-A"), c("SSC-H", "SSC-A")))
ssc.mean <- c("SSC-H"=2e4, "SSC-A"=2.5e4)
ssc.ellipsoid <- ellipsoidGate(filterId="SSC.Ellipse", .gate=ssc.cov, mean=ssc.mean)
ssc.filter <- filter(P3.singlet, ssc.ellipsoid)
ssc.gates <- lapply(ssc.filter, 
                    function(res) filterDetails(res, "SSC.Ellipse")$filter)

# # gate on side scatter
# ggcyto(P3.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.ellipsoid)

# subset all singlet cells before all other analyses
P3.ssc <- Subset(P3.singlet, ssc.gates)

extract <- numeric(0)
for(x in seq_len(length(P3.ssc))){
  events <- nrow(P3.ssc[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P3.ssc <- P3.ssc[-(extract),]}

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P3.ssc)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("B3-A", "B4-A")
P3.warping <- warpSet(P3.ssc, stains=setdiff(pars, second.warp))
P3.warping <- warpSet(P3.warping, stains=second.warp, peakNr=2)
#P3.warping <- warpSet(P3.warping, stains=third.warp, bwFac = 0.75, nbreaks=20)

# densityplot(x=~`R1-A`,
#             data=P3.warping,
#             channel=c("R1-A"))

# set up gating scheme using a GatingSet object
P3.gs <- GatingSet(P3.warping)

rm(list=c("P3.warping", "P3.flow.data", "P3.fs_filt",
          "P3.singlet", "P3.ssc", "P3.cd3"))
gc()

####################################################
# First separate the CD3 and TCR gamma-delta cells #
####################################################
gdtcr.points <- c(0, 6.5, 6.5, 0, 0)
cd3.points <- c(7.5, 7.5, 12, 12, 7.5)
abtecell.gate <- polygonGate(filterId="CD3", .gate=cbind("V1-A"=cd3.points, "B2-A"=gdtcr.points))

flowWorkspace::add(P3.gs, abtecell.gate, parent="root")
recompute(P3.gs)

gdtcr.points <- c(6.5, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 10.5, 10, 9.5, 9, 8.5, 8, 7.5, 6.5)
cd3.points <- c(7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 11.5, 10)
gdtecell.gate <- polygonGate(filterId="gdT", .gate=cbind("V1-A"=cd3.points, "B2-A"=gdtcr.points))

flowWorkspace::add(P3.gs, gdtecell.gate, parent="root")
recompute(P3.gs)

# ggcyto(P3.gs[1:6], aes(x=`B2-A`, y=`V1-A`), subset="root") +
#   geom_hex(bins=256) +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(abtecell.gate) +
#   geom_gate(gdtecell.gate)

# # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P3.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="CD3") +
#   geom_hex(bins=256) +
#   lims(x=c(0, 12), y=c(0, 12))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7.5, 12, 12, 7.5, 7.5)
cd4.points <- c(0, 0, 7.5, 7.5, 0)
cd8.gate <- polygonGate(filterId="CD8", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7.5, 7.5, 0, 0, 7.5)
cd4.points <- c(7.5, 11, 11, 7.5, 7.5)
cd4.gate <- polygonGate(filterId="CD4", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(0, 7.5, 7.5, 0, 0)
cd4.points <- c(0, 0, 7.5, 7.5, 0)
DN.gate <- polygonGate(filterId="DN", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(8, 12, 12, 8, 8)
cd4.points <- c(7.5, 7.5, 11, 11, 7.5)
DP.gate <- polygonGate(filterId="DP", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# # # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P3.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="root") +
#   geom_hex(bins=256) +
#   labs(x="CD4", y="CD8") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cd8.gate) +
#   geom_gate(cd4.gate) +
#   geom_gate(DN.gate) +
#   geom_gate(DP.gate)

# add the T cell gates
flowWorkspace::add(P3.gs, cd4.gate, parent="CD3")
recompute(P3.gs)

flowWorkspace::add(P3.gs, cd8.gate, parent="CD3")
recompute(P3.gs)

flowWorkspace::add(P3.gs, DN.gate, parent="CD3")
recompute(P3.gs)

flowWorkspace::add(P3.gs, DP.gate, parent="CD3")
recompute(P3.gs)

# T cell groups defined by this population:
# gamma-delta T cells
# NKT cells
# MAIT cells

###############
# NKT TCRV24+ #
###############
# NKT cells co-express TCR V24 and CD161
v24.points <- c(7.75, 7.75, 11, 11, 7.75)
cd161.points <- c(7.5, 11, 11, 7.5, 7.5)
nkt.gate <- polygonGate(filterId="NKT", .gate=cbind("B3-A"=cd161.points, "R1-A"=v24.points))

flowWorkspace::add(P3.gs, nkt.gate, parent="root")
recompute(P3.gs)

##################
# MAIT TCR V7.2+ #
##################
# MAIT cells co-express TCR V.2 and CD161
v72.points <-  c(7.75, 7.75, 11, 11, 7.75)
cd161.points <- c(7.5, 11, 11, 7.5, 7.5)
mait.gate <- polygonGate(filterId="MAIT", .gate=cbind("B1-A"=v72.points, "R1-A"=cd161.points))

flowWorkspace::add(P3.gs, mait.gate, parent="/root")
recompute(P3.gs)

# ggcyto(P3.gs[1:6], aes(x=`B1-A`, y=`R1-A`), subset="root") +
#   geom_hex(bins=128) +
#   labs(x="CD161", y="TCRV7.2") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(mait.gate)

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P3.gs)[3:length(getNodes(P3.gs))]
samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P3.gs, y=node)
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
P3.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P3.df <- P3.df[!duplicated(P3.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD3", "HLADR", "TCRv7.2", 
                     "TCRgd", "CD161", "CD8B", "TCRv24", "CD4")
P3.df <- P3.df[, parm.goi]
colnames(P3.df) <- names(parm.goi)

celltype <- unlist(lapply(strsplit(rownames(P3.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P3.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P3.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P3.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P3.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
rm(list=c("P3.gs"))
gc()

individuals <- levels(P3.meta$Individual)
celltypes <- levels(P3.meta$CellType)

print("Performing bootstraps")
cell.summary.list <- list()
n.boots <- 100
boot.prop <- 0.66

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P3.df[grepl(rownames(P3.df), pattern=celltype), ]
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
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel3C_bootstrap.txt",
            quote=FALSE, row.names=FALSE, sep="\t")