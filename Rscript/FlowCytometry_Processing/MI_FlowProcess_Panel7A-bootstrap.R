#! /usr/bin/env Rscript

## MI flow cytometry processing
############################
### Panel 7 - PMN subsets ##
############################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 7A")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P7.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P07")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P7.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P7.flow.data) <- sample.ids

# subset samples into group A
P7.flow.data <- P7.flow.data[1:324]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P7.flow.data))){
  events <- nrow(P7.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P7.fs_filt <- P7.flow.data[seq(along=P7.flow.data)]
if(length(extract)){P7.fs_filt <- P7.fs_filt[-(extract),]}

# remove channels that aren't interesting
param.desc <-  pData(parameters(P7.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P7.flow.data[[1]])
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

P7.fs_filt <- transform(P7.fs_filt, transform1)
P7.fs_filt <- transform(P7.fs_filt, transform2)
P7.fs_filt <- transform(P7.fs_filt, transform3)
P7.fs_filt <- transform(P7.fs_filt, transform4)
P7.fs_filt <- transform(P7.fs_filt, transform5)
P7.fs_filt <- transform(P7.fs_filt, transform6)
P7.fs_filt <- transform(P7.fs_filt, transform7)
P7.fs_filt <- transform(P7.fs_filt, transform8)

# gate on leukocytes from FSC and SSC
# ggcyto(P7.fs_filt[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5))

# set up a singlets gate
# setup forward and side scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))
singlet.filter <- filter(P7.fs_filt, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# ggcyto(P7.fs_filt[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(singlet.gate)

P7.singlet <- Subset(P7.fs_filt, singlet.gates)

# forward and side scater gate
area.points <- c(0, 0, 2.5e5, 2.5e5, 2.5e4, 0)
height.points <- c(0, 2e4, 2.25e5, 1.5e5, 0, 0)
ssc.gate <- polygonGate(filterId="SSC", .gate=cbind("SSC-A"=area.points, "SSC-H"=height.points))
ssc.filter <- filter(P7.singlet, ssc.gate)
ssc.gates <- lapply(ssc.filter, 
                    function(res) filterDetails(res, "SSC")$filter)

# ggcyto(P7.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.gate)

P7.ssc <- Subset(P7.singlet, ssc.gates)
extract <- numeric(0)
for(x in seq_len(length(P7.ssc))){
  events <- nrow(P7.ssc[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P7.ssc <- P7.ssc[-(extract),]}

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P7.ssc)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("V2-A", "R1-A")
third.warp <- c("R2-A")
P7.warping <- warpSet(P7.ssc, stains=setdiff(pars, c(second.warp, third.warp)))
P7.warping <- warpSet(P7.warping, stains=second.warp, bwFac = 0.75, nbreaks=15)
P7.warping <- warpSet(P7.warping, stains=third.warp, bwFac = 2, nbreaks=20, peakNr=2)

densityplot(x=~`R2-A`,
            data=P7.warping,
            channel=c("R2-A"))

# set up gating scheme using a GatingSet object
P7.gs <- GatingSet(P7.warping)

rm(list=c("P7.warping", "P7.flow.data", "P7.fs_filt",
          "P7.singlet", "P7.ssc"))
gc()

###############
## Set up gates for PMNs
##############
# PMN defined by this panel
# Neutrophils - CD16hiCDw125-
# Eosinophils - SSChiCDw125+
# Basophils - FCER1hiCDw125-

######
# Neutrophils
######
# Gate on the big CD16+ cells first
cd16.points <- c(10, 10, 9.75, 9.5, 10, 12, 12, 10.5, 10)
fsc.points <- c(7.5e4, 1e5, 1.25e5, 1.5e5, 2.5e5, 2.5e5, 7.5e4, 6e4, 7.5e4)
cd16.gate <- polygonGate(filterId="CD16", .gate=cbind("B3-A"=cd16.points, "FSC-A"=fsc.points))

cd16.points <- c(0, 0, 7.5, 7.5, 0)
fsc.points <- c(2.5e4, 2.5e5, 2.5e5, 0, 0)
notcd16.gate <- polygonGate(filterId="notCD16", .gate=cbind("B3-A"=cd16.points, "FSC-A"=fsc.points))

flowWorkspace::add(P7.gs, cd16.gate, parent="root")
recompute(P7.gs)

flowWorkspace::add(P7.gs, notcd16.gate, parent="root")
recompute(P7.gs)

# get the CDw125- (CD125) CD16+ neutrophils
cd16.points <- c(9.5, 9.5, 12, 12, 9.5)
cd125.points <- c(0, 7, 7, 0, 0)
neut.gate <- polygonGate(filterId="Neutrophil", .gate=cbind("B3-A"=cd16.points, "B2-A"=cd125.points))

flowWorkspace::add(P7.gs, neut.gate, parent="CD16")
recompute(P7.gs)

######
# Eosinophils
######
# Gate on the big CD16+ cells first
cd32.points <- c(7, 6.75, 6.75, 6.75, 8, 8.75, 9, 9, 9, 9, 9, 9, 9)
cd125.points <- c(6, 7, 7.75, 8, 8, 7.75, 7.75, 7.75, 7.5, 7.25, 7, 6.75, 6.5)
cd125.gate <- polygonGate(filterId="CD125", .gate=cbind("B4-A"=cd32.points, "B2-A"=cd125.points))

flowWorkspace::add(P7.gs, cd125.gate, parent="notCD16")
recompute(P7.gs)

# get the big granular eosinophils
fsc.points <- c(5e4, 5e4, 2.5e5, 2.5e5, 5e4)
ssc.points <- c(1e5, 2.5e5, 2.5e5, 1e5, 1e5)
eosin.gate <- polygonGate(filterId="Eosinophil", .gate=cbind("FSC-A"=fsc.points, "SSC-A"=ssc.points))

flowWorkspace::add(P7.gs, eosin.gate, parent="CD125")
recompute(P7.gs)

######
# Basophils
######
# Gate on CD203+ FCER1+
fcer.points <- c(8, 11, 11, 8, 8)
cd203.points <- c(3, 3, 10, 10, 3)
fcer.gate <- polygonGate(filterId="FCER", .gate=cbind("B1-A"=fcer.points, "R1-A"=cd203.points))

flowWorkspace::add(P7.gs, fcer.gate, parent="notCD16")
recompute(P7.gs)

# Gate on CD203+ FCER1+
fcer.points <- c(7.5, 7.5, 11, 11, 7.5)
cd125.points <- c(2, 7.5, 7.5, 2, 2)
baso.gate <- polygonGate(filterId="Basophil", .gate=cbind("B1-A"=fcer.points, "B2-A"=cd125.points))

flowWorkspace::add(P7.gs, baso.gate, parent="FCER")
recompute(P7.gs)


# ggcyto(P7.gs[1:6], aes(x=`B2-A`, y=`B1-A`), subset="FCER") +
#   geom_hex(bins=128) +
#   labs(y="FCER1", x="CD125") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(baso.gate)


###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P7.gs)[4:length(getNodes(P7.gs))]

# remove superfluous intermediate nodes
all.nodes <- all.nodes[!grepl(all.nodes, pattern="CD125$")]
all.nodes <- all.nodes[!grepl(all.nodes, pattern="FCER$")]

samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P7.gs, y=node)
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
P7.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P7.df <- P7.df[!duplicated(P7.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD62L", "Viability", "FcERIa",
                     "CD125", "CD16", "CD32", "CD203c", "BLANK")
P7.df <- P7.df[, parm.goi]
colnames(P7.df) <- names(parm.goi)

celltype <- unlist(lapply(strsplit(rownames(P7.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P7.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P7.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P7.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P7.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
rm(list=c("P7.gs"))
gc()

individuals <- levels(P7.meta$Individual)
celltypes <- levels(P7.meta$CellType)

print("Performing bootstraps")
cell.summary.list <- list()
n.boots <- 100
boot.prop <- 0.66

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P7.df[grepl(rownames(P7.df), pattern=celltype), ]
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
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel7A_bootstrap.txt",
            quote=FALSE, row.names=FALSE, sep="\t")