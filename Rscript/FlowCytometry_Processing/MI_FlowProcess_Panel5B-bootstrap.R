#! /usr/bin/env Rscript

## MI flow cytometry processing
#####################################
### Panel 5 - Major lineages panel ##
#####################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 5B")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P5.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P05")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P5.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P5.flow.data) <- sample.ids

# subset samples into group B
P5.flow.data <- P5.flow.data[321:640]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P5.flow.data))){
  events <- nrow(P5.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P5.fs_filt <- P5.flow.data[seq(along=P5.flow.data)]
if(length(extract)){P5.fs_filt <- P5.fs_filt[-(extract),]}

# remove channels that aren't interesting
param.desc <-  pData(parameters(P5.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P5.flow.data[[1]])
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

P5.fs_filt <- transform(P5.fs_filt, transform1)
P5.fs_filt <- transform(P5.fs_filt, transform2)
P5.fs_filt <- transform(P5.fs_filt, transform3)
P5.fs_filt <- transform(P5.fs_filt, transform4)
P5.fs_filt <- transform(P5.fs_filt, transform5)
P5.fs_filt <- transform(P5.fs_filt, transform6)
P5.fs_filt <- transform(P5.fs_filt, transform7)
P5.fs_filt <- transform(P5.fs_filt, transform8)

# gate on all CD45+ cells
# ggcyto(P5.fs_filt[1:6], aes(x=`FSC-A`, y=`B2-A`)) +
#   geom_hex(bins=128) +
#   labs(x="FSC-A", y="CD45") +
#   lims(x=c(0, 2.5e5), y=c(0, 12))

# set up a CD45+ gate
cd45.points <- c(8, 8, 12.5, 12.5, 8)
fsc.points <- c(0, 2.5e5, 2.5e5, 0, 0)

haem.gate <- polygonGate(filterId="CD45", .gate=cbind("B2-A"=cd45.points, "FSC-A"=fsc.points))
haem.filter <- filter(P5.fs_filt, haem.gate)
haem.gates <- lapply(haem.filter, 
                     function(res) filterDetails(res, "CD45")$filter)

# ggcyto(P5.fs_filt[1:6], aes(x=`FSC-A`, y=`B2-A`)) +
#   geom_hex(bins=128) +
#   labs(x="FSC-A", y="CD45") +
#   lims(x=c(0, 2.5e5), y=c(0, 12)) +
#   geom_gate(haem.gate) 

# subset all CD3+ cells before all other analyses
P5.cd45 <- Subset(P5.fs_filt, haem.gates)


# set up a singlets gate
# setup forward and side scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))
singlet.filter <- filter(P5.cd45, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# ggcyto(P5.cd45[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   labs(x="FSC-H", y="FSC-W") +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(singlet.gate)

P5.singlet <- Subset(P5.cd45, singlet.gates)

# set up a side scatter gate requries a covariance matrix
# setup covariance gate
area.points <- c(0, 0, 2.5e5, 2.5e5, 1e4, 0)
height.points <- c(0, 1e4, 2.25e5, 1.75e5, 0, 0)
ssc.gate <- polygonGate(filterId="SSC", .gate=cbind("SSC-H"=height.points, "SSC-A"=area.points))
ssc.filter <- filter(P5.singlet, ssc.gate)
ssc.gates <- lapply(ssc.filter, 
                    function(res) filterDetails(res, "SSC")$filter)
# 
# ggcyto(P5.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   labs(x="SSC-H", y="SSC-A") +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.gate)

P5.ssc <- Subset(P5.singlet, ssc.gates)
# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P5.ssc))){
  events <- nrow(P5.ssc[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P5.ssc <- P5.ssc[-(extract),]}

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P5.ssc)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("B4-A", "R1-A")
third.warp <- c("V2-A", "R2-A")
P5.warping <- warpSet(P5.ssc, stains=setdiff(pars, c(second.warp, third.warp)))
P5.warping <- warpSet(P5.warping, stains=second.warp, peakNr=2)
P5.warping <- warpSet(P5.warping, stains=third.warp, bwFac = 1, nbreaks=20, peakNr=2)

# densityplot(x=~`V2-A`,
#             data=P5.warping,
#             channel=c("V2-A"))

# set up gating scheme using a GatingSet object
P5.gs <- GatingSet(P5.warping)

rm(list=c("P5.warping", "P5.flow.data", "P5.fs_filt",
          "P5.ssc"))
gc()

## First define CD19+ B cells and non-Bcells
cd19.points <- c(8, 12, 12, 8, 8)
cd16.points <- c(0, 0, 9, 9, 0)
cd19.gate <- polygonGate(filterId="Bcells", .gate=cbind("R1-A"=cd19.points, "B3-A"=cd16.points))

cd19.points <- c(0, 8, 8, 12, 12, 0, 0)
cd16.points <- c(0, 0, 9, 9, 12, 12, 0)
notcd19.gate <- polygonGate(filterId="notBcells", .gate=cbind("R1-A"=cd19.points, "B3-A"=cd16.points))

# add the B cell gates
flowWorkspace::add(P5.gs, cd19.gate, parent="root")
recompute(P5.gs)

# add the B cell gates
flowWorkspace::add(P5.gs, notcd19.gate, parent="root")
recompute(P5.gs)

## Define T cells CD3+
cd3.points <- c(8.5, 8.5, 12.5, 12.5, 8.5)
cd45.points <- c(10, 12.5, 12.5, 10, 10)
cd3.gate <- polygonGate(filterId="CD3", .gate=cbind("V1-A"=cd3.points, "B2-A"=cd45.points))

cd3.points <- c(0, 0, 7.5, 7.5, 7.5, 0)
cd45.points <- c(7, 12.5, 12.5, 12.5, 7, 7)
notcd3.gate <- polygonGate(filterId="notCD3", .gate=cbind("V1-A"=cd3.points, "B2-A"=cd45.points))

# add the B cell gates
flowWorkspace::add(P5.gs, cd3.gate, parent="notBcells")
recompute(P5.gs)

# add the B cell gates
flowWorkspace::add(P5.gs, notcd3.gate, parent="notBcells")
recompute(P5.gs)

##########################
## Major T cell subsets ##
##########################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(8, 12, 12, 8, 8)
cd4.points <- c(0, 0, 7.5, 7.5, 0)
cd8.gate <- polygonGate(filterId="CD8", .gate=cbind("B4-A"=cd8.points, "R2-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(8, 8, 0, 0, 8)
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
# ggcyto(P5.gs[1:6], aes(x=`R2-A`, y=`B4-A`), subset="CD3") +
#   geom_hex(bins=256) +
#   labs(x="CD4", y="CD8") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cd8.gate) +
#   geom_gate(cd4.gate) +
#   geom_gate(DN.gate) +
#   geom_gate(DP.gate)

# add the T cell gates
flowWorkspace::add(P5.gs, cd4.gate, parent="CD3")
recompute(P5.gs)

flowWorkspace::add(P5.gs, cd8.gate, parent="CD3")
recompute(P5.gs)

flowWorkspace::add(P5.gs, DN.gate, parent="CD3")
recompute(P5.gs)

flowWorkspace::add(P5.gs, DP.gate, parent="CD3")
recompute(P5.gs)

############################
# NK cells on CD65 & CD16 ##
############################
## NK cells CD56+CD14-
cd56.points <- c(7.5, 7.5, 9, 12, 12, 7.5)
cd14.points <- c(0, 7.25, 7.75, 7.75, 0, 0)
nkcell.gate <- polygonGate(filterId="NKcell", .gate=cbind("V2-A"=cd14.points, "B1-A"=cd56.points))

## not NK cells CD56-
cd56.points <- c(0, 0, 12, 12, 9, 7.5, 7.5, 0)
cd14.points <- c(0, 12, 12, 7.75, 7.75, 7.25, 0, 0)
notnkcell.gate <- polygonGate(filterId="notNKcell", .gate=cbind("V2-A"=cd14.points, "B1-A"=cd56.points))

# ggcyto(P5.gs[1:6], aes(x=`V2-A`, y=`B1-A`), subset="notCD3") +
#   geom_hex(bins=256) +
#   labs(x="CD14", y="CD56") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(nkcell.gate) +
#   geom_gate(notnkcell.gate)

# add the CD16 bright gate
flowWorkspace::add(P5.gs, nkcell.gate, parent="notCD3")
recompute(P5.gs)

flowWorkspace::add(P5.gs, notnkcell.gate, parent="notCD3")
recompute(P5.gs)

## within NK cell subsets
## CD16+CD56hi
## CD16hiCD56+
cd56.points <- c(0, 0, 10, 9, 7.5, 0)
cd16.points <- c(9, 12, 12, 10, 9, 9)
cd16bright.gate <- polygonGate(filterId="CD16hi", .gate=cbind("B3-A"=cd16.points, "B1-A"=cd56.points))

## CD16+CD56hi
cd56.points <- c(7, 7, 7.5, 9, 10, 11, 11, 7)
cd16.points <- c(0, 9, 9, 10, 12, 12, 0, 0)
cd56bright.gate <- polygonGate(filterId="CD56hi", .gate=cbind("B3-A"=cd16.points, "B1-A"=cd56.points))

# ggcyto(P5.gs[1:6], aes(x=`B3-A`, y=`B1-A`), subset="NKcell") +
#   geom_hex(bins=256) +
#   labs(x="CD16", y="CD56") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cd16bright.gate) +
#   geom_gate(cd56bright.gate)

# add the CD16 bright gate
flowWorkspace::add(P5.gs, cd16bright.gate, parent="NKcell")
recompute(P5.gs)

# add the CD16 bright gate
flowWorkspace::add(P5.gs, cd56bright.gate, parent="NKcell")
recompute(P5.gs)

##############################
## monocytes and neutrophils #
##############################
# small-medium sized CD14+ monocytes
cd16.points <- c(0, 2, 5, 7.5, 9, 10.5, 10.5, 10.5)
ssc.points <- c(0, 1.5e5, 1.75e5, 1.75e5, 1.5e5, 1e5, 7.5e4, 0, 0)
monocyte.gate <- polygonGate(filterId="Mono", .gate=cbind("SSC-A"=ssc.points, "B3-A"=cd16.points))

# large CD14- granulocytes
# small-medium sized CD14+ monocytes
cd16.points <- c(10.5, 9, 9, 12, 12, 10.5)
ssc.points <- c(1e5, 1.5e5, 2.5e5, 2.5e5, 1e5, 1e5)
gran.gate <- polygonGate(filterId="Granulocyte", .gate=cbind("SSC-A"=ssc.points, "B3-A"=cd16.points))

# ggcyto(P5.gs[1:6], aes(x=`SSC-A`, y=`B3-A`), subset="notNKcell") +
#   geom_hex(bins=256) +
#   labs(x="SSC-A", y="CD16") +
#   lims(x=c(0, 2.5e5), y=c(0, 12)) +
#   geom_gate(monocyte.gate) +
#   geom_gate(gran.gate)

# add the monocyte & granulocyte gates
flowWorkspace::add(P5.gs, monocyte.gate, parent="notNKcell")
recompute(P5.gs)

flowWorkspace::add(P5.gs, gran.gate, parent="notNKcell")
recompute(P5.gs)

#####################
## Monocyte subsets #
#####################
# split by CD14/CD16 into classical and non-classical (CD16+)
# classical CD14+CD16lo
cd14.points <- c(8, 12, 12, 10, 8, 8)
cd16.points <- c(0, 0, 9, 9, 7.5, 0)
classical.gate <- polygonGate(filterId="ClassicMono", .gate=cbind("B3-A"=cd16.points, "V2-A"=cd14.points))

# non-classical CD14loCD16+
cd14.points <- c(4, 8, 9.5, 9.5, 4, 4)
cd16.points <- c(8, 8, 10, 12, 12, 8)
nonclassic.gate <- polygonGate(filterId="CD16Mono", .gate=cbind("B3-A"=cd16.points, "V2-A"=cd14.points))

# 
# ggcyto(P5.gs[1:6], aes(x=`V2-A`, y=`B3-A`), subset="Mono") +
#   geom_hex(bins=256) +
#   labs(x="CD14", y="CD16") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(classical.gate) +
#   geom_gate(nonclassic.gate)

flowWorkspace::add(P5.gs, classical.gate, parent="Mono")
recompute(P5.gs)

flowWorkspace::add(P5.gs, nonclassic.gate, parent="Mono")
recompute(P5.gs)

## finally, neutrophils as large CD16+ non-monocytes 
forward.points <- c(7.5e4, 2.5e5, 2.5e5, 7.5e4, 7.5e4)
cd16.points <- c(9, 9, 12, 12, 9, 9)
neutrophil.gate <- polygonGate(filterId="Neutrophil", .gate=cbind("B3-A"=cd16.points, "FSC-A"=forward.points))

# ggcyto(P5.gs[1:6], aes(x=`FSC-A`, y=`B3-A`), subset="Granulocyte") +
#   geom_hex(bins=256) +
#   labs(y="CD16", x="FSC-A") +
#   lims(x=c(0, 2.5e5), y=c(0, 12)) +
#   geom_gate(neutrophil.gate)

flowWorkspace::add(P5.gs, neutrophil.gate, parent="Granulocyte")
recompute(P5.gs)

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ and notNodes nodes
all.nodes <- getNodes(P5.gs)[2:length(getNodes(P5.gs))]

# remove the "not" mid-gates
all.nodes <- all.nodes[!grepl(all.nodes, pattern="CD3$")]
all.nodes <- all.nodes[!grepl(all.nodes, pattern="notBcells$")]
all.nodes <- all.nodes[!grepl(all.nodes, pattern="NKcell$")]
all.nodes <- all.nodes[!grepl(all.nodes, pattern="Granulocyte$")]
all.nodes <- all.nodes[!grepl(all.nodes, pattern="/Mono$")]

samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P5.gs, y=node)
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
P5.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P5.df <- P5.df[!duplicated(P5.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD3", "CD14", "CD56",
                     "CD45", "CD16", "CD8B", "CD19", "CD4")
P5.df <- P5.df[, parm.goi]
colnames(P5.df) <- names(parm.goi)

print("Making meta-data")
celltype <- unlist(lapply(strsplit(rownames(P5.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P5.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P5.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P5.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P5.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
rm(list=c("P5.gs"))
gc()

individuals <- levels(P5.meta$Individual)
celltypes <- levels(P5.meta$CellType)

print("Performing bootstraps")
cell.summary.list <- list()
n.boots <- 100
boot.prop <- 0.66

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P5.df[grepl(rownames(P5.df), pattern=celltype), ]
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
      gc()
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
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel5B_bootstrap.txt",
            quote=FALSE, row.names=FALSE, sep="\t")