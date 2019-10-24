#! /usr/bin/env Rscript

## MI flow cytometry processing
###########################
### Panel 10 - ILC panel ##
###########################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 10A")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P10.flow.data <- read.flowSet(path=path.to.files,
                              transformation=FALSE, pattern="P10")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P10.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P10.flow.data) <- sample.ids

# subset samples into group A
P10.flow.data <- P10.flow.data[1:324]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P10.flow.data))){
  events <- nrow(P10.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P10.fs_filt <- P10.flow.data[seq(along=P10.flow.data)]
if(length(extract)){P10.fs_filt <- P10.fs_filt[-(extract),]}

# remove channels that aren't interesting
param.desc <-  pData(parameters(P10.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P10.flow.data[[1]])
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

P10.fs_filt <- transform(P10.fs_filt, transform1)
P10.fs_filt <- transform(P10.fs_filt, transform2)
P10.fs_filt <- transform(P10.fs_filt, transform3)
P10.fs_filt <- transform(P10.fs_filt, transform4)
P10.fs_filt <- transform(P10.fs_filt, transform5)
P10.fs_filt <- transform(P10.fs_filt, transform6)
P10.fs_filt <- transform(P10.fs_filt, transform7)
P10.fs_filt <- transform(P10.fs_filt, transform8)

# # gate on all small live cells
side.points <- c(0, 7e4, 7e4, 0, 0)
viable.points <- c(0, 0, 8.5, 8.5, 0)
live.gate <- polygonGate(filterId="Live", .gate=cbind("V2-A"=viable.points, "SSC-A"=side.points))
live.filter <- filter(P10.fs_filt, live.gate)
live.gates <- lapply(live.filter, 
                     function(res) filterDetails(res, "Live")$filter)

# ggcyto(P10.fs_filt[1:6], aes(x=`SSC-A`, y=`V2-A`)) +
#   geom_hex(bins=128) +
#   labs(x="SSC-A", y="Viability") +
#   lims(x=c(0, 2.5e5), y=c(0, 12)) +
#   geom_gate(live.gate)

P10.live <- Subset(P10.fs_filt, live.gates)

# set up a singlets gates
# setup height and width forward scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))
singlet.filter <- filter(P10.live, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

#ggcyto(P10.live[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#  geom_hex(bins=128) +
#  labs(x="FSC-H", y="FSC-W") +
#  lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#  geom_gate(singlet.gate)

P10.singlet <- Subset(P10.live, singlet.gates)

# set up a side scatter gate requries a covariance matrix
# setup covariance gate
area.points <- c(0, 0, 1.5e5, 1.5e5, 5e4, 0)
height.points <- c(0, 1e4, 1.5e5, 1e5, 0, 0)
ssc.gate <- polygonGate(filterId="SSC", .gate=cbind("SSC-H"=height.points, "SSC-A"=area.points))
ssc.filter <- filter(P10.singlet, ssc.gate)
ssc.gates <- lapply(ssc.filter, 
                    function(res) filterDetails(res, "SSC")$filter)

# ggcyto(P10.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   labs(x="SSC-H", y="SSC-A") +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.gate)

P10.ssc <- Subset(P10.singlet, ssc.gates)
extract <- numeric(0)
for(x in seq_len(length(P10.ssc))){
  events <- nrow(P10.ssc[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P10.ssc <- P10.ssc[-(extract),]}

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P10.ssc)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("B1-A")
P10.warping <- warpSet(P10.ssc, stains=setdiff(pars, second.warp))
P10.warping <- warpSet(P10.warping, stains=second.warp, bwFac = 1.3, nbreaks=20)

# densityplot(x=~`B1-A`,
#             data=P10.warping,
#             channel=c("B1-A"))

# set up gating scheme using a GatingSet object
P10.gs <- GatingSet(P10.warping)

rm(list=c("P10.flow.data", "P10.fs_filt", "P10.live",
          "P10.ssc", "P10.singlet"))
gc()

#################
# Remove all Lineage+ cells
################
# remove cells from the dump+ channels and check viability
line.gate <- rangeGate(x=P10.warping[1:6], filterId="Lineage", stain="V1-A",
                       positive=FALSE, refLine=6)

# ggcyto(P10.gs[1:6], aes(x=`V1-A`), subset="root") +
#   geom_density() +
#   labs(x="Lineage") +
#   lims(x=c(0, 12)) +
#   geom_gate(line.gate)

flowWorkspace::add(P10.gs, line.gate, parent="root")
recompute(P10.gs)

##########################
## Defining ILC subsets ##
##########################
# ILC types defined by this panel
# ILC1 - CRTh2-CD56-CD117-
# ILC2 - CD161+CRTh2+
# ILC3 - CD56-CD117+CRTh2-
# CD56+ ILC - CD56+CRTh2-

## First define Lin-CD127+ ILCs
lineage.points <- c(0, 0, 7, 7, 0)
cd127.points <- c(6, 12, 12, 6, 6)
ilc.gate <- polygonGate(filterId="ILCs", .gate=cbind("V1-A"=lineage.points, "B1-A"=cd127.points))

flowWorkspace::add(P10.gs, ilc.gate, parent="Lineage")
recompute(P10.gs)

##################
### ILC2 - CRTh2+
crth2.points <- c(7, 12, 12, 7, 7)
cd161.points <- c(7.25, 7.25, 12, 12, 7.25)
ilc2.gate <- polygonGate(filterId="ILC2", .gate=cbind("B2-A"=crth2.points, "B3-A"=cd161.points))

crth2.points <- c(0, 7, 7, 0, 0)
cd161.points <- c(0, 0, 12, 12, 0)
crth2neg.gate <- polygonGate(filterId="CRTh2neg", .gate=cbind("B2-A"=crth2.points, "B3-A"=cd161.points))

flowWorkspace::add(P10.gs, ilc2.gate, parent="ILCs")
recompute(P10.gs)

flowWorkspace::add(P10.gs, crth2neg.gate, parent="ILCs")
recompute(P10.gs)

##################
### ILC1 - CD56-CD117-
cd56.points <- c(0, 0, 6.75, 6.75, 0)
cd117.points <- c(0, 6.5, 6.5, 0, 0)
ilc1.gate <- polygonGate(filterId="ILC1", .gate=cbind("R2-A"=cd56.points, "B4-A"=cd117.points))

##################
### ILC3 - CD56-CD117+
cd56.points <- c(0, 0, 6.75, 6.75, 0)
cd117.points <- c(7, 12, 12, 7, 7)
ilc3.gate <- polygonGate(filterId="ILC3", .gate=cbind("R2-A"=cd56.points, "B4-A"=cd117.points))

##################
### CD56+ILC - CD56+
cd56.points <- c(7, 7, 12, 12, 7)
cd117.points <- c(0, 12, 12, 0, 0)
cd56.gate <- polygonGate(filterId="CD56posILC", .gate=cbind("R2-A"=cd56.points, "B4-A"=cd117.points))

# ggcyto(P10.gs[1:6], aes(x=`B4-A`, y=`R2-A`), subset="CRTh2neg") +
#   geom_hex(bins=128) +
#   labs(x="CD117", y="CD56") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(ilc1.gate) +
#   geom_gate(ilc3.gate) +
#   geom_gate(cd56.gate)


flowWorkspace::add(P10.gs, ilc1.gate, parent="CRTh2neg")
recompute(P10.gs)

flowWorkspace::add(P10.gs, ilc3.gate, parent="CRTh2neg")
recompute(P10.gs)

flowWorkspace::add(P10.gs, cd56.gate, parent="CRTh2neg")
recompute(P10.gs)

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ and notNodes nodes
all.nodes <- getNodes(P10.gs)[4:length(getNodes(P10.gs))]

# remove the "not" mid-gates
all.nodes <- all.nodes[!grepl(all.nodes, pattern="CRTh2neg$")]

samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P10.gs, y=node)
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
P10.df <- do.call(rbind.data.frame,
                  samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P10.df <- P10.df[!duplicated(P10.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "Lin", "Viability", "CD127",
                     "CRTh2", "CD161", "CD117", "NKp44", "CD56")
P10.df <- P10.df[, parm.goi]
colnames(P10.df) <- names(parm.goi)

print("Making meta-data")
celltype <- unlist(lapply(strsplit(rownames(P10.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P10.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P10.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P10.meta <- do.call(cbind.data.frame,
                    list("Sample"=rownames(P10.df),
                         "CellType"=celltype,
                         "Individual"=indiv,
                         "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
rm(list=c("P10.gs"))
gc()

individuals <- levels(P10.meta$Individual)
celltypes <- levels(P10.meta$CellType)

print("Performing bootstraps")
cell.summary.list <- list()
n.boots <- 100
boot.prop <- 0.66

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P10.df[grepl(rownames(P10.df), pattern=celltype), ]
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
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel10A_bootstrap.txt",
            quote=FALSE, row.names=FALSE, sep="\t")