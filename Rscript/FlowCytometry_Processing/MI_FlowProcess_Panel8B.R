#! /usr/bin/env Rscript

## MI flow cytometry processing
######################################
### Panel 8 - dendritic cell  panel ##
######################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 8B")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P8.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P08")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P8.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P8.flow.data) <- sample.ids

# subset samples into group B
P8.flow.data <- P8.flow.data[324:646]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P8.flow.data))){
  events <- nrow(P8.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P8.fs_filt <- P8.flow.data[seq(along=P8.flow.data)]
if(length(extract)){P8.fs_filt <- P8.fs_filt[-(extract),]}

# remove channels that aren't interesting
param.desc <-  pData(parameters(P8.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P8.flow.data[[1]])
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

P8.fs_filt <- transform(P8.fs_filt, transform1)
P8.fs_filt <- transform(P8.fs_filt, transform2)
P8.fs_filt <- transform(P8.fs_filt, transform3)
P8.fs_filt <- transform(P8.fs_filt, transform4)
P8.fs_filt <- transform(P8.fs_filt, transform5)
P8.fs_filt <- transform(P8.fs_filt, transform6)
P8.fs_filt <- transform(P8.fs_filt, transform7)
P8.fs_filt <- transform(P8.fs_filt, transform8)

# set up a singlets gates
# setup height and width forward scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))
singlet.filter <- filter(P8.fs_filt, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

P8.singlet <- Subset(P8.fs_filt, singlet.gates)

# set up a side scatter gate requries a covariance matrix
# setup covariance gate
area.points <- c(0, 0, 2.5e5, 2.5e5, 5e4, 0)
height.points <- c(0, 1e4, 2e5, 1.5e5, 0, 0)
ssc.gate <- polygonGate(filterId="SSC", .gate=cbind("SSC-H"=height.points, "SSC-A"=area.points))
ssc.filter <- filter(P8.singlet, ssc.gate)
ssc.gates <- lapply(ssc.filter, 
                    function(res) filterDetails(res, "SSC")$filter)

# ggcyto(P8.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   labs(x="SSC-H", y="SSC-A") +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.gate)

P8.ssc <- Subset(P8.singlet, ssc.gate)

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P8.ssc)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("V1-A", "V2-A")
third.warp <- c("B3-A")
P8.warping <- warpSet(P8.ssc, stains=setdiff(pars, c(second.warp, third.warp)))
P8.warping <- warpSet(P8.warping, stains=second.warp, bwFac = 1.1, nbreaks=20)
P8.warping <- warpSet(P8.warping, stains=third.warp, bwFac = 1.75, nbreaks=20, peakNr=2)

# densityplot(x=~`B3-A`,
#             data=P8.warping,
#             channel=c("B3-A"))

# set up gating scheme using a GatingSet object
P8.gs <- GatingSet(P8.warping)

rm(list=c("P8.flow.data", "P8.fs_filt",
          "P8.ssc", "P8.singlet"))
gc()

# remove cells from the dump+ channels and check viability
live.gate <- rangeGate(x=P8.warping[1:6], filterId="Live", stain="V2-A",
                       positive=FALSE, refLine=4)

flowWorkspace::add(P8.gs, live.gate, parent="root")
recompute(P8.gs)

rm(list=c("P8.warping"))
gc()

# ggcyto(P8.gs[1:6], aes(x=`V2-A`), subset="root") +
#   geom_density() +
#   #geom_hex(bins=128) +
#   #labs(x="CD14", y="HLA-DR") +
#   lims(x=c(0, 12)) +
#   geom_gate(live.gate)

##################################
# # gate on all CD14 and HLA-DR
# monocytes
cd14.points <- c(9.5, 12, 12, 9.5, 9.5)
hla.points <- c(4, 4, 12, 12, 4)
mono.gate <- polygonGate(filterId="Monocytes", .gate=cbind("V1-A"=cd14.points, "B3-A"=hla.points))

# DCs
cd14.points <- c(0, 6.5, 7.5, 8.5, 8.5, 0, 0)
hla.points <- c(7, 7, 7.5, 9, 12, 12, 7.25)
dc.gate <- polygonGate(filterId="DCs", .gate=cbind("V1-A"=cd14.points, "B3-A"=hla.points))

# ggcyto(P8.gs[1:6], aes(x=`V1-A`, y=`B3-A`), subset="Live") +
#   geom_hex(bins=128) +
#   labs(x="CD14", y="HLA-DR") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(mono.gate) +
#   geom_gate(dc.gate)

# add the gates
flowWorkspace::add(P8.gs, mono.gate, parent="Live")
recompute(P8.gs)

# add the gates
flowWorkspace::add(P8.gs, dc.gate, parent="Live")
recompute(P8.gs)

##############################
## Defining DC cell subsets ##
##############################
## DC subsets defined by this panel
# pDC - BDCA2+BDCA4+-
# cDC1 - BDCA1+BDCA3-
# cDC3 - BDCA1-BDCA3+

## plasmacytoid dendritic cells
bdca2.points <- c(7.5, 11, 11, 7.5, 7.5)
bdca4.points <- c(6, 6, 11, 11, 6)
pdc.gate <- polygonGate(filterId="pDC", .gate=cbind("B2-A"=bdca2.points, "R1-A"=bdca4.points))

## select non-pDCs
bdca2.points <- c(0, 7, 7, 0, 0)
bdca4.points <- c(0, 0, 7.5, 7.5, 0)
notpdc.gate <- polygonGate(filterId="NotpDC", .gate=cbind("B2-A"=bdca2.points, "R1-A"=bdca4.points))

# ggcyto(P8.gs[1:6], aes(x=`B2-A`, y=`R1-A`), subset="DCs") +
#   geom_hex(bins=128) +
#   labs(x="BDCA-2", y="BDCA-4") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(pdc.gate) + 
#   geom_gate(notpdc.gate)

# add the gates
flowWorkspace::add(P8.gs, pdc.gate, parent="DCs")
recompute(P8.gs)

flowWorkspace::add(P8.gs, notpdc.gate, parent="DCs")
recompute(P8.gs)

###############################
## cDC subsets
# cDC3
bdca1.points <- c(4.5, 7.5, 7.5, 4.5, 4.5)
bdca3.points <- c(7.75, 7.75, 11, 11, 7.75)
cdc3.gate <- polygonGate(filterId="cDC3", .gate=cbind("B1-A"=bdca1.points, "R2-A"=bdca3.points))

# cDC1
bdca1.points <- c(7.5, 12.5, 12.5, 7.5, 7.5)
bdca3.points <- c(0, 0, 7.5, 7.5, 0)
cdc1.gate <- polygonGate(filterId="cDC1", .gate=cbind("B1-A"=bdca1.points, "R2-A"=bdca3.points))

flowWorkspace::add(P8.gs, cdc3.gate, parent="NotpDC")
recompute(P8.gs)

flowWorkspace::add(P8.gs, cdc1.gate, parent="NotpDC")
recompute(P8.gs)

# ggcyto(P8.gs[1:6], aes(x=`B1-A`, y=`R2-A`), subset="Monocytes") +
#   geom_hex(bins=64) +
#   labs(x="BDCA-1", y="BDCA-3") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cdc3.gate) +
#   geom_gate(cdc1.gate)

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ and notNodes nodes
all.nodes <- getNodes(P8.gs)[3:length(getNodes(P8.gs))]

# remove the "not" mid-gates
all.nodes <- all.nodes[!grepl(all.nodes, pattern="NotpDC$")]
all.nodes <- all.nodes[!grepl(all.nodes, pattern="DCs$")]

samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P8.gs, y=node)
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
P8.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P8.df <- P8.df[!duplicated(P8.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD14", "DUMP", "BDCA1",
                     "BDCA2", "HLADR", "CD86", "BDCA4", "BDCA3")
P8.df <- P8.df[, parm.goi]
colnames(P8.df) <- names(parm.goi)

print("Making meta-data")
celltype <- unlist(lapply(strsplit(rownames(P8.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P8.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P8.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P8.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P8.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
print("Calculating summary statistics")
individuals <- levels(P8.meta$Individual)
celltypes <- levels(P8.meta$CellType)

cell.summary.list <- list()

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P8.df[grepl(rownames(P8.df), pattern=celltype), ]
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
print("Removing records with < 100 cells")
all.cell.df <- all.cell.df[all.cell.df$NCells > 100, ]

write.table(all.cell.df,
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel8B_summary.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

