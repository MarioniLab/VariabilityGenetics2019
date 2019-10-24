#! /usr/bin/env Rscript

## MI flow cytometry processing
################################
### Panel 4 - NK cell subsets ##
################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 4C")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P4.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P04")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P4.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P4.flow.data) <- sample.ids

# subset samples into group C
P4.flow.data <- P4.flow.data[644:967]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P4.flow.data))){
  events <- nrow(P4.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P4.fs_filt <- P4.flow.data[seq(along=P4.flow.data)]
if(length(extract)){P4.fs_filt <- P4.fs_filt[-(extract),]}

# remove channels that aren't interesting
param.desc <-  pData(parameters(P4.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P4.flow.data[[1]])
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

P4.fs_filt <- transform(P4.fs_filt, transform1)
P4.fs_filt <- transform(P4.fs_filt, transform2)
P4.fs_filt <- transform(P4.fs_filt, transform3)
P4.fs_filt <- transform(P4.fs_filt, transform4)
P4.fs_filt <- transform(P4.fs_filt, transform5)
P4.fs_filt <- transform(P4.fs_filt, transform6)
P4.fs_filt <- transform(P4.fs_filt, transform7)
P4.fs_filt <- transform(P4.fs_filt, transform8)

# # gate on NK cells separate from other lineages, i.e. CD3/CD14
# ggcyto(P4.fs_filt[1:6], aes(x=`V1-A`, y=`B2-A`)) +
#   geom_hex(bins=128) +
#   labs(x="CD3/CD14", y="NKp46") +
#   lims(x=c(0, 12), y=c(0, 12))

# setup NK cell selection CD3/CD14- NKp46+
lin.points <- c(0, 6, 7, 7.5, 7.75, 7.75, 7.75, 5, 5, 0)
nkp46.points <- c(7, 7, 7.25, 7.5, 7.75, 8, 11, 11, 11, 11)

nkcell.gate <- polygonGate(filterId="NKCell", .gate=cbind("V1-A"=lin.points, "B2-A"=nkp46.points))
nkcell.filter <- filter(P4.fs_filt, nkcell.gate)
nkcell.gates <- lapply(nkcell.filter, 
                       function(res) filterDetails(res, "NKCell")$filter)

# # gate on NK cells separate from other lineages, i.e. CD3/CD14
# ggcyto(P4.fs_filt[1:6], aes(x=`V1-A`, y=`B2-A`)) +
#   geom_hex(bins=128) +
#   labs(x="CD3/CD14", y="NKp46") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(nkcell.gate)

# subset all CD3+ cells before all other analyses
P4.nkcells <- Subset(P4.fs_filt, nkcell.gates)

# set up a live cell gate
viability.points <- c(0, 0, 8.5, 8.5, 0)
side.points <- c(0, 7.5e4, 7.5e4, 0, 0)
live.gate <- polygonGate(filterId="Live", .gate=cbind("V2-A"=viability.points, "SSC-A"=side.points))
live.filter <- filter(P4.nkcells, live.gate)
live.gates <- lapply(live.filter, 
                     function(res) filterDetails(res, "Live")$filter)

# ggcyto(P4.nkcells[1:6], aes(x=`SSC-A`, y=`V2-A`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 12)) +
#   geom_gate(live.gate)
P4.live <- Subset(P4.nkcells, live.gates)

# set up a singlet gate
height.points <- c(0, 2e5, 2e5, 0, 0)
width.points <- c(1.5e5, 1.5e5, 2.5e5, 2.5e5, 1.5e5)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-W"=width.points, "FSC-H"=height.points))
singlet.filter <- filter(P4.live, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)
# 
# ggcyto(P4.live[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(singlet.gate)

P4.singlet <- Subset(P4.live, singlet.gates)
extract <- numeric(0)
for(x in seq_len(length(P4.singlet))){
  events <- nrow(P4.singlet[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

if(length(extract)){P4.singlet <- P4.singlet[-(extract),]}
#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P4.singlet)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("R1-A")
third.warp <- c("B4-A")
P4.warping <- warpSet(P4.singlet, stains=setdiff(pars, c(second.warp, third.warp)), bwFac=1.5)
P4.warping <- warpSet(P4.warping, stains=second.warp, bwFac = 2.5)
P4.warping <- warpSet(P4.warping, stains=third.warp, bwFac = 1.5, nbreaks=20, peakNr=2)

densityplot(x=~`B4-A`,
            data=P4.warping,
            channel=c("B4-A"))

# set up gating scheme using a GatingSet object
P4.gs <- GatingSet(P4.warping)

rm(list=c("P4.warping", "P4.flow.data", "P4.fs_filt",
          "P4.singlet", "P4.live", "P4.nkcells"))
gc()

########################################
# CD16hi and CD56 brights are defined ##
########################################
## CD16hiCD56+
cd56.points <- c(0, 0, 10, 9, 7.5, 0)
cd16.points <- c(9, 12, 12, 10, 9, 9)
cd16bright.gate <- polygonGate(filterId="CD16hi", .gate=cbind("R1-A"=cd16.points, "R2-A"=cd56.points))

# add the CD16 bright gate
flowWorkspace::add(P4.gs, cd16bright.gate, parent="root")
recompute(P4.gs)

## CD16+CD56hi
cd56.points <- c(7, 7, 7.5, 9, 10, 11, 11, 7)
cd16.points <- c(0, 9, 9, 10, 12, 12, 0, 0)
cd56bright.gate <- polygonGate(filterId="CD56hi", .gate=cbind("R1-A"=cd16.points, "R2-A"=cd56.points))

# add the CD56 bright gate
flowWorkspace::add(P4.gs, cd56bright.gate, parent="root")
recompute(P4.gs)

# ggcyto(P4.gs[5:10], aes(x=`R1-A`, y=`R2-A`), subset="root") +
#   geom_hex(bins=256) +
#   labs(x="CD16", y="CD56") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(cd16bright.gate) +
#   geom_gate(cd56bright.gate)


###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P4.gs)[2:length(getNodes(P4.gs))]
samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P4.gs, y=node)
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
P4.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P4.df <- P4.df[!duplicated(P4.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "Lin", "Viability", "CD69",
                     "NKp46", "HLADR", "CD8A", "CD16", "CD56")
P4.df <- P4.df[, parm.goi]
colnames(P4.df) <- names(parm.goi)

print("Making meta-data")
celltype <- unlist(lapply(strsplit(rownames(P4.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P4.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P4.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P4.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P4.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
print("Calculating summary statistics")
individuals <- levels(P4.meta$Individual)
celltypes <- levels(P4.meta$CellType)

cell.summary.list <- list()

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P4.df[grepl(rownames(P4.df), pattern=celltype), ]
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
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel4C_summary.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

