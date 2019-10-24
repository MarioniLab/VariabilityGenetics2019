#! /usr/bin/env Rscript

## TwinsUK flow cytometry processing
#################################
### Panel 7 - Dendritic cells  ##
#################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/fcs_files/"
P7.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P7")

twin.ids <- unlist(lapply(strsplit(x=sampleNames(P7.flow.data), split=" ", fixed=TRUE),
                          FUN=function(X) paste("Twin", X[1], X[2], sep="_")))
sampleNames(P7.flow.data) <- twin.ids

# apply compensation
comp.matrix <- keyword(P7.flow.data[[1]])$`$SPILLOVER`
P7.flow.data <- compensate(P7.flow.data, spillover <- comp.matrix)
P7.flow.data <- P7.flow.data[303:646]

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
transform1 <- transformList("B515-A", biexpTrans)
transform2 <- transformList("B710-A", biexpTrans)
transform3 <- transformList("V450-A", biexpTrans)
transform4 <- transformList("V545-A", biexpTrans)
transform5 <- transformList("V565-A", biexpTrans)
transform6 <- transformList("V585-A", biexpTrans)
transform7 <- transformList("V605-A", biexpTrans)
transform8 <- transformList("V655-A", biexpTrans)
transform9 <- transformList("V705-A", biexpTrans)
transform10 <- transformList("V800-A", biexpTrans)
transform11 <- transformList("R660-A", biexpTrans)
transform12 <- transformList("R710-A", biexpTrans)
transform13 <- transformList("R780-A", biexpTrans)
transform14 <- transformList("G560-A", biexpTrans)
transform15 <- transformList("G610-A", biexpTrans)
transform16 <- transformList("G660-A", biexpTrans)
transform17 <- transformList("G710-A", biexpTrans)
transform18 <- transformList("G780-A", biexpTrans)
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
P7.fs_filt <- transform(P7.fs_filt, transform9)
P7.fs_filt <- transform(P7.fs_filt, transform10)
P7.fs_filt <- transform(P7.fs_filt, transform11)
P7.fs_filt <- transform(P7.fs_filt, transform12)
P7.fs_filt <- transform(P7.fs_filt, transform13)
P7.fs_filt <- transform(P7.fs_filt, transform14)
P7.fs_filt <- transform(P7.fs_filt, transform15)
P7.fs_filt <- transform(P7.fs_filt, transform16)
P7.fs_filt <- transform(P7.fs_filt, transform17)
P7.fs_filt <- transform(P7.fs_filt, transform18)

# P7.fs_filt <- transform(P7.fs_filt, transformForwardH)
# P7.fs_filt <- transform(P7.fs_filt, transformSideA)
# P7.fs_filt <- transform(P7.fs_filt, transformForwardA)

# # gate on singlets
# ggcyto(P7.fs_filt[1:6], aes(x=`FSC-A`, y=`FSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(35000, 2.7e5), y=c(35000, 2.7e5))

# set up a singlets gate
# setup forward and side scatter gates
area.points <- c(3e4, 5e4, 2.3e5, 2.3e5, 3.5e4, 3.2e4)
height.points <- c(3.2e4, 3.2e4, 2e5, 2.5e5, 5e4, 3.2e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-A"=area.points, "FSC-H"=height.points))
singlet.filter <- filter(P7.fs_filt, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# ggcyto(P7.fs_filt[1:6], aes(x=`FSC-A`, y=`FSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(30000, 2.7e5), y=c(30000, 2.7e5)) +
#   geom_gate(singlet.gates)

# subset all singlet cells before all other analyses
P7.singlet <- Subset(P7.fs_filt, singlet.gate)

# select live cells on Side scatter and AqBlue
live.range <- rectangleGate(filterId="Live", list("V545-A"=c(0, 7.5)))
live.filter <- filter(P7.singlet, live.range)
live.gates <- lapply(live.filter, 
                     function(res) filterDetails(res, "Live")$filter)

# ggcyto(P7.singlet[1:6], aes(x=`V545-A`, y=`SSC-A`)) +
#   geom_hex(bins=128)  +
#   geom_gate(live.gates)

# subset all live cells before all other analyses
P7.filter <- Subset(P7.singlet, live.range)

# apply side and forward scatter transformations
P7.filter <- transform(P7.filter, transformForwardH)
P7.filter <- transform(P7.filter, transformSideA)
P7.filter <- transform(P7.filter, transformForwardA)

# ggcyto(P7.filter[1:6], aes(x=`SSC-A`, y=`FSC-A`)) +
#   geom_hex(bins=128) 

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I may need to warp some parameters separately to others
#densityplot(x=~`V605-A`,
#            data=P7.filter[4:9],
#            channel=c("V605-A"))

pars <- colnames(P7.filter)[4:21]
second.warp <- c("V605-A", "G660-A", "B515-A")
first.warp <- pars[!pars %in% second.warp]
P7.warping <- warpSet(P7.filter, stains=first.warp)
P7.warping <- warpSet(P7.warping, stains=second.warp, peakNr=3)

#densityplot(x=~`G710-A`,
#            data=P7.warping[1:9],
#            channel=c("G710-A"))

# set up gating scheme using a GatingSet object
P7.gs <- GatingSet(P7.warping)

rm(list=c("P7.warping", "P7.flow.data", "P7.filter", "P7.singlet",
          "P7.fs_filt"))
gc()

# plot HLA-DR vs lineage markers to select APCs
lin.points <- c(0, 0, 8.5, 8.5, 0)
hla.points <- c(7, 12, 12, 7, 7)
dc.gate <- polygonGate(filterId="DC", .gate=cbind("R710-A"=hla.points, "V605-A"=lin.points))

flowWorkspace::add(P7.gs, dc.gate, parent="root")
recompute(P7.gs)

# dendritic cell groups defined by this panel:
# CD4 myeloid DCs - CD11c+CD123loCD1c+CD141-
# CD8 myeloid DCs - CD11c+CD123loCD1c+CD141+
# Inflammatory CD1c- DC - CD11c+Cd123loCD1c-CD16+
# CD16-CD1c- DC - CD11c+Cd123loCD1c-CD16-
# PDC - CD11c-CD123+
# CD123+CD11c+ DC - CD11c+CD123+

#######
## Remove CD14+ monocytes
#######
# plot HLA-DR vs lineage markers to select APCs
cd14.points <- c(0, 0, 7, 7, 0)
hla.points <- c(7, 12, 12, 7, 7)
mono.gate <- polygonGate(filterId="Neg", .gate=cbind("R710-A"=hla.points, "V800-A"=cd14.points))

flowWorkspace::add(P7.gs, mono.gate, parent="DC")
recompute(P7.gs)

#############################
# myeloid DCs - CD11c+CD123lo
#############################
cd11c.points <- c(10.5, 10.5, 8.5, 10.5, 12, 12)
cd123.points <- c(0, 5.0, 7.5, 10.5, 10.5, 0)
MDC.gate <- polygonGate(filterId="MDC", .gate=cbind("G660-A"=cd11c.points, "G560-A"=cd123.points))

flowWorkspace::add(P7.gs, MDC.gate, parent="Neg")

######################
# PDC - CD11c-CD123+ #
######################
cd11c.points <- c(8.5, 10.5, 8, 7.0, 8.5)
cd123.points <- c(7.5, 10.5, 10.5, 8.5, 7.5)
PDC.gate <- polygonGate(filterId="PDC", .gate=cbind("G660-A"=cd11c.points, "G560-A"=cd123.points))

flowWorkspace::add(P7.gs, PDC.gate, parent="Neg")

##################################
# CD123+CD11c+ DC - CD123+CD11c+ #
##################################
cd11c.points <- c(2, 7.5, 8.5, 2.5)
cd123.points <- c(10, 10, 12, 12)
cd123DC.gate <- polygonGate(filterId="CD123", .gate=cbind("G660-A"=cd11c.points, "G560-A"=cd123.points))

flowWorkspace::add(P7.gs, cd123DC.gate, parent="Neg")
recompute(P7.gs)

##############################################
# CD4 myeloid DCs - CD11c+CD123loCD1c+CD141- #
##############################################
cd1c.points <- c(7.5, 12, 12, 7.5)
cd141.points <- c(0, 0, 6.5, 6.5)
CD4.gate <- polygonGate(filterId="CD4", .gate=cbind("R780-A"=cd1c.points, "B515-A"=cd141.points))

flowWorkspace::add(P7.gs, CD4.gate, parent="MDC")
recompute(P7.gs)

##############################################
# CD8 myeloid DCs - CD11c+CD123loCD1c+CD141+ #
##############################################
cd1c.points <- c(0, 12, 12, 0)
cd141.points <- c(6.5, 6.5, 12, 12)
CD8.gate <- polygonGate(filterId="CD8", .gate=cbind("R780-A"=cd1c.points, "B515-A"=cd141.points))

flowWorkspace::add(P7.gs, CD8.gate, parent="MDC")
recompute(P7.gs)

#################
# Double neg gate
#################
cd1c.points <- c(0, 7, 7, 0)
cd141.points <- c(0, 0, 6.5, 6.5)
Dneg.gate <- polygonGate(filterId="DoubleNeg", .gate=cbind("R780-A"=cd1c.points, "B515-A"=cd141.points))

flowWorkspace::add(P7.gs, Dneg.gate, parent="MDC")
recompute(P7.gs)

###################################################
# Inflammatory CD1c- DC - CD11c+Cd123loCD1c-CD16+ #
###################################################
cd16.points <- c(7.5, 12, 12, 7.5)
cd141.points <- c(0, 0, 6.5, 6.5)
inflam.gate <- polygonGate(filterId="Inflammatory", .gate=cbind("G710-A"=cd16.points, "B515-A"=cd141.points))

flowWorkspace::add(P7.gs, inflam.gate, parent="DoubleNeg")
recompute(P7.gs)

###################################################
# Inflammatory CD1c- DC - CD11c+Cd123loCD1c-CD16+ #
###################################################
cd16.points <- c(0, 7, 7, 0)
cd141.points <- c(0, 0, 6.5, 6.5)
cd16neg.gate <- polygonGate(filterId="CD16neg", .gate=cbind("G710-A"=cd16.points, "B515-A"=cd141.points))

flowWorkspace::add(P7.gs, cd16neg.gate, parent="DoubleNeg")
recompute(P7.gs)


# # plot Cd3 vs forward scatter
# ggcyto(P7.gs[4:9], aes(x=`G710-A`, y=`B515-A`), subset="DoubleNeg") +
#   geom_hex(bins=128) +
#   lims(y=c(0, 12), x=c(0, 12)) +
#   geom_gate(inflam.gate) +
#   geom_gate(cd16neg.gate)


###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD3+, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P7.gs)[2:length(getNodes(P7.gs))]

# also need to ignore the just CCR4-CXCR3- nodes as these are intermediates
all.nodes <- all.nodes[!grepl(all.nodes, pattern="(/DC$)|(Neg$)|(MDC$)")]
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

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "B515-A", "V450-A", "V545-A", "V585-A",
              "V605-A", "V800-A", "R660-A", "R710-A", "R780-A", "G560-A",
              "G610-A", "G660-A", "G710-A", "G780-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD141", "CD83", "AqBlu", "CD8", "Lin",
                     "CD14", "CD32", "HLA.DR", "CD1c", "CD123", "CD64", "CD11c",
                     "CD16", "CD274")
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

# test.df <- P7.df[grepl(rownames(P7.df), pattern="Twin_TS10_264\\._CD3_CD8_Exhausted"), ]
# # forward scatter represents the cell size, so FSC-A^3 should be the volume
# # most parameters are positively correlated with cell size.
# apply(test.df, 2, FUN=function(X) cor(test.df[, "FSC.A"]^3, X))
# apply(test.df, 2, FUN=function(X) cor(test.df[, "SSC.A"]^3, X))
# apply(test.df, 2, FUN=function(X) cor(test.df[, "FSC.H"]^3, X))

individuals <- levels(P7.meta$Individual)
celltypes <- levels(P7.meta$CellType)

cell.summary.list <- list()

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
	file="/mnt/scratcha/jmlab/morgan02/FACS/summary_files/Panel7b_summary.txt",
	quote=FALSE, row.names=FALSE, sep="\t")

# ggplot(all.cell.df, aes(x=Mean, y=Variance, fill=Individual)) +
#   geom_point(shape=21) +
#   theme_mike() + 
#   scale_fill_Publication() +
#   guides(fill=guide_legend(ncol=3))
# 
# ggplot(all.cell.df, aes(x=Mean, y=CV2, fill=Parameters)) +
#   geom_point(shape=21) +
#   theme_mike() + 
#   scale_fill_Publication() +
#   scale_y_log10() +
#   guides(fill=guide_legend(ncol=3))

# there are parameter-specific mean-CV2 relationships
# is this due to the different antibodies and fluourophores used?
# It would be useful to check this for proteins measured in multiple
# panels with different fluoruophores, and for flourescence between
# different panels for smilarly expressed proteins.  My suspicion 
# is that there are flourescence-specific effects as AqBu is a consistent 
# control staining across panels.  This will also be useful for estimating 
# any batch effects between panels and experiments.

# # remove scatter paramters
# ggplot(all.cell.df[!all.cell.df$Parameters %in% c("FSC.A", "SSC.A", "FSC.H"), ], 
#        aes(x=Mean, y=CV2, fill=Parameters)) +
#   geom_point(shape=21, size=3) +
#   theme_mike() + 
#   scale_fill_Publication() +
#   scale_y_log10() +
#   guides(fill=guide_legend(ncol=3, size=4)) +
#   facet_wrap(~CellType, ncol=4, scales="free")





