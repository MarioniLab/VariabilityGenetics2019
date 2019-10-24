#! /usr/bin/env Rscript

## TwinsUK flow cytometry processing
######################################
### Panel 3 - Helper T cell subsets ##
######################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/fcs_files/"
P3.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P3")

twin.ids <- unlist(lapply(strsplit(x=sampleNames(P3.flow.data), split=" ", fixed=TRUE),
                          FUN=function(X) paste("Twin", X[1], X[2], sep="_")))
sampleNames(P3.flow.data) <- twin.ids

# apply compensation
comp.matrix <- keyword(P3.flow.data[[1]])$`$SPILLOVER`
P3.flow.data <- compensate(P3.flow.data, spillover <- comp.matrix)

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

P3.fs_filt <- transform(P3.fs_filt, transform1)
P3.fs_filt <- transform(P3.fs_filt, transform2)
P3.fs_filt <- transform(P3.fs_filt, transform3)
P3.fs_filt <- transform(P3.fs_filt, transform4)
P3.fs_filt <- transform(P3.fs_filt, transform5)
P3.fs_filt <- transform(P3.fs_filt, transform6)
P3.fs_filt <- transform(P3.fs_filt, transform7)
P3.fs_filt <- transform(P3.fs_filt, transform8)
P3.fs_filt <- transform(P3.fs_filt, transform9)
P3.fs_filt <- transform(P3.fs_filt, transform10)
P3.fs_filt <- transform(P3.fs_filt, transform11)
P3.fs_filt <- transform(P3.fs_filt, transform12)
P3.fs_filt <- transform(P3.fs_filt, transform13)
P3.fs_filt <- transform(P3.fs_filt, transform14)
P3.fs_filt <- transform(P3.fs_filt, transform15)
P3.fs_filt <- transform(P3.fs_filt, transform16)
P3.fs_filt <- transform(P3.fs_filt, transform17)
P3.fs_filt <- transform(P3.fs_filt, transform18)

# P3.fs_filt <- transform(P3.fs_filt, transformForwardH)
# P3.fs_filt <- transform(P3.fs_filt, transformSideA)
# P3.fs_filt <- transform(P3.fs_filt, transformForwardA)

# # gate on singlets
# ggcyto(P3.fs_filt[1:6], aes(x=`FSC-A`, y=`FSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(35000, 2.7e5), y=c(35000, 2.7e5))

# set up a singlets gate
# setup forward and side scatter gates
area.points <- c(3e4, 5e4, 2.3e5, 2.3e5, 3.5e4, 3.2e4)
height.points <- c(3.2e4, 3.2e4, 2e5, 2.5e5, 5e4, 3.2e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-A"=area.points, "FSC-H"=height.points))
singlet.filter <- filter(P3.fs_filt, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# ggcyto(P3.fs_filt[1:6], aes(x=`FSC-A`, y=`FSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(30000, 2.7e5), y=c(30000, 2.7e5)) +
#   geom_gate(singlet.gates)

# subset all singlet cells before all other analyses
P3.singlet <- Subset(P3.fs_filt, singlet.gate)

# select live cells on Side scatter and AqBlue
live.range <- rectangleGate(filterId="Live", list("V545-A"=c(0, 7.5)))
live.filter <- filter(P3.singlet, live.range)
live.gates <- lapply(live.filter, 
                     function(res) filterDetails(res, "Live")$filter)

# ggcyto(P3.singlet[1:6], aes(x=`V545-A`, y=`SSC-A`)) +
#   geom_hex(bins=128)  +
#   geom_gate(live.gates)

# subset all live cells before all other analyses
P3.filter <- Subset(P3.singlet, live.range)

# apply side and forward scatter transformations
P3.filter <- transform(P3.filter, transformForwardH)
P3.filter <- transform(P3.filter, transformSideA)
P3.filter <- transform(P3.filter, transformForwardA)

# ggcyto(P3.filter[1:6], aes(x=`SSC-A`, y=`FSC-A`)) +
#   geom_hex(bins=128) 

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I may need to warp some parameters separately to others
# densityplot(x=~`G780-A`,
#             data=P3.filter[4:9],
#             channel=c("G780-A"))

pars <- colnames(P3.filter)[4:21]
P3.warping <- warpSet(P3.filter, stains=pars)
#P3.warping <- warpSet(P3.warping, stains=second.warp, peakNr=2)

# set up gating scheme using a GatingSet object
P3.gs <- GatingSet(P3.warping)

rm(list=c("P3.warping", "P3.flow.data", "P3.filter", "P3.singlet", 
          "P3.fs_filt"))
gc()

# # plot Cd3 vs side scatter
# ggcyto(P3.gs[1:6], aes(x=`SSC-A`, y=`G610-A`), subset="root") +
#   geom_hex(bins=256) +
#   lims(y=c(0, 12))

# plot CD3 vs side scater to select T cells
tcell.gate <- rectangleGate(filterId="CD3", list("R780-A"=c(7.4, 100)))
flowWorkspace::add(P3.gs, tcell.gate, parent="root")
recompute(P3.gs)

# # # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P3.gs[4:9], aes(x=`V585-A`, y=`V800-A`), subset="CD3") +
#   geom_hex(bins=256) +
#   lims(x=c(0, 12), y=c(0, 12))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7, 12, 12, 7, 7)
cd4.points <- c(2, 2, 9, 7.25, 2)
cd8.gate <- polygonGate(filterId="CD8", .gate=cbind("V585-A"=cd8.points, "V800-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7, 7, 0, 0, 7)
cd4.points <- c(6.5, 11, 11, 6.5, 6.5)
cd4.gate <- polygonGate(filterId="CD4", .gate=cbind("V585-A"=cd8.points, "V800-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(0, 7, 7, 0, 0)
cd4.points <- c(0, 0, 6.5, 6.5, 0)
DN.gate <- polygonGate(filterId="DN", .gate=cbind("V585-A"=cd8.points, "V800-A"=cd4.points))

# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd8.points <- c(7, 12, 12, 7, 7)
cd4.points <- c(7.25, 9, 11, 11, 7.25)
DP.gate <- polygonGate(filterId="DP", .gate=cbind("V585-A"=cd8.points, "V800-A"=cd4.points))

# # plot CD4 Vs CD8 to get major T cell subsets
# ggcyto(P3.gs[4:9], aes(x=`V585-A`, y=`V800-A`), subset="CD3") +
#   geom_hex(bins=256) +
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
# split each in to naive and memory
# TFH - CD4+CXCR5+
# TH1 - CD4+CXCR3+
# TH2 - CD4+CCR4+
# TH9 - CD4+CCR4-CXCR3-CCR6-
# TH17 - CD4+CCR6+CD161+
# TH21 - ???
# TH22 - CD4+CCR10+CCR4+CCR6+


##########################################
# Follicular helper T cells - CD4+CXCR5+ #
##########################################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
cxcr5.points <- c(5, 11, 11, 5, 5)
cd4.points <- c(5, 5, 11, 11, 5)
TFH.gate <- polygonGate(filterId="Tfh", .gate=cbind("R660-A"=cxcr5.points, "V800-A"=cd4.points))

flowWorkspace::add(P3.gs, TFH.gate, parent="/CD3/CD4")
recompute(P3.gs)

############################
# Th1 T cells - CD4+CXCR3+ #
############################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
cxcr3.points <- c(5, 11, 11, 5, 5)
cd4.points <- c(4, 4, 11, 11, 4)
TH1.gate <- polygonGate(filterId="Th1", .gate=cbind("G660-A"=cxcr3.points, "V800-A"=cd4.points))

flowWorkspace::add(P3.gs, TH1.gate, parent="/CD3/CD4")
recompute(P3.gs)

############################
# Th2 T cells - CD4+CCR4+ #
############################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
ccr4.points <- c(8, 13, 13, 8, 8)
cd4.points <- c(0, 0, 10, 10, 0)
TH2.gate <- polygonGate(filterId="Th2", .gate=cbind("G780-A"=ccr4.points, "V800-A"=cd4.points))

flowWorkspace::add(P3.gs, TH2.gate, parent="/CD3/CD4")
recompute(P3.gs)

######################################
# Th9 T cells - CD4+CCR4-CXCR3-CCR6- #
######################################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
ccr4.points <- c(0, 8, 8, 0, 0)
cxcr3.points <- c(0, 0, 5, 5, 0)
ccr4.cxcr3.gate <- polygonGate(filterId="CCR4negCXCR3neg", .gate=cbind("G780-A"=ccr4.points, "G660-A"=cxcr3.points))
flowWorkspace::add(P3.gs, ccr4.cxcr3.gate, parent="/CD3/CD4")
recompute(P3.gs)

ccr4.points <- c(0, 8, 8, 0, 0)
ccr6.points <- c(0, 0, 5, 5, 0)
TH9.gate <- polygonGate(filterId="Th9", .gate=cbind("G780-A"=ccr4.points, "V605-A"=ccr6.points))
flowWorkspace::add(P3.gs, TH9.gate, parent="/CD3/CD4/CCR4negCXCR3neg")
recompute(P3.gs)

##################################
# Th17 T cells - CD4+CCR6+CD161+ #
##################################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
ccr6.points <- c(7.75, 11, 11, 7.75, 7.75)
cd161.points <- c(6.5, 6.5, 11, 11, 6.5)
TH17.gate <- polygonGate(filterId="Th17", .gate=cbind("V605-A"=ccr6.points, "B515-A"=cd161.points))

flowWorkspace::add(P3.gs, TH17.gate, parent="/CD3/CD4")
recompute(P3.gs)

#######################################
# Th22 T cells - CD4+CCR10+CCR4+CCR6+ #
#######################################
# create CCR4+CCR6+ gate
# might need 4 different gates to get the two SP populations, then the DP and DN cells
ccr4.points <- c(8.5, 8.5, 12, 12, 8.5)
ccr10.points <- c(5.5, 11, 11, 5.5, 5.5)
TH22.gate <- polygonGate(filterId="Th22", .gate=cbind("G780-A"=ccr4.points, "G560-A"=ccr10.points))

flowWorkspace::add(P3.gs, TH22.gate, parent="/CD3/CD4")
recompute(P3.gs)


# #plot the remaining T cell markers
# ggcyto(P3.gs[4:9], aes(x=`G560-A`, y=`G780-A`), subset="/CD3/CD4") +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(TH22.gate)
# 
# 
# densityplot(x=~`G560-A`,
#             data=getData(P3.gs[4:9], "/CD3/CD4"),
#             channel=c("G560-A"))



###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD3+, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P3.gs)[3:length(getNodes(P3.gs))]

# also need to ignore the just CCR4-CXCR3- nodes as these are intermediates
#all.nodes <- all.nodes[!grepl(all.nodes, pattern="CCR4negCXCR3neg$")]
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

gc()

# turn this into one massive data set?
P3.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P3.df <- P3.df[!duplicated(P3.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "B515-A", "V450-A", "V545-A", "V585-A",
              "V605-A", "V655-A", "V800-A", "R660-A", "R710-A", "G560-A", "G610-A",
              "G660-A", "G780-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD161", "PD1", "AqBlu", "CD8", "CCR6",
                     "CD45RA", "CD4", "CXCR5", "CCR7", "CCR10", "CD3", "CXCR3", "CCR4")
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

# test.df <- P3.df[grepl(rownames(P3.df), pattern="Twin_TS10_264\\._CD3_CD8_Exhausted"), ]
# # forward scatter represents the cell size, so FSC-A^3 should be the volume
# # most parameters are positively correlated with cell size.
# apply(test.df, 2, FUN=function(X) cor(test.df[, "FSC.A"]^3, X))
# apply(test.df, 2, FUN=function(X) cor(test.df[, "SSC.A"]^3, X))
# apply(test.df, 2, FUN=function(X) cor(test.df[, "FSC.H"]^3, X))

individuals <- levels(P3.meta$Individual)
celltypes <- levels(P3.meta$CellType)

cell.summary.list <- list()

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
# for some reason this filters out all of the different subsets
# I'll have to change it to 50 cells for this panel only.
#all.cell.df <- all.cell.df[all.cell.df$NCells > 50, ]

write.table(all.cell.df,
	file="/mnt/scratcha/jmlab/morgan02/FACS/summary_files/Panel3_summary.txt",
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





