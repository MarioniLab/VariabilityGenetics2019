#! /usr/bin/env Rscript

## TwinsUK flow cytometry processing
######################################
### Panel 4 - Helper T cell subsets ##
######################################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/fcs_files/"
P4.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P4")

twin.ids <- unlist(lapply(strsplit(x=sampleNames(P4.flow.data), split=" ", fixed=TRUE),
                          FUN=function(X) paste("Twin", X[1], X[2], sep="_")))
sampleNames(P4.flow.data) <- twin.ids

# apply compensation
comp.matrix <- keyword(P4.flow.data[[1]])$`$SPILLOVER`
P4.flow.data <- compensate(P4.flow.data, spillover <- comp.matrix)

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

P4.fs_filt <- transform(P4.fs_filt, transform1)
P4.fs_filt <- transform(P4.fs_filt, transform2)
P4.fs_filt <- transform(P4.fs_filt, transform3)
P4.fs_filt <- transform(P4.fs_filt, transform4)
P4.fs_filt <- transform(P4.fs_filt, transform5)
P4.fs_filt <- transform(P4.fs_filt, transform6)
P4.fs_filt <- transform(P4.fs_filt, transform7)
P4.fs_filt <- transform(P4.fs_filt, transform8)
P4.fs_filt <- transform(P4.fs_filt, transform9)
P4.fs_filt <- transform(P4.fs_filt, transform10)
P4.fs_filt <- transform(P4.fs_filt, transform11)
P4.fs_filt <- transform(P4.fs_filt, transform12)
P4.fs_filt <- transform(P4.fs_filt, transform13)
P4.fs_filt <- transform(P4.fs_filt, transform14)
P4.fs_filt <- transform(P4.fs_filt, transform15)
P4.fs_filt <- transform(P4.fs_filt, transform16)
P4.fs_filt <- transform(P4.fs_filt, transform17)
P4.fs_filt <- transform(P4.fs_filt, transform18)

# P4.fs_filt <- transform(P4.fs_filt, transformForwardH)
# P4.fs_filt <- transform(P4.fs_filt, transformSideA)
# P4.fs_filt <- transform(P4.fs_filt, transformForwardA)

# # gate on singlets
# ggcyto(P4.fs_filt[1:6], aes(x=`FSC-A`, y=`FSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(35000, 2.7e5), y=c(35000, 2.7e5))

# set up a singlets gate
# setup forward and side scatter gates
area.points <- c(3e4, 5e4, 2.3e5, 2.3e5, 3.5e4, 3.2e4)
height.points <- c(3.2e4, 3.2e4, 2e5, 2.5e5, 5e4, 3.2e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-A"=area.points, "FSC-H"=height.points))
singlet.filter <- filter(P4.fs_filt, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# ggcyto(P4.fs_filt[1:6], aes(x=`FSC-A`, y=`FSC-H`)) +
#   geom_hex(bins=128) +
#   lims(x=c(30000, 2.7e5), y=c(30000, 2.7e5)) +
#   geom_gate(singlet.gates)

# subset all singlet cells before all other analyses
P4.singlet <- Subset(P4.fs_filt, singlet.gate)

# select live cells on Side scatter and AqBlue
live.range <- rectangleGate(filterId="Live", list("V545-A"=c(0, 7.5)))
live.filter <- filter(P4.singlet, live.range)
live.gates <- lapply(live.filter, 
                     function(res) filterDetails(res, "Live")$filter)

# ggcyto(P4.singlet[1:6], aes(x=`V545-A`, y=`SSC-A`)) +
#   geom_hex(bins=128)  +
#   geom_gate(live.gates)

# subset all live cells before all other analyses
P4.filter <- Subset(P4.singlet, live.range)

# apply side and forward scatter transformations
P4.filter <- transform(P4.filter, transformForwardH)
P4.filter <- transform(P4.filter, transformSideA)
P4.filter <- transform(P4.filter, transformForwardA)

# ggcyto(P4.filter[1:6], aes(x=`SSC-A`, y=`FSC-A`)) +
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
#             data=P4.filter[4:9],
#             channel=c("G780-A"))

pars <- colnames(P4.filter)[4:21]
P4.warping <- warpSet(P4.filter, stains=pars)
#P4.warping <- warpSet(P4.warping, stains=second.warp, peakNr=2)

# set up gating scheme using a GatingSet object
P4.gs <- GatingSet(P4.warping)

rm(list=c("P4.warping", "P4.flow.data", "P4.filter", "P4.singlet", 
          "P4.fs_filt"))
gc()

# # # plot Cd3 vs side scatter
# ggcyto(P4.gs[4:9], aes(x=`SSC-A`, y=`R780-A`), subset="root") +
#   geom_hex(bins=256) +
#   lims(y=c(0, 12))

# plot CD3 vs side scater to select T cells
nk.gate <- rectangleGate(filterId="CD3neg", list("R780-A"=c(0, 7.5)))
flowWorkspace::add(P4.gs, nk.gate, parent="root")
recompute(P4.gs)

# NK cell groups defined by this population:
# split each in to naive and memory
# Early NK - CD56+CD16-
# Mature NK - CD56+CD16+
# Terminal NK - CD56-CD16+


############################
# Early NK cells CD56+CD16-#
############################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd56.points <- c(7.5, 12, 12, 7.5, 7.5)
cd16.points <- c(0, 0, 6.5, 6.5, 0)
Early.gate <- polygonGate(filterId="EarlyNK", .gate=cbind("V450-A"=cd16.points, "V585-A"=cd56.points))

flowWorkspace::add(P4.gs, Early.gate, parent="/CD3neg")
recompute(P4.gs)

#############################
# Mature NK cells CD56+CD16+#
#############################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd56.points <- c(6, 8, 9, 9.5, 8.5, 7.75, 6, 5, 6)
cd16.points <- c(7.5, 7.5, 9, 11.5, 11.5, 11, 9.5, 8)
Mature.gate <- polygonGate(filterId="MatureNK", .gate=cbind("V450-A"=cd16.points, "V585-A"=cd56.points))

flowWorkspace::add(P4.gs, Mature.gate, parent="/CD3neg")
recompute(P4.gs)

#################################
# Terminal NK cells CD56dimCD16+#
#################################
# might need 4 different gates to get the two SP populations, then the DP and DN cells
cd56.points <- c(5, 6, 8, 4.5, 2.5, 2.5, 5)
cd16.points <- c(8, 9.5, 11.5, 11.5, 11.5, 8, 8)
Terminal.gate <- polygonGate(filterId="TerminalNK", .gate=cbind("V450-A"=cd16.points, "V585-A"=cd56.points))

flowWorkspace::add(P4.gs, Terminal.gate, parent="/CD3neg")
recompute(P4.gs)


# #plot the remaining T cell markers
# ggcyto(P4.gs[4:7], aes(x=`V585-A`, y=`V450-A`), subset="/CD3neg") +
#   geom_hex(bins=128) +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(Early.gate) +
#   geom_gate(Mature.gate) +
#   geom_gate(Terminal.gate)
# 
# 
# densityplot(x=~`G560-A`,
#             data=getData(P4.gs[4:9], "/CD3/CD4"),
#             channel=c("G560-A"))

###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD3+, all CD4+ and all CD8+ nodes
all.nodes <- getNodes(P4.gs)[2:length(getNodes(P4.gs))]

# also need to ignore the just CCR4-CXCR3- nodes as these are intermediates
all.nodes <- all.nodes[!grepl(all.nodes, pattern="CD3neg$")]
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

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "B515-A", "V450-A", "V545-A", "V585-A",
              "V605-A", "R660-A", "R710-A", "R780-A", "G560-A", "G610-A", "G660-A",
              "G710-A", "G780-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "CD158a.h", "CD16", "AqBlu", "CD56", "CD4",
                     "CD337", "CD62L", "CD3", "CD158b", "CCR7", "CD335", "CD2", "CD314")
P4.df <- P4.df[, parm.goi]
colnames(P4.df) <- names(parm.goi)

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

# test.df <- P4.df[grepl(rownames(P4.df), pattern="Twin_TS10_264\\._CD3_CD8_Exhausted"), ]
# # forward scatter represents the cell size, so FSC-A^3 should be the volume
# # most parameters are positively correlated with cell size.
# apply(test.df, 2, FUN=function(X) cor(test.df[, "FSC.A"]^3, X))
# apply(test.df, 2, FUN=function(X) cor(test.df[, "SSC.A"]^3, X))
# apply(test.df, 2, FUN=function(X) cor(test.df[, "FSC.H"]^3, X))

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
all.cell.df <- all.cell.df[all.cell.df$NCells > 100, ]

write.table(all.cell.df,
	file="/mnt/scratcha/jmlab/morgan02/FACS/summary_files/Panel4_summary.txt",
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





