#! /usr/bin/env Rscript

## MI flow cytometry processing
###############################
### Panel 6 - B cells  panel ##
###############################
library(flowCore)
library(flowWorkspace)
library(flowStats)
library(openCyto)
library(flowViz)
library(ggcyto)
library(ggplot2)

print("Panel 6A")
path.to.files <- "/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/fcs_files"
#path.to.files <- "~/Dropbox/Noise_FACS/Replication_cohorts/Milieu_Interieur/fcs_files/"
P6.flow.data <- read.flowSet(path=path.to.files,
                             transformation=FALSE, pattern="P06")
sample.ids <- unlist(lapply(strsplit(x=sampleNames(P6.flow.data), split="-", fixed=TRUE),
                            FUN=function(X) paste(X[1], X[2], sep="_")))
sampleNames(P6.flow.data) <- sample.ids

# subset samples into group A
P6.flow.data <- P6.flow.data[1:317]

# remove samples with < 1000 events recorded
extract <- numeric(0)
for(x in seq_len(length(P6.flow.data))){
  events <- nrow(P6.flow.data[[x]])
  if(events < 1000){
    extract <- c(extract, x)
  }
}

P6.fs_filt <- P6.flow.data[seq(along=P6.flow.data)]
if(length(extract)){P6.fs_filt <- P6.fs_filt[-(extract),]}

# remove channels that aren't interesting
param.desc <-  pData(parameters(P6.flow.data[[1]]))[, "desc"]
names(param.desc) <- colnames(P6.flow.data[[1]])
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

P6.fs_filt <- transform(P6.fs_filt, transform1)
P6.fs_filt <- transform(P6.fs_filt, transform2)
P6.fs_filt <- transform(P6.fs_filt, transform3)
P6.fs_filt <- transform(P6.fs_filt, transform4)
P6.fs_filt <- transform(P6.fs_filt, transform5)
P6.fs_filt <- transform(P6.fs_filt, transform6)
P6.fs_filt <- transform(P6.fs_filt, transform7)
P6.fs_filt <- transform(P6.fs_filt, transform8)

# # gate on all CD19+ cells
# ggcyto(P6.fs_filt[1:6], aes(x=`R1-A`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   labs(x="CD19", y="SSC-A") +
#   lims(x=c(0, 12), y=c(0, 2.5e5))

# set up a CD19+ gate
cd19.points <- c(8, 8, 12.5, 12.5, 8)
ssc.points <- c(0, 1.5e5, 1.5e5, 0, 0)

cd19.gate <- polygonGate(filterId="CD19", .gate=cbind("R1-A"=cd19.points, "SSC-A"=ssc.points))
cd19.filter <- filter(P6.fs_filt, cd19.gate)
cd19.gates <- lapply(cd19.filter, 
                     function(res) filterDetails(res, "CD19")$filter)

# # gate on all CD19+ cells
# ggcyto(P6.fs_filt[1:6], aes(x=`R1-A`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   labs(x="CD19", y="SSC-A") +
#   lims(x=c(0, 12), y=c(0, 2.5e5)) +
#   geom_gate(cd19.gate) 

# subset all CD19+ cells before all other analyses
P6.cd19 <- Subset(P6.fs_filt, cd19.gates)

# set up a singlets gates
# setup height and width forward scatter gates
width.points <- c(1.5e5, 1.5e5, 2.2e5, 2.2e5, 1.5e5)
height.points <- c(1e4, 2e5, 2e5, 1e4, 1e4)
singlet.gate <- polygonGate(filterId="Singlet", .gate=cbind("FSC-H"=height.points, "FSC-W"=width.points))
singlet.filter <- filter(P6.cd19, singlet.gate)
singlet.gates <- lapply(singlet.filter, 
                        function(res) filterDetails(res, "Singlet")$filter)

# ggcyto(P6.cd19[1:6], aes(x=`FSC-H`, y=`FSC-W`)) +
#   geom_hex(bins=128) +
#   labs(x="FSC-H", y="FSC-W") +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(singlet.gate)

P6.singlet <- Subset(P6.cd19, singlet.gates)

# set up a side scatter gate requries a covariance matrix
# setup covariance gate
area.points <- c(0, 0, 1.5e5, 1.5e5, 5e4, 0)
height.points <- c(0, 1e4, 1.5e5, 1e5, 0, 0)
ssc.gate <- polygonGate(filterId="SSC", .gate=cbind("SSC-H"=height.points, "SSC-A"=area.points))
ssc.filter <- filter(P6.singlet, ssc.gate)
ssc.gates <- lapply(ssc.filter, 
                    function(res) filterDetails(res, "SSC")$filter)

# ggcyto(P6.singlet[1:6], aes(x=`SSC-H`, y=`SSC-A`)) +
#   geom_hex(bins=128) +
#   labs(x="SSC-H", y="SSC-A") +
#   lims(x=c(0, 2.5e5), y=c(0, 2.5e5)) +
#   geom_gate(ssc.gate)

P6.ssc <- Subset(P6.singlet, ssc.gates)

# double check CD21 & CD19 staining
ggcyto(P6.singlet[1:6], aes(x=`R1-A`, y=`B4-A`)) +
  geom_hex(bins=128) +
  labs(x="CD19", y="CD21") +
  lims(x=c(0, 12), y=c(0, 12))

#################################################
## need to align main parameters before gating ##
#################################################
### some parameters aren't being registered properly
# might need to only select cells with values > 0

# might need to split the parameters to get those with 2 peaks
# properly aligned compare to those with shoulder peaks.
# I need to warp some parameters separately to others
pars <- colnames(P6.ssc)[c(8, 11, 14, 17, 20, 23, 26, 29)]
second.warp <- c("V2-A")
third.warp <- c("B3-A")
P6.warping <- warpSet(P6.ssc, stains=setdiff(pars, second.warp))
P6.warping <- warpSet(P6.warping, stains=second.warp, bwFac = 1.1, nbreaks=20)
P6.warping <- warpSet(P6.warping, stains=third.warp, bwFac = 1.5, nbreaks=20)

# densityplot(x=~`B3-A`,
#             data=P6.warping,
#             channel=c("B3-A"))

# set up gating scheme using a GatingSet object
P6.gs <- GatingSet(P6.warping)

rm(list=c("P6.warping", "P6.flow.data", "P6.fs_filt",
          "P6.ssc", "P6.cd19", "P6.singlet"))
gc()

#############################
## Defining B cell subsets ##
#############################
# B cell types defined by this panel
# CD27+ IgM+
# CD24hi memory B cells
# CD24int memory B cells
# CD24low memory B cells
# Plasmocytes (CD38hi)
# Germinal center (CD38int)
# CD27-IgD- memory
# Naive (CD24+CD38-)
# Transitional (CD24+CD38+)
# Founder (CD24-CD38-)

## First define CD27+ and IgD+ cells
igd.points <- c(7, 12, 12, 7, 7)
cd27.points <- c(8.25, 8.25, 12, 12, 8.25)
doublepos.gate <- polygonGate(filterId="IgDposCD27pos", .gate=cbind("V2-A"=igd.points, "B3-A"=cd27.points))

igd.points <- c(0, 7, 7, 0, 0)
cd27.points <- c(0, 0, 7.5, 7.5, 0)
doubleneg.gate <- polygonGate(filterId="IgDnegCD27neg", .gate=cbind("V2-A"=igd.points, "B3-A"=cd27.points))

igd.points <- c(7, 12, 12, 7, 7)
cd27.points <- c(0, 0, 8.25, 8.25, 0)
igdpos.gate <- polygonGate(filterId="IgDposCD27neg", .gate=cbind("V2-A"=igd.points, "B3-A"=cd27.points))

igd.points <- c(0, 7, 7, 0, 0)
cd27.points <- c(7.5, 7.5, 12, 12, 7.5)
cd27pos.gate <- polygonGate(filterId="IgDnegCD27pos", .gate=cbind("V2-A"=igd.points, "B3-A"=cd27.points))

# add the B cell gates
flowWorkspace::add(P6.gs, doublepos.gate, parent="root")
recompute(P6.gs)

# add the B cell gates
flowWorkspace::add(P6.gs, doubleneg.gate, parent="root")
recompute(P6.gs)

# add the B cell gates
flowWorkspace::add(P6.gs, igdpos.gate, parent="root")
recompute(P6.gs)

# add the B cell gates
flowWorkspace::add(P6.gs, cd27pos.gate, parent="root")
recompute(P6.gs)

# # plot CD27 Vs IgD to get major B cell subsets
# ggcyto(P6.gs[1:6], aes(x=`V2-A`, y=`B3-A`), subset="root") +
#   geom_hex(bins=256) +
#   labs(x="IgD", y="CD27") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(doublepos.gate) +
#   geom_gate(doubleneg.gate) +
#   geom_gate(igdpos.gate) +
#   geom_gate(cd27pos.gate)

####
## Define memory subsets, GC and plasmocytes
####
# CD24 & CD38
# memory B cells CD24+ CD38lo
cd24.points <- c(0, 0, 12, 12, 0)
cd38.points <- c(0, 7.5, 7.5, 0, 0)
memory.gate <- polygonGate(filterId="Memory", .gate=cbind("B2-A"=cd38.points, "R2-A"=cd24.points))

# germinal center B cells CD38int
cd24.points <- c(0, 0, 12, 12, 0)
cd38.points <- c(7.5, 9.5, 9.5, 7.5, 7.5)
gc.gate <- polygonGate(filterId="GC", .gate=cbind("B2-A"=cd38.points, "R2-A"=cd24.points))

# plasmocytes B cells CD38hi
cd24.points <- c(0, 0, 12, 12, 0)
cd38.points <- c(9.5, 12.5, 12.5, 9.5, 9.5)
plasmo.gate <- polygonGate(filterId="Plasmocyte", .gate=cbind("B2-A"=cd38.points, "R2-A"=cd24.points))

flowWorkspace::add(P6.gs, memory.gate, parent="IgDnegCD27pos")
recompute(P6.gs)

flowWorkspace::add(P6.gs, gc.gate, parent="IgDnegCD27pos")
recompute(P6.gs)

flowWorkspace::add(P6.gs, plasmo.gate, parent="IgDnegCD27pos")
recompute(P6.gs)

# ggcyto(P6.gs[1:6], aes(x=`B2-A`, y=`R2-A`), subset="IgDnegCD27pos") +
#   geom_hex(bins=128) +
#   labs(x="CD38", y="CD24") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(memory.gate) +
#   geom_gate(gc.gate) +
#   geom_gate(plasmo.gate)

###
# Memory B cell subsets
###
cd24.points <- c(0, 0, 6, 6, 0)
cd21.points <- c(0, 12, 12, 0, 0)
cd24lo.gate <- polygonGate(filterId="CD24lo", .gate=cbind("B4-A"=cd21.points, "R2-A"=cd24.points))

cd24.points <- c(6, 6, 8, 8, 6)
cd21.points <- c(0, 12, 12, 0, 0)
cd24int.gate <- polygonGate(filterId="CD24int", .gate=cbind("B4-A"=cd21.points, "R2-A"=cd24.points))

cd24.points <- c(8, 8, 12, 12, 8)
cd21.points <- c(0, 12, 12, 0, 0)
cd24hi.gate <- polygonGate(filterId="CD24hi", .gate=cbind("B4-A"=cd21.points, "R2-A"=cd24.points))

flowWorkspace::add(P6.gs, cd24lo.gate, parent="Memory")
recompute(P6.gs)

flowWorkspace::add(P6.gs, cd24int.gate, parent="Memory")
recompute(P6.gs)

flowWorkspace::add(P6.gs, cd24hi.gate, parent="Memory")
recompute(P6.gs)

 # ggcyto(P6.gs[1:6], aes(x=`B4-A`, y=`R2-A`), subset="Memory") +
 #  geom_hex(bins=128) +
 #  labs(x="CD21", y="CD24") +
 #  lims(x=c(0, 12), y=c(0, 12)) +
 #   geom_gate(cd24lo.gate) +
 #   geom_gate(cd24int.gate) +
 #   geom_gate(cd24hi.gate)

######
# Define naive, transitional & founder
######
cd24.points <- c(0, 0, 12, 12, 0)
cd38.points <- c(0, 9.5, 9.5, 0, 0)
naive.gate <- polygonGate(filterId="Naive", .gate=cbind("B2-A"=cd38.points, "R2-A"=cd24.points))

cd24.points <- c(6.5, 9, 12, 12, 6.5)
cd38.points <- c(10, 12, 12, 10, 10)
trans.gate <- polygonGate(filterId="Transitional", .gate=cbind("B2-A"=cd38.points, "R2-A"=cd24.points))

cd24.points <- c(2.5, 2.5, 8.75, 6.25, 2.5)
cd38.points <- c(10, 12, 12, 10, 10)
found.gate <- polygonGate(filterId="Founder", .gate=cbind("B2-A"=cd38.points, "R2-A"=cd24.points))

flowWorkspace::add(P6.gs, naive.gate, parent="IgDposCD27neg")
recompute(P6.gs)

flowWorkspace::add(P6.gs, trans.gate, parent="IgDposCD27neg")
recompute(P6.gs)

flowWorkspace::add(P6.gs, found.gate, parent="IgDposCD27neg")
recompute(P6.gs)

# ggcyto(P6.gs[1:6], aes(x=`B2-A`, y=`R2-A`), subset="IgDposCD27neg") +
#   geom_hex(bins=128) +
#   labs(x="CD38", y="CD24") +
#   lims(x=c(0, 12), y=c(0, 12)) +
#   geom_gate(naive.gate) +
#   geom_gate(trans.gate) +
#   geom_gate(found.gate)


###########################################################################
## Extract single cell expression values for each sample for all markers ##
###########################################################################
# iterate over the nodes and extract expression data for each individual for all
# single cells
# ignore root, all CD4+ and all CD8+ and notNodes nodes
all.nodes <- getNodes(P6.gs)[2:length(getNodes(P6.gs))]

# remove the "not" mid-gates
all.nodes <- all.nodes[!grepl(all.nodes, pattern="IgDposCD27neg$")]
all.nodes <- all.nodes[!grepl(all.nodes, pattern="IgDnegCD27pos$")]

samp.list <- list()
for(i in seq_along(all.nodes)){
  node <- all.nodes[i]
  print(node)
  
  node.flow <- getData(P6.gs, y=node)
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
P6.df <- do.call(rbind.data.frame,
                 samp.list)
# nearly 1.2 million cells!!  Hopefully I've done this so that cells aren't in overlapping categories
# I'll double check for duplicates and remove them
P6.df <- P6.df[!duplicated(P6.df), ]

parm.goi <- c("FSC-A", "SSC-A", "FSC-H", "V1-A", "V2-A", "B1-A", "B2-A",
              "B3-A", "B4-A", "R1-A", "R2-A")
names(parm.goi) <- c("FSC.A", "SSC.A", "FSC.H", "IgM", "IgD", "IgG",
                     "CD38", "CD27", "CD21", "CD19", "CD24")
P6.df <- P6.df[, parm.goi]
colnames(P6.df) <- names(parm.goi)

print("Making meta-data")
celltype <- unlist(lapply(strsplit(rownames(P6.df), split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))
indiv <- unlist(lapply(strsplit(rownames(P6.df), split=".", fixed=TRUE),
                       FUN=function(X) paste0(X[2])))
cell.id <- unlist(lapply(strsplit(rownames(P6.df), split=".", fixed=TRUE),
                         FUN=function(X) paste(X[2:4], collapse=".")))

P6.meta <- do.call(cbind.data.frame,
                   list("Sample"=rownames(P6.df),
                        "CellType"=celltype,
                        "Individual"=indiv,
                        "CellID"=cell.id))

##################################
## Calculate summary statistics ##
##################################
print("Calculating summary statistics")
individuals <- levels(P6.meta$Individual)
celltypes <- levels(P6.meta$CellType)

cell.summary.list <- list()

for(i in seq_along(celltypes)){
  celltype <- celltypes[i]
  celltype.df <- P6.df[grepl(rownames(P6.df), pattern=celltype), ]
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
            file="/mnt/scratcha/jmlab/morgan02/FACS/Milieu_Interieur/summary_files/Panel6A_summary.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

