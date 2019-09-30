#! /usr/bin/env Rscript

# This is the script for calculating the bootstrap SE for the variance components
# I need to use a bootstrap because the variance components analysis uses a Cholesky decomposition to estimate the 
# components, so it is deterministic.

library(ggplot2)
library(ggrepel)
library(scales)
library(reshape2)
library(pheatmap)
library(ICC)
library(umx)
library(viridis)
library(stringr)
source("GGMike/theme_mike.R")
mxOption(NULL,"Default optimizer","SLSQP")

panel.dir <- "TwinsUK/summary/"
panel.files <- list.files(panel.dir, pattern="_residuals")

panel.list <- list()
for(i in seq_along(panel.files)){
  pan <- panel.files[i]
  pan.df <- read.table(paste0(panel.dir, pan),
                       sep="\t", head=TRUE, stringsAsFactors=FALSE)
  pan.df <- pan.df[, c("Variance", "Mean", "CV2", "CV", "Fano", "MeanSize", "SizeCor", "Individual", "NCells",
                       "Parameters", "CellType", "Z.Noise")]
  pan.df$Panel <- unlist(strsplit(pan, split="_", fixed=TRUE))[1]
  panel.list[[pan]] <- pan.df
}

panel.df <- do.call(rbind.data.frame,
                    panel.list)
panel.df$FlowJo.Sample.ID <- unlist(lapply(strsplit(panel.df$Individual, split="_", fixed=TRUE),
                                           FUN=function(X) paste0(X[3])))
panel.df$Parameters[panel.df$Parameters %in% c("HLA.DR")] <- "HLADR"

# read in Twins demographic data
twin.meta <- read.table("TwinsUK/ID_table_to_E919.txt",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)

panel.twins <- merge(panel.df, twin.meta[twin.meta$dataset != "Disc/Repl", ], by='FlowJo.Sample.ID')
panel.twins$Measurement <- paste(panel.twins$Parameters, panel.twins$CellType, sep=":")

# add a variable that is the inverse of the number of cells to adjust in the GWAS
panel.twins$InvNCells <- 1/panel.twins$NCells


#### read in FCGR2A genotypes to adjust as covariates
fcgr2.snps <- read.table("TwinsUK/FCGR2A_snps.raw",
                         header=TRUE, stringsAsFactors=FALSE, sep="\t")

boot.dir <- panel.dir
boot.files <- list.files(boot.dir, pattern="bootstrap-residuals")

boot.list <- list()
for(i in seq_along(boot.files)){
  pan <- boot.files[i]
  pan.df <- read.table(paste0(boot.dir, pan),
                       sep="\t", head=TRUE, stringsAsFactors=FALSE)
  pan.df$boot <- unlist(strsplit(pan, split="_", fixed=TRUE))[1]
  
  boot.list[[pan]] <- pan.df
}

boot.df <- do.call(rbind.data.frame,
                   boot.list)
boot.df$Parameter[boot.df$Parameter %in% c("HLA.DR")] <- "HLADR"
colnames(boot.df) <- c("ResidVar", "CV2Var", "ZVar", "Parameters", "CellType", "Individual", "Panel")

panel.boot <- merge(panel.twins, boot.df, by=c('Parameters', 'CellType', 'Individual', 'Panel'))

umx_set_auto_plot(FALSE)

# for each parameter in each panel and cell type I need to split into MZ and DZ twins
noise_sem.list <- list()
mean_sem.list <- list()

panels <- unique(panel.boot$Panel)
counter <- 1
n.boot <- 10
for(x in seq_along(panels)){
  p.twins <- panel.boot[panel.boot$Panel %in% panels[x], ]
  
  # split by parameter and cell type
  panel.proteins <- unique(p.twins$Parameters)
  panel.cells <- unique(p.twins$CellType)
  
  for(i in seq_along(panel.proteins)){
    prot <- panel.proteins[i]
    
    for(j in seq_along(panel.cells)){
      cell <- panel.cells[j]
      
      noise.data <- p.twins[(p.twins$CellType %in% cell) & (p.twins$Parameters %in% prot), ]
      #noise.data$Z.Noise <- as.numeric(scale(noise.data$SimpleResiduals, center=TRUE))
      noise.data$Z.Mean <- as.numeric(scale(noise.data$Mean, center=TRUE))
      
      # split into MZ and DZ twins, then pair them up based on family ID
      # select all twins that are proper pairs
      # to get standardized coefficients, I need to standardize the noise and mean
      
      ### MZ ###
      mz.data <- noise.data[noise.data$Zygosity %in% c("MZ"), ]
      # mz.data$Z.Noise <- as.numeric(scale(mz.data$SimpleResiduals, center=TRUE))
      # mz.data$Z.Mean <- as.numeric(scale(mz.data$Mean, center=TRUE))
      
      mz.data <- mz.data[duplicated(mz.data$PublicID_Family) | duplicated(mz.data$PublicID_Family, fromLast=TRUE), 
                         c("Individual", "PublicID_Family", "Z.Mean", "Z.Noise", "Age", "dataset")]
      
      ### DZ ###
      dz.data <- noise.data[noise.data$Zygosity %in% c("DZ"), ]
      # dz.data$Z.Noise <- as.numeric(scale(dz.data$SimpleResiduals, center=TRUE))
      # dz.data$Z.Mean <- as.numeric(scale(dz.data$Mean, center=TRUE))
      
      dz.data <- dz.data[duplicated(dz.data$PublicID_Family) | duplicated(dz.data$PublicID_Family, fromLast=TRUE), 
                         c("Individual", "PublicID_Family", "Z.Mean", "Z.Noise", "Age", "dataset")]
      
      # check both MZ and DZ data are available
      # there also needs to be some variance in the trait!!
      if(dim(mz.data)[1] >= 8 & dim(dz.data)[1] >= 8 & all(var(noise.data$Mean) != 0)){
        ### MZ ###
        m1.noise <- mz.data[duplicated(mz.data$PublicID_Family), ]
        colnames(m1.noise) <- c("ID.1", "PublicID_Family", "Mean1", "Noise1", "Age1", "dataset1")
        m2.noise <- mz.data[duplicated(mz.data$PublicID_Family, fromLast=TRUE), ]
        colnames(m2.noise) <- c("ID.2", "PublicID_Family", "Mean2", "Noise2", "Age2", "dataset2")
        
        mz.matched <- merge(m1.noise, m2.noise, by='PublicID_Family')
        mz.matched$dataset1 <- factor(mz.matched$dataset1,
                                      levels=c("Repl", "Disc"))
        mz.matched$dataset2 <- factor(mz.matched$dataset2,
                                      levels=c("Repl", "Disc"))
        
        ### DZ ###
        d1.noise <- dz.data[duplicated(dz.data$PublicID_Family), ]
        colnames(d1.noise) <- c("ID.1", "PublicID_Family", "Mean1", "Noise1", "Age1", "dataset1")
        d2.noise <- dz.data[duplicated(dz.data$PublicID_Family, fromLast=TRUE), ]
        colnames(d2.noise) <- c("ID.2", "PublicID_Family", "Mean2", "Noise2", "Age2", "dataset2")
        
        dz.matched <- merge(d1.noise, d2.noise, by='PublicID_Family')
        dz.matched$dataset1 <- factor(dz.matched$dataset1,
                                      levels=c("Repl", "Disc"))
        dz.matched$dataset2 <- factor(dz.matched$dataset2,
                                      levels=c("Repl", "Disc"))
        
        # we can estimate the components for mean and noise simultaneously, and therefore get their genetic correlation as well
        # umx gives a stupid amount of output, suppress with sink()
        # iterate 10 or 100 times? I need to sample pairs
        
        mean.boot.mat <- matrix(nrow=3, ncol=n.boot)
        noise.boot.mat <- matrix(nrow=3, ncol=n.boot)
        
        for(k in seq_along(c(1:n.boot))){
          mz.sample <- sample(c(1:nrow(mz.matched)), size=ceiling(nrow(mz.matched)*0.75))
          dz.sample <- sample(c(1:nrow(dz.matched)), size=ceiling(nrow(dz.matched)*0.75))
          
          sink("aux")
          twin.ace <- umxACE(selDVs=c("Mean", "Noise"), dzData=dz.matched[dz.sample, ], mzData=mz.matched[mz.sample, ], sep="")
          mean.ace <- umxACE(selDVs=c("Mean"), dzData=dz.matched[dz.sample, ], mzData=mz.matched[mz.sample, ], sep="")
          noise.ace <- umxACE(selDVs=c("Noise"), dzData=dz.matched[dz.sample, ], mzData=mz.matched[mz.sample, ], sep="")
          sink(NULL)
          
          mean.estimates <- parameters(mean.ace, digits=4)
          noise.estimates <- parameters(noise.ace, digits=4)
          twin.estimates <- parameters(twin.ace, digits=4)
          
          # get the variances interms of %variance
          mean.var <- twin.estimates[c(3, 6, 9), 2]/sum(twin.estimates[c(3, 6, 9), 2])
          mean.boot.mat[, k] <- mean.var
          
          noise.var <- twin.estimates[c(5, 8, 11), 2]/sum(twin.estimates[c(5, 8, 11), 2])
          noise.boot.mat[, k] <- noise.var
        }
        meas <- paste(unique(noise.data$Measurement), panels[x], sep="_")
        mean.se <- apply(mean.boot.mat, 1, sd)
        noise.se <- apply(noise.boot.mat, 1, sd)
        
        mean_sem.list[[counter]] <- list("Panel"=panels[x], "Parameter"=prot, "CellType"=cell,
                                         "Mean.A"=mean.var[1], "Mean.A.SE"=mean.se[1],
                                         "Mean.C"=mean.var[2], "Mean.C.SE"=mean.se[2],
                                         "Mean.E"=mean.var[3], "Mean.E.SE"=mean.se[3])
        
        noise_sem.list[[counter]] <- list("Panel"=panels[x], "Parameter"=prot, "CellType"=cell,
                                          "Noise.A"=noise.var[1], "Noise.A.SE"=noise.se[1],
                                          "Noise.C"=noise.var[2], "Noise.C.SE"=noise.se[2],
                                          "Noise.E"=noise.var[3], "Noise.E.SE"=noise.se[3])
        counter <- counter + 1
      }
    }
  }
}


sem.h2.noise <- do.call(rbind.data.frame,
                        noise_sem.list)
colnames(sem.h2.noise) <- c("Panel", "Parameter", "CellType", "Noise.A", "Noise.A.SE", "Noise.C", "Noise.C.SE", "Noise.E", "Noise.E.SE")

sem.h2.mean <- do.call(rbind.data.frame,
                       mean_sem.list)
colnames(sem.h2.mean) <- c("Panel", "Parameter", "CellType", "Mean.A", "Mean.A.SE", "Mean.C", "Mean.C.SE", "Mean.E", "Mean.E.SE")

sem.h2.df <- merge(sem.h2.noise, sem.h2.mean, by=c('Panel', 'Parameter', 'CellType'))

write.table(sem.h2.df,
            file="~/Dropbox/Noise_FACS/TwinsUK/heritability/SEM_rCV2-ACE-bootstrap_SE.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

sem.h2.melt <- melt(sem.h2.df, id.vars=c("Panel", "Parameter", "CellType"))

varcomp.cols <- c("#386cb0", "#fdb462", "#7fc97f")
names(varcomp.cols) <- c("E", "C", "A")

sem.h2.melt$Estimate <- gsub(sem.h2.melt$variable, pattern=".[A|C|E]", replacement="")
sem.h2.melt$VarComp <- gsub(sem.h2.melt$variable, pattern="(Mean.)|(Noise.)", replacement="")
sem.h2.melt$VarComp <- gsub(sem.h2.melt$VarComp, pattern=".SE", replacement="")
sem.h2.melt$VarComp <- factor(sem.h2.melt$VarComp,
                              labels=c("A", "C", "E"),
                              levels=c("A", "C", "E"))
sem.h2.melt$Summary <- "Point"
sem.h2.melt$Summary[grep(sem.h2.melt$variable, pattern="SE")] <- "SE"

# combine cell type and parameter into a single variable for ease
sem.h2.melt$Measurement <- paste(sem.h2.melt$Parameter, sem.h2.melt$CellType, sep=":")

h2.histo <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A") & 
                                 sem.h2.melt$Summary %in% c("SE") , ], 
                   aes(x=value, fill=variable)) +
  geom_histogram(colour='black') +
  theme_mike() +
  facet_wrap(~variable, ncol=3, scales="free") +
  scale_fill_Publication() +
  labs(x=expression(paste("H"^2, " estimate")), y=expression("Frequency"))

# plot the distributions of SE - the errors on the variability traits are larger than the mean
# the errors between components are most similar for the variability than the mean though.
se.box.p <- ggplot(sem.h2.melt[sem.h2.melt$Summary %in% c("SE") &
                                 !sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "Lin", "SSC.A", "AqBlu"), ], 
       aes(x=VarComp, y=value, fill=Estimate)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_Publication(labels=c("Mean", "Variability")) +
  labs(x="Variance Component", y="Standard Error") +
  #guides(fill=guide_legend(label=c("Mean", "Variability"))) +
  geom_jitter(shape=21, position=position_jitterdodge(jitter.width = 0.1))

ggsave(se.box.p,
       filename="plot.dir/StandardErrors_varianceComp-boxplot.png",
       height=3.25, width=4.25, dpi=300)


se.by.cell <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "Lin", "SSC.A", "AqBlu") & 
                     sem.h2.melt$Summary %in% c("SE"), ], 
       aes(x=VarComp, y=value, fill=Estimate)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_Publication(labels=c("Mean", "Variability")) +
  labs(x="Variance Component", y="Standard Error") +
  #guides(fill=guide_legend(label=c("Mean", "Variability"))) +
  facet_wrap(~Parameter) +
  geom_jitter(shape=21, position=position_jitterdodge(jitter.width = 0.1))

ggsave(se.by.cell,
       filename="plot.dir/VarCompSE_byCell-boxplot.png",
       height=7.25, width=7.25, dpi=300)

h2.by.cell <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "Lin", "SSC.A", "AqBlu") & 
                                   sem.h2.melt$Summary %in% c("Point"), ], 
                     aes(x=VarComp, y=value*100, fill=Estimate)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_Publication(labels=c("Mean", "Variability")) +
  labs(x="Variance Component", y="% Phenotypic Variance") +
  facet_wrap(~Parameter) +
  geom_jitter(shape=21, position=position_jitterdodge(jitter.width = 0.1))

ggsave(h2.by.cell,
       filename="plot.dir/VarComp_byCell-boxplot.png",
       height=9.25, width=9.25, dpi=300)



# which h2+SE estimates are in [0, 1]
gt.zero <- sem.h2.df[sem.h2.df$Noise.A > 0 &
                       !sem.h2.df$Parameter %in% c("AqBlu", "Lin", "FSC.A", "FSC.H", "SSC.A"), ]
#lt.one <- gt.zero[gt.zero$Noise.A - gt.zero$Noise.A.SE > 0, ]
#lt.one <- lt.one[lt.one$Noise.A + lt.one$Noise.A.SE < 1, ]
gt.zero <- gt.zero[order(gt.zero$Noise.A, decreasing=TRUE), ]
gt.zero$Index <- c(1:nrow(gt.zero))
gt.zero$Bounds <- "1"
gt.zero$Bounds[gt.zero$Noise.A + gt.zero$Noise.A.SE > 1] <- "0"
gt.zero$Bounds[gt.zero$Noise.A - gt.zero$Noise.A.SE < 0] <- "0"

h2.ord.p <- ggplot(gt.zero) +
  geom_rect(data=data.frame(xmin=-Inf, xmax=Inf,
                            ymin=0, ymax=1), alpha=0.5,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            colour='grey80', fill='grey80') +
  geom_errorbar(data=gt.zero, aes(x=Index, ymin=Noise.A - Noise.A.SE, ymax=Noise.A + Noise.A.SE)) +
  geom_point(data=gt.zero,
             aes(x=Index, y=Noise.A, fill=Bounds),
             shape=21, size=3) +
  scale_y_continuous(limits=c(-0.5, 1.5), oob=squish) +
  theme_mike() +
  scale_fill_manual(values=c("grey80", "orange")) +
  labs(x="Variability Trait", y=expression(bold(paste("h"^2)))) +
  theme(axis.ticks.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank()) +
  guides(fill=guide_legend(title="Within [0, 1]")) +
  NULL


ggsave(h2.ord.p,
       filename="plot.dir/H2_ordered_errorbar.png",
       height=3.95, width=5.95, dpi=300)



