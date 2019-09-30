#! /usr/bin/env Rscript

# This is the script for calculating the permutation null distributions for the variance components analysis
# will permute the family labels within each zyogisty class to generate the null distribution.

library(ggplot2)
library(ggrepel)
library(scales)
library(reshape2)
library(pheatmap)
library(ICC)
library(umx)
library(viridis)
library(stringr)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")
mxOption(NULL,"Default optimizer","SLSQP")

panel.dir <- "~/Dropbox/Noise_FACS/TwinsUK/summary/"
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

twin.meta <- read.table("~/Dropbox/Noise_FACS/TwinsUK/ID_table_to_E919.txt",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)

panel.twins <- merge(panel.df, twin.meta[twin.meta$dataset != "Disc/Repl", ], by='FlowJo.Sample.ID')
panel.twins$Measurement <- paste(panel.twins$Parameters, panel.twins$CellType, sep=":")

# add a variable that is the inverse of the number of cells to adjust in the GWAS
panel.twins$InvNCells <- 1/panel.twins$NCells


#### read in FCGR2A genotypes to adjust as covariates
fcgr2.snps <- read.table("~/CI_filesystem/mnt/nas-data/jmlab/group_folders/morgan02/Noise_genetics/TwinsUK/FCGR2A_snps.raw",
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
n.boot <- 100
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
      
      # permute the MZ and DZ labels at this point
      mean.boot.mat <- matrix(nrow=3, ncol=n.boot)
      noise.boot.mat <- matrix(nrow=3, ncol=n.boot)
      
      for(k in seq_along(1:n.boot)){
        # For the permutation I think I should randomly swap zygosity labels, but do I need to preserve family structure?
        # Therefore I should shuffle the family IDs within the zygosity classes.
        
        ### MZ ###
        mz.data <- noise.data[noise.data$Zygosity %in% c("MZ"), ]

        mz.data <- mz.data[duplicated(mz.data$PublicID_Family) | duplicated(mz.data$PublicID_Family, fromLast=TRUE), 
                           c("Individual", "PublicID_Family", "Z.Mean", "Z.Noise", "Age", "dataset")]
        # shuffle family labels
        mz.data$PublicID_Family <- sample(mz.data$PublicID_Family)
        
        ### DZ ###
        dz.data <- noise.data[noise.data$Zygosity %in% c("DZ"), ]
        # dz.data$Z.Noise <- as.numeric(scale(dz.data$SimpleResiduals, center=TRUE))
        # dz.data$Z.Mean <- as.numeric(scale(dz.data$Mean, center=TRUE))
        
        dz.data <- dz.data[duplicated(dz.data$PublicID_Family) | duplicated(dz.data$PublicID_Family, fromLast=TRUE), 
                           c("Individual", "PublicID_Family", "Z.Mean", "Z.Noise", "Age", "dataset")]
        dz.data$PublicID_Family <- sample(dz.data$PublicID_Family)
        
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
          
          sink("aux")
          twin.ace <- umxACE(selDVs=c("Mean", "Noise"), dzData=dz.matched, mzData=mz.matched, sep="")
          mean.ace <- umxACE(selDVs=c("Mean"), dzData=dz.matched, mzData=mz.matched, sep="")
          noise.ace <- umxACE(selDVs=c("Noise"), dzData=dz.matched, mzData=mz.matched, sep="")
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
      }
      
      meas <- paste(unique(noise.data$Measurement), panels[x], sep="_")
        
      mean.null.df <- as.data.frame(t(mean.boot.mat))
      colnames(mean.null.df) <- c("Mean.A", "Mean.C", "Mean.E")
      mean.null.df$Panel <- panels[x]
      mean.null.df$Parameter <- prot
      mean.null.df$CellType <- cell
        
      noise.null.df <- as.data.frame(t(noise.boot.mat))
      colnames(noise.null.df) <- c("Noise.A", "Noise.C", "Noise.E")
      noise.null.df$Panel <- panels[x]
      noise.null.df$Parameter <- prot
      noise.null.df$CellType <- cell
        
      mean_sem.list[[counter]] <- mean.null.df
      noise_sem.list[[counter]] <- noise.null.df
      counter <- counter + 1
    }
  }
}


null.h2.noise <- do.call(rbind.data.frame,
                         noise_sem.list)

null.h2.mean <- do.call(rbind.data.frame,
                       mean_sem.list)

null.h2.df <- merge(null.h2.noise, null.h2.mean, by=c('Panel', 'Parameter', 'CellType'))

write.table(null.h2.df,
            file="~/Dropbox/Noise_FACS/TwinsUK/heritability/NULL_rCV2-ACE-bootstrap.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

null.h2.df <- read.table("~/Dropbox/Noise_FACS/TwinsUK/heritability/NULL_rCV2-ACE-bootstrap.tsv",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)
sem.h2.df <- read.table("~/Dropbox/Noise_FACS/TwinsUK/heritability/SEM_rCV2-ACE-bootstrap_SE.tsv",
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)

# loop over each trait and record the proportion of permuted variance component values which are <
# the observed variance component value.
noise.pval.list <- list()
mean.pval.list <- list()
ncounter <- 1

for(x in seq_along(panels)){
  p.twins <- panel.boot[panel.boot$Panel %in% panels[x], ]
  
  # split by parameter and cell type
  panel.proteins <- unique(p.twins$Parameters)
  panel.cells <- unique(p.twins$CellType)
  
  for(i in seq_along(panel.proteins)){
    prot <- panel.proteins[i]
    
    for(j in seq_along(panel.cells)){
      cell <- panel.cells[j]
      
      trait <- paste(prot, cell, panels[x], sep="_")
      trait.vars <- sem.h2.df[sem.h2.df$Panel %in% panels[x] &
                                sem.h2.df$Parameter %in% prot &
                                sem.h2.df$CellType %in% cell, ]
      
      null.vars  <- null.h2.df[null.h2.df$Panel %in% panels[x] &
                                 null.h2.df$Parameter %in% prot &
                                 null.h2.df$CellType %in% cell, ]
      # get null and test values
      noise.a <- trait.vars$Noise.A
      null.noise.a <- null.vars$Noise.A
      # add correction as per Phipson et al
      noise.a.pval <- 1 - (sum(noise.a > null.noise.a) + 1)/(length(null.noise.a) + 1)
      
      noise.c <- trait.vars$Noise.C
      null.noise.c <- null.vars$Noise.C
      noise.c.pval <- 1 - (sum(noise.c > null.noise.c) + 1)/(length(null.noise.c) + 1)
      
      noise.e <- trait.vars$Noise.E
      null.noise.e <- null.vars$Noise.E
      noise.e.pval <- 1 - (sum(noise.e > null.noise.e) + 1)/(length(null.noise.e) + 1)
      
      mean.a <- trait.vars$Mean.A
      null.mean.a <- null.vars$Mean.A
      mean.a.pval <- 1 - (sum(mean.a > null.mean.a) + 1)/(length(null.mean.a) + 1)
      
      mean.c <- trait.vars$Mean.C
      null.mean.c <- null.vars$Mean.C
      mean.c.pval <- 1 - (sum(mean.c > null.mean.c) + 1)/(length(null.mean.c) + 1)
      
      mean.e <- trait.vars$Mean.E
      null.mean.e <- null.vars$Mean.E
      mean.e.pval <- 1 - (sum(mean.e > null.mean.e) + 1)/(length(null.mean.e) + 1)
      
      noise.pval.list[[paste0(ncounter)]] <- list("Panel"=panels[x], "Parameter"=prot, "CellType"=cell,
                                                  "Trait"=trait, "Noise.A.P"=noise.a.pval, "Noise.C.P"=noise.c.pval, "Noise.E.P"=noise.e.pval,
                                                  "Mean.A.P"=mean.a.pval, "Mean.C.P"=mean.c.pval, "Mean.E.P"=mean.e.pval)
      ncounter <- ncounter + 1
      
    }
  }
}

pval.df <- do.call(rbind.data.frame,
                   noise.pval.list)

# merge with all variance components
varcomp.df <- merge(sem.h2.df, pval.df, by=c("Panel", "Parameter", "CellType"))
varcomp.df$Trait <- as.character(varcomp.df$Trait)

h2.pval.hist <- ggplot(varcomp.df[!grepl(varcomp.df, pattern="(FSC)|(SSC)|(AqBlu)|(Lin)"), ], aes(x=Noise.A.P)) +
  geom_histogram(colour='black', fill='grey80', bins=50) +
  theme_mike() +
  labs(x=expression(paste("h"^2, " permutation p-value")), y="Count")

ggsave(h2.pval.hist,
       filename="~/Dropbox/Noise_FACS/TwinsUK/plot.dir/H2_permutePval_histogram.pdf",
       height=2.75, width=4.45, useDingbats=FALSE)


sem.h2.melt <- melt(varcomp.df, id.vars=c("Panel", "Parameter", "CellType", "Trait"))

varcomp.cols <- c("#386cb0", "#fdb462", "#7fc97f")
names(varcomp.cols) <- c("E", "C", "A")

sem.h2.melt$Estimate <- gsub(sem.h2.melt$variable, pattern=".[A|C|E]", replacement="")
sem.h2.melt$VarComp <- gsub(sem.h2.melt$variable, pattern="(Mean.)|(Noise.)", replacement="")
sem.h2.melt$VarComp <- gsub(sem.h2.melt$VarComp, pattern=".SE", replacement="")
sem.h2.melt$VarComp <- gsub(sem.h2.melt$VarComp, pattern=".P", replacement="")
sem.h2.melt$VarComp <- factor(sem.h2.melt$VarComp,
                              labels=c("A", "C", "E"),
                              levels=c("A", "C", "E"))
sem.h2.melt$Summary <- "Point"
sem.h2.melt$Summary[grep(sem.h2.melt$variable, pattern="SE")] <- "SE"
sem.h2.melt$Summary[grep(sem.h2.melt$variable, pattern="P")] <- "P"
sem.h2.melt$Trait <- as.character(sem.h2.melt$Trait)

non.zero.traits <- unique(c(varcomp.df$Trait[varcomp.df$Noise.A > 0], varcomp.df$Trait[varcomp.df$Mean.A > 0]))
non.zero.traits <- non.zero.traits[!grepl(non.zero.traits, pattern="(^SSC)|(^FSC)|(^AqBlu)|(^Lin)")]

### plotting
nonzero.h2.box <- ggplot(sem.h2.melt[(sem.h2.melt$value > 0) & 
                     (sem.h2.melt$VarComp %in% c("A")) &
                     (sem.h2.melt$Summary %in% c("Point")), ],
       aes(x=Estimate, y=value)) +
  geom_boxplot(fill='grey80') +
  theme_mike() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24)) +
  labs(x="Trait", y=expression(paste("h"^2)))

ggsave(nonzero.h2.box,
       filename="~/Dropbox/Noise_FACS/TwinsUK/plot.dir/H2_nonzero-boxplot.pdf",
       height=3.25, width=3.25, useDingbats=FALSE)

nonzero.h2.box


