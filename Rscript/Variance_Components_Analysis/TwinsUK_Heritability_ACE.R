#! /usr/bin/env Rscript

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

# data containing processed flow cytometry data
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

# demographic data from TwinsUK
twin.meta <- read.table("TwinsUK/ID_table_to_E919.txt",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)

panel.twins <- merge(panel.df, twin.meta[twin.meta$dataset != "Disc/Repl", ], by='FlowJo.Sample.ID')
panel.twins$Measurement <- paste(panel.twins$Parameters, panel.twins$CellType, sep=":")

# file of Twins with phenotype data for downstream analysis
write.table(panel.twins[!duplicated(panel.twins[, c("PublicID_Family", "PublicID")]), c("PublicID", "PublicID")],
            file="TwinsUK_keep.txt", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# add a variable that is the inverse of the number of cells to adjust in the GWAS
panel.twins$InvNCells <- 1/panel.twins$NCells

# write out a covariates file for all possible covariates
write.table(panel.twins[!duplicated(panel.twins[, c("PublicID_Family", "PublicID")]), c("PublicID", "PublicID", "Age", "InvNCells", "MeanSize", "SizeCor", "dataset")],
            file="TwinsUK_covar.txt", sep="\t",
            row.names=FALSE, col.names=FALSE, quote=FALSE)


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

for(x in seq_along(panels)){
  p.twins <- panel.boot[panel.boot$Panel %in% panels[x], ]
  
  # split by parameter and cell type
  panel.proteins <- unique(p.twins$Parameters)
  panel.proteins <- panel.proteins[!grepl(panel.proteins, pattern="(AqBlu)|(SSC)|(FSC)")]
  panel.cells <- unique(p.twins$CellType)
  
  for(i in seq_along(panel.proteins)){
    prot <- panel.proteins[i]
    
    for(j in seq_along(panel.cells)){
      cell <- panel.cells[j]
      
      noise.data <- p.twins[(p.twins$CellType %in% cell) & (p.twins$Parameters %in% prot), ]
      #noise.data$Z.Noise <- as.numeric(scale(noise.data$SimpleResiduals, center=TRUE))
      
      # the fluorescence has been transformed onto a log scale to account for cell volume
      # transform this back onto the original scale for QTL mapping
      noise.data$OrigMean <- 10**noise.data$Mean
      noise.data$Z.Mean <- as.numeric(scale(noise.data$OrigMean, center=TRUE))
      
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
        sink("aux")
        twin.ace <- umxACE(selDVs=c("Mean", "Noise"), dzData=dz.matched, mzData=mz.matched, sep="", intervals=TRUE, addCI=TRUE)
        mean.ace <- umxACE(selDVs=c("Mean"), dzData=dz.matched, mzData=mz.matched, sep="", intervals=TRUE, addCI=TRUE)
        noise.ace <- umxACE(selDVs=c("Noise"), dzData=dz.matched, mzData=mz.matched, sep="", intervals=TRUE, addCI=TRUE)
        sink(NULL)
        
        mean.estimates <- parameters(mean.ace, digits=4)
        noise.estimates <- parameters(noise.ace, digits=4)
        twin.estimates <- parameters(twin.ace, digits=4)
        
        # get the variances interms of %variance
        mean.var <- twin.estimates[c(3, 6, 9), 2]/sum(twin.estimates[c(3, 6, 9), 2])
        noise.var <- twin.estimates[c(5, 8, 11), 2]/sum(twin.estimates[c(5, 8, 11), 2])
        
        meas <- paste(unique(noise.data$Measurement), panels[x], sep="_")
        if(noise.var[1] > 0){
          noise.data$DataSet <- as.numeric(noise.data$dataset %in% c("Disc"))
          write.table(noise.data[, c("PublicID", "PublicID", "Z.Noise")],
                      file=paste0("TwinsUK/trait_files/Noise_traits-", meas, ".pheno"),
                      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
          
          ###  write out FCGR2A variants as continuous covariates
          noise.fcgr2 <- merge(noise.data, fcgr2.snps, by.x="PublicID", by.y="IID")
          write.table(noise.fcgr2[, c("PublicID", "PublicID", "Age", "DataSet", "rs4657041_T")],
                      file=paste0("TwinsUK/FCGR2A_covar_files/Noise_traits-", meas, ".covar"),
                      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
          
          # write out the matching covariate files
          noise.data$logResidVar <- log10(noise.data$ResidVar)
          write.table(noise.data[, c("PublicID", "PublicID", "Age", "DataSet")],
                      file=paste0("TwinsUK/covar_files/Noise_traits-", meas, ".covar"),
                      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
        }
        
        if(mean.var[1] > 0){
          # transform the mean back onto the original scale - it has been normalised by cell volume, but on a log10 scale
          # write out the matching covariate files
          noise.data$DataSet <- as.numeric(noise.data$dataset %in% c("Disc"))
          noise.fcgr2 <- merge(noise.data, fcgr2.snps, by.x="PublicID", by.y="IID")
          
          write.table(noise.data[, c("PublicID", "PublicID", "Z.Mean")],
                      file=paste0("TwinsUK/trait_files/Mean_traits-", meas, ".pheno"),
                      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
          
          # write out the matching covariate files
          write.table(noise.data[, c("PublicID", "PublicID", "Age", "DataSet")],
                      file=paste0("TwinsUK/covar_files/Mean_traits-", meas, ".covar"),
                      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
          
          ###  write out FCGR2A variants as continuous covariates
          write.table(noise.fcgr2[, c("PublicID", "PublicID", "Age", "DataSet", "rs4657041_T")],
                      file=paste0("TwinsUK/FCGR2A_covar_files/Mean_traits-", meas, ".covar"),
                      quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
          
        }
        
        # calculate the genetic correlations
        if(twin.estimates[3, 2] != 0 & twin.estimates[5, 2] != 0){
          twin.rg <- twin.estimates[4, 2]/sqrt(twin.estimates[3, 2] * twin.estimates[5, 2])
        }
        else{
          twin.rg <- 0
        }
        
        mean_sem.list[[counter]] <- list("Panel"=panels[x], "Parameter"=prot, "CellType"=cell,
                                         "Mean.A"=mean.var[1],
                                         "Mean.C"=mean.var[2],
                                         "Mean.E"=mean.var[3],
                                         "Rg"=twin.rg)
        
        noise_sem.list[[counter]] <- list("Panel"=panels[x], "Parameter"=prot, "CellType"=cell,
                                          "Noise.A"=noise.var[1],
                                          "Noise.C"=noise.var[2],
                                          "Noise.E"=noise.var[3],
                                          "Rg"=twin.rg)
        counter <- counter + 1
      }
    }
  }
}

sem.h2.noise <- do.call(rbind.data.frame,
                        noise_sem.list)
colnames(sem.h2.noise) <- c("Panel", "Parameter", "CellType", "Noise.A", "Noise.C", "Noise.E", "Rg")

sem.h2.mean <- do.call(rbind.data.frame,
                       mean_sem.list)
colnames(sem.h2.mean) <- c("Panel", "Parameter", "CellType", "Mean.A", "Mean.C", "Mean.E", "Rg")

sem.h2.df <- merge(sem.h2.noise, sem.h2.mean, by=c('Panel', 'Parameter', 'CellType'))

write.table(sem.h2.df,
            file="TwinsUK/heritability/SEM_rCV2-ACE.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)


sem.h2.melt <- melt(sem.h2.df, id.vars=c("Panel", "Parameter", "CellType"))

h2.histo <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A") & 
                                 (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ], 
                   aes(x=value, fill=variable)) +
  geom_histogram(colour='black') +
  theme_mike() +
  facet_wrap(~variable, ncol=3, scales="free") +
  scale_fill_Publication() +
  labs(x=expression(paste("H"^2, " estimate")), y=expression("Frequency"))

varcomp.cols <- c("#386cb0", "#fdb462", "#7fc97f")
names(varcomp.cols) <- c("E", "C", "A")

sem.h2.melt$Estimate <- gsub(sem.h2.melt$variable, pattern=".[A|C|E]", replacement="")
sem.h2.melt$VarComp <- gsub(sem.h2.melt$variable, pattern="(Mean.)|(Noise.)", replacement="")
sem.h2.melt$VarComp <- factor(sem.h2.melt$VarComp,
                              labels=c("A", "C", "E"),
                              levels=c("A", "C", "E"))

# combine cell type and parameter into a single variable for ease
sem.h2.melt$Measurement <- paste(sem.h2.melt$Parameter, sem.h2.melt$CellType, sep=":")

h2.panel1.p <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") & 
                                    sem.h2.melt$Panel %in% c("Panel1") &
                                    (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                      aes(x=Measurement, fill=VarComp, y=value)) +
  geom_bar(stat='identity') +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=11),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=11)) +
  facet_grid(~Estimate) +
  coord_flip() +
  scale_fill_manual(values=varcomp.cols)

ggsave(h2.panel1.p,
       filename="plot.dir/Panel1_rCV2_ACE-bar.png",
       height=11.75, width=11.75, dpi=300)

h2.panel2.p <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") & 
                                    sem.h2.melt$Panel %in% c("Panel2") &
                                    (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                      aes(x=Measurement, fill=VarComp, y=value)) +
  geom_bar(stat='identity') +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=11),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=11)) +
  facet_grid(~Estimate) +
  coord_flip() +
  scale_fill_manual(values=varcomp.cols)

ggsave(h2.panel2.p,
       filename="plot.dir/Panel2_rCV2_ACE-bar.png",
       height=11.75, width=11.75, dpi=300)

#### Panel 3
h2.panel3.p <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") & 
                                    sem.h2.melt$Panel %in% c("Panel3") &
                                    (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                      aes(x=Measurement, fill=VarComp, y=value)) +
  geom_bar(stat='identity') +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=11),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=11)) +
  facet_grid(~Estimate) +
  coord_flip() +
  scale_fill_manual(values=varcomp.cols)

ggsave(h2.panel3.p,
       filename="plot.dir/Panel3_rCV2_ACE-bar.png",
       height=11.75, width=11.75, dpi=300)

#### Panel 4
h2.panel4.p <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") & 
                                    sem.h2.melt$Panel %in% c("Panel4") &
                                    (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                      aes(x=Measurement, fill=VarComp, y=value)) +
  geom_bar(stat='identity') +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=11),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=11)) +
  facet_grid(~Estimate) +
  coord_flip() +
  scale_fill_manual(values=varcomp.cols)

ggsave(h2.panel4.p,
       filename="plot.dir/Panel4_rCV2_ACE-bar.png",
       height=11.75, width=11.75, dpi=300)

#### Panel 5
h2.panel5.p <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") & sem.h2.melt$Panel %in% c("Panel5") &
                                    (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                      aes(x=Measurement, fill=VarComp, y=value)) +
  geom_bar(stat='identity') +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=11),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=11)) +
  facet_grid(~Estimate) +
  coord_flip() +
  scale_fill_manual(values=varcomp.cols)

ggsave(h2.panel5.p,
       filename="plot.dir/Panel5_rCV2_ACE-bar.png",
       height=11.75, width=11.75, dpi=300)

#### Panel 6
h2.panel6.p <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") & 
                                    sem.h2.melt$Panel %in% c("Panel6") &
                                    (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                      aes(x=Measurement, fill=VarComp, y=value)) +
  geom_bar(stat='identity') +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=11),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=11)) +
  facet_grid(~Estimate) +
  coord_flip() +
  scale_fill_manual(values=varcomp.cols)

ggsave(h2.panel6.p,
       filename="plot.dir/Panel6_rCV2_ACE-bar.png",
       height=11.75, width=11.75, dpi=300)

#### Panel 7
h2.panel7.p <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") & 
                                    sem.h2.melt$Panel %in% c("Panel7") &
                                    (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                      aes(x=Measurement, fill=VarComp, y=value)) +
  geom_bar(stat='identity') +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=11),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=11)) +
  facet_grid(~Estimate) +
  coord_flip() +
  scale_fill_manual(values=varcomp.cols)

ggsave(h2.panel7.p,
       filename="plot.dir/Panel7_rCV2_ACE-bar.png",
       height=11.75, width=11.75, dpi=300)

### H^2 estimate plotting --------------------------------------------------

h2.compare <- ggplot(sem.h2.melt[!sem.h2.melt$Parameter %in% c("FSC.A", "FSC.H", "SSC.A", "AqBlu", "Lin") &
                                   (!sem.h2.melt$variable %in% c("Rg.x", "Rg.y")), ],
                     aes(x=Estimate, fill=VarComp, y=value*100)) +
  geom_boxplot() +
  #geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge(jitter.width=0.5)) +
  theme_mike() +
  theme(axis.text.x=element_text(vjust=0.5, size=14),
        axis.text.y=element_text(vjust=0.5, hjust=0, size=14)) +
  scale_fill_manual(values=varcomp.cols) +
  scale_x_discrete(labels=c("Mean", "Variability")) +
  labs(x="Trait", y="% Phenotypic Variance")

ggsave(h2.compare,
       filename="plot.dir/H2_rCV2_compare-boxplot.pdf",
       height=4.75, width=4.25*1.618, useDingbats=FALSE)

non.zero.h2 <- sem.h2.melt[sem.h2.melt$VarComp %in% c("A") & sem.h2.melt$value > 0, ]
non.zero.c <- sem.h2.melt[sem.h2.melt$VarComp %in% c("C") & sem.h2.melt$value > 0, ]
non.zero.e <- sem.h2.melt[sem.h2.melt$VarComp %in% c("E") & sem.h2.melt$value > 0, ]

h2.nonzero <- ggplot(non.zero.h2[!non.zero.h2$Parameter %in% c("FSC.A", "SSC.A", "FSC.H", "Lin", "AqBlu"), ], 
                     aes(x=Estimate)) +
  geom_bar(colour='black') +
  labs(x="Trait", y="Count") + 
  scale_x_discrete(labels=c("Mean", "Variability")) +
  theme_mike()

ggsave(h2.nonzero,
       filename="plot.dir/Nonzero_H2_rCV2-bar.pdf",
       height=3.25, width=4.75, useDingbats=FALSE)

# plot h2 estimates for h2>0 traits
non.zero.h2$Estimate[non.zero.h2$Estimate %in% c("Noise")] <- "Variability"
h2.nonzero.box <- ggplot(non.zero.h2[non.zero.h2$VarComp %in% c("A") &
                                       !non.zero.h2$Parameter %in% c("FSC.A", "SSC.A", "FSC.H", "Lin", "AqBlu"),],
                         aes(x=Estimate, y=value)) +
  geom_boxplot(fill='grey80') +
  theme_mike() +
  scale_y_continuous(limits=c(0, 1)) +
  labs(x="Measure", y=expression(paste(bold("h"^2))))

ggsave(h2.nonzero.box,
       filename="plot.dir/Nonzero_H2_rCV2-boxplot.pdf",
       height=3.25, width=3.75, useDingbats=FALSE)


# these following plots illustrate the low genetic correlaiton between mean and variability traits which is a useful sanity check
# that our QTL mapping is unlikely to be biased by any residual mean-variance relationship

h2.rg.hist <- ggplot(sem.h2.df[!sem.h2.df$Parameter %in% c("FSC.A", "FSC.H", "SSC.A") & 
                                 (sem.h2.df$Mean.A != 0) & (sem.h2.df$Noise.A != 0), ],
                     aes(Rg.x)) +
  geom_histogram(colour='black') +
  theme_mike() +
  scale_x_continuous(limits=c(-1, 1), oob=squish) +
  labs(x="Genetic correlation", y="Frequency")

ggsave(h2.rg.hist,
       filename="plot.dir/H2_rCV2-Rg_histogram.png",
       height=3.25, width=4.75, dpi=300)


rg.df <- sem.h2.df[!sem.h2.df$Parameter %in% c("FSC.A", "FSC.H", "SSC.A") & (sem.h2.df$Mean.A != 0) & (sem.h2.df$Noise.A != 0), ]
high.rg <- rg.df[abs(rg.df$Rg.x) > 0.8, ]
low.rg <- rg.df[abs(rg.df$Rg.x) >= 0 & abs(rg.df$Rg.x) <= 0.2, ]

high.rg.p <- ggplot(high.rg[abs(high.rg$Rg.x) <= 1, ],
                    aes(x=Noise.A, y=Mean.A, fill=Rg.x)) +
            geom_point(shape=21, size=3) +
            theme_mike() +
            #guides(fill=FALSE) +
            scale_x_continuous(limits=c(0, 1), oob=squish) +
            scale_fill_viridis() +
            labs(x=expression(paste("Additive H"^2, " Noise")), y=expression(paste("Additive H"^2, " Mean Expression")))

ggsave(high.rg.p,
       filename="plot.dir/Rg_high_rCV2-scatter.png",
       height=4.25, width=4.75, dpi=300)

# This is the concordance between mean and noise $H^2$ for high Rg traits.  Interestingly most of the genetic correlations are < 1, which would imply that selection
# would act in opposite directions, i.e. it might increase noise and decrease mean expression, and _vice versa_.  There are two traits which buck this trend, and
# which have a high positive genetic correlation.

low.rg.p <- ggplot(low.rg[abs(low.rg$Rg.x) <= 1, ],
aes(x=Noise.A, y=Mean.A, fill=Rg.x)) +
geom_point(shape=21, size=3) +
theme_mike() +
scale_x_continuous(limits=c(0, 1), oob=squish) +
#guides(fill=FALSE) +
scale_fill_viridis() +
labs(x=expression(paste("Additive H"^2, " Noise")), y=expression(paste("Additive H"^2, " Mean Expression")))

ggsave(low.rg.p,
       filename="plot.dir/Rg_low_rCV2-scatter.png",
       height=4.25, width=4.75, dpi=300)

# This is the concordance between mean and noise $H^2$ for low Rg traits.  There is a lot more variability in genetic correlations amongst the low correlation traits, 
# perhaps that is to be expected.

# There ~200 proteins with a non-zero additive genetic variance component, and >350 with mean expression non-zero heritability.  That gives me the potential to perform 
# $>$ 550 eQTL analyses!
