##########################################################
## Determine enrichment/depletion of lead SNPs in 
## genomic annotations compared to a background 
## set of SNPs genome-wide that are matched based on 
## the same minor allele frequency distribution
##########################################################
library(biomaRt)
library(ggplot2)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

# link to biomaRt and extract gene information
biomaRt.connection <- useMart("ensembl", "hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")

allgene.df <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "transcript_length",
                                   "start_position", "end_position", "strand", "gene_biotype"),
                    mart=biomaRt.connection)

# filter out non-protein coding genes
allgene.df <- allgene.df[allgene.df$gene_biotype %in% c("protein_coding"), ]
allgene.df <- allgene.df[!duplicated(allgene.df$ensembl_gene_id), ]

# input SNPs and allele frequencies
lead.snps <- read.table("~/Dropbox/Noise_FACS/varQTLs-clump_results.txt",
                        sep="\t", header=TRUE, stringsAsFactors=FALSE)
lead.snps <- lead.snps[lead.snps$SigTier %in% c("Top", "Second"), ]
var.snps <- lead.snps[lead.snps$Measure %in% c("Variability"), ]


# retrieve SNP allele frequences (calculated in 1K genomes)
# only keep SNPs that appear in GeneAtlas results
maf.dir <- "~/CI_filesystem/mnt/nas-data/jmlab/group_folders/morgan02/Noise_genetics/MAF.dir/"
maf.files <- list.files(maf.dir, pattern="frq.gz")
maf.files <- maf.files[!grepl(maf.files, pattern="(SGDP)|(phase_1)|(MT)")]

# I'll also randomly sample a set of background SNPs with the same allele frequency
# I think ~1% per chromosome should suffice
random.snps <- data.frame()

for(f in seq_along(maf.files)){
  # only test the nuclear genome
  if(!grepl(maf.files[f], pattern="MT")){
    f.maf <- read.table(paste0(maf.dir, maf.files[f]), header=TRUE, stringsAsFactors=FALSE)
    
    # randomly sample snps not in the gwas lead set
    n.snps <- nrow(f.maf)
    samp.size <- ceiling(n.snps * 0.01)
    rando.snps <- sample(f.maf$SNP, size=samp.size)
    
    if(length(rando.snps) > 0){
      random.snps <- rbind.data.frame(random.snps, f.maf[f.maf$SNP %in% rando.snps, ])
    }
    
    sink(file="/dev/null")
    rm(list=c("f.maf"))
    gc()
    sink(file=NULL)
    message(paste0("Random SNP MAF samples:", nrow(random.snps)))  
  }
}

# need to sample a set of SNPs that match on MAF from both the GWAS set and genome-wide set

############################################
## Sample a random genome-wide SNP in the same 
## allele frequency range
############################################

varsnps <- unique(var.snps$SNP)
background.match.snps.list <- c()

limit <- 20
for(i in seq_along(varsnps)){
  i.snp <- varsnps[i]
  i.maf <- var.snps$MAF[var.snps$SNP %in% i.snp]
  
  if(!is.na(i.maf)){
    # randomly sample a GWAS hit with the same MAF +/- 10%
    eps.maf <- i.maf * 0.25
    min.maf <- i.maf - eps.maf
    max.maf <- i.maf + eps.maf
    if(min.maf < 0){
      min.maf <- eps.maf
    }
    if(max.maf > 0.5){
      max.maf <- 0.5
    }

    # get genome-wide SNPs with same MAF profile
    bg.max.snps <- random.snps$SNP[random.snps$MAF >= min.maf]
    bg.min.snps <- random.snps$SNP[random.snps$MAF <= max.maf]
    bg.match.snps <- intersect(bg.max.snps, bg.min.snps)
    bg.match.snps <- bg.match.snps[!is.na(bg.match.snps)]
    
    # if there are no SNPs in the range then
    # bump it up to MAF +/- 20%
    if(length(gwas.match.snps) > 0 | length(bg.match.snps)){
      eps.maf <- i.maf * 0.5
      min.maf <- i.maf - eps.maf
      max.maf <- i.maf + eps.maf
      if(min.maf < 0){
        min.maf <- eps.maf
      }
      if(max.maf > 0.5){
        max.maf <- 0.5
      }

      # get genome-wide SNPs with same MAF profile
      bg.max.snps <- random.snps$SNP[random.snps$MAF >= min.maf]
      bg.min.snps <- random.snps$SNP[random.snps$MAF <= max.maf]
      bg.match.snps <- intersect(bg.max.snps, bg.min.snps)
      bg.match.snps <- bg.match.snps[!is.na(bg.match.snps)]
      
      bg.match.snps <- bg.match.snps[!is.na(bg.match.snps)]
    }
    
    bg.new.snp <- TRUE
    # if the SNP is already included, then randomly sample a new one until a new SNP is included
    # if the list doesn't grow after 20 random samples then switch the MAF bin size
    iters <- 0
    while(bg.new.snp){
      if(iters > limit){
        eps.maf <- i.maf * 0.5
        min.maf <- i.maf - eps.maf
        max.maf <- i.maf + eps.maf
        if(min.maf < 0){
          min.maf <- eps.maf
        }
        if(max.maf > 0.5){
          max.maf <- 0.5
        }
        
        # get genome-wide SNPs with same MAF profile
        bg.max.snps <- random.snps$SNP[random.snps$MAF >= min.maf]
        bg.min.snps <- random.snps$SNP[random.snps$MAF <= max.maf]
        bg.match.snps <- intersect(bg.max.snps, bg.min.snps)
        bg.match.snps <- bg.match.snps[!is.na(bg.match.snps)]
        
        bg.match.snps <- bg.match.snps[!is.na(bg.match.snps)]
      }
      
      bg.rand.snp <- sample(bg.match.snps, size=1)
      
      if(sum(bg.rand.snp %in% gwas.match.snps.list) == 0){
        background.match.snps.list <- unique(c(background.match.snps.list, bg.rand.snp))  
        bg.new.snp <- FALSE
      } else if(sum(bg.match.snps %in% background.match.snps.list) == length(bg.match.snps)){
        break
      }
      
      iters <- iters + 1
    }
  }
}

# I should now have two vectors of SNP IDs with the same MAF profile as the var-QTLs
# One for the GWAS lead SNPs
# One for the genome-wide random SNPs
# plot.new()
# hist(var.snps$MAF, breaks=20, col='green')
# hist(random.snps$MAF[random.snps$SNP %in% background.match.snps.list], breaks=20, col='orange', add=TRUE)

# I can now pull out the annotations for these SNP sets
# I'll extract the overlaps for genes, regulome, upstream, downstream
# all annotation files are tabix queried so I can just extract them on
# coordinates - if I had them!!! I can ge these from biomaRt though

# link to biomaRt and extract gene information
biomaRt.snp.connection <- useMart("ENSEMBL_MART_SNP", "hsapiens_snp", host="grch37.ensembl.org", path="/biomart/martservice")

allsnp.df <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end",
                                  "variation_names", "distance_to_transcript", "ensembl_gene_stable_id"),
                   filters="snp_filter", values=unique(c(background.match.snps.list, gwas.match.snps.list, var.snps$SNP)),
                   mart=biomaRt.snp.connection)
# keep snp gene ID unique combinations
allsnp.df <- allsnp.df[!duplicated(allsnp.df[, c("refsnp_id", "ensembl_gene_stable_id")]), ]

annot.dir <- "~/CI_filesystem/mnt/nas-data/jmlab/group_folders/morgan02/Noise_genetics/Enrichment/snp_intersect_bed/"
random.annot.files <- list.files(annot.dir, pattern="ALL.chr")
random.annot.files <- random.annot.files[grepl(random.annot.files, pattern="bgz$")]

# loop over the random SNPs, extract the overlapping regions and count
random.annots <- list()
for(i in seq_along(background.match.snps.list)){
  i.snp <- background.match.snps.list[i]
  i.chr <- unique(allsnp.df$chr_name[allsnp.df$refsnp_id %in% i.snp])
  i.start <- unique(allsnp.df$chrom_start[allsnp.df$refsnp_id %in% i.snp])
  chr.files <- random.annot.files[grepl(random.annot.files, pattern=paste0("chr", i.chr))]
  snp.list <- list()
  # loop over files
  for(x in seq_along(chr.files)){
    chr.file <- chr.files[x]
    annot.type <- unlist(lapply(strsplit(chr.file, split="-", fixed=TRUE),
                                FUN=function(X) paste0(X[-1])))
    annot.type <- gsub(paste(annot.type, collapse="_"), pattern="\\.bed\\.bgz", replacement="")
    
    tabix.query <- paste0("tabix ", paste0(annot.dir, chr.file), " ", i.chr, ":", i.start, "-", i.start+1)
    tabix.call <- system(command=tabix.query, intern=TRUE)
    tabix.out <- do.call(rbind.data.frame,
                         sapply(tabix.call, FUN=function(X) strsplit(X, split="\t", fixed=TRUE)))
    if(nrow(tabix.out) > 0){
      if(ncol(tabix.out) == 9){
        colnames(tabix.out) <- c("SNP.CHR", "SNP.START", "SNP.END", "SNP", 
                                 "ANNOT.CHR", "ANNOT.START", "ANNOT.END", "ANNOT",
                                 "OVERLAP")
      } else if(ncol(tabix.out) == 11){
        colnames(tabix.out) <- c("SNP.CHR", "SNP.START", "SNP.END", "SNP", 
                                 "ANNOT.CHR", "ANNOT.START", "ANNOT.END", "ANNOT", 
                                 "SCORE", "STRAND", "OVERLAP")
      }
      # take the unique lines
      tabix.out <- tabix.out[!duplicated(tabix.out[, c("SNP", "ANNOT", "ANNOT.START")]), ]
      tabix.tab <- as.data.frame(table(tabix.out$ANNOT))
      tabix.tab$Set <- annot.type
      snp.list[[annot.type]] <- tabix.tab
    }
    
    if(length(snp.list)){
      annot.df <- do.call(rbind.data.frame,
                          snp.list)
      colnames(annot.df) <- c("ANNOT", "COUNT", "SET")
      annot.df$SNP <- i.snp
      random.annots[[i.snp]] <- annot.df
    }
  }
}

random.annot.df <- do.call(rbind.data.frame,
                           random.annots)
# drop gene as gene and transcript are redundant annotations
# that represent introns
random.annot.df$ANNOT <- as.character(random.annot.df$ANNOT)

random.annot.df <- random.annot.df[!random.annot.df$ANNOT %in% c("gene"), ]
random.annot.df$Annot <- random.annot.df$ANNOT
random.annot.df$Annot[random.annot.df$SET %in% c("GRCh37_100kb")] <- "100kb.TSS"
random.annot.df$Annot[random.annot.df$SET %in% c("GRCh37_upstream")] <- "Upstream"
random.annot.df$Annot[random.annot.df$SET %in% c("GRCh37_downstream")] <- "Downstream"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("three_prime_utr", "five_prime_utr")] <- "UTR"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("transcript")] <- "Transcript"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("CTCF_binding_site")] <- "CTCF.site"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("promoter")] <- "Promoter"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("promoter_flanking_region")] <- "Upstream"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("exon", "CDS")] <- "Genic"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("open_chromatin_region")] <- "Open.chromatin"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("enhancer")] <- "Enhancer"
random.annot.df$Annot[random.annot.df$ANNOT %in% c("TF_binding_site")] <- "TF.site"

# deduplicate the upstream/downstream/100kb + SNP combo
random.annot.df  <- random.annot.df[!duplicated(random.annot.df[, c("SNP", "Annot")]), ]
random.freq.df <- as.data.frame(table(random.annot.df$Annot)/length(background.match.snps.list))


############################################
## do the same for the variability SNPs
############################################
varqtl.annot.files <- list.files(annot.dir, pattern="ALL_lead")
varqtl.annot.files <- varqtl.annot.files[grepl(varqtl.annot.files, pattern="bgz$")]

# loop over the varqtl SNPs, extract the overlapping regions and count
varqtl.annots <- list()
for(i in seq_along(varsnps)){
  i.snp <- varsnps[i]
  i.chr <- unique(allsnp.df$chr_name[allsnp.df$refsnp_id %in% i.snp])
  i.start <- unique(allsnp.df$chrom_start[allsnp.df$refsnp_id %in% i.snp])
  snp.list <- list()
  # loop over files
  for(x in seq_along(varqtl.annot.files)){
    chr.file <- varqtl.annot.files[x]
    annot.type <- unlist(lapply(strsplit(chr.file, split="-", fixed=TRUE),
                                FUN=function(X) paste0(X[-1])))
    annot.type <- gsub(paste(annot.type, collapse="_"), pattern="\\.bed\\.bgz", replacement="")
    
    tabix.query <- paste0("tabix ", paste0(annot.dir, chr.file), " ", i.chr, ":", i.start, "-", i.start+1)
    tabix.call <- system(command=tabix.query, intern=TRUE)
    tabix.out <- do.call(rbind.data.frame,
                         sapply(tabix.call, FUN=function(X) strsplit(X, split="\t", fixed=TRUE)))
    if(nrow(tabix.out) > 0){
      if(ncol(tabix.out) == 9){
        colnames(tabix.out) <- c("SNP.CHR", "SNP.START", "SNP.END", "SNP", 
                                 "ANNOT.CHR", "ANNOT.START", "ANNOT.END", "ANNOT",
                                 "OVERLAP")
      } else if(ncol(tabix.out) == 11){
        colnames(tabix.out) <- c("SNP.CHR", "SNP.START", "SNP.END", "SNP", 
                                 "ANNOT.CHR", "ANNOT.START", "ANNOT.END", "ANNOT", 
                                 "SCORE", "STRAND", "OVERLAP")
      }
      # take the unique lines
      tabix.out <- tabix.out[!duplicated(tabix.out[, c("SNP", "ANNOT", "ANNOT.START")]), ]
      tabix.tab <- as.data.frame(table(tabix.out$ANNOT))
      tabix.tab$Set <- annot.type
      snp.list[[annot.type]] <- tabix.tab
    }
    
    if(length(snp.list)){
      annot.df <- do.call(rbind.data.frame,
                          snp.list)
      colnames(annot.df) <- c("ANNOT", "COUNT", "SET")
      annot.df$SNP <- i.snp
      varqtl.annots[[i.snp]] <- annot.df
    }
  }
}

varqtl.annot.df <- do.call(rbind.data.frame,
                           varqtl.annots)

# drop gene as gene and transcript are redundant annotations
# that represent introns
varqtl.annot.df$ANNOT <- as.character(varqtl.annot.df$ANNOT)

varqtl.annot.df <- varqtl.annot.df[!varqtl.annot.df$ANNOT %in% c("gene"), ]
varqtl.annot.df$Annot <- varqtl.annot.df$ANNOT
varqtl.annot.df$Annot[varqtl.annot.df$SET %in% c("GRCh37_100kb")] <- "100kb.TSS"
varqtl.annot.df$Annot[varqtl.annot.df$SET %in% c("GRCh37_upstream")] <- "Upstream"
varqtl.annot.df$Annot[varqtl.annot.df$SET %in% c("GRCh37_downstream")] <- "Downstream"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("three_prime_utr", "five_prime_utr")] <- "UTR"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("transcript")] <- "Transcript"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("CTCF_binding_site")] <- "CTCF.site"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("promoter")] <- "Promoter"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("promoter_flanking_region")] <- "Upstream"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("exon", "CDS")] <- "Genic"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("open_chromatin_region")] <- "Open.chromatin"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("enhancer")] <- "Enhancer"
varqtl.annot.df$Annot[varqtl.annot.df$ANNOT %in% c("TF_binding_site")] <- "TF.site"

# deduplicate the upstream/downstream/100kb + SNP combo
varqtl.annot.df  <- varqtl.annot.df[!duplicated(varqtl.annot.df[, c("SNP", "Annot")]), ]
varqtl.freq.df <- as.data.frame(table(varqtl.annot.df$Annot)/length(varsnps))

## aggregate results and compare
varqtl.freq.df$Set <- "VarSNP"
random.freq.df$Set <- "Background"

comp.freq.df <- do.call(rbind.data.frame,
                        list("var"=varqtl.freq.df,
                             "bg"=random.freq.df))

## merge UTR annotations and make names pretty
comp.freq.df$Var1 <- as.character(comp.freq.df$Var1)

vsnp.ensembl <- ggplot(comp.freq.df[comp.freq.df$Set %in% c("VarSNP"),], 
                       aes(x=reorder(Var1, Freq), y=Freq)) +
  geom_bar(position='dodge', stat='identity',
           fill='grey80', colour='black') +
  theme_mike() +
  labs(y="Proportion of VarSNPs", x="Ensembl Annotation") +
  #scale_fill_Publication() +
  guides(fill=FALSE) +
  scale_y_continuous(limit=c(0, 0.4)) +
  coord_flip() 

ggsave(vsnp.ensembl,
       filename="~/Dropbox/Noise_FACS/plot.dir/VarSNP-Ensembl_location-bar.pdf",
       height=3.25, width=4.95, useDingbats=FALSE)

vsnp.ensembl

### fun a fishers exact test between varSNPs and background set
annots <- unique(comp.freq.df$Var1)
fisher.res.list <- list()

bg.annot.tab <- as.data.frame(table(random.annot.df$Annot))
bg.annot.tab$Var1 <- as.character(bg.annot.tab$Var1)

var.annot.tab <- as.data.frame(table(varqtl.annot.df$Annot))
var.annot.tab$Var1 <- as.character(var.annot.tab$Var1)

for(x in seq_along(annots)){
  x.annot <- annots[x]
  bg.counts <- bg.annot.tab[bg.annot.tab$Var1 %in% x.annot, ]
  vq.counts <- var.annot.tab[var.annot.tab$Var1 %in% x.annot, ]

  # setup the contingency table
  # a=var in annot
  # b=vat Not in annot
  # c=bg in annot
  # d=bg Not in annot
  a.var <- vq.counts$Freq
  b.var <- length(varsnps) - a.var

  c <- bg.counts$Freq
  if(length(c) == 0){
    c <- 0
  }
  d <- length(background.match.snps.list) - c
  # need to add 1 correction for 0-cells
  bg.fish.mat <- matrix(c(a.var, c, b.var, d), ncol=2,
                        dimnames=list(c("VarSNP", "BgSNP"), 
                                      c(x.annot, paste0("Not", x.annot)))) + 0.5
  bg.fish.res <- fisher.test(bg.fish.mat, or=1)
  bg.fish.sum <- list("Annot"=x.annot, "OR"=bg.fish.res$estimate, "P"=bg.fish.res$p.value,
                      "L95"=bg.fish.res$conf.int[1], "U95"=bg.fish.res$conf.int[2])

  fisher.res.list[[x.annot]] <- bg.fish.sum
}

# VarSNPs compared to background
bg.enrich.df <- do.call(rbind.data.frame,
                        fisher.res.list)
bg.enrich.df$logOR <- log(bg.enrich.df$OR)
bg.enrich.df$logL95 <- log(bg.enrich.df$L95)
bg.enrich.df$logU95 <- log(bg.enrich.df$U95)

bg.enrich.df$Sig <- as.character(as.numeric(bg.enrich.df$P <= 0.05))

v.enrich <- ggplot(bg.enrich.df, aes(x=reorder(Annot, OR), y=OR, fill=Sig)) +
  #geom_errorbar(aes(ymin=logL95, ymax=logU95)) +
  geom_segment(aes(x=reorder(Annot, OR), xend=reorder(Annot, OR),
                   y=0, yend=OR)) +
  geom_point(shape=21, size=5, colour='black') +
  theme_mike() +
  scale_fill_manual(values=c("0"='grey80', "1"='purple')) +
  labs(y="Odds Ratio", x="Ensembl Annotation") +
  coord_flip() 



### check distance to TSS
# loop over lead SNPs and get the distance to TSS and TES
varqtl.dist.list <- list()

for(i in seq_along(varsnps)){
  varqtl.snp.i <- varsnps[i]
  snp.chr <- allsnp.df[allsnp.df$refsnp_id %in% varqtl.snp.i, ]$chr_name
  snp.bp <- allsnp.df[allsnp.df$refsnp_id %in% varqtl.snp.i, ]$chrom_start
  
  chr.genes <- allgene.df[allgene.df$chromosome_name %in% snp.chr, ]
  positions <- c(abs(chr.genes$start_position - snp.bp))
  end.positions <- c(abs(chr.genes$end_position - snp.bp))
  
  nearest.start <- chr.genes$start_position[which(abs(chr.genes$start_position - snp.bp) == min(positions))]
  nearest.end <- chr.genes$end_position[which(abs(chr.genes$end_position - snp.bp) == min(end.positions))]
  
  varqtl.dist.list[[varqtl.snp.i]] <- list("TSSDist"=unique(snp.bp - nearest.start),
                                           "TESDist"=unique(snp.bp - nearest.end))
}

varqtl.tss.dist <- as.data.frame(do.call(rbind.data.frame,
                                         varqtl.dist.list))
varqtl.tss.dist$SNP <- gsub(rownames(varqtl.tss.dist), pattern="\\.[0-9]", replacement="")

tss.p <- ggplot(varqtl.tss.dist,
                aes(x=TSSDist)) +
  geom_density(bw=1e5, lwd=0.5, alpha=0.75, fill='grey80') +
  geom_vline(mapping=aes(xintercept=-1e5), lty=2, col='purple') +
  geom_vline(mapping=aes(xintercept=1e5), lty=2, col='purple') +
  geom_text(x=-3e5, y=2.2e-6, label="-100kb", fontface='plain', size=5) +
  geom_text(x=3e5, y=2.2e-6, label="100kb", fontface='plain', size=5) +
  scale_x_continuous(limits=c(-1e6, 1e6), oob=squish,
                     breaks=c(-1e6, -5e5, 0, 5e5, 1e6),
                     labels=c("-1Mb","-500kb", "0", "500kb", "1Mb")) +
  #scale_y_continuous(breaks=c(0, 5e-7, 1e-6, 2e-6)) +
  theme_mike()
  

tes.p <- ggplot(varqtl.tss.dist,
                aes(x=TESDist)) +
  geom_density(bw=1e5, lwd=0.5, alpha=0.75, fill='grey80') +
  geom_vline(mapping=aes(xintercept=-1e5), lty=2, col='purple') +
  geom_vline(mapping=aes(xintercept=1e5), lty=2, col='purple') +
  geom_text(x=-3e5, y=2.2e-6, label="-100kb", fontface='plain', size=5) +
  geom_text(x=3e5, y=2.2e-6, label="100kb", fontface='plain', size=5) +
  scale_x_continuous(limits=c(-1e6, 1e6), oob=squish,
                     breaks=c(-1e6, -5e5, 0, 5e5, 1e6),
                     labels=c("-1Mb","-500kb", "0", "500kb", "1Mb")) +
  #scale_y_continuous(breaks=c(0, 5e-7, 1e-6, 2e-6)) +
  theme_mike()


