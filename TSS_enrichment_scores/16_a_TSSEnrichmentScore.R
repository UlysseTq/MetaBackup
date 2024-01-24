projPath <- "/DATA/WORK/INRA/BFP/PROJETS/2023_ATACseq02_BGI"
setwd(projPath)

#Libraires
library(GenomicFeatures)
library(GenomicAlignments)
library(ggplot2)
library(dplyr)

#options
options(width=160)

#functions
source("/DATA/LIB/Rfun/RleList2matrix.r")
source("/DATA/LIB/Rfun/GetAverageProfileWithCI.r")
# TSS enrichment score
# described here: https://www.encodeproject.org/data-standards/terms/#enrichment


#function to get default ggplot2 colors:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#function to make TSS/TES plots:
#~ myplotfunGR <- function(mydf, mycolors) {
#~   p <- ggplot(mydf, aes(x=Position, y=Coverage, group=GeneGroup))+
#~     geom_ribbon(aes(ymin = Lower,
#~                     ymax = Upper,
#~                     fill = GeneGroup),linetype=0,alpha=0.2)+
#~     geom_line(aes(x = Position,
#~                   y = Coverage,
#~                   color = GeneGroup),lwd=1) +
#~     scale_colour_manual(values = mycolors) +
#~     scale_fill_manual(values = mycolors) +
#~     ylab("Coverage") +
#~     theme_bw() +
#~     theme(text=element_text(size=18, color = "black"),
#~           panel.grid.minor.x=element_blank())

#~   return(p)
#~ }

#Get the list of 10339 well expressed genes:
WEG <- readRDS("/DATA/WORK/INRA/BFP/PROJETS/2023_RNAseq_Pericarp/Aligned/STAR_SL4.0_MtPt/RData/WellExpressedGenes_GR.rds")

##-------------
## Function to obtain the normalized TSS profiles
##-------------

GetNormTSSprof <- function(tssrle, GOF) {
    # Convert TSS rle to matrix and extract genes in GOF
    tssmat <- RleList2matrix(tssrle[names(tssrle) %in% names(WEG)])
    # get the internal 4001bp (RLE should be centered on the TSS)
    tssmat2kb <- tssmat[,(((ncol(tssmat) - 1) / 2) - 2000 + 1):(4001 + (((ncol(tssmat) - 1) / 2) - 2000))]
    # Get the normalization factor
    NormFactor <- mean(tssmat2kb[,c(1:100,3902:4001)] , na.rm = TRUE)
    #Normalize the profiles by this normalization factor
    tss2kb_norm <- tssmat2kb / NormFactor
    #Get average profile
    normprof_tss2kb <- GetAverageProfileWithCI(tss2kb_norm, 
                                               conftype="basic",
                                               normtype = "none",
                                               filterOutliers = FALSE)
    colnames(normprof_tss2kb)[2] <- "NormCoverage"
    #Return the average profile
    return(normprof_tss2kb)
}

##-------------
## Get the normalized profiles for all files
##-------------

###------
### Profiles with all reads
###------

#Profiles path
fnALL <- dir(file.path(projPath, "Aligned", "bowtie2", "SL4.0_Mt_Pt", "cutadapt", "RDataQ2_noChrCM"),
             pattern = "TSS4kb_cutadaptfiltrQ2_normcovr_.+.rds",
             full = TRUE)
#sample names
snALL <- gsub("TSS4kb_cutadaptfiltrQ2_normcovr_|.rds", "", basename(fnALL))

#Get average profiles
profALL <- list()
for (i in 1:length(fnALL)) {
    profALL[[snALL[i]]] <- GetNormTSSprof(readRDS(fnALL[i]), WEG)
}

#Convert to data frame
profALLdf <- reshape2::melt(profALL, id.vars=1:4)
profALLdf <- dplyr::rename(profALLdf, "Sample" = "L1")

#plot
pp <- profALLdf |>
    ggplot(aes(x=Position, y=NormCoverage, group = Sample))+
    geom_ribbon(aes(ymin = Lower,
                    ymax = Upper,
                    group = Sample),
                fill = "blue", linetype=0,alpha=0.2)+
    geom_line(aes(x = Position,
                  y = NormCoverage, 
                  group = Sample),
              color = "blue", lwd=1) +
    geom_vline(xintercept=c(1801,2201), linetype = "dashed", color = "red") +
    ylab("Normalized coverage") +
    scale_x_continuous(breaks = c(1, 1001, 2001, 3001, 4001),
                       labels = c("-2kb", "-1kb", "TSS", "+1kb", "+2kb")) +
    theme_bw() +
    theme(text=element_text(size=18, color = "black"),
          panel.grid.minor.x=element_blank()) +
    facet_wrap(vars(Sample))
    NULL

myres=150
png(file.path(projPath, "Aligned", "bowtie2", "SL4.0_Mt_Pt", "cutadapt", "RDataQ2_noChrCM", "Plots", "TSSenrichmentScore_ATACseq02_BGI.png"),
    width = 12 * myres,
    height = 12 * myres,
    res = myres)
print(pp)
dev.off()


#Extract the max values within +/- 200bp of the TSS => this is the TSS enrichment score
TSS_ES_ALL <- unlist(bplapply(profALL, function(x) {filter(x, Position > 1800, Position < 2202) %>% pull(NormCoverage) %>% max}))
TSS_ES_ALL
#~    M11S2    M11S4    M12S1    M12S2    M12S3    M12S4    M12S5    M12S6    M12S7    M12S8    M12S9     M8S3 
#~ 1.247406 1.275686 1.375647 1.677606 1.826764 1.472155 1.674340 1.561388 1.336292 1.253725 1.226449 1.470624 

###------
### Profiles with fragments <102bp
###------

#Profiles path
fnshort <- dir(file.path(projPath, "Aligned", "bowtie2", "SL4.0_Mt_Pt", "cutadapt", "RDataQ2_MaxFragSize102"),
             pattern = "TSS4kb_filtQ2_normcovr_smallFrags_.+.rds",
             full = TRUE)
#sample names
snshort <- gsub("TSS4kb_filtQ2_normcovr_smallFrags_|.rds", "", basename(fnshort))

#Get average profiles
profShort <- list()
for (i in 1:length(fnshort)) {
    profShort[[snshort[i]]] <- GetNormTSSprof(readRDS(fnshort[i]), WEG)
}

#Convert to data frame
profShortdf <- reshape2::melt(profShort, id.vars=1:4)
profShortdf <- dplyr::rename(profShortdf, "Sample" = "L1")

#plot
pps <- profShortdf |>
    ggplot(aes(x=Position, y=NormCoverage, group = Sample))+
    geom_ribbon(aes(ymin = Lower,
                    ymax = Upper,
                    group = Sample),
                fill = "blue", linetype=0,alpha=0.2)+
    geom_line(aes(x = Position,
                  y = NormCoverage, 
                  group = Sample),
              color = "blue", lwd=1) +
    geom_vline(xintercept=c(1801,2201), linetype = "dashed", color = "red") +
    ylab("Normalized coverage") +
    scale_x_continuous(breaks = c(1, 1001, 2001, 3001, 4001),
                       labels = c("-2kb", "-1kb", "TSS", "+1kb", "+2kb")) +
    theme_bw() +
    theme(text=element_text(size=18, color = "black"),
          panel.grid.minor.x=element_blank()) +
    facet_wrap(vars(Sample))
    NULL

myres=150
png(file.path(projPath, "Aligned", "bowtie2", "SL4.0_Mt_Pt", "cutadapt", "RDataQ2_MaxFragSize102", "Plots", "TSSenrichmentScore_ATACseq02_smallFrags.png"),
    width = 12 * myres,
    height = 12 * myres,
    res = myres)
print(pps)
dev.off()

#Extract the max values within +/- 200bp of the TSS => this is the TSS enrichment score
TSS_ES_Short <- unlist(bplapply(profShort, function(x) {filter(x, Position > 1800, Position < 2202) %>% pull(NormCoverage) %>% max}))
TSS_ES_Short
#~    M11S2    M11S4    M12S1    M12S2    M12S3    M12S4    M12S5    M12S6    M12S7    M12S8    M12S9     M8S3 
#~ 1.766228 1.639397 1.631518 2.426816 2.724471 1.955810 2.309583 2.024621 2.039357 1.797171 1.647638 2.603228 


### Create and save a data frame with the enrichment scores
TSSES_df <- data.frame(Experiment = "ATACseq02_BGI",
                       SampleName = names(TSS_ES_ALL),
                       TSS_EnrichmentScore = TSS_ES_ALL,
                       TSS_EnrichmentScore_smallFrags = TSS_ES_Short)

saveRDS(TSSES_df, file.path(projPath, "Aligned", "bowtie2", "SL4.0_Mt_Pt", "cutadapt", "TSS_EnrichmentScores_ATACseq02_BGI.rds"))
write.table(TSSES_df, file.path(projPath, "Aligned", "bowtie2", "SL4.0_Mt_Pt", "cutadapt", "TSS_EnrichmentScores_ATACseq02_BGI.tsv"),
            sep="\t", dec=".", quote = FALSE, row.names = FALSE)


