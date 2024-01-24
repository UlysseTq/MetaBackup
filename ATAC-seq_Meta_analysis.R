setwd("~/Ulysse_Tuquoi/NGS")

### Libraries used ####
library(tidyverse)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(gridExtra)
library(mixOmics)


### Different ways of importing our data here ####
DATA <- read.csv("Graphs.csv", header =  TRUE, sep = ",")

DF<- head(DATA[c(1:36)], 36)


### Basic graphs while learning ggplot2 ####

#Graph of filtered reads(ratio_subSURmono)
ggplot(DF, mapping = aes(x = ratio_subSURmono, y = filtered_reads)) + 
  geom_jitter(mapping=aes(color=batch)) 

#Four Graphs of unimap(ratio_subSURmono) in different incubation times
ggplot(DF) +
  geom_point(mapping = aes(x=ratio_subSURmono, y=unimap, color=batch, shape=buffer, size=tn5)) +
  facet_wrap(~ incubation_time, nrow=2, labeller = "label_both")

#Four Graphs of unimap(chPt) in different incubation times
ggplot(DF) +
  geom_point(mapping = aes(x=chPt, y=unimap, color=batch, shape=buffer, size=filtered_reads)) +
  facet_wrap(~ incubation_time, nrow=2, labeller = "label_both")

#Effect of incubation time on filtered reads in our libraries
ggplot(DF) +
  geom_bar(mapping = aes(x=name, y=filtered_reads, fill='incubation_time'), stat="identity") +
  coord_flip()

#Chromosome00 values color-coded by unimape percentage
ggplot(DF) +
  geom_bar(mapping = aes(x=name, y=SL4.0ch00, fill=unimap), stat="identity")

#Chromosome00 values color-coded by ChPt values
ggplot(DF) +
  geom_bar(mapping = aes(x=name, y=SL4.0ch00, fill=chPt), stat="identity")

#Variance of the ratio_subSURmono between our batches
ggplot(DF) + 
  stat_summary(
    mapping = aes(x = batch, y = ratio_subSURmono),
    fun.min = min,
    fun.max = max,
    fun = median
  )

#Percentage of unimap and its spread when using either Gif or sysmex
ggplot(DF) +
  geom_boxplot(mapping = aes(x=buffer,y=unimap))


### Trying to produce PCAs ####

#Two ways of removing qualitative datas from the dataframe
#DF.PCA <- head(DATA[c(3:34)], 36)
DF.PCA <- DATA %>% dplyr::select(-name, -buffer)

#Calculating the values to choose the number of dimensions
res.pca = PCA(X = DF.PCA, scale.unit = TRUE, graph = FALSE)
eig.val = get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100))

#Correlations circle
fviz_pca_var(res.pca, axes = c(1,2), 
  col.var = "cos2",
  gradient.cols = c("blue", "orange", "red"), 
  repel = TRUE # Evite le chevauchement de texte
)

#Calculate each variable's contribution and describe the PC
var = get_pca_var(res.pca) 
corrplot(var$contrib, is.corr=FALSE)

ind = get_pca_ind(res.pca)
corrplot(ind$contrib, is.corr=FALSE)

dimdesc(res.pca, proba = 0.05)

#Bare individuals graph
fviz_pca_ind(res.pca, axes = c(1,2),
  repel = TRUE # Evite le chevauchement de texte 
)

#Individuals graph grouped by buffer
fviz_pca_ind(res.pca, axes = c(1,2),
  col.ind = DF$buffer, # variable qualitative à représenter
  legend.title = "Groupes" 
)

#Individuals graph grouped by batch
fviz_pca_ind(res.pca, axes = c(1,2),
             col.ind = DF$batch, # variable qualitative à représenter
             legend.title = "Groupes" 
)

#Individuals graph grouped by buffer with confidence ellipses
fviz_pca_ind(res.pca, axes = c(1,2),
  geom.ind = "point", # Montre les points seulement
  col.ind = DF$buffer, # variable qualitative à représenter
  addEllipses = TRUE, # Ellipses de concentration 
  legend.title = "Groupes" 
)

#PCA biplot with too much data
PCA <- fviz_pca_biplot(res.pca, axes = c(1,2),
                geom.ind = c("point", "text"), # Add text to show labels
                col.ind = DF$buffer, # variable qualitative à représenter
                legend.title = "Groupes",
                repel = TRUE
)

table_data <- data.frame(Names = rownames(DF), DF$name)
table_plot <- tableGrob(table_data)

grid.arrange(PCA, table_plot, ncol = 2)

### Reducing my dataframe by creating the mean of my chromosomes values (excluding 00,Mt and Pt) ####

# Assuming 'DF' is your data frame and 'cols_to_remove' is a vector of column positions
cols_to_remove <- 21:32 

# Get the column names for these positions
col_names_to_remove <- colnames(DF)[cols_to_remove]

DF2 <- DF %>%
  mutate("mean_Chr1-12" = rowMeans(dplyr::select(., all_of(col_names_to_remove)), na.rm = TRUE)) %>% dplyr::select(-all_of(col_names_to_remove))


#Removing qualitative datas from the dataframe
DF2.PCA <- DF2 %>% dplyr::select(-name, -buffer, -batch, -tween, -digitonin, -ploidy)

#Calculating the values to choose the number of dimensions
res.pca2 = PCA(X = DF2.PCA, scale.unit = TRUE, graph = FALSE)
eig.val2 = get_eigenvalue(res.pca2)
fviz_eig(res.pca2, addlabels = TRUE, ylim = c(0, 100))

#Correlations circle
fviz_pca_var(res.pca2, axes = c(1,2), 
             col.var = "cos2",
             gradient.cols = c("blue", "orange", "red"), 
             repel = TRUE # Evite le chevauchement de texte
)

fviz_pca_var(res.pca2, axes = c(2,3), 
             col.var = "cos2",
             gradient.cols = c("blue", "orange", "red"), 
             repel = TRUE # Evite le chevauchement de texte
)

fviz_pca_var(res.pca2, axes = c(1,3), 
             col.var = "cos2",
             gradient.cols = c("blue", "orange", "red"), 
             repel = TRUE # Evite le chevauchement de texte
)

fviz_pca_var(res.pca2, axes = c(2,5), 
             col.var = "cos2",
             gradient.cols = c("blue", "orange", "red"), 
             repel = TRUE # Evite le chevauchement de texte
)

#Calculate each variable's contribution and describe the PC
var = get_pca_var(res.pca2) 
corrplot(var$contrib, is.corr=FALSE)

ind = get_pca_ind(res.pca2)
corrplot(ind$contrib, is.corr=FALSE)

dimdesc(res.pca2, proba = 0.05)

#Bare individuals graph
fviz_pca_ind(res.pca2, axes = c(1,2),
             repel = TRUE # Evite le chevauchement de texte 
)

fviz_pca_ind(res.pca2, axes = c(2,3),
             repel = TRUE # Evite le chevauchement de texte 
)

fviz_pca_ind(res.pca2, axes = c(1,3),
             repel = TRUE # Evite le chevauchement de texte 
)

#Individuals graph grouped by buffer
fviz_pca_ind(res.pca2, axes = c(1,2),
             col.ind = DF$buffer, # variable qualitative à représenter
             legend.title = "Groupes" 
)

#Individuals graph grouped by TSS
fviz_pca_ind(res.pca2, axes = c(1,2),
             gradient.cols = c("blue", "orange", "red"), 
             col.ind = DF$TSS_EnrichmentScore, # variable qualitative à représenter
             legend.title = "Groupes" 
)

fviz_pca_ind(res.pca2, axes = c(2,3),
             gradient.cols = c("blue", "orange", "red"), 
             col.ind = DF$TSS_EnrichmentScore, # variable qualitative à représenter
             legend.title = "Groupes" 
)

fviz_pca_ind(res.pca2, axes = c(1,3),
             gradient.cols = c("blue", "orange", "red"), 
             col.ind = DF$TSS_EnrichmentScore, # variable qualitative à représenter
             legend.title = "Groupes" 
)

#Individuals graph grouped by agitation
fviz_pca_ind(res.pca2, axes = c(1,2),
             geom.ind = "point", # Montre les points seulement
             gradient.cols = c("blue", "orange", "red"), 
             col.ind = DF$incubation_agitation, # variable qualitative à représenter
             legend.title = "Agitation" 
)

fviz_pca_ind(res.pca2, axes = c(1,2),
             geom.ind = "point", # Montre les points seulement
             gradient.cols = c("blue", "orange", "red"), 
             col.ind = DF$tn5, # variable qualitative à représenter
             legend.title = "Tn5" 
)

#PCA biplot with too much data
fviz_pca_biplot(res.pca2, axes = c(1,2),
              geom.ind = c("point", "text"), # Add text to show labels
              gradient.cols = c("blue", "orange", "red"), 
              col.ind = DF2$TSS_EnrichmentScore_smallFrags, #variable qualitative à représenter
              legend.title = "TSS_EnrichmentScore_smallFrags",
              repel = TRUE
)

fviz_pca_biplot(res.pca2, axes = c(2,3),
                geom.ind = c("point", "text"), # Add text to show labels
                gradient.cols = c("blue", "orange", "red"), 
                col.ind = DF2$TSS_EnrichmentScore_smallFrags, #variable qualitative à représenter
                legend.title = "TSS_EnrichmentScore_smallFrags",
                repel = TRUE
)

fviz_pca_biplot(res.pca2, axes = c(1,3),
                geom.ind = c("point", "text"), # Add text to show labels
                gradient.cols = c("blue", "orange", "red"), 
                col.ind = DF2$TSS_EnrichmentScore_smallFrags, #variable qualitative à représenter
                legend.title = "TSS_EnrichmentScore_smallFrags",
                repel = TRUE
)

fviz_pca_biplot(res.pca2, axes = c(2,5),
                geom.ind = c("point", "text"), # Add text to show labels
                gradient.cols = c("blue", "orange", "red"), 
                col.ind = DF2$TSS_EnrichmentScore_smallFrags, #variable qualitative à représenter
                legend.title = "TSS_EnrichmentScore_smallFrags",
                repel = TRUE
)


### Trying to cluster our values through CAH ####

DF3.cluster <- DF2.PCA

DF.reduced <- DF2.PCA %>% dplyr::select(-age, -tn5, -incubation_time, -incubation_agitation, -pcr_cycles)

# Distance matrix
distance = dist(scale(DF.reduced))

#Ward.D2 CAH
model = hclust(distance, method = "ward.D2")

#Show the dendro
plot(model, hang=-1)

#Add the number of wanted clusters (k) to the dendro
plot(model, hang=-1)
rect.hclust(model, k = 3)

#Cuts the dendrogram
classe = cutree(model, k = 3)

#Adds the classes to the DATA
DF3.cluster$classeCAH = as.factor(classe)


### Trying to cluster our values through k-means ####

#Setting a seed as a reference
set.seed(123)

#Creating the k-means dataset on DF.reduced with scaling
res.km <- kmeans(scale(DF.reduced), 3, nstart = 100)

#Showing the graph of the clusters on the PCA dimensions
fviz_cluster(res.km, data = DF.reduced, 
             axes = c(1,2),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom.ind = c("point", "text"), 
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

fviz_cluster(res.km, data = DF.reduced, 
             axes = c(2,3),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom.ind = c("point", "text"), 
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

fviz_cluster(res.km, data = DF.reduced, 
             axes = c(1,3),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom.ind = c("point", "text"), 
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#Adds the classes to the DATA
DF3.cluster$classeKm = as.factor(res.km$cluster)


### Partial Least Square analysis ####

##PLS

#Setting a seed as a reference
set.seed(123)

#
X <- DF2.PCA %>% dplyr::select(-TSS_EnrichmentScore_smallFrags, -TSS_EnrichmentScore) # use the gene expression data as the X matrix
Y <- DF2.PCA$TSS_EnrichmentScore_smallFrags+DF2.PCA$TSS_EnrichmentScore # use the clinical data as the Y matrix

pls.result <- pls(X, Y)

plotIndiv(pls.result)
x11() # plot the samples
plotVar(pls.result)     # plot the variables

##sPLS

spls.result <- spls(X, Y, keepX = c(8, 17), keepY = c(1, 1))  # run the method
plotIndiv(spls.result) # plot the samples
plotVar(spls.result)   # plot the variables
selectVar(spls.result, comp = 1)$X$name
plotLoadings(spls.result, method = 'mean', contrib = 'max') # depict weight assigned to each of these variables
