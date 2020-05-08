# Just the Heatmaps from sequencing results (TopHat counts)
#
# by Spela Konjar and Pedro L Fernandes (@pfern)
# Version 1.0
# Paço de Arcos, April 30th 2020

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
  

# Working directory
#
# For Linux:
setwd ("~/Desktop/Spela_Heatmaps")

Sample_Condition<-c("CD8_mem","CD8_mem","CD8_mem","WT","WT","WT")
  
Genes_of_Interest <- read.table(file="Genes_of_Interest.txt", header = TRUE)
Names_of_Genes_of_Interest <-  Genes_of_Interest[,1]
IDs_of_Genes_of_Interest <- Genes_of_Interest[,2]
  
  
# Read the counts from the trimmed tab delimited files from Tophat output
lane_1 <- read.table(file="lane1_SK5_CD8_mem_TAGCTT_L001_R1_trimmed.tophat-htseq-count.tab", header = FALSE)
lane_2 <- read.table(file="lane1_SK4_CD8_mem_GATCAG_L001_R1_trimmed.tophat-htseq-count.tab", header = FALSE)
lane_3 <- read.table(file="lane1_SK6_CD8_mem_GGCTAC_L001_R1_trimmed.tophat-htseq-count.tab", header = FALSE)
lane_4 <- read.table(file="lane5_WT_SK1_ATCACG_L005_R1_trimmed-htseq-count.tab", header = FALSE)
lane_5 <- read.table(file="lane5_WT_SK2_TTAGGC_L005_R1_trimmed-htseq-count.tab", header = FALSE)
lane_6 <- read.table(file="lane5_WT_SK3_ACTTGA_L005_R1_trimmed-htseq-count.tab", header = FALSE)
  
# Get the list of IDs from the data files (Col 1)
IDs <- lane_1[,1]
  
# Juxtapose the SIX lanes in single a dataframe
six_lanes <- cbind(lane_1[,2],lane_2[,2],lane_3[,2],lane_4[,2], lane_5[,2], lane_6[,2])
# head(six_lanes)
  
# Create a selector to pick the rows of interest
selector <- match(IDs_of_Genes_of_Interest, IDs)
# selector
  
# Select the rows for the heatmap
selected_six_lanes <- six_lanes[selector, c(1:6)]
# IDs[selector]
# selected_six_lanes
  
# Fill the expression matrix, to prepare for the Heatmap

# Logarithic transformation  
# Counts +1  to avouid the log of ZERO
ssl_log <- log2(selected_six_lanes + 1)
  
expression <-  as.matrix(ssl_log)
expression

plot.new()
  
# Draw the Heatmap
# Colour Palette selection 
  
# hmcol controls the colour palette
# dark freen to dark red
# not for the colour blind
#

hmcol  <- colorRampPalette(c("darkgreen", "#7BE67F", "white", "#E84F4F", "darkred"))(100)

# Heatmap itself!
  
heatmap.2( expression, col = hmcol, 
            labCol = Sample_Condition, labRow = Names_of_Genes_of_Interest, 
            scale="none", srtRow = 0, srtCol =35, trace="none", 
            # No Row and No Column order, just as in the files
            Rowv = FALSE, Colv = FALSE, dendrogram="none", 
            # Key and Layout
            keysize = 1, 
            #( "bottom.magin", "left.margin", "top.margin", "right.margin" )
            key.par=list(mar=c(4, 0, 1.5, 35)),
            # lmat -- added 2 lattice sections (5 and 6) for padding
            lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(3, 9), lwid=c(1, 7, 1),
            key.title=NA, key.xlab=NA, key.ylab = NA,  key.xtickfun = NA, key.ytickfun = NA)
  
# And the values for inspection
expression


  


  