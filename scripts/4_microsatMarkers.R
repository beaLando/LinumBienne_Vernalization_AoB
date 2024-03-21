########################################################################
# LIBRARIES                                                            #
########################################################################

library(tidyverse)
library(reshape2)
library(lattice)
library(patchwork)

library(adegenet) #Ho, Hs (unbiased)
library(hierfstat) #Hs (biased), pairwise Fst, alleles summaries
library(poppr) #MLG, DAPC (and other clustering methods)
library(pegas) #Hardy-Weinberg Eq., Fst-Fis-Fit
library(pophelper) #check online instructions for installation https://www.royfrancis.com/pophelper/#:~:text=pophelper%20is%20an%20R%20package,such%20as%20ADMIXTURE%20or%20fastSTRUCTURE.

library(flextable)
library(leaflet)
library(leaflet.minicharts)



########################################################################
# REMINDER                                                             #
########################################################################
# AA = x    # Aa = y    # aa = z    x + y + z = N (sample size)

# f(AA) = x / N       f(Aa) = y / N       f(aa) = z / N

# f(A) = (2x + y) / 2N          f(a) = (2z + y) / 2N

# let p = f(A), q = f(a)    p & q are allele frequencies

## SO..
##  p + q = 1     p = 1 - q    q = 1 - p

##  (p + q)2  =  p2 + 2pq + q2  =  1

##  (1 - q)2 + 2(1 - q)(q) + q2 = 1

###   Ho = f(Aa) = observed heterozygosity 
###   He = 2pq   = expected heterozygosity (for two alleles)

###   He = 2pq + 2pr + 2qr = 1 - (p2 + q2 + r2)    for three alleles ...and so on...  He = 1 - SUM(qi)2, for i to n alleles



########################################################################
# DATA                                                                 #
########################################################################

# Raw (only dominant markers analysed in poppr)
ssr <- read.genalex("./data/genotype/Microsats_dominant_conversion_poppr.csv", ploidy = 1)
ssr.info <- read.csv("./data/genotype/Microsats_dominant_conversion_poppr_info.csv")
ssr.info$Sample.no. <- as.character(ssr.info$Sample.no.)

other(ssr)$xy <- ssr.info[, c(4, 5)]

str(ssr)
ssr

# Structure Output
path.str.dom <- "./data/genotype/structure_output/dominant" #this is best, since our markers show signs of duplications. The ancestor of Lb-Lus underwent a few rounds of WGD in ancient times.
path.str.codom <- "./data/genotype/structure_output/codominant"

## Dominant
sfiles <- list.files(path = path.str.dom, pattern = "\\_f$", full.names=T)
slist <- readQ(files = sfiles, filetype="structure")

names(attributes(slist[[1]]))

tr1 <- tabulateQ(qlist = slist)
sr1 <- summariseQ(tr1)
summariseQ(tr1, writetable = TRUE, exportpath = path.str.dom)

str.info <- ssr.info
colnames(str.info)[3] <- "V1"

if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist, "rownames<-", str.info$V1)

## Co-dominant
sfiles.co <- list.files(path = path.str.codom, pattern = "\\_f$", full.names=T)
slist.co <- readQ(files = sfiles.co, filetype="structure")

names(attributes(slist.co[[1]]))

tr1.co <- tabulateQ(qlist = slist.co)
sr1.co <- summariseQ(tr1.co)
summariseQ(tr1.co, writetable = TRUE, exportpath = path.str.codom)

str.info.co <- ssr.info
colnames(str.info.co)[3] <- "V1"

if(length(unique(sapply(slist.co,nrow)))==1) slist.co <- lapply(slist.co, "rownames<-", str.info.co$V1)


########################################################################
# ANALYSIS                                                             #
########################################################################

# POPPR ANALYSIS of SSR as DOMINANT MARKERS
## Stats & filtering
gac <- genotype_curve(ssr, sample = 1000, quiet = TRUE)
gac

info_table(ssr, type = "missing", plot = TRUE, plotlab = FALSE)

ssr <- missingno(ssr, type = "geno", cutoff = 0.25) #Removing 2 genotypes: VIL_19 and LLA_14

mlg.table(ssr)
informloci(ssr)

# cutoff value: 1.2987012987013 % ( 2 samples ).
# MAF         : 0.01
# 
# Found 38 uninformative loci 
# ============================ 
#   38 loci found with a cutoff of 2 samples :
#   ssr11_1_193, ssr10_1_119, ssr10_1_121, ssr10_1_125, ssr10_1_123, ssr10_1_127, ssr6_1_161, ssr6_1_157, ssr6_1_163, ssr6_11_167,
# ssr2_1_215, ssr2_1_230, ssr2_1_227, ssr2_1_224, ssr2_1_233, ssr11_2_243, ssr11_2_246, ssr11_2_249, ssr11_2_264, ssr11_2_261,
# ssr11_2_255, ssr4_2_203, ssr4_2_206, ssr4_2_200, ssr2a_2_280, ssr2a_2_282, ssr2b_2_335, ssr12_3_222, ssr12_4_190, ssr3_4_246,
# ssr3_4_250, ssr3_4_230, ssr3_4_248, ssr3_4_252, ssr3_4_242, ssr3_4_234, ssr3_4_232, ssr3_4_226 
# 8 loci found with MAF < 0.01 :
#   ssr6_11_167, ssr2_1_224, ssr2_1_233, ssr2b_2_335, ssr12_3_222, ssr12_4_190, ssr3_4_242, ssr3_4_226
# 
# This is a genclone object
# -------------------------
#   Genotype information:
#   
#   75 original multilocus genotypes 
# 154 haploid individuals
# 26 dominant loci
# 
# Population information:
#   
#   1 stratum - Pop
# 6 populations defined - 11, 6, IOW2, LLA, SUT, VIL


## DAPC
### k-means
maxK <- 6 #there are 6 populations
myMat <- matrix(nrow = 10, ncol = maxK)

colnames(myMat) <- 1:ncol(myMat)

for(i in 1:nrow(myMat)){
  grp <- find.clusters(ssr, n.pca = 20, choose.n.clust = FALSE,  max.n.clust = maxK) #20 pcs for ~90% variance explained
  myMat[i,] <- grp$Kstat
  }

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1 #2 seems the best number of clusters

ggsave("./output/ssr_kmeans_bic.pdf", p1)

### DAPC
my_k <- 2:6
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(666)
  grp_l[[i]] <- find.clusters(ssr, n.pca = 20, n.clust = my_k[i])
  ssr@strata$cluster <- grp_l[[i]]$grp
  setPop(ssr) <- ~ cluster
  dapc_l[[i]] <- dapc(ssr, grp = grp_l[[i]]$grp, n.pca = 20, n.da = my_k[i], scale = TRUE)
  }

#### Structure-like plot
my_df <- NULL

for(i in 1:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Family <- rownames(tmp)
  tmp <- melt(tmp, id = c("Family", "K"))
  names(tmp)[3:4] <- c("Cluster", "Posterior")
  tmp$Pop <- ssr@strata$Pop
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

str(my_df)
my_df$Cluster <- factor(my_df$Cluster, levels = c("1", "2", "3", "4", "5", "6", "7"))
my_df$Pop <- factor(my_df$Pop, levels = c("6", "11", "LLA", "VIL", "IOW2", "SUT"))

ps <- ggplot(my_df, aes(x = Family, y = Posterior, fill = Cluster)) + 
  geom_bar(stat = "identity") + 
  facet_grid(K ~ Pop, scales = "free", space = "free", 
             labeller = labeller(K = grp.labs)) + 
  theme_bw() + 
  ylab("Posterior membership probability") + 
  theme(legend.position='none') + 
  scale_fill_manual(values = c("slateblue4", "slategray3", "seashell2", "tan3", "navajowhite2", "indianred4", "gray")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))

ps

#### PCA-like plot
my_df <- NULL
my_df <- as.data.frame(dapc_l[[2]]$ind.coord)
my_df$Cluster <- dapc_l[[2]]$grp
my_df$Pop <- ssr@strata$Pop

str(my_df)
my_df$Cluster <- factor(my_df$Cluster, levels = c("1", "2", "3"))
my_df$Pop <- factor(my_df$Pop, levels = c("6", "11", "LLA", "VIL", "IOW2", "SUT"))

pp <- ggplot(my_df, aes(x = LD1, y = LD2, color = Cluster, fill = Cluster)) + 
  geom_point(size = 6, aes(shape = Pop), alpha = 0.5) + 
  theme_classic(base_size = 20) + 
  theme(legend.position = "top") +
  guides(shape = guide_legend(nrow = 1)) +
  scale_color_manual(values = c("#003C67FF", "#8A919799", "#CD534CFF", "#709AE199", "moccasin", "wheat4")) +
  scale_fill_manual(values = c("#003C67FF", "#8A919799", "#CD534CFF", "#709AE199", "moccasin", "wheat4"))

pp



# PROCESSING of STRUCTURE OUTPUT
p <- evannoMethodStructure(data = sr1, exportplot = TRUE, writetable = TRUE, basesize = 12, linesize = 0.7, exportpath = path.str.dom)

## SSR as dominant markers
slist1 <-alignK(slist)
labs1 <- str.info[, 2, drop = FALSE]

### All runs
p1 <- plotQ(slist1[1:10], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p2 <- plotQ(slist1[11:20], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p3 <- plotQ(slist1[21:30], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p4 <- plotQ(slist1[31:40], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p5 <- plotQ(slist1[41:50], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p6 <- plotQ(slist1[51:60], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))

grid.arrange(p1$plot[[1]], p2$plot[[1]], p3$plot[[1]], p4$plot[[1]], p5$plot[[1]], p6$plot[[1]])

### Subset
psd <- plotQ(slist1[c(11, 21, 31)], #I pick a few representative runs #c(11, 21, 31, 41, 51)
            imgoutput = "join", 
            clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"), basesize = 18, showyaxis = T, showticks = T,
            #showindlab = T, useindlab = T, indlabcol = "black", indlabangle = 45,indlabvjust = 1, indlabsize = 5, 
            sortind = "all", sharedindlab = FALSE, grplab = labs1, subsetgrp = c("6", "11", "LLA", "VIL", "IOW2", "SUT"), grplabsize = 4, linesize = 0.8, pointsize = 4,
            splab = c("K = 2", "K = 3", "K = 4"), #c("K = 2", "K = 3", "K = 4", "K = 5", "K = 6")
            returnplot = T, exportplot = F, 
            showlegend = TRUE, legendpos = "right", legendkeysize = 15, legendtextsize = 15, legendspacing = 5, legendrow = 1)

grid.arrange(psd$plot[[1]])


psd.unordered <- plotQ(slist1[c(11, 21, 31)], 
             imgoutput = "join", 
             clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"), basesize = 11, showyaxis = T, #showticks = T,
             #showindlab = T, useindlab = T, sharedindlab = TRUE, indlabcol = "black", indlabangle = 45,indlabvjust = 1, indlabsize = 5, 
             grplab = labs1, subsetgrp = c("6", "11", "LLA", "VIL", "IOW2", "SUT"), grplabsize = 4, linesize = 0.8, pointsize = 4,
             splab = c("K = 2", "K = 3", "K = 4"), 
             returnplot = T, exportplot = F, 
             showlegend = TRUE, legendpos = "right", legendkeysize = 15, legendtextsize = 15, legendspacing = 5, legendrow = 1)

grid.arrange(psd.unordered$plot[[1]])


## SSR as codominant markers
slist1.co <-alignK(slist.co) ## for comparison between STRUCTURE output using SSR markers as codominant and dominant I used CLUMPAK online (https://tau.evolseq.net/clumpak/compareDifferentPrograms.html), option "Compare"

### All runs
p1 <- plotQ(slist1.co[1:10], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p2 <- plotQ(slist1.co[11:20], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p3 <- plotQ(slist1.co[21:30], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p4 <- plotQ(slist1.co[31:40], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p5 <- plotQ(slist1.co[41:50], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))
p6 <- plotQ(slist1.co[51:60], imgoutput="join", returnplot=T,exportplot=F,basesize=11, clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"))

grid.arrange(p1$plot[[1]], p2$plot[[1]], p3$plot[[1]], p4$plot[[1]], p5$plot[[1]], p6$plot[[1]])

### Subset
pscd <- plotQ(slist1.co[c(11, 21, 31)], #I pick a few representative runs #c(11, 21, 31, 41, 51)
             imgoutput = "join", 
             clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"), basesize = 18, showyaxis = T, showticks = T,
             #showindlab = T, useindlab = T, indlabcol = "black", indlabangle = 45,indlabvjust = 1, indlabsize = 5, 
             sortind = "all", sharedindlab = FALSE, grplab = labs1, subsetgrp = c("6", "11", "LLA", "VIL", "IOW2", "SUT"), grplabsize = 4, linesize = 0.8, pointsize = 4,
             splab = c("K = 2", "K = 3", "K = 4"), #c("K = 2", "K = 3", "K = 4", "K = 5", "K = 6")
             returnplot = T, exportplot = F, 
             showlegend = TRUE, legendpos = "right", legendkeysize = 15, legendtextsize = 15, legendspacing = 5, legendrow = 1)

grid.arrange(pscd$plot[[1]])


pscd.unordered <- plotQ(slist1.co[c(11, 21, 31)], 
                       imgoutput = "join", 
                       clustercol = c("#8A919799", "#CD534CFF", "#003C67FF", "#709AE199", "moccasin", "wheat4"), basesize = 11, showyaxis = T, #showticks = T,
                       #showindlab = T, useindlab = T, sharedindlab = TRUE, indlabcol = "black", indlabangle = 45,indlabvjust = 1, indlabsize = 5, 
                       grplab = labs1, subsetgrp = c("6", "11", "LLA", "VIL", "IOW2", "SUT"), grplabsize = 4, linesize = 0.8, pointsize = 4,
                       splab = c("K = 2", "K = 3", "K = 4"), 
                       returnplot = T, exportplot = F, 
                       showlegend = TRUE, legendpos = "right", legendkeysize = 15, legendtextsize = 15, legendspacing = 5, legendrow = 1)

grid.arrange(pscd.unordered$plot[[1]])



# PLOT DAPC & STRUCTURE OUTPUTS TOGETHER (only for dominant)
## Scatter + Bars
pt <- (wrap_elements(grid.arrange(psd$plot[[1]])) / pp)
#ggsave("./output/ssr.pdf", pt, width = 20, height = 15, units = "in")