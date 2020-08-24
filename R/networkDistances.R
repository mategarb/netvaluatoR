library("NetworkDistance")
library("corrplot")
library("VisuNet")
library("reshape2")
library("ggplot2")
library("ggpubr")

### balance rules ### ASD groups
rls_min <- min(table(newModel$decision))
Md_l1 <- list()
Md_l2 <- list()
Md_l3 <- list()
Md_l4 <- list()
Md_l5 <- list()
Md_l6 <- list()

for(i in 1:20){
  set.seed(i)
pddnos_n <- sample(which(as.character(newModel$decision)=="PDD-NOS"), rls_min)
asperg_n <- sample(which(as.character(newModel$decision)=="ASPERGER'S DISORDER"), rls_min)
autism_n <- sample(which(as.character(newModel$decision)=="AUTISM"), rls_min)
control_n <- sample(which(as.character(newModel$decision)=="CONTROL"), rls_min)
newModel2 <- newModel[c(pddnos_n, asperg_n, autism_n, control_n), ]

net0 <- network_gen(newModel2, decision="PDD-NOS", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 
net1 <- network_gen(newModel2, decision="ASPERGER'S DISORDER", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 
net2 <- network_gen(newModel2, decision="AUTISM", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 
net3 <- network_gen(newModel2, decision="CONTROL", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 

vis_out <- list(net0[[1]], net1[[1]], net2[[1]], net3[[1]])
names(vis_out) <- c("PDD-NOS", "ASPERGER'S DISORDER", "AUTISM", "CONTROL")
#vis_out <- visunet(newModel2, type = "RDF")

#### create adjacency matrices
ams <- createAdjMat(vis_out)

Md1 <- nd.centrality(ams, mode="Degree", directed = FALSE)
Md2 <- nd.centrality(ams, mode="Close", directed = FALSE)
Md3 <- nd.centrality(ams, mode="Between", directed = FALSE)
Md4 <- nd.hamming(ams, out.dist = TRUE)
Md5 <- nd.gdd(ams, out.dist = TRUE)
Md6 <- nd.dsd(ams, out.dist = TRUE, type="Adj")

Md1 <- as.matrix(Md1$D)
Md2 <- as.matrix(Md2$D)
Md3 <- as.matrix(Md3$D)
Md4 <- as.matrix(Md4$D)
Md5 <- as.matrix(Md5$D)
Md6 <- as.matrix(Md6$D)

Md_l1[[i]] <- Md1
Md_l2[[i]] <- Md2
Md_l3[[i]] <- Md3
Md_l4[[i]] <- Md4
Md_l5[[i]] <- Md5
Md_l6[[i]] <- Md6

print(i)
}
Md_fin1 <- Reduce("+", Md_l1)/length(Md_l1)
#Md_fin2 <- Reduce("+", Md_l2)/length(Md_l2)
Md_fin3 <- Reduce("+", Md_l3)/length(Md_l3)
#Md_fin4 <- Reduce("+", Md_l4)/length(Md_l4)
#Md_fin5 <- Reduce("+", Md_l5)/length(Md_l5)
#Md_fin6 <- Reduce("+", Md_l6)/length(Md_l6)

Md_plot <- Md_fin1

Md_plot[upper.tri(Md_plot)] <- Inf
nams0 <- c("PDD-NOS","AS","autism","control")
colnames(Md_plot) <- nams0 #gsub(" ","_",tolower(names(ams)))
rownames(Md_plot) <- nams0 #gsub(" ","_",tolower(names(ams)))

melted_cormat <- melt(Md_plot)
head(melted_cormat)

vcn <- as.numeric(Md_plot)
vcn[vcn==0] <- NA

melted_cormat$value[melted_cormat$value==0] <- NA
melted_cormat2 <- melted_cormat[-which(is.infinite(melted_cormat$value)),]
#melted_cormat2$value <- as.numeric(scale(melted_cormat2$value))

ggplot(data = melted_cormat2 , aes(x=Var1, y=Var2, fill=value), label = ifelse(Padj < 0.001, "***", "ns")) + 
geom_tile(color = "lightgray") +
  scale_fill_gradient2(low = "chartreuse3", mid="snow", high = "orchid3", midpoint = (min(melted_cormat2$value, na.rm = T)+max(melted_cormat2$value, na.rm = T))/2, space = "Lab", na.value="lightgray",
                        name="distance") +
  coord_fixed() +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.x = element_text(angle = 0, vjust = 1, size = 14, hjust = 0.5),
        axis.text.y = element_text(angle = 0, vjust = 0.5, size = 14, hjust = 1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

#### PERMUTATION TEST #####
newModel0 <- newModel
Md_fin1_perm <- list()
Md_fin2_perm <- list()
Md_fin3_perm <- list()
Md_fin4_perm <- list()
Md_fin5_perm <- list()
Md_fin6_perm <- list()

#permutations
for(j in 1:500){
Md_l1_perm <- list()
#Md_l2_perm <- list()
Md_l3_perm <- list()
#Md_l4_perm <- list()
#Md_l5_perm <- list()
#Md_l6_perm <- list()

set.seed(j)
#newModel0[,3:8] <- newModel0[sample(1:length(newModel0$decision)),3:8] # need to try with only decision
newModel0[,3] <- newModel0[sample(1:length(newModel0$decision)),3] 
for(i in 1:20){
  pddnos_n <- sample(which(as.character(newModel0$decision)=="PDD-NOS"), rls_min)
  asperg_n <- sample(which(as.character(newModel0$decision)=="ASPERGER'S DISORDER"), rls_min)
  autism_n <- sample(which(as.character(newModel0$decision)=="AUTISM"), rls_min)
  control_n <- sample(which(as.character(newModel0$decision)=="CONTROL"), rls_min)
  newModel2 <- newModel0[c(pddnos_n, asperg_n, autism_n, control_n), ]
  
  net0 <- network_gen(newModel2, decision="PDD-NOS", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 
  net1 <- network_gen(newModel2, decision="ASPERGER'S DISORDER", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 
  net2 <- network_gen(newModel2, decision="AUTISM", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 
  net3 <- network_gen(newModel2, decision="CONTROL", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0, minAcc=0.25, type = 'RDF', topN = 100) 
  
  vis_out <- list(net0[[1]], net1[[1]], net2[[1]], net3[[1]])
  names(vis_out) <- c("PDD-NOS", "ASPERGER'S DISORDER", "AUTISM", "CONTROL")
  #vis_out <- visunet(newModel2, type = "RDF")
  
  #### create adjacency matrices
  ams <- createAdjMat(vis_out)
  
  Md1_perm <- nd.centrality(ams, mode="Degree", directed = FALSE)
  #Md2_perm <- nd.centrality(ams, mode="Close", directed = FALSE)
  Md3_perm <- nd.centrality(ams, mode="Between", directed = FALSE)
  #Md4_perm <- nd.hamming(ams, out.dist = TRUE)
  #Md5_perm <- nd.gdd(ams, out.dist = TRUE)
  #Md6_perm <- nd.dsd(ams, out.dist = TRUE, type="Adj")
  
  Md1_perm <- as.matrix(Md1_perm$D)
  #Md2_perm <- as.matrix(Md2_perm$D)
  Md3_perm <- as.matrix(Md3_perm$D)
  #Md4_perm <- as.matrix(Md4_perm$D)
  #Md5_perm <- as.matrix(Md5_perm$D)
  #Md6_perm <- as.matrix(Md6_perm$D)
  
  Md_l1_perm[[i]] <- Md1_perm
  #Md_l2_perm[[i]] <- Md2_perm
  Md_l3_perm[[i]] <- Md3_perm
 # Md_l4_perm[[i]] <- Md4_perm
  #Md_l5_perm[[i]] <- Md5_perm
  #Md_l6_perm[[i]] <- Md6_perm
  
}
#Md2_perm[[j]] <- Reduce("+", Md_l2)/length(Md_l2)

Md_fin1_perm[[j]] <- Reduce("+", Md_l1_perm)/length(Md_l1_perm)
#Md_fin2_perm[[j]] <- Reduce("+", Md_l2_perm)/length(Md_l2_perm)
Md_fin3_perm[[j]] <- Reduce("+", Md_l3_perm)/length(Md_l3_perm)
#Md_fin4_perm[[j]] <- Reduce("+", Md_l4_perm)/length(Md_l4_perm)
#Md_fin5_perm[[j]] <- Reduce("+", Md_l5_perm)/length(Md_l5_perm)
#Md_fin6_perm[[j]] <- Reduce("+", Md_l6_perm)/length(Md_l6_perm)

print(j)
}
Md_fin_perm <- Md_fin1_perm
Md_fin <- Md_fin1
###########################
saveRDS(Md_fin1,"/Users/mateuszgarbulowski/Desktop/ASD_subtypes_permutation_test/ori_cent_degree.rds")
saveRDS(Md_fin3,"/Users/mateuszgarbulowski/Desktop/ASD_subtypes_permutation_test/ori_cent_between.rds")

saveRDS(Md_fin1_perm,"/Users/mateuszgarbulowski/Desktop/ASD_subtypes_permutation_test/perm_cent_degree.rds")
saveRDS(Md_fin3_perm,"/Users/mateuszgarbulowski/Desktop/ASD_subtypes_permutation_test/perm_cent_between.rds")

Md_fin <- readRDS("/Users/mateuszgarbulowski/Desktop/ASD_subtypes_permutation_test/subtype/ori_cent_between.rds")
Md_fin_perm <- readRDS("/Users/mateuszgarbulowski/Desktop/ASD_subtypes_permutation_test/subtype/perm_cent_between.rds")
nams0 <- c("PDD-NOS","AS","autism","control")
plotNetPerm <- function(crds){
  nms2 <- nams0 
  nms3 <- paste0(nms2[crds], collapse = " vs ")
  crds2 <- rev(crds)
  vals <- data.frame(values = unlist(lapply(Md_fin_perm, function(x) x[crds[1],crds[2]])))
  thrt <- Md_fin[crds2[1],crds2[2]]
  colnames(vals) <- "values"
  # We add a 1 to the numerator and denominator to account for misestimation of the p-value 
  # (for more details see Phipson and Smyth, Permutation P-values should never be zero).
  
  #the proportion of permutations with larger difference
  # two tailed p-value:
  pval <- round((sum(vals <= thrt)+1/500)/(dim(vals)[1]+1/500), digits = 3)
  pval <- format.pval(pval, digits = 3, eps = 0.001, nsmall = 2)
  ggplot(vals, aes(x = values)) + 
    geom_histogram(aes(y=..density..), bins = 25, col="gray35", fill="gray70")+
    geom_density(alpha=.3, fill="tomato", col="tomato3") + 
    geom_vline(aes(xintercept=thrt), color="gray30", linetype="dashed", size=1) +
    ggtitle(nms3)+
    xlab("scaled distance")+
    theme_classic()+
    labs(subtitle = paste0("p-value: ",pval))
}

p1 <- plotNetPerm(c(1,2))
p2 <- plotNetPerm(c(1,3))
p3 <- plotNetPerm(c(2,3))
p4 <- plotNetPerm(c(3,4))
p5 <- plotNetPerm(c(1,4))
p6 <- plotNetPerm(c(2,4))

ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, labels = c("a", "b","c","d","e","f"),
          common.legend = TRUE, legend = "bottom")
