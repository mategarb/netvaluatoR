#johnson
ads_j <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Adults_Results2020-09-11/AllGenes/Johnson/NetworksAllGenes-Johnson2020-09-11.txt")
lietal_j <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Lietal_Results2020-09-11/AllGenes/Johnson/NetworksAllGenes-Johnson2020-09-11.txt")
ped_j <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Pediatric_Results2020-09-11/AllGenes/Johnson/NetworksAllGenes-Johnson2020-09-11.txt")
target_j <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_TARGET_Results2020-09-11/AllGenes/Johnson/NetworksAllGenes-Johnson2020-09-11.txt")

colnames(ads_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
colnames(lietal_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
colnames(ped_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
colnames(target_j) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")

ads_j <- ads_j[grep(",", ads_j$features),]
lietal_j <- lietal_j[grep(",", lietal_j$features),]
ped_j <- ped_j[grep(",", ped_j$features),]
target_j <- target_j[grep(",", target_j$features),]

#genetic
ads_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Adults_Results2020-09-11/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-11.txt")
lietal_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Lietal_Results2020-09-11/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-11.txt")
ped_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_Pediatric_Results2020-09-11/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-11.txt")
target_g <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/Merged Models Features/Merged_TARGET_Results2020-09-11/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-11.txt")
colnames(ads_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
colnames(lietal_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
colnames(ped_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")
colnames(target_g) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")

ads_g <- ads_g[grep(",", ads_g$features),]
lietal_g <- lietal_g[grep(",", lietal_g$features),]
ped_g <- ped_g[grep(",", ped_g$features),]
target_g <- target_g[grep(",", target_g$features),]

ads_g <- ads_g[which(quantile(ads_g$coverageRHS, probs = c(0.05, 0.95))[2] < ads_g$coverageRHS),]
lietal_g <- lietal_g[which(quantile(lietal_g$coverageRHS, probs = c(0.05, 0.95))[2] < lietal_g$coverageRHS),]
ped_g <- ped_g[which(quantile(ped_g$coverageRHS, probs = c(0.05, 0.95))[2] < ped_g$coverageRHS),]
target_g  <- target_g[which(quantile(target_g$coverageRHS, probs = c(0.05, 0.95))[2]<target_g$coverageRHS),]

ads_jg <- rbind(ads_j, ads_g) #ads_j #
lietal_jg <- rbind(lietal_j, lietal_g) #lietal_j #
ped_jg <- rbind(ped_j, ped_g) #ped_j #
target_jg <- rbind(target_j, target_g) #target_j #

local <- rbind(ads_jg, ped_jg) #ads_jg 
li_target <- rbind(lietal_jg, target_jg) #lietal_jg

## aggregate the same
local <- aggregate(.~features+decision, FUN=mean, data = local, na.action = na.pass)
local$supportRHS <- round(local$supportRHS)
li_target <- aggregate(.~features+decision, FUN=mean, data = li_target, na.action = na.pass)
li_target$supportRHS <- round(li_target$supportRHS)

nams <- c("Local P+A", "Li+TARGET")
#nams <- c("Local", "TARGET")
local$decision <- gsub("1","",local$decision)
li_target$decision <- gsub("1","",li_target$decision)


#nams <- c("Local", "Li et al")
#nams <- c("Local", "TARGET")

### balance rules ### ASD groups
rls_min <- min(table(c(paste0(as.character(local$decision),"_local"), paste0(as.character(li_target$decision),"_lietal"))))
Md_l1 <- list()
Md_l2 <- list()
Md_l3 <- list()
Md_l4 <- list()
Md_l5 <- list()
Md_l6 <- list()
source("R/VNtoAdjMat.R")
for(i in 1:30){
  set.seed(i)
  adsd_n <- sample(which(as.character(local$decision)=="Diagnosis"), rls_min)
  adsr_n <- sample(which(as.character(local$decision)=="Relapse"), rls_min)
  lietald_n <- sample(which(as.character(li_target$decision)=="Diagnosis"), rls_min)
  lietalr_n <- sample(which(as.character(li_target$decision)=="Relapse"), rls_min)
  local2 <- local[c(adsd_n, adsr_n), ]
  li_target2 <- li_target[c(lietald_n, lietalr_n), ]
  
  #net0 <- network_gen(local2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
  #net1 <- network_gen(local2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
  #net2 <- network_gen(li_target2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
  #net3 <- network_gen(li_target2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
  
  net0 <- network_gen(local2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
  net1 <- network_gen(local2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
  net2 <- network_gen(li_target2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
  net3 <- network_gen(li_target2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
  
  vis_out <- list(net0[[1]], net1[[1]], net2[[1]], net3[[1]])
  names(vis_out) <- c("Local adults\ndiagnosis", "Local adults\nrelapse", "Li et al.\ndiagnosis", "Li et al.\nrelapse")
  #vis_out <- visunet(newModel2, type = "RDF")
  
  #### create adjacency matrices
  ams <- VNtoAdjMat(vis_out)
  
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

Md_plot <- Md_fin3
Md_plot <- Md3


Md_plot[upper.tri(Md_plot)] <- Inf
nams0 <- c("Local_Diagnosis", "Local_Relapse", "Li_Diagnosis", "Li_Relapse")
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
  local[,2] <- local[sample(1:length(local$decision)),2] 
  li_target[,2] <- li_target[sample(1:length(li_target$decision)),2] 
  for(i in 1:30){
    set.seed(i)
    adsd_n <- sample(which(as.character(local$decision)=="Diagnosis"), rls_min)
    adsr_n <- sample(which(as.character(local$decision)=="Relapse"), rls_min)
    lietald_n <- sample(which(as.character(li_target$decision)=="Diagnosis"), rls_min)
    lietalr_n <- sample(which(as.character(li_target$decision)=="Relapse"), rls_min)
    local2 <- local[c(adsd_n, adsr_n), ]
    li_target2 <- li_target[c(lietald_n, lietalr_n), ]
    
    #net0 <- network_gen(local2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
    #net1 <- network_gen(local2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
    #net2 <- network_gen(li_target2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
    #net3 <- network_gen(li_target2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.25, minAcc=0.75, type = 'L', topN = 100)
    
    net0 <- network_gen(local2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
    net1 <- network_gen(local2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
    net2 <- network_gen(li_target2, decision="Diagnosis", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
    net3 <- network_gen(li_target2, decision="Relapse", ruleSet = 'own', NodeColor = 'DL', Param = 'Min Cov', minValue=0.1, minAcc=0.5, type = 'L', topN = 50)
    
    vis_out <- list(net0[[1]], net1[[1]], net2[[1]], net3[[1]])
    names(vis_out) <- c("Local adults\ndiagnosis", "Local adults\nrelapse", "Li et al.\ndiagnosis", "Li et al.\nrelapse")
    #vis_out <- visunet(newModel2, type = "RDF")
    
    #### create adjacency matrices
    ams <- VNtoAdjMat(vis_out)
    
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
nams0 <- c("Local_Diagnosis", "Local_Relapse", "Li_Diagnosis", "Li_Relapse")
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
