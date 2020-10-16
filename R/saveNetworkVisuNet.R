#rules <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/WithoutIntegration/MergePar-Linda_Adults_PCBeforeOldMCFS-5MinRead2020-09-24/AllGenes/Genetic/NetworksAllGenes-Genetic2020-09-24.txt")
#rules <- read.table("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/WithoutIntegration/MergePar-Linda_Adults_PCBeforeOldMCFS-5MinRead2020-09-24/AllGenes/Johnson/NetworksAllGenes-Johnson2020-09-24.txt")

# GENETIC ADULTS
rules <- read.csv("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/WithoutIntegration/MergePar-Linda_Adults_PCBeforeOldMCFS-5MinRead2020-09-24/AllGenes/Genetic/RulesAllGenes-Genetic2020-09-24.csv")
rules <- rules[,-1]

# GENETIC PEDIATRICS
rules <- read.csv("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/WithoutIntegration/MergePar-Linda_Pediatric_PCBefore-5MinRead-ParameterG2020-09-24/AllGenes/Genetic/RulesAllGenes-Genetic2020-09-24.csv")
rules <- rules[,-1]

# JOHNSON ADULTS
rules <- read.csv("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/WithoutIntegration/MergePar-Linda_Adults_PCBeforeOldMCFS-5MinRead2020-09-24/AllGenes/Johnson/RulesAllGenes-Johnson2020-09-24.csv")
rules <- rules[,-1]

# JOHNSON PEDIATRICS
rules <- read.csv("/Users/mateuszgarbulowski/Desktop/AML-RNA-seq/WithoutIntegration/MergePar-Linda_Pediatric_PCBefore-5MinRead-ParameterG2020-09-24/AllGenes/Johnson/RulesAllGenes-Johnson2020-09-24.csv")
rules <- rules[,-1]

hist(rules[which(p.adjust(rules$pValue, method="fdr") <= 0.05),]$pValue)


#colnames(rules) <- c("features","decision","accuracyRHS","supportRHS","coverageRHS","pValue")

source("R/network_gen.R")

dataset <- 'Local adults'
ruleSet = 'ten'
decision = 'all'
height ='1000px'
width = '1500px'

# GENETIC ADULTS
net <- network_gen(rules[which(p.adjust(rules$pValue, method="fdr") <= 0.05),], decision=decision, ruleSet = ruleSet, NodeColor = 'DL', Param = 'Min Coverage', minValue=0.2, minAcc=0.5, fract = 1, type = 'RDF', topN = 20) 

# GENETIC PEDIATRICS
net <- network_gen(rules[which(p.adjust(rules$pValue, method="fdr") <= 0.05),], decision=decision, ruleSet = ruleSet, NodeColor = 'DL', Param = 'Min Coverage', minValue=0.2, minAcc=0.5, fract = 0.5, type = 'RDF', topN = 20) 


# JOHNSON ADULTS
net <- network_gen(rules, decision=decision, ruleSet = ruleSet, NodeColor = 'DL', Param = 'Min Coverage', minValue=0.2, minAcc=0.5, fract = 0.05, type = 'RDF', topN = 20) 

# JOHNSON PEDIATRICS
net <- network_gen(rules, decision=decision, ruleSet = ruleSet, NodeColor = 'DL', Param = 'Min Coverage', minValue=0.2, minAcc=0.5, fract = 1, type = 'RDF', topN = 20) 

# GENETIC ADULTS
net$font.size <- 20
visNetwork(nodes = net[[1]]$nodes, edges = net[[1]]$edges, width = width, main = decision, height = height)%>% 
  visLayout(randomSeed = 123) %>%
  visInteraction(hover = TRUE, navigationButtons = TRUE)  %>% 
  visOptions(selectedBy = list(variable = "group" , multiple = TRUE, main = "Select by decision", style = 'width: 200px; height: 30px;
                               padding-left: 80px;
                               font-size: 15px;
                               color: black;
                               border:none;
                               outline:none;')) %>% 
  visPhysics(solver = "repulsion", repulsion = list(nodeDistance=150)) %>% 
  visSave(file = '/Users/mateuszgarbulowski/Desktop/local_pediatric_top20_johnson.html')

