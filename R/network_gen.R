network_gen <- function(rules, decision = 'all', type , topN = 0, NodeColor = 'DL', ruleSet = 'ten', Param = 'Min Support', minValue=0, minAcc =0, fract = 0.1 ,CustObjectNodes=list(), CustObjectEdges=list()){
  source('/Users/mateuszgarbulowski/Desktop/VisuNet/R/generate_object.R')
  source('/Users/mateuszgarbulowski/Desktop/VisuNet/R/generateNet.R')
  source('/Users/mateuszgarbulowski/Desktop/VisuNet/R/filtration_rules.R')
  source('/Users/mateuszgarbulowski/Desktop/VisuNet/R/filtration_rules_10per.R')
  source('/Users/mateuszgarbulowski/Desktop/VisuNet/R/data_input.R')
  source('/Users/mateuszgarbulowski/Desktop/VisuNet/R/RDF_columns_test.R')
  source('/Users/mateuszgarbulowski/Desktop/VisuNet/R/Viewrules_type_L.R')
  # rules <- ros$main
  rules <- data_input( rules, type = type)
  decs = unique(as.matrix(rules$decision))
  
  if(ruleSet == 'ten'){
    rules_10per_param <-  filtration_rules_10per(rules, fraction = fract)
    minAcc = rules_10per_param$minAcc
    if(Param == 'Min Support'){
      minValue = rules_10per_param$minSupp
    }else{
      minValue = rules_10per_param$minDecisionCoverage
    }
  }else if(ruleSet == 'own'){
    minValue = minValue
    minAcc = minAcc
  }
  if("pValue" %in% colnames(rules) == FALSE){ rules$pValue = 0.05}
  RulesFiltr =  filtration_rules(rules, minAcc, Param, minValue)
  data_input=generate_object(decs, RulesFiltr,type, topN, Param, NodeColor,  CustObjectNodes, CustObjectEdges)
  out <- list(data <- data_input[[decision]], 
              Acc <- minAcc,
              Val <-  minValue)
  return(out)
}