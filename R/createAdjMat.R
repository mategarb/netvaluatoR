createAdjMat <- function(vis_out, opt="binary"){
  
  outs <- names(vis_out)
  if(length(which(outs=="all"))!=0){
    outs <- outs[-which(outs=="all")]
  }

  vis_out <- vis_out[outs]
  newLevs <- unique(c(as.character(unlist(lapply(vis_out, function(x) x$edges$from))), as.character(unlist(lapply(vis_out, function(x) x$edges$to)))))
  lt <- list()
  for(i in 1:length(outs)){
    
    fromf <- as.factor(as.character(vis_out[[outs[i]]]$edges$from))
    tof <- as.factor(as.character(vis_out[[outs[i]]]$edges$to))
    
    levels(fromf) <- newLevs
    levels(tof) <- newLevs
    
    tbl <- table(fromf, tof)
    tbl[tbl > 1] <- 1
    tbl <- as.matrix(tbl)
    lt[[i]] <- tbl 
  }
  names(lt) <- outs
  return(lt)
}
