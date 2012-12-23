##################################################################################
### Purpose: Generate latex code for the stdudy-specific measure table
### Input: the object "multipletables"
### Output:  Latex code
### Author: Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data: 7/13/2012
##################################################################################
xtable.multipletables<-function(x,caption = NULL, label = NULL, align = NULL,
    digits = NULL, display = NULL,...) {
  if (!inherits(x, "multipletables"))
    stop("Use only with 'multiple' objects.\n")
  measure <- x$measure
  model <- x$model
  reports <- study_specifc(x)
  tables <- as.data.frame(reports)
  table.latex <- xtable(tables,label=paste("Study-specific",measure),...)
  return(table.latex)
}
