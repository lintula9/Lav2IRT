# Libraries
package_names <- c("lavaan", "mirt", "pbapply", "RColorBrewer", "ggplot2", "StanHeaders", "Rcpp", "BH", "RcppEigen", "RcppParallel",
                   "devtools")

for (i in package_names){
  if ( !requireNamespace( i, 
                          quietly = T )) {
    message( "\n\n\t\tNo package ", i, ". Installing package ", i, ".\n\n")
    install.packages( i )
  }
  library( i, character.only = TRUE )}

