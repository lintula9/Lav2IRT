# Libraries
package_names <- c("lavaan", "mirt", "pbapply", "RColorBrewer", "ggplot2", "StanHeaders", "Rcpp", "BH", "RcppEigen", "RcppParallel")

if( !requireNamespace("StanHeaders", quietly = F) ) install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

for (i in package_names){
  if ( !requireNamespace( i, 
                          quietly = T )) {
    message( "\n\n\t\tNo package ", i, ". Installing package ", i, ".\n\n")
    install.packages( i )
  }
  library( i, character.only = TRUE )}

# Graphical settings
cols <- brewer.pal(n = 8, name = "Dark2")[1:3]
names(cols) <- c("Total", "Male", "Female")
theme_set(theme_bw())
options(ggplot2.discrete.colour= cols)
par(family = "serif")


