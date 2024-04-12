# Gradient computations test.

specific_corr = factorcorr[ -grep(dimname, colnames(factorcorr)) , -grep(dimname, rownames(factorcorr)) ]
itemloadings_specific = output$lambda[ which( rownames( output$lambda ) == varname ),
                                       which( colnames( output$lambda ) != dimname ) ]
MarginalizingConstant = sqrt( 1 + as.numeric( t(itemloadings_specific) %*% specific_corr %*% itemloadings_specific ) )

# Item loading is changed to marginalized item loading, and item thresholds are changed similarly.
itemloading = itemloading / MarginalizingConstant
itemthresholds = itemthresholds / MarginalizingConstant



# Toy
C = matrix(c(1,.5,.5,
             .5,1.,.5,
             .5,.5,1), ncol = 3)
lambda = 2
lambda_s = c(1.5,1,0.5)
p = length(lambda_s)
theta = c(lambda, lambda_s, as.vector(C))

MargLamda <- function(theta) {
  
  return( theta[1] / sqrt(1 + (  (theta[2:(1+p)]) %m% t(theta[2:(1+p)])  %m% theta[(1+1+p):(1+p+p^2)] ) )
   )
  }
MargLamdaGrad = makeGradFunc(MargLamda)


MargLamdaGrad(x = c(lambda,lambda_s,as.vector(C)))
