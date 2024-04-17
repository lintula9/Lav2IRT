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
C_lt = C[lower.tri(C, diag = T)]
lambda = 2
lambda_s = c(1.5,1,0.5)
p = length(lambda_s)
theta = c( lambda, lambda_s, C_lt )

MargLamda <- function(theta) {
  S = matrix(0, ncol = p, nrow = p)
  S[ lower.tri(S,diag = T) ] = theta[(1+1+p):(p*(p+1)/2 + 1 + p)]
  S[ upper.tri(S) ] = S[ lower.tri(S) ]
  return( theta[1] / sqrt(1 + (  t(as.vector((theta[2:(1+p)]) %m% t(theta[2:(1+p)])))  %m% as.vector(S) ) ) )
   
  }
MargLamdaGrad = makeGradFunc(MargLamda)


MargLamdaGrad( x = c( lambda, lambda_s, C_lt ) )





numDeriv::grad(function (theta)  { S = matrix(0, ncol = p, nrow = p)
               S[ lower.tri(S,diag = T) ] = theta[(1+1+p):(p*(p+1)/2 + 1 + p)]
               S[ upper.tri(S) ] = S[ lower.tri(S) ]
               return( theta[1] / sqrt(1 + (  t(as.vector((theta[2:(1+p)]) %m% t(theta[2:(1+p)])))  %m% as.vector(S) ) ) )}, 
               x = theta)
               
# Get the vcov of parameter estimates
varnames = lavaan::lavNames(LavaanResult_ideal, type = "ov")
LambdaNames = paste("p","=~", varnames, sep = "")
sfactornames = paste(lavaan::lavNames(LavaanResult_ideal, type = "lv")); sfactornames[ ]

ScovNames = paste(sfactornames, "~~", sfactornames, sep = "")

vcov(LavaanResult_ideal)[ LambdaNames , LambdaNames ]


str(LavaanResult_ideal)
