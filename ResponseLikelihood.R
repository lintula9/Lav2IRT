# Likelihood curves for response pattern - given a WLSMV estimator result.
# Note, that these are only meaningful for marginalized or unidimensional Item Factor Analysis parameters.
# Also note that the calculations are based on the normal ogive model. Logistic approx are not provided, but should be practically equivalent.
# All observed variables that are provided are assumed to be loaded on the dimension of interest. All residual variances are assumed 1 as default in WLSMV.
# Responses are given so that 1 is the first category - integer valued.

# This is very much in development 22.11.2023.

library(lavaan)

# Define the basic -function for normal ogive model, WLSMV estimator:
basic <- function( tau, lambda, theta) {tau-lambda*theta}



# Define theta arbitrarily, with the assigned name:
theta_range = seq(-6,6,.01)
dimname = "p" # change to param.
theta_name = dimname

# Retrieve tau, lambda
output = inspect( LavaanResult2nuisance, "est" )
lambdas = output$lambda[ , which( colnames( out$lambda ) == dimname ) ]
lambda_specific = matrix(output$lambda[ , which( colnames( output$lambda ) != dimname ) ], ncol = ncol(output$lambda) - 1)
thrshlds = output$tau
varnames = rownames(output$lambda)

# Marginalize?
margC = sqrt( 1 +  diag(lambda_specific %*% t(lambda_specific))  )
lambdas = lambdas / margC
nTau = colSums(sapply(names(lambdas),function(i)grepl(i,rownames(thrshlds)))) # This finds number of categories.
thrshlds = thrshlds / rep(margC, times = nTau)
lambdas = rep(lambdas, times = nTau) # Extend lambda to match tau

# Define the likelihood function: densit of standard normal at a point:
matrix(rep(thrshlds, times = length(theta_range)), nrow = length(thrshlds)) - (lambdas %*% t(theta_range) )
p_f = pnorm( matrix(rep(thrshlds, times = length(theta_range)), nrow = length(thrshlds)) - (lambdas %*% t(theta_range) ) )[ , ]



# This needs refining. The problem is that the likelihoods get very small and become non-comparable (unstandardized).
plot(theta_probability, x = theta_range, type = "l", ylab = paste("Probability of ", dimname ), xlab = dimname)
for(i in 1:10){
  
  response = as.integer(X_2nuisance[i,]) # change to param.
  
theta_probability = 1
for ( i in 1:length(response)) {
  resp = response[ i ]
  varplaces = grep( pattern = paste("^",varnames[ i ],"\\b",sep=""), x = rownames( output$tau ) ) 
  if( resp == 1) {p_r =  p_f[ varplaces , ][ 1, ] }
  if( resp > 1 & resp <= length(varplaces) ) {p_r =  p_f[ varplaces , ][ resp, ] - p_f[ varplaces , ][ resp - 1, ] }
  if( resp == length(varplaces) ) {p_r =  1 - p_f[ varplaces , ][ resp - 1, ] }
  
  theta_probability = theta_probability * p_r 
}
print(pracma::trapz(theta_range,theta_probability))

lines(theta_probability, x = theta_range, type = "l", ylab = paste("Probability of ", dimname ), xlab = dimname)
}
