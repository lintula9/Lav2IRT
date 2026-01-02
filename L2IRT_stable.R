# Lavaan to IRT stable.

# Libraries
package_names <- c("lavaan", "mirt", "pbapply", "RColorBrewer")

for (i in package_names){
  if ( !requireNamespace( i, 
                          quietly = T )) {
    message( "\n\n\t\tNo package ", i, ". Installing package ", i, ".\n\n")
    install.packages( i )
  }
  library( i, character.only = TRUE )}


# Lavaan to IRT 20.3.2024.

LavaanIRTProbabilities <- function( lavaanfit, # Probability of some latent factor level, given the model and observing some value of some item.
                                    varname, 
                                    dimname,
                                    dimmin = -6, 
                                    dimmax = 6, 
                                    marginalize = F,
                                    silent = F,
                                    std = F) {
  
  if( std ) {output = inspect( object = lavaanfit, what = "std" )
  } else { output = inspect( object = lavaanfit, what = "est" ) }
  
  # Error messages.
  if ( !(varname %in% rownames( output$lambda )) ) stop( paste( varname, "not found in lavaan object." ) )
  if ( lavaanfit@Options$mimic != "Mplus") stop( "Please use mimic = 'Mplus' argument in lavaan." )
  if ( lavaanfit@Options$parameterization != "theta") stop( "Please use parameterization = 'theta' argument in lavaan." )
  if ( lavaanfit@Options$std.lv != T ) stop( "Please use std.lv = TRUE argument in lavaan." )
  
  # Pick parameters.
  itemloading = output$lambda[ which( rownames( output$lambda ) == varname ), dimname ]
  itemthresholds = output$tau[ grep( pattern = paste("^",varname,"\\b",sep=""), x = rownames( output$tau ) ) ]
  itemtheta = output$theta[ varname, varname ]
  factorcorr = output$psi
  
  itemloc = which( lavaanfit@Data@ov.names[[1]] == varname )
  itemlevels = as.character( 1 : ( length( itemthresholds ) + 1 ))
  nCat <- length( itemlevels ) 
  
  factormean = output$alpha[ which( rownames( output$alpha ) == dimname ) ]
  factorvar = output$psi[ dimname, dimname ]
  
  factorX = seq(dimmin, dimmax, .01)
  
  if ( marginalize & !std ) {
    specific_corr = factorcorr[ -grep(dimname, colnames(factorcorr)) , -grep(dimname, rownames(factorcorr)) ]
    itemloadings_specific = output$lambda[ which( rownames( output$lambda ) == varname ),
                                           which( colnames( output$lambda ) != dimname ) ]
    MarginalizingConstant = sqrt( 1 + as.numeric( t(itemloadings_specific) %*% specific_corr %*% itemloadings_specific ) )
    
    # Item loading is changed to marginalized item loading, and itemthresholds are changed similarly.
    itemloading = itemloading / MarginalizingConstant
    itemthresholds = itemthresholds / MarginalizingConstant
    
  } 
  if( marginalize & std  ) {
    warning(" ------------------- Marginalization is not currently supported for a standardized solution. ------------")
    
    itemloadings_specific = output$lambda[ which( rownames( output$lambda ) == varname ),
                                           which( colnames( output$lambda ) != dimname ) ]
    MarginalizingConstant = sqrt( itemtheta + as.numeric( t(itemloadings_specific) %*% specific_corr %*% itemloadings_specific ) )
    
    # Item loading is changed to marginalized item loading. itemthresholds are changed similarly.
    itemloading = itemloading / MarginalizingConstant
    itemthresholds = itemthresholds / MarginalizingConstant
  }
  
  if( !marginalize & !silent) {message("Unmarginalized loadings are used.")}
  
  # Item Y indicates P( latentVariable = at some level | X is some category ):
  ProbTheta <- matrix(ncol = nCat, nrow = length(factorX))
  ProbTheta[ , 1                ] <- pnorm( q =  ( itemthresholds[ 1 ] - itemloading * factorX ) / sqrt( itemtheta ) ) # First category probability.
  ProbTheta[ , 2 : ( nCat - 1 ) ] <- sapply( 2 : ( nCat - 1 ), FUN = function( j ) { # Middle category probabilities.
    pnorm( q =( itemthresholds[ j ] - itemloading * factorX ) / sqrt( itemtheta ) )  - pnorm( q = ( itemthresholds[ j - 1 ] - itemloading * factorX ) / sqrt( itemtheta ) ) } ) 
  ProbTheta[ , nCat             ] <- ( 1 - pnorm( q =  ( itemthresholds[ nCat - 1 ] - itemloading * factorX ) / sqrt( itemtheta ) ) ) # Last category probability.
  
  attr(ProbTheta, "Variable name") <- varname
  attr(ProbTheta, "Thresholds") <- itemthresholds
  attr(ProbTheta, "Loading") <- itemloading
  attr(ProbTheta, "lvarname") <- dimname
  attr(ProbTheta, "Dimension interval") <- range(factorX)
  attr(ProbTheta, "ThetaParameter") <- output$theta[ varname, varname ]
  attr(ProbTheta, "Plotting method") <- "ICC"
  attr(ProbTheta, "IsProbMat") <- T
  
  class(ProbTheta) <- "Lav2IRT"
  
  return( ProbTheta )
}

L2IRTPointP <- function( lavaanfit, # Probability of some latent factor level, given the model and observing some value of some item.
                         varname, 
                         dimname,
                         points = 0,
                         marginalize = F,
                         silent = F,
                         std = F ) {
  
  if( std ) { output = inspect( object = lavaanfit, what = "std" )
  } else { output = inspect( object = lavaanfit, what = "est" ) }
  
  if ( !(varname %in% rownames( output$lambda )) ) stop( paste( varname, "not found in lavaan object." ) )
  if ( lavaanfit@Options$mimic != "Mplus") stop( "Please use mimic = 'Mplus' argument in lavaan." )
  if ( lavaanfit@Options$parameterization != "theta") stop( "Please use parameterization = 'theta' argument in lavaan." )
  if ( lavaanfit@Options$std.lv != T ) stop( "Please use std.lv = TRUE argument in lavaan." )
  
  itemloading = output$lambda[ which( rownames( output$lambda ) == varname ), dimname ]
  itemthresholds = output$tau[ grep( pattern = paste("^",varname,"\\b",sep=""), x = rownames( output$tau ) ) ]
  itemtheta = output$theta[ varname, varname ]
  factorcorr = output$psi
  
  itemloc = which( lavaanfit@Data@ov.names[[1]] == varname )
  itemlevels = as.character( 1 : ( length( itemthresholds ) + 1 ))
  nCat <- length( itemlevels ) 
  
  factormean = output$alpha[ which( rownames( output$alpha ) == dimname ) ]
  factorvar = output$psi[ dimname, dimname ]
  
  if ( marginalize & !std ) {
    
    specific_corr = factorcorr[ -grep(dimname, colnames(factorcorr)) , -grep(dimname, rownames(factorcorr)) ]
    itemloadings_specific = output$lambda[ which( rownames( output$lambda ) == varname ),
                                           which( colnames( output$lambda ) != dimname ) ]
    MarginalizingConstant = sqrt( 1 + as.numeric( t(itemloadings_specific) %*% specific_corr %*% itemloadings_specific ) )
    
    # Item loading is changed to marginalized item loading, and item thresholds are changed similarly.
    itemloading = itemloading / MarginalizingConstant
    itemthresholds = itemthresholds / MarginalizingConstant
    
  } else if( marginalize & std  ) {
    warning(" ------------------- Marginalization is not currently supported for a standardized solution. ------------")
    
    
    itemloadings_specific = output$lambda[ which( rownames( output$lambda ) == varname ),
                                           which( colnames( output$lambda ) != dimname ) ]
    MarginalizingConstant = sqrt( itemtheta + as.numeric( t(itemloadings_specific) %*% specific_corr %*% itemloadings_specific ) )
    
    # Item loading is changed to marginalized item loading. itemthresholds are changed similarly.
    itemloading = itemloading / MarginalizingConstant
    itemthresholds = itemthresholds / MarginalizingConstant
  }
  
  # Item Y indicates P( latentVariable = at some level | X is some category ):
  ProbTheta <- matrix(ncol = nCat, nrow = length(points))
  ProbTheta[ , 1                ] <- pnorm( q =  ( itemthresholds[ 1 ] - itemloading * points ) / sqrt( itemtheta ) ) # First category probability.
  ProbTheta[ , 2 : ( nCat - 1 ) ] <- sapply( 2 : ( nCat - 1 ), FUN = function( j ) { # Middle category probabilities.
    pnorm( q =( itemthresholds[ j ] - itemloading * points ) / sqrt( itemtheta ) )  - pnorm( q = ( itemthresholds[ j - 1 ] - itemloading * points ) / sqrt( itemtheta ) ) } ) 
  ProbTheta[ , nCat             ] <- ( 1 - pnorm( q =  ( itemthresholds[ nCat - 1 ] - itemloading * points ) / sqrt( itemtheta ) ) ) # Last category probability.
  
  attr(ProbTheta, "Variable name") <- varname
  attr(ProbTheta, "Thresholds") <- itemthresholds
  attr(ProbTheta, "Loading") <- itemloading
  attr(ProbTheta, "lvarname") <- dimname
  attr(ProbTheta, "Dimension interval") <- range(points)
  attr(ProbTheta, "ThetaParameter") <- output$theta[ varname, varname ]
  attr(ProbTheta, "Plotting method") <- "ICC"
  attr(ProbTheta, "IsProbMat") <- T
  
  class(ProbTheta) <- c("Lav2IRT")
  
  return( ProbTheta )
}

ItemInformation <- function( ProbabilityMatrix ) {
  if(class(ProbabilityMatrix) != "Lav2IRT") stop( paste( "Provided object is not of class Lav2IRT." ))
  if(!(attr(ProbabilityMatrix, "IsProbMat"))) stop( paste( "Provided object is not a Lav2IRT matrix of probabilities." ))
  
  nLevels <- ncol( ProbabilityMatrix )
  lambda <- attr( ProbabilityMatrix , "Loading" )
  
  # Note: theta (i.e., unique factor variance, 'error', is set to 1 with theta parameterization by default.)
  # Hence, including it does not make a difference. Only included currently for future improvements.
  
  theta <- attr( ProbabilityMatrix , "ThetaParameter" ) 
  
  Item_Q <- matrix(ncol = nLevels + 1 , nrow = nrow(ProbabilityMatrix)) # Number of columns is nlevels + 1 because we need '0' level/category.
  Item_Q[ , 1 ] <- 0 # When level is 0.
  Item_Q[ , 2 ] <- ProbabilityMatrix[ , 1] # When level is 1.
  Item_Q[ , 3 : ( nLevels ) ] <-
    sapply( 2:(nLevels - 1), FUN = function( j ) { # For 1 to nlevels - 1. Last column is spared for the last level.
      apply( ProbabilityMatrix[ , 1 : j, drop = F ], MARGIN = 1, sum )
    } )
  Item_Q[ , nLevels + 1 ] <- 1 # When level is last category.
  
  IIC <- (
    
    (3.29 * (lambda^2) / theta) * 
      rowSums(
        matrix( sapply( 2:( ncol( Item_Q )), FUN = function( r ) {
          ( ( Item_Q[ , r, drop = F  ] * ( 1 - Item_Q[ , r, drop = F  ] ) ) - ( Item_Q[ , r -1, drop = F  ] * (( 1 - Item_Q[ , r - 1, drop = F  ] )) ) )^2 / 
            (  ProbabilityMatrix[ , r - 1, drop = F  ] )
        }, simplify = "matrix"), ncol = ncol(Item_Q) - 1 ) ) )
  
  attr(IIC, "Plotting method") <- "IIC"
  attr(IIC, "Dimension interval") <- attr( ProbabilityMatrix, "Dimension interval" )
  attr(IIC, "lvarname") <- attr( ProbabilityMatrix, "lvarname" )
  attr(IIC, "xsequence") <- seq(attr(ProbabilityMatrix, "Dimension interval")[1], 
                                attr(ProbabilityMatrix, "Dimension interval")[2], .1)
  class(IIC) <- "Lav2IRT"
  
  return( IIC )
  
}

SetInformation <- function( lavaanfit, itemset, points = 0, dimname, marginalize = F, std = F) {
  
  info = 0
  for(i in itemset) {
    info = info + ItemInformation(L2IRTPointP(lavaanfit = lavaanfit , varname = i, points = points, 
                                              dimname = dimname, marginalize = marginalize, silent = T, std = std))
  }
  if( any( is.nan(info) ) ) message("Setting NaN values to zero.")
  info[ is.nan(info) ] <- 0 # Set NaN to zero.
  return(info)
}

TestInformation <- function( lavaanfit, dimname ) { # Test information is the sum of all item informations.
  
  output = inspect( object = lavaanfit, what = "est" )
  
  if ( lavaanfit@Options$mimic != "Mplus") stop( "Please use mimic = 'Mplus' argument in lavaan." )
  if ( lavaanfit@Options$parameterization != "theta") stop( "Please use parameterization = 'theta' argument in lavaan." )
  
  testInfo <- rowSums( sapply(lavaanfit@Data@ov.names[[1]], FUN = function( x ) ItemInformation( LavaanIRTProbabilities( lavaanfit, dimname, varname = x ) ) ) )
  attr(testInfo, "Plotting method") <- "testInfo"
  class(testInfo) <- "Lav2IRT"
  
  return(testInfo)
}

