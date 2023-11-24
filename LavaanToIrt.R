# Lavaan to IRT 16.11.2023.
# Recently added dimname and plotting options.
# Recently added marginalizing procedure.

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

  if ( !(varname %in% rownames( output$lambda )) ) stop( paste( varname, "not found in lavaan object." ) )
  if ( lavaanfit@Options$mimic != "Mplus") stop( "Please use mimic = 'Mplus' argument in lavaan." )
  if ( lavaanfit@Options$parameterization != "theta") stop( "Please use parameterization = 'theta' argument in lavaan." )
  if ( lavaanfit@Options$std.lv != T ) stop( "Please use std.lv = TRUE argument in lavaan." )
  
  itemloading = output$lambda[ which( rownames( output$lambda ) == varname ), dimname ]
  itemthresholds = output$tau[ grep( pattern = paste("^",varname,"\\b",sep=""), x = rownames( output$tau ) ) ]
  itemtheta = output$theta[ varname, varname ]
  
  itemloc = which( lavaanfit@Data@ov.names[[1]] == varname )
  itemlevels = as.character( 1 : ( length( itemthresholds ) + 1 ))
  nCat <- length( itemlevels ) 
  
  factormean = output$alpha[ which( rownames( output$alpha ) == dimname ) ]
  factorvar = output$psi[ dimname, dimname ]
  
  factorX = seq(dimmin, dimmax, .01)
  
  if ( marginalize & !std ) {

    
    itemloadings_specific = output$lambda[ which( rownames( output$lambda ) == varname ),
                                   which( colnames( output$lambda ) != dimname ) ]
    MarginalizingConstant = sqrt( 1 + as.numeric( itemloadings_specific %*% itemloadings_specific ) )
    
    # Item loading is changed to marginalized item loading, and itemthresholds are changed similarly.
    itemloading = itemloading / MarginalizingConstant
    itemthresholds = itemthresholds / MarginalizingConstant
    
  } else if( marginalize & std  ) {
    warning(" ------------------- Marginalization is not currently supported for a standardized solution. ------------")

    
    itemloadings_specific = output$lambda[ which( rownames( output$lambda ) == varname ),
                                           which( colnames( output$lambda ) != dimname ) ]
    MarginalizingConstant = sqrt( itemtheta + as.numeric( itemloadings_specific %*% itemloadings_specific ) )
    
    # Item loading is changed to marginalized item loading. itemthresholds are changed similarly.
    itemloading = itemloading / MarginalizingConstant
    itemthresholds = itemthresholds / MarginalizingConstant
  }
  
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

ItemInformation <- function( PorbabilityMatrix ) {
  if(class(PorbabilityMatrix) != "Lav2IRT") stop( paste( "Provided object is not of class Lav2IRT." ))
  if(!(attr(PorbabilityMatrix, "IsProbMat"))) stop( paste( "Provided object is not a Lav2IRT matrix of probabilities." ))
  
  nLevels <- ncol( PorbabilityMatrix )
  lambda <- attr( PorbabilityMatrix , "Loading" )
  
  # Note: theta (i.e., unique factor variance, 'error', is set to 1 with theta parameterization by default.)
  # Hence, including it does not make a difference. Only included currently for future improvements.
  
  theta <- attr( PorbabilityMatrix , "ThetaParameter" ) 
  
  Item_Q <- matrix(ncol = nLevels + 1 , nrow = nrow(PorbabilityMatrix)) # Number of columns is nlevels + 1 because we need '0' level/category.
  Item_Q[ , 1 ] <- 0 # When level is 0.
  Item_Q[ , 2 ] <- PorbabilityMatrix[ , 1] # When level is 1.
  Item_Q[ , 3 : ( nLevels ) ] <-
    sapply( 2:(nLevels - 1), FUN = function( j ) { # For 1 to nlevels - 1. Last column is spared for the last level.
      apply( PorbabilityMatrix[ , 1 : j], MARGIN = 1, sum )
    } )
  Item_Q[ , nLevels + 1 ] <- 1 # When level is last category.
  
  IIC <- (
    
    (3.29 * (lambda^2) / theta) * 
      rowSums(
        sapply( 2:( ncol( Item_Q )), FUN = function( r ) {
          ( ( Item_Q[ , r ] * ( 1 - Item_Q[ , r ] ) ) - ( Item_Q[ , r -1 ] * (( 1 - Item_Q[ , r - 1 ] )) ) )^2 / 
            (  PorbabilityMatrix[ , r - 1 ] )
        }, simplify = "matrix")) )
  
  attr(IIC, "Plotting method") <- "IIC"
  attr(IIC, "Dimension interval") <- attr( PorbabilityMatrix, "Dimension interval" )
  attr(IIC, "lvarname") <- attr( PorbabilityMatrix, "lvarname" )
  attr(IIC, "xsequence") <- seq(attr(PorbabilityMatrix, "Dimension interval")[1], 
                                attr(PorbabilityMatrix, "Dimension interval")[2], .1)
  class(IIC) <- "Lav2IRT"
  
  return( IIC )
  
}

TestInformation <- function( lavaanfit ) { # Test information is the sum of all item informations.
  
  output = inspect( object = lavaanfit, what = "est" )
  
  if ( dim( output$lambda )[ 2 ] > 1 ) stop( "Plots only given for one factor models." )
  if ( lavaanfit@Options$mimic != "Mplus") stop( "Please use mimic = 'Mplus' argument in lavaan." )
  if ( lavaanfit@Options$parameterization != "theta") stop( "Please use parameterization = 'theta' argument in lavaan." )
  
  testInfo <- rowSums( sapply(lavaanfit@Data@ov.names[[1]], FUN = function( x ) ItemInformation( LavaanIRTProbabilities( lavaanfit, varname = x ) ) ) )
  attr(testInfo, "Plotting method") <- "testInfo"
  class("testInfo") <- "Lav2IRT"
  
  return(testInfo)
}

# Plotting methods -----------
setClass("Lav2IRT", contains = "matrix")
setMethod("plot", "Lav2IRT", 
          definition = function( x, ... ) {
            
            if(!(class( x ) == "Lav2IRT")) stop( "This method is for Lav2IRT classes only." )
            
            else if( attr( x, "Plotting method") == "ICC" ) {
              par( mar = c(5, 4, 4, 8), xpd = TRUE)
              matplot( 
                x = seq(attr( x, "Dimension interval")[ 1 ], attr( x, "Dimension interval")[ 2 ], .01),
                y = x, 
                type = "l", main = "Item Characteristic Curves",
                ylab = expression( paste( "P(", theta, ")" ) ), ... )
              legend("right",
                     xjust = 0,
                     inset = c( -0.1, 0 ),
                     legend = as.character( 1 : ncol( x ) ), 
                     col = 1 : ncol( x ), 
                     lty = 1 : ncol( x ),
                     bty = "n", 
                     title = "Item") }
            
            else if( attr(x, "Plotting method") == "IIC" ){
              plot( 
                x = seq(attr(x, "Dimension interval")[ 1 ], attr( x, "Dimension interval")[ 2 ], .01),
                y = x, 
                ylab = expression( paste( "I(", theta, ")" ) ),
                ...)}
            
            
            par( mar = c( 5.1, 4.1, 4.1, 2.1 ), xpd = F)
          })
