# Lavaan to irt simulated data and lavaan fit 
# This needs to be updated to give the user confidence that it gives results identical to mirt package.
# 11.11.2023 - Sakari
library(lavaan);library(mirt)
par(family = "serif")

IdealModel <- "
    
    # Big general factor
    p =~ .6 * X1 + .6 * X2 + .6 * X3 + .6 * X4 + .6 * X5 + .6 * X6

    X1 | -1*t1 + 0*t2 + 1*t3
    X2 | -1.5*t1 + 0*t2 + 1.5*t3
    X3 | -2*t1 + 0*t2 + 2*t3
    X4 | -1*t1 + 0*t2 + 1*t3
    X5 | -1*t1 + 0*t2 + 1*t3
    X6 | -1*t1 + 0*t2 + 1*t3

    p ~~ 1*p "

TrueModel <- "
    
    # Big general factor
    p =~ .6 * X1 + .6 * X2 + .6 * X3 + .6 * X4 + .6 * X5 + .6 * X6

    # Nuisance factor
    nuisance =~ .8 * X1 + .3 * X2 + .1 * X3

    X1 | -1*t1 + 0*t2 + 1*t3
    X2 | -1.5*t1 + 0*t2 + 1.5*t3
    X3 | -2*t1 + 0*t2 + 2*t3
    X4 | -1*t1 + 0*t2 + 1*t3
    X5 | -1*t1 + 0*t2 + 1*t3
    X6 | -1*t1 + 0*t2 + 1*t3

    # Ensure orthonormality of factors.
    p ~~ 0*nuisance
    nuisance ~~ 1*nuisance
    p ~~ 1*p "

X <- tibble(simulateData(TrueModel, model.type = F, 
                         int.lv.free = F, 
                         std.lv = T, 
                         sample.nobs = 10000))
X_ideal <- tibble(simulateData(IdealModel, model.type = F, 
                         int.lv.free = F, 
                         std.lv = T, 
                         sample.nobs = 10000))
HypoModel <-  "
    
    p =~ X1 + X2 + X3 + X4 + X5 + X6

    # Nuisance factor
    nuisance =~ X1 + X2 + X3

    # Ensure orthonormality of factors.
    p ~~ 0*nuisance
    nuisance ~~ 1*nuisance
    p ~~ 1*p "

HypoModel_ideal <-  "
    
    p =~ X1 + X2 + X3 + X4 + X5 + X6

    # Ensure orthonormality of factors.
    p ~~ 1*p "

LavaanResult <- sem(HypoModel, 
                    data = X, 
                    ordered = T, 
                    estimator = "WLSMV", 
                    mimic = "Mplus", 
                    parameterization = "theta", 
                    std.lv = T)
LavaanResult_ideal <- sem(HypoModel_ideal, 
                          data = X_ideal, 
                          ordered = T, 
                          estimator = "WLSMV", 
                          mimic = "Mplus", 
                          parameterization = "theta", 
                          std.lv = T)
latentVar <- seq(-6,6,.01)


# Example of marginalization procedure:
out = inspect(LavaanResult, "est")
out_ideal = inspect(LavaanResult_ideal, "est")
lambda_ideal = out_ideal$lambda
thrshlds_ideal = out_ideal$tau

lambda_start = out$lambda[ , which( colnames( out$lambda ) == "p" ) ]
lambda_specific = out$lambda[ , which( colnames( out$lambda ) != "p" ) ]
thrshlds_start = out$tau
margC = sqrt( 1 + as.numeric( lambda_specific^2 ) )
lambda = lambda_start / margC
thrshlds = thrshlds_start / rep(margC, each = 3)

simures1 <- data.frame(round(cbind(cbind(thrshlds_start, thrshlds),thrshlds_ideal), 3))
names(simures1) <- c("Conditional", "Marginal", "Ideal")
simures2 <- data.frame(round(cbind(cbind(lambda_start, lambda), lambda_ideal), 3))
names(simures2) <- c("Conditional", "Marginal", "Ideal")
rownames(simures2) <- paste("lambda", 1:nrow(simures2))

# The central observation here is that conditional estimates are brought closer to the ideal 'true' estimates 
# this simulation:
simures1;simures2

# The above procedure is done within Lav2IRT functions - in a more general format.


# The specified models are basic 'restricted' 'bi-item-factor' models (the jargon is so thick man).
##### Run tests.#

testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X3", dimname = "p") # 
testRes2 <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X3", dimname = "p", 
                                    marginalize = T) # 

# Tests are made in the future.











# Plots (not updated 21.11.2023) ------
# Assumed scenario:
plot(x=latentVar,as.numeric(ItemInformation(testRes)), ylim = c(0,1), type = "l", 
     ylab = "Information",  main = "A", cex.lab = 1.2, xlab = "General psychopathology dimension", yaxt = "n", xaxt = "n")
lines(x=latentVar + .25,ItemInformation(testRes) / 1.4, lty = 2)
lines(x=latentVar -.25,ItemInformation(testRes) / 1.4, lty = 3)

# Low information scenario
plot(x=latentVar,as.numeric(ItemInformation(testRes)) / 2, ylim = c(0,1), type = "l", 
     ylab = "Information",  main = "B", cex.lab = 1.2, xlab = "General psychopathology dimension", yaxt = "n", xaxt = "n")
lines(x=latentVar + .25,ItemInformation(testRes) / 2.8, lty = 2)
lines(x=latentVar -.25,ItemInformation(testRes) / 2.8, lty = 3)

# High information, varying locations scenario
plot(x=latentVar,as.numeric(ItemInformation(testRes)), ylim = c(0,1), type = "l", 
     ylab = "Information",  main = "C", cex.lab = 1.2, xlab = "General psychopathology dimension", yaxt = "n", xaxt = "n")
lines(x=latentVar + 2.5,ItemInformation(testRes) / 1.4, lty = 2)
lines(x=latentVar -2.5,ItemInformation(testRes) / 1.4, lty = 3)

# Wrong information scenario
plot(x=latentVar,as.numeric(ItemInformation(testRes)) / 2.8, ylim = c(0,1), type = "l", 
     ylab = "Information",  main = "D", cex.lab = 1.2, xlab = "General psychopathology dimension", yaxt = "n", xaxt = "n")
lines(x=latentVar + .25,ItemInformation(testRes) , lty = 2)
lines(x=latentVar -.25,ItemInformation(testRes) , lty = 3)
dev.off()


# mirt -package comparisons for reference. Note, that the estimation method is not exactly the same.
comparisonModel <- mirt(X, model = 1, itemtype = "graded" )
plot(comparisonModel, main = "", type = "info")
itemplot(comparisonModel, item = "X3", main = "", type = "info")
