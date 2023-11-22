# Lavaan to irt simulated data tests and demonstrations for the effect of marginalization procedure implemented in Lav2IRT.



# 1st section: evaluation in scenario with 1 nuisance variable, comparison to mirt package. -----

library(lavaan);library(mirt);par(family = "serif")

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
                         sample.nobs = 10000, orthogonal = T, 
                         standardized = T))
X_ideal <- tibble(simulateData(IdealModel, model.type = F, 
                         int.lv.free = F, 
                         std.lv = T, 
                         sample.nobs = 10000, 
                         standardized = T))
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
LavaanResult_wrong <- sem(HypoModel_ideal, 
                         data = X, 
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
out_wrong = inspect(LavaanResult_wrong, "est")

lambda_start = out$lambda[ , which( colnames( out$lambda ) == "p" ) ]
lambda_specific = out$lambda[ , which( colnames( out$lambda ) != "p" ) ]
thrshlds_start = out$tau
margC = sqrt( 1 + as.numeric( lambda_specific^2 ) )
lambda = lambda_start / margC
thrshlds = thrshlds_start / rep(margC, each = 3)


simures1 <- data.frame(round(cbind(cbind(cbind(thrshlds_start, thrshlds),thrshlds_ideal), out_wrong$tau), 3))
names(simures1) <- c("Conditional", "Marginal", "Ideal", "Wrong")
simures2 <- data.frame(round(cbind(cbind(cbind(lambda_start, lambda), lambda_ideal), out_wrong$lambda), 3))
names(simures2) <- c("Conditional", "Marginal", "Ideal", "Wrong")
rownames(simures2) <- paste("lambda", 1:nrow(simures2))

# The central observation here is that conditional estimates are brought closer to the ideal 'true' estimates 
# this simulation:
simures1;simures2

# The above procedure is done within Lav2IRT functions - in a more general format.


# The specified models are basic 'restricted' 'bi-item-factor' models (the jargon is so thick man).
##### Run tests.#

if(FALSE) {dev.new(noRStudioGD = T);par(mfrow = c(3,3))
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X1", dimname = "p");plot(ItemInformation(testRes), type = "l", main = "Conditional")
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X1", dimname = "p", std = T);plot(ItemInformation(testRes), type = "l", main = "Conditional, standardized")
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X1", dimname = "p", std = T, marginalize = T);plot(ItemInformation(testRes), type = "l", main = "Marginalized, standardized")
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X1", dimname = "p",marginalize = T);plot(ItemInformation(testRes), type = "l", main = "Marginalized")
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult_ideal,  varname = "X1", dimname = "p");plot(ItemInformation(testRes), type = "l", main = "Ideal unidimensional")
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult_wrong,  varname = "X1", dimname = "p");plot(ItemInformation(testRes), type = "l", main = "Wrong unidimensional naive model")
comparisonModel <- mirt(X_ideal, model = 1, itemtype = "graded" );plot(iteminfo(extract.item(comparisonModel,item = "X1"), Theta = latentVar), type = "l", main = "Mirt ideal")
comparisonModel <- mirt(X, model = 1, itemtype = "graded" );plot(iteminfo(extract.item(comparisonModel,item = "X1"), Theta = latentVar), type = "l", main = "Mirt wrong naive model")
par(mfrow = c(1,1))}


if(FALSE){
# Note that mirt package uses different estimator.
# Yet the difference between the marginal information and ideal scenario mirt -package estimated information is negligible:
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X1", dimname = "p",marginalize = T) # 
comparisonModel <- mirt(X_ideal, model = 1, itemtype = "graded" );
plot(ItemInformation(testRes) - iteminfo(extract.item(comparisonModel,item = "X1"), Theta = latentVar),ylim = c(-1,1),lwd=2, type = "l", main = "Difference in information.")
mean(ItemInformation(testRes) - iteminfo(extract.item(comparisonModel,item = "X1"), Theta = latentVar))
# The difference is comparable (due to sampling error, can be even smaller) to 
# difference when using the ideal scenario WLSMV estimator instead of the marginalized WLSMV result.
testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult_ideal,  varname = "X1", dimname = "p");
lines(ItemInformation(testRes) - iteminfo(extract.item(comparisonModel,item = "X1"), Theta = latentVar),x=latentVar, ylim = c(-1,1), type = "l", lty=2,lwd=2)
mean(ItemInformation(testRes) - iteminfo(extract.item(comparisonModel,item = "X1"), Theta = latentVar))
}


# 2nd section, evaluation in 2 nuisance latent variable scenario ------
# Additional test with 2 nuisance latent variables in a complex setting with multiple factors loading on the same item.
# X4 endogenous variable is set to be most complex; loads on 3 factors. 

TrueModel2Nuisance <- "
    # Big general factor
    p =~ .6 * X1 + .6 * X2 + .6 * X3 + .6 * X4 + .6 * X5 + .6 * X6 + .6 * X7 + .6 * X8 + .6*X9

    # Nuisance factor
    nuisance =~ .8*X1 + .3*X2 + .1*X3 + .4*X4
    nuisance2 =~ .4*X4 + .4*X5 + .7*X6 + .3*X7

    X1 | -1*t1 + 0*t2 + 1*t3
    X2 | -1.5*t1 + 0*t2 + 1.5*t3
    X3 | -2*t1 + 0*t2 + 2*t3
    X4 | -1*t1 + 0*t2 + 1*t3
    X5 | -1*t1 + 0*t2 + 1*t3
    X6 | -2*t1 + 0*t2 + 2*t3
    X7 | -1*t1 + 0*t2 + 1*t3
    X8 | -1*t1 + 0*t2 + 1*t3
    X9 | -1*t1 + 0*t2 + 1*t3

    # Ensure orthonormality of factors.
    p ~~ 0*nuisance
    p ~~ 0*nuisance2
    nuisance ~~ 0*nuisance2
    nuisance2 ~~ 1*nuisance2
    nuisance ~~ 1*nuisance
    p ~~ 1*p "
HypoModel2nuisance <-  "
    
    p =~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9

    # Nuisance factor
    nuisance =~ X1 + X2 + X3 + X4
    nuisance2 =~ X4 + X5 + X6 + X7 + X8
    # Ensure orthonormality of factors.
    p ~~ 0*nuisance
    nuisance ~~ 1*nuisance
    nuisance2 ~~ 0*p + 0*nuisance
    nuisance2 ~~ 1*nuisance2
    p ~~ 1*p "


X_2nuisance <- tibble(simulateData(TrueModel2Nuisance, 
                         model.type = F, 
                         int.lv.free = F, 
                         std.lv = T, 
                         sample.nobs = 10000, 
                         orthogonal = T, 
                         standardized = T))

LavaanResult2nuisance <- sem(HypoModel2nuisance, 
                    data = X_2nuisance, 
                    ordered = T, 
                    estimator = "WLSMV", 
                    mimic = "Mplus", 
                    parameterization = "theta", 
                    std.lv = T)

# Calculate lambdas etc.
# Example of marginalization procedure:
out = inspect(LavaanResult2nuisance, "est")
out_ideal = inspect(LavaanResult_ideal, "est")
lambda_ideal = out_ideal$lambda
thrshlds_ideal = out_ideal$tau
out_wrong = inspect(LavaanResult_wrong, "est")

lambda_start = out$lambda[ , which( colnames( out$lambda ) == "p" ) ]
lambda_specific = out$lambda[ , which( colnames( out$lambda ) != "p" ) ]
thrshlds_start = out$tau
margC = sqrt( 1 + as.numeric( rowSums(lambda_specific^2) ) )
lambda = lambda_start / margC
thrshlds = thrshlds_start / rep(margC, each = 3)


simures1 <- data.frame(round(cbind(thrshlds_start, thrshlds), 3))
names(simures1) <- c("Conditional", "Marginal")
simures2 <- data.frame(round(cbind(lambda_start, lambda), 3))
names(simures2) <- c("Conditional", "Marginal")
rownames(simures2) <- paste("lambda", 1:nrow(simures2))

# The central observation here is that conditional estimates are brought closer to the ideal 'true' estimates 
# this simulation:
simures1;simures2
# Note, that unstandardized estimates for the true 0.6 simulation estimate of lambdas are correct, since the simulation model uses different observed variance than
# what is used in estimation for endogenous variables.
# The most pathological variable, X4, which loads onto 3 different orthogonal factors gets a marginal lambda of ~.73, which is acceptable.
# Also the marginal tau estimates are working as intended and resemble those obtained from ideal models.
# Also note, that the differneces between marginal and conditional parameters are LARGE. Conditional parameters are not recommended for interpretation.


# More tests are made in the future.


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
comparisonModel <- mirt(X_ideal, model = 1, itemtype = "graded" )
plot(comparisonModel, main = "", type = "info")
itemplot(comparisonModel, item = "X3", main = "", type = "info")
