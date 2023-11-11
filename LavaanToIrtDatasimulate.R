# Lavaan to irt simulated data and lavaan fit 

par(family = "serif")

TrueModel <- "
    
    p =~ .6 * X1 + .6 * X2 + .6 * X3

    X1 | -1*t1 + 0*t2 + 1*t3
    X2 | -1*t1 + 0*t2 + 1*t3
    X3 | -1*t1 + 0*t2 + 1*t3

    p ~ 0
    p ~~ 1*p "
X <- tibble(simulateData(TrueModel, model.type = F, 
                         int.lv.free = F, 
                         std.lv = T, sample.nobs = 10000))

HypoModel <-  "
    
    p =~ beta_1 * X1 + beta_2 * X2 + beta_3 * X3

    X1 | t1 + t2
    X2 | t1 + t2
    X3 | t1 + t2

    p ~ 0
    p ~~ 1*p "

LavaanResult <- sem(HypoModel, 
                    data = X, 
                    ordered = T, 
                    estimator = "WLSMV", 
                    mimic = "Mplus", 
                    parameterization = "theta", 
                    std.lv = T)
latentVar <- seq(-6,6,.01)

####            # Run tests.#


testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X3") # Only works manually....
testRes2 <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X3") # Only works manually....

tiff("C:/Users/lintu/OneDrive - University of Helsinki/Väitöskirja/General/Application figures/ApplicationFigure3.tiff",
     family = "serif", width = 6,height = 6,units = "in",res = 480, pointsize = 9)
par(mfrow = c(2,2))
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

comparisonModel <- mirt(X, model = 1, itemtype = "graded" )
plot(comparisonModel, main = "", type = "info")
itemplot(comparisonModel, item = "X3", main = "", type = "info")
