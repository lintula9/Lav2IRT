# Lavaan to irt simulated data and lavaan fit 

par()$font
library(lavaan)

TrueModel <- "
    
    F1 =~ .5 * X1 + .8 * X2 + .9 * X3

    X1 | .5*t1 + 1*t2 + 1.5*t3
    X2 | 0*t1 + 0.5*t2 
    X3 | -1*t1 + 0*t2 + .5*t3 + .7*t4

    F1 ~ 0
    F1 ~~ 1*F1 "
X <- tibble(simulateData(TrueModel, model.type = F, 
                         int.lv.free = F, 
                         std.lv = T, sample.nobs = 10000))

HypoModel <-  "
    
    F1 =~ beta_1 * X1 + beta_2 * X2 + beta_3 * X3

    X1 | t1 + t2
    X2 | t1 + t2
    X3 | t1 + t2

    F1 ~ 0
    F1 ~~ 1*F1 "

LavaanResult <- sem(HypoModel, 
                    data = X, 
                    ordered = T, 
                    estimator = "WLSMV", 
                    mimic = "Mplus", 
                    parameterization = "theta", 
                    std.lv = T)

inspect(LavaanResult, "est")

# Run tests.

library(mirt)

comparisonModel <- mirt(X, model = 1, itemtype = "graded" )

testRes <- LavaanIRTProbabilities( lavaanfit = LavaanResult,  varname = "X3")
testRes <- new("Lav2IRT")

itemplot( comparisonModel, item = "X3", main = "" )

plot(testRes)

itemplot(comparisonModel, item = "X3", main = "", type = "info")
plot(ItemInformation(testRes))

plot(comparisonModel, main = "", type = "info")

coef(comparisonModel)
inspect(LavaanResult, "est")
