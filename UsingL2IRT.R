# Lavaan to Irt example.

source("Libraries.R"); source("LavaanToIRT.R")

# Example 1, bifactor model simulated data.

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

X <- data.frame(lavaan::simulateData(TrueModel, model.type = F, 
                         int.lv.free = F, 
                         std.lv = T, 
                         sample.nobs = 10000, orthogonal = T, 
                         standardized = T))


lavaanresult = lavaan::efa(X,
                           estimator = "WLSMV",
                           parameterization = "theta", # Required
                           rotation = "bigeomin",
                           std.lv = T, # Required
                           ordered = T,
                           mimic = "mplus", # Required
                           nfactors = 2, 
                           rotation.args = list(orthogonal = T))$nf2


# One item

Probs = LavaanIRTProbabilities(lavaanfit = lavaanresult, 
                       varname = "X1", 
                       dimname = "f2", 
                       dimmin = -9, 
                       dimmax = 9, 
                       marginalize = T, 
                       std = F #Default, F, is required for marginalization.
                       )
ItemInfos = ItemInformation(Probs)

plot(ItemInfos, type = "l") 
ggplot() + 
  geom_line(aes(y = as.numeric(ItemInfos), # manual switch to numeric class is currently required.
                x = seq(-9, 9, length.out = length(ItemInfos))))

# Multiple items

setinfos = SetInformation(lavaanresult, 
               itemset = c("X1", "X2"), points = seq(-9,9,by = 0.01), dimname = "f2", marginalize = T, std = F)
setinfos
