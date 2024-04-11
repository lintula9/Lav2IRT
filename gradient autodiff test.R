# Gradient computations, test

# R into C++
if(!exists("package_names")) source("Libraries.R")

## update PKG_CXXFLAGS and PKG_LIBS
Sys.setenv(PKG_CXXFLAGS = StanHeaders:::CxxFlags(as_character = TRUE))
SH <- system.file(ifelse(.Platform$OS.type == "windows", "libs", "lib"), .Platform$r_arch, package = "StanHeaders", mustWork = TRUE)
Sys.setenv(PKG_LIBS = paste0(StanHeaders:::LdFlags(as_character = TRUE), " -L", shQuote(SH), " -lStanHeaders"))

Rcpp::sourceCpp("LamdaMargGradient.cpp" )
