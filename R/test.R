

library(NMGS)
dyn.load("/usr/lib/libgsl.so", local = FALSE, now = FALSE)
dyn.load("/usr/lib/libgslcblas.so", local = FALSE, now = FALSE)
dyn.load("pkg/src/NMGS.so")


#testf(2)
#.Call("safeexp", as.double(10), PACKAGE = "NMGS")
#.Call("Sum", as.integer(1), as.integer(2), PACKAGE = "NMGS")
#.C("Sum")


