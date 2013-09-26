

library(NMGS)

dyn.load("/usr/lib/libgsl.so", local = FALSE, now = FALSE)
dyn.load("/usr/lib/libgslcblas.so", local = FALSE, now = FALSE)
dyn.load("pkg/src/NMGS.so")

testf(1,2)
#.Call("compare_doubles", 1, 2, PACKAGE = "NMGS")

