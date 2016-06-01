library(devtools)
document("pkg")

# Tricks needed to avoid problems with dynlib
#library(NMGS)
#document("pkg", reload = FALSE, clean = TRUE)
#document("pkg", clean = TRUE)

