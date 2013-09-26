# Copyright (C) 2013 Leo Lahti, Christopher Quince, Keith Harris et al.
# Contact: Leo Lahti <leo.lahti@iki.fi>
#
# Licence: FreeBSD


#' testf
#' Test function
#' 
#' @param a A number 
#' @param b Another Number
#'
#' @export
#'
#' @details TBA
#'
#' @return TBA
#' @seealso TBA
#'
#' @references See citation("NMGS") 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples #  TBA
#' @keywords utilities


testf <- function( a, b ){

  qOFz <- .Call("compare_doubles", a, b, PACKAGE = "NMGS")

  qOFz
  
}


