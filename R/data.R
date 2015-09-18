#' Pantheon data
#'
#' Dataframe of globally known people according to MIT's Pantheon 1.0.
#' Dataset includes number of globally known people for a sample of 10 countries
#' that are associated to 53 different occupations. 
#' The dataframe provided with diveR package is in shape of edges.
#' The full dataset is described in [Yu et al., 2015].
#'
#' @format A dataframe with variables:
#' \describe{
#' \item{Country}{Name of the country}
#' \item{Occupation}{Occupation according to Pantheon's taxonomy}
#' \item{Value}{Quantity of globally known that waws born in that country} 
#' }
#'
#' @source \url{http://pantheon.media.mit.edu/}
#' @references 
#' A. Z. Yu, S. Ronen, K. Hu, T. Lu, and C. A. Hidalgo, ‘Pantheon: A Dataset for the Study of Global Cultural Production’, arXiv:1502.07310 [physics], Feb. 2015.
#' @examples
#' str(pantheon)
#' summary(pantheon)
#' pantheon[pantheon$Country=='Chile',]
"pantheon"
