#' Pantheon data
#'
#' Dataframe of globally known people according to MIT's Pantheon
#' Dataset includes number of people globally know for a sample of 10 countries.
#' People is assigned by occupation. This is an example of a dataframe
#' in edges shape.
#'
#' @format A data frame with variables:
#' \describe{
#' \item{Country}{Name of the country}
#' \item{Occupation}{Occupation according to Pantheon's taxonomy}
#' \item{Export}{Quantity of globally known people exported born in that country}
#' }
#'
#' @source \url{http://pantheon.media.mit.edu/}
#' @examples
#' str(pantheon)
#' summary(pantheon)
#' pantheon[pantheon$Country=='Chile',]
"pantheon"
