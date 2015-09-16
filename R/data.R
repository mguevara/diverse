#' Pantheon data
#'
#' Dataframe of globally known people according to MIT's Pantheon  [!Please Add reference to Amy's Nature dataset!]
#' Dataset includes number of globally known people for a sample of 10 countries.
#' People are associated to [XXX] different occupations. 
#' The dataframe in a edge shape 
#'
#' @format A dataframe with variables:
#' \describe{
#' \item{Country}{Name of the country}
#' \item{Occupation}{Occupation according to Pantheon's taxonomy}
#' \item{Export}{Quantity of globally known people exported born in that country}  #!!!Don't write export people!!!
#' }
#'
#' @source \url{http://pantheon.media.mit.edu/}
#' @examples
#' str(pantheon)
#' summary(pantheon)
#' pantheon[pantheon$Country=='Chile',]
"pantheon"
