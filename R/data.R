#' Pantheon dataset
#'
#' Dataframe of globally known people according to MIT's Pantheon 1.0.
#' Dataset includes number of globally known people for a \strong{sample} 
#' of 10 countries and 53 different occupations. 
#' The \strong{complete} dataset is described in [Yu et al., 2015].
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
#' A. Z. Yu, S. Ronen, K. Hu, T. Lu, and C. A. Hidalgo, 'Pantheon: A Dataset for the Study of Global Cultural Production', arXiv:1502.07310 [physics], Feb. 2015.
#' @keywords dataset
#' @examples
#' data(pantheon)
#' str(pantheon)
#' summary(pantheon)
#' pantheon[pantheon$Country=="Chile",]
"pantheon"

#' Geese dataset
#'
#' A matrix of species of geese.
#' Dataset includes the quantity of 4 species of geese observed by year in the Netherlands.
#' Data was curated by the Dutch bird protection organisation Sovon.
#'
#' @format A matrix with variables:
#' \describe{
#' \item{Columns}{Year of observation}
#' \item{Rows}{Specie}
#' }
#'
#' @source \url{https://www.sovon.nl/en}
#' @source \url{http://www.compass-project.eu/applets/3/index_EN.html}
#' @keywords dataset
#' @examples
#' data(geese)
#' str(geese)
#' summary(geese)
#' geese[,"2000"]
#' geese["Mute Swan",]
"geese"