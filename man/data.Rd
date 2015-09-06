\name{data}
\alias{data}
\docType{data}
\title{CNT 62}
\description{
A historical SITC data set of production per country across categories for 1962. Data from The MIT Observatory of Economic Complexity (http://atlas.media.mit.edu) based on The Center for International Data from Robert Feenstra. The production is measured in thousands of US dollars.  
}
\usage{load("~/MyDiversity/data/cnt_62.Rda")}
\format{
  A data frame with 35151 observations on the following 3 variables.
  \describe{
    \item{\code{country}}{151 different countries coded using 5 chars: __ (continent), ___ (country). Ex.: asind, indicates Asia - India}
    \item{\code{category}}{689 different production categories coded using 4 chars according to the SITC (Standard International Trade Classification) four digit level standard: _ (group), _ (sector), __ ()}
    \item{\code{value}}{an integer value representing the production of the country in the category in thousands of US dollars.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Data description: https://atlas.media.mit.edu/about/data/sources/
}
\references{
%%
}
\examples{
load("~/MyDiversity/data/cnt_62.Rda")
diversity(data, type="gini")
}
\keyword{datasets}
