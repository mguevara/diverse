# diverse
##Description
- Title: Easy to Use Diversity Measures
- Version: 0.1.0
- Authors: Miguel R. Guevara <miguel.guevara@upla.cl>, Dominik Hartmann
        <d.hartmann@uni-hohenheim.de>, Marcelo Mendoza
        <marcelo.mendoza@usm.cl>
- Description: *diverse* package computes the most common diversity measures used in social and natural sciences.
- Depends: proxy, reshape2, foreign
- URL: https://github.com/mguevara/diverse
- Repository: CRAN
- Date: 2015-10-19
- BugReports: https://github.com/mguevara/diverse/issues
- License: CC BY-NC 3.0
- LazyData: true
- NeedsCompilation: no

##Examples
```R
 diversity(pantheon)
 diversity(pantheon, type='variety')
 diversity(geese, type='berger-parker', entity_col=TRUE)
 #reading csv data matrix
 path_to_file <- system.file("extdata", "PantheonMatrix.csv", package = "diver")
 X <- read_data(path = path_to_file)
 diversity(data=X, type="gini")
 diversity(data=X, type="rao-stirling", method="cosine")
 diversity(data=X, type="all", method="jaccard")

 #reading csv dataframe
 path_to_file <- system.file("extdata", "PantheonEdges.csv", package = "diver")
 X <- read_data(path = path_to_file)
 #true diversity
 diversity(data=X, type="td", q=1)
 #rao stirling with differente parameters
 diversity(data=X, type="rao-stirling", method="euclidean", alpha=0, beta=1)
```
