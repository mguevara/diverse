#' @title  \strong{Main} function to compute diversity measures
#' @description \strong{Main} function of the package. The diversity function computes diversity measures for a dataset with entities, categories and values.
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the parameter "entity_col" has to be set to TRUE. The matrix must include names for the rows and the columns. The parameter "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param type A mnemonic string referencing to the available diversity measures. The available measures are: "variety", (Shannon) "entropy", "gini-simpson", "simpson", "true-diversity", "herfindahl-hirschman", "berger-parker", "renyi", (Shannon) "evenness", "rao", "rao-stirling". A list of short mnemonics for each measure: "v", "e", "gs", "s", "td", "hh", "bp", "re", "ev", "r", and "rs". The default for type is "all" which computes all available formulas. More information for each measure can be found in the sections on 'Details and Examples'.
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @param dis Optional square matrix of distances or dissimilarities between categories. It allows the user to provide her own matrix of dissimilarities between categories. The category names have to be both in the rows and in the columns, and these must be the exact same names used by the categories in the parameter "data". Only the upper triangle will be used. If  the parameter "dis" is not defined, and the user requires a measure that uses disparities (e.g. Rao), then a matrix of disparities is computed internally using the method defined by the parameter 'method'. The default value is NULL.
#' @param method The "rao-stirling" and "rao"-diversity indices use a disparity function to measure the distance between objects. If the user does not provide a matrix with disparities by using the parameter 'dis', then a matrix of disparities is computed using the method especified in this parameter (method). Possible values for this parameter are distance or dissimilarity methods available in "proxy" package as for example "Euclidean", "Kullback" or "Canberra". This parameter also accepts a similarity method available in the "proxy" package, as for example: "cosine", "correlation" or "Jaccard" among others. In the latter case, a correspondent transformation to a dissimilarity measure will be retrieved. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). The default value is Euclidean distance.
#' @param q The parameter used for the true diversity index. This parameter is also used for the Renyi entropy. The default value is 0.
#' @param alpha Parameter for Rao-Stirling diversity. The default value is 1.
#' @param beta Parameter for Rao-Stirling diversity. The default value is 1.
#' @details  
#' Notation used in the following formulas: \eqn{N}, category count; \eqn{p_i}, proportion of entity comprises category \eqn{i}; \eqn{d_{ij}}, disparity between \eqn{i} and \eqn{j};  \eqn{q},\eqn{\alpha} and \eqn{\beta}, parameters.
#' 
#' The available diversity measures included in the package are listed above. The titles of the formulas are the possible mnemonic values that the parameter "type" might take to compute that formula (i.e. diversity(data, type='variety') or diversity(data, type='v'):
#' 
#' 
#' \strong{variety, v:}
#'  N, category counts per entity [MacArthur 1965]
#' 
#' 
#' \strong{entropy, e:}
#' Shannon entropy per entity [Shannon 1948] \deqn{- \sum_i(p_i \log p_i)}
#' 
#' \strong{Herfindahl-Hirschman, hh, hhi:} The Herfindahl-Hirschman Index used in economy to measure the concentration of markets. \deqn{\sum_i(p_i^2)}
#' 
#' 
#' \strong{gini-simpson, gs:}
#' Gini-Simpson index per object [Gini 1912].  This measure is also known as the Gibbs-Martin index or the Blau index in sociology, psychology and management studies.  \deqn{1 - \sum_i(p_i^2)}
#' 
#' 
#' \strong{simpson, s:}
#' Simpson index per entity [Simpson 1949].   \deqn{ D = \sum_i n_i(n_i-1) / N(N-1)} 
#' When this measure is required, then also associated variations Simpson's Index of Diversity \eqn{1-D} and the Reciprocal Simpson \eqn{1/D} will be computed.
#' 
#' 
#' \strong{true-diversity, td:}
#' True diversity index per entity [Hill 1973]. This measure is \eqn{q} parameterized. Default for \eqn{q} is 0.  \deqn{(\sum_ip_{i}^q)^{1/(1-q)}}
#' 
#' 
#' \strong{berger-parker, bp:}
#' Berger-Parker index is equals to the maximum \eqn{p_i} value in the entity, i.e. the proportional abundance of the most abundant type. When this measure is required, the reciprocal measure is also computed.
#'  
#'  
#' \strong{renyi, re:}
#'  Renyi entropy per object. This measure is a generalization of the Shannon entropy parameterized by \eqn{q}. It corresponds to the logarithm of the true diversity index. The default value for \eqn{q} is 0. \deqn{(1-q)^{-1} \log(\sum_i p_i^q)}
#' 
#' 
#' \strong{evenness, ev:}
#'  Shannon evenness per object across categories [Pielou, 1969] \deqn{-\sum_i(p_i \log p_i)/\log{N} }
#' 
#' 
#' \strong{rao:}
#' Rao diversity. \deqn{\sum_{ij}d_{ij} p_i p_j }
#' 
#' 
#' \strong{rao-stirling, rs:}
#'  Rao-Stirling diversity per object across categories [Stirling, 2007]. Default values are \eqn{\alpha=1} and \eqn{\beta=1}.
#' For the pairwise disparities the measure allows to consider the Jaccard Index, Euclidean distances, Cosine Similarity among others. \deqn{\sum_{ij}{d_{ij}}^\alpha {(p_i p_j )}^\beta}
#'  
#'  
#'  
#' @return A data frame with diversity measures as columns for each entity.
#' @references
#' Gini, C. (1912). "Italian: Variabilita e mutabilita" 'Variability and Mutability', Memorie di metodologica statistica.
#' 
#' Hill, M. (1973). "Diversity and evenness: a unifying notation and its consequences". Ecology 54: 427-432.
#' 
#' MacArthur, R. (1965). "Patterns of Species Diversity". Biology Reviews 40: 510-533.  
#' 
#' Pielou, E. (1969). "An Introduction to Mathematical Ecology". Wiley.
#' 
#' Shannon, C. (1948). "A Mathematical Theory of Communication". Bell entity Technical Journal 27 (3): 379-423.
#' 
#' Simpson, A. (1949). "Measurement of Diversity". Nature 163: 41-48.
#' 
#' Stirling, A. (2007). "A General Framework for Analysing Diversity in Science, Technology and Society". Journal of the Royal Society Interface 4: 707-719.
#' @examples
#' diversity(pantheon)
#' diversity(pantheon, type='variety')
#' diversity(geese, type='berger-parker', entity_col=TRUE)
#' #reading csv data matrix
#' path_to_file <- system.file("extdata", "PantheonMatrix.csv", package = "diverse")
#' X <- read_data(path = path_to_file)
#' diversity(data=X, type="gini")
#' diversity(data=X, type="rao-stirling", method="cosine")
#' diversity(data=X, type="all", method="jaccard")
#' 
#' #reading csv dataframe
#' path_to_file <- system.file("extdata", "PantheonEdges.csv", package = "diverse")
#' X <- read_data(path = path_to_file)
#' #true diversity
#' diversity(data=X, type="td", q=1)
#' #rao stirling with differente parameters
#' diversity(data=X, type="rao-stirling", method="euclidean", alpha=0, beta=1)
#' @export
diversity <- function(data, type="all", entity_col=FALSE, dis=NULL, method='euclidean', q=0, alpha=1, beta=1){
  X <- get_data(data, entity_col)
	diversity <- data.frame(row.names=rownames(X))
	
  if (type == 'variety' || type =='v' || type== 'all') {
    m_d <- as.data.frame(rowSums(X>0, na.rm=TRUE))
    colnames(m_d) <- c('variety')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
	sumsX <- rowSums(X, na.rm=TRUE)
  propX <- X / rowSums(X, na.rm=TRUE)
  if(type == 'entropy' || type=='e' || type == 'all') { 
    m_d <- as.data.frame(-1 * rowSums(propX * log(propX), na.rm=TRUE))
    colnames(m_d) <- c('entropy')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
	
	if(type == 'all'  || type == 'herfindahl-hirschman' || type == 'hh' || type == 'hhi'){
		m_d <- as.data.frame(rowSums(propX ^ 2, na.rm=TRUE))
		colnames(m_d) <- c('HHI')
		diversity <- merge(diversity,m_d, by=0, all=TRUE)
		rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	}
  if(type == 'gini-simpson' || type=='gs' || type == 'all') {
    m_d <- as.data.frame(1 - rowSums(propX ^ 2, na.rm=TRUE))
  	colnames(m_d) <- c('gini.simpson')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
    X_simp <- X
  	X_simp[X_simp==0] <- NA
  	m_d <- as.data.frame(rowSums((X_simp*(X_simp-1))/matrix(sumsX*(sumsX-1), ncol=ncol(X_simp), nrow=nrow(X_simp)), na.rm=TRUE)) 
  	colnames(m_d) <- c('simpson.D')
  	m_d['simpson.I'] <- 1-m_d$simpson.D
  	m_d['simpson.R'] <- 1/m_d$simpson.D
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'true' || type=='td' || type == 'all') {
    p <- 1/(1-q)
    m_d <- as.data.frame((rowSums(propX ^ q, na.rm=TRUE)) ^ p)
    colnames(m_d) <- c('true.diversity')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'berger-parker' || type=='bp' || type == 'all') {
    m_d <- as.data.frame(apply(propX, 1, max))
    colnames(m_d) <- c('berger.parker.D')
  	m_d['berger.parker.I'] <- 1/m_d$berger.parker
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'renyi' || type=='re' || type == 'all') {
    p <- 1/(1-q)
    m_d <- as.data.frame(log(rowSums(propX ^ q, na.rm=TRUE) ^ p))
    colnames(m_d) <- c('renyi-entropy')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'evenness' || type=='ev' || type == 'all') {
    m_d <- as.data.frame(-1 * rowSums(propX * log(propX), na.rm=TRUE)/log(rowSums(X>0, na.rm=TRUE))) # if N == 1 -> NaN
    colnames(m_d) <- c('evenness')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'rao-stirling' || type=='rs' || type == 'all' || type=='rao' || type=='r' || type=='disparity' || type=='d'){
  	if(is.null(dis))
  	{
  		disX <- u_distances(X, entity_col = entity_col, method=method) #compute distances first		
  	}
  	else
  	{
  		disX <- dis[,colnames(propX)] #reordering cols of matrix distances
  		disX <- disX[colnames(propX),] #reordering rows of matrix distances 
  	}
  	
  	disX_mask <- disX
  	disX_mask[ (!is.na(disX_mask))] <- 1
  	disX_mask[lower.tri(disX_mask)] <- 0
  	diag(disX_mask) <- 0
  	disX_mask[is.na(disX_mask)] <- 0
  	
  	#print(propX)
  	N <- ncol(propX)
  	m_d <- data.frame(row.names =  rownames(propX))
  	#str(m_d)
  	if(type == 'rao-stirling' || type=='rs' || type=='all')
  	{
  		m_d[,'rao.stirling'] <- NA; ms_label <- 'rao.stirling'	
  	}
  	if(type=='rao' || type=='r')
  	{
  		m_d[,'rao'] <- NA	 ;  ms_label <- 'rao'
  	}
  	if(type=='disparity' || type=='d')
  	{
  		m_d[,'disparity.sum'] <- NA	
  		m_d[,'disparity.mean'] <- NA	
  		ms_label <- 'disparity.sum'
  	}
  	for(entity in row.names(propX)) #go into each entity
  	  {
  	  	entity_data <- propX[entity,]
  	  	prop_i <- matrix(entity_data, nrow=N,ncol=N, byrow=TRUE)
  	  	prop_j <- matrix(entity_data, nrow=N,ncol=N, byrow=FALSE)
  	  	  	  	
  	  	if(ms_label=='disparity.sum')
  	  	{
  	  		prop_i[prop_i>0] <- 1 #binarizing proportion
  	  		prop_j[prop_j>0] <- 1 #binarizing proportion
  	  		alpha <- 1; beta<-1
  	  		entity_variety <- length(entity_data[entity_data>0])
  	  		#print(paste("entity", entity, 'ent data', entity_data, 'ent variety', entity_variety))
  	  	}
  	  	
  	  	p_ij <- prop_i * prop_j
				p_ij_mask <- p_ij
  			p_ij_mask[upper.tri(p_ij_mask)] <- 1
  			
  		  diag(p_ij_mask) <- 0
  		  p_ij_mask[lower.tri(p_ij_mask)] <- 0
  	  	rs_entity <- ((disX^alpha)*disX_mask) * ((p_ij^beta)*p_ij_mask) #masks ensures for proportions to use only upper triangle matrix. For distances, to use only existant distances and upper triangle matrix.
  		  m_d[entity, ms_label] <- sum(rs_entity)
  		  
  		  if(ms_label=='disparity.sum')
  		  { 
  		  	m_d[entity, 'disparity.mean'] <- sum(rs_entity)/entity_variety 
  		  }
  		
  	  }
  
  	
      diversity <- merge(diversity,m_d, by=0, all=TRUE)
      rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  
  return(diversity)
}



#' @title Variety or Richness
#' @description It computes the variety (number of distinct types) or simple diversity of an entity. It is also known as richness. 
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the parameter "entity_col" has to be set to TRUE. The matrix must include names for the rows and the columns. The parameter "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param sort Indicates whether results should be ordered or not. Define it to FALSE to avoid ordering.
#' @param decreasing If parameter "sort" is set to TRUE, this parameter indicates descending order. The default value is TRUE. 
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @examples 
#' dim_variety(data=pantheon)
#' dim_variety(data=pantheon, sort=FALSE)
#' @return A dataframe with values of variety for each entity.
#' @export
dim_variety <- function(data, sort=TRUE, decreasing=TRUE, entity_col=FALSE)
{
	vari <- diversity(data, type='v', entity_col=entity_col)
	if(sort != FALSE)
	{
		vari['category'] <- row.names(vari)
		vari <- vari[order(vari$variety, decreasing = decreasing), ]
		vari['category'] <- NULL	
	}
	
	return(vari)
}

#' @title Ubiquity of categories across entities
#' @description It computes the ubiquity or the rareness of the categories
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the parameter "entity_col" has to be set to TRUE. The matrix must include names for the rows and the columns. The parameter "data", also accepts a dataframe with three columns in the following order: entity, category and value.  
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @examples 
#' ub <- u_ubiquity(data=pantheon)
#' @return A dataframe with values of number of entities where the category is present. Ordered in decreasing order.
#' @export
u_ubiquity <- function(data, entity_col = FALSE)
{
	ubiq <- diversity(data, type='v', method='euclidean' , entity_col= (!entity_col))
	colnames(ubiq) <- 'ubiquity'
	ubiq['category'] <- row.names(ubiq)
	ubiq <- ubiq[order(ubiq$ubiquity, decreasing = TRUE), ]
	ubiq['category'] <- NULL
	return(ubiq)
}

#' @title A procedure to read data of a data file in formats csv, dta or spss
#' @description This function reads a file with data shaped as a matrix or as edges list. Several types of formats are allowed.
#' @param path A string representing the path to data file. If the data contained in the file is shaped as a matrix, the first column must include the names of the categories. If the data is shaped as edges list, it must contain three columns: Entity, category and value. 
#' @param sep Separator character used in the file to separate columns. Only for CSV file. Default value is comma.
#' @param type It indicates the type of data to be read. This parameter facilitates the input of diverse types of data files, such as spss or stata. Possible options are the names of the mentioned software. The default value is csv.
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @return A data frame with three columns (entity, category, value).
#' @examples 
#' #reading an edges list or panel shape, source data must include three columns
#' path <-   system.file("extdata", "PantheonEdges.csv", package = "diverse")
#' sep <- ','
#' data <- read_data(path)
#' #reading a table
#' path  <- system.file("extdata", "PantheonMatrix.csv", package = "diverse")
#' sep <- ','
#' data <- read_data(path)
#' #reading a table which includes the entities in the columns
#' path <- system.file("extdata", "Geese.csv", package = "diverse")
#' data <- read_data(path, entity_col=TRUE)
#' @export
#' @importFrom reshape2 melt
#' @importFrom foreign read.spss read.dta
#' @importFrom utils read.csv
read_data <- function(path, type='csv',sep=',', entity_col=FALSE){

	if(type=='csv')
	{
		data_temp <- read.csv(path, sep=sep, nrow=1)
		if(ncol(data_temp)>3)
		{#matrix shape
			data <- read.csv(path,sep=sep, check.names=FALSE)
			data <- melt(data, na.rm = TRUE)
		}
		else
		{
			col_classes <- c('factor', 'factor', NA)
			data <- read.csv(path, sep=sep, colClasses=col_classes)		
		}
		
	}
		if(type=='spss')
	{
			data <- read.spss(path, use.value.labels = FALSE)
	}
	if(type=='stata')
	{
		data <- read.dta(path)
	}
	if(ncol(data)> 3)
	{
		if(entity_col==TRUE)
		{
			data <- t(data)
		}
			
		data <- melt(data, na.rm=TRUE)
	}
	row.names(data) <- NULL
  return(data)
}

#' @title A procedure to create a disparity matrix from a data frame or a matrix
#' @description It takes a data frame or a matrix to create a disparity matrix
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the parameter "entity_col" has to be set to TRUE. The matrix must include names for the rows and the columns. The parameter "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param method A distance or dissimilarity method available in "proxy" package as for example "Euclidean", "Kullback" or "Canberra". This parameter also accepts a similarity method available in the "proxy" package, as for example: "cosine", "correlation" or "Jaccard" among others. In the latter case, a correspondent transformation to a dissimilarity measure will be retrieved. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). The default value is Euclidean distance.
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @return A distance or dissimilarity square matrix
#' @examples 
#' Xdis <- u_distances(pantheon)
#' Xdis <- u_distances(pantheon, method="jaccard", entity_col=TRUE)
#' Xdis <- u_distances(pantheon, method="cosine", entity_col=TRUE)
#' @export
#' @importFrom proxy dist
u_distances <- function(data, method='euclidean', entity_col=FALSE){
    X <- get_data(data=data, entity_col=entity_col)
	  disX <- as.matrix(dist(t(X), method=method), diag=1) 
  	return(disX)
}

#' @title A procedure to compute the sum and average of disparities of entities
#' @description Computes the sum and the average of disparities between the categories present in each entity.
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the parameter "entity_col" has to be set to TRUE. The matrix must include names for the rows and the columns. The parameter "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param method A distance or dissimilarity method available in "proxy" package as for example "Euclidean", "Kullback" or "Canberra". This parameter also accepts a similarity method available in the "proxy" package, as for example: "cosine", "correlation" or "Jaccard" among others. In the latter case, a correspondent transformation to a dissimilarity measure will be retrieved. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). The default value is Euclidean distance.
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @return A data frame with disparity measures for each entity in the dataset. Both the sum of disparities and the average of disparities are computed.
#' @examples 
#' dim_disparity(pantheon)
#' dim_disparity(data = pantheon, method='Canberra')
#' @export
dim_disparity <- function(data, method='euclidean', entity_col=FALSE) {
  disparity <- diversity(data=data, method=method, type='disparity')
  return(disparity)
}

#' @title Main measures of balance
#' @description  A procedure to compute several measures associated with the balance or evenness of categories.
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the parameter "entity_col" has to be set to TRUE. The matrix must include names for the rows and the columns. The parameter "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @return A data frame that includes the measures of balance: Shannon entropy, Herfindahl-Hirschman Index (HHI), Gini-Simpson (and its derivated measures I, index of diversity and R, reciprocal Simpson) and Shannon evenness. 
#' @examples 
#' dim_balance(pantheon)
#' @export
dim_balance <- function(data, entity_col=FALSE )
{
	balance <- diversity(data, type='entropy', entity_col=entity_col) #first balance measure
	measures <- c( 'hh', 'gini-simpson','evenness' )
	for(measure in measures)
	{
		m_b <- diversity(data, type=measure, entity_col=entity_col)
		balance <- merge(balance,m_b, by=0, all=TRUE)
		rownames(balance) <- balance$Row.names; balance$Row.names <- NULL
	}
	
	return(balance)
}

#' @title Most common measures used in Ecology to analyze biodiversity
#' @description A procedure to compute the most common measures used to analyze the biodiversity in an ecosystem, such as Berger-Parker, Shannon Entropy and Simpson with their variations.
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the parameter "entity_col" has to be set to TRUE. The matrix must include names for the rows and the columns. The parameter "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param entity_col A flag to indicate that entities are in the columns. The analysis assumes that the entities are in the rows of the matrix. If the entities are in the columns and the categories in the rows, then the parameter "entity_col" has to be set to TRUE. The default value is FALSE.
#' @return A dataframe with common measures of diversity in ecosystems.
#' @examples 
#' str(geese)
#' diver_bio(geese)
#' @export
diver_bio <- function(data, entity_col=FALSE)
{
	biodiv <- diversity(data, type='entropy', entity_col=entity_col) #first balance measure
	measures <- c( 'berger-parker','simpson')
	for(measure in measures)
	{
		m_b <- diversity(data, type=measure, entity_col=entity_col)
		biodiv <- merge(biodiv,m_b, by=0, all=TRUE)
		rownames(biodiv) <- biodiv$Row.names; biodiv$Row.names <- NULL
	}
	
	return(biodiv)
}
