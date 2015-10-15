#' @title  \strong{Main} function to compute diversity measures
#' @description \strong{Main} function of the package. It receives an object with data especifying entitie, categories and values of abundance, then, it computes the required diversity measure.
#' @param data A numeric matrix with: entities as rows, categories as columns and cells as value of abundance. It could also be a matrix with categories in rows and entities in columns, but in that case, the paramenter "entity_col" should be set to TRUE. The matrix must include proper names for rows and columns. The parameter data, can also be a dataframe with three columns in this order: entities, categories, value of abundance.
#' @param type A mnemonic string referencing the diversity measure. List of available measures: "variety", "entropy", "gini-simpson", "simpson", "true-diversity", "herfindah-hirschman", "berger-parker", "renyi", "evenness", "rao", "rao-stirling". A list of short mnemonics for each measure: "v", "e", "gs", "s", "td", "hh", "bp", "re", "ev", "r",and "rs". The default for type is "all" which computes all available formulas. More information for each measure in details and examples.
#' @param entity_col Entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @param dis Square matrix of distances or disimilarities between categories. It must include in the rownames the exact name used for each category in the dataset. Only the upper triangle will be used. If this parameter is not defined, and the user requieres a mesure that uses disparities (e.g. Rao), then a matrix of disparities is computed internally using the method defined by the parameter 'method'.
#' @param method "rao-stirling" and "rao" measures, use a disparity function to measure the distance between objects. If the user does not provide a matrix with disparities by using the paramenter 'dis', then a matrix of disparities is computed using the method especified in this parameter (method). For example: "cosine", "jaccard", "euclidean". The default method is cosine. The user can choose, one of the disparity measures availables in package proxy.
#' @param q parameter used for true diversity index. This parameter is also used for the Renyi entropy. Default is 0.
#' @param alpha Parameter for Rao-Stirling diversity. Default is 1.
#' @param beta Parameter for Rao-Stirling diversity. Default is 1.
#' @details  
#' Notation used in the following formulas: \eqn{N},category count; \eqn{p_i}, proportion of entity comprises category \eqn{i}; \eqn{d_{ij}}, disparity between \eqn{i} and \eqn{j};  \eqn{q},\eqn{\alpha} and \eqn{\beta}, parameters.
#' 
#' The available diversity measures included in the package, are listed above. The titles of the formulas, are the possible values that the parameter "type" might take to compute that formula:
#' 
#' 
#' \strong{variety v:}
#'  N, category counts per entity [MacArthur 1965]
#' 
#' 
#' \strong{entropy e:}
#' Shannon entropy per entity [Shannon 1948] \deqn{- \sum_i(p_i \log p_i)}
#' 
#' 
#' \strong{gini-simpson gs:}
#' Gini-Simpson index per object [Gini 1912].  It is also known as the Gibbs-Martin index or the Blau index in sociology, psychology and management studies. This measure is also known as the Herfindahl-Hirschman Index in economy \deqn{1 - \sum_i(p_i^2)}
#' 
#' 
#' \strong{simpson s:}
#' Simpson index per entity [Simpson 1949].   \deqn{ D = \sum_i n_i(n_i-1) / N(N-1)} 
#' When this measure is required, other associated measures are also retrieved, as Simpson's Index of Diversity \eqn{1-D} and the Reciprocal Simpson \eqn{1/D}.
#' 
#' 
#' \strong{true-diversity td:}
#' True diversity index per entity [Hill 1973]. This measure is \eqn{q} parameterized. Default for \eqn{q} is 0.  \deqn{(\sum_ip_{i}^q)^{1/(1-q)}}
#' 
#' 
#' \strong{berger-parker bp:}
#' It is equals to the maximum \eqn{p_i} value in the entity, i.e. the proportional abundance of the most abundant type. 
#'  
#'  
#' \strong{renyi re}
#'  Renyi entropy per object. It is a generalization of the Shannon entropy parameterized by \eqn{q}. It corresponds to the logarithm of the true diversity. Default for \eqn{q} is 0. \deqn{(1-q)^{-1} \log(\sum_i p_i^q)}
#' 
#' 
#' \strong{evenness ev:}
#'  Shannon evenness per object across categories [Pielou, 1969] \deqn{-\sum_i(p_i \log p_i)/\log{N} }
#' 
#' 
#' \strong{rao rao:}
#' Rao diversity. \deqn{\sum_{ij}d_{ij} p_i p_j }
#' 
#' 
#' \strong{rao-stirling rs:}
#'  Rao-Stirling diversity per object across categories [Stirling, 2007]. Default values are \eqn{\alpha=1} and \eqn{\beta=1}.
#' As pairwise disparities (d_ij) the measure considers Jaccard, Euclidean, Cosine or others measures available in the package proxy. \deqn{\sum_{ij}{d_{ij}}^\alpha {(p_i p_j )}^\beta}
#'  
#'  
#'  
#' @return A data frame with diversity measures as columns for each object of data
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
#' path_to_file <- system.file("extdata", "PantheonMatrix.csv", package = "diver")
#' X <- read_data(path = path_to_file)
#' diversity(data=X, type="gini")
#' diversity(data=X, type="rao-stirling", method="cosine")
#' diversity(data=X, type="all", method="jaccard")
#' 
#' #reading csv dataframe
#' path_to_file <- system.file("extdata", "PantheonEdges.csv", package = "diver")
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
  if(type == 'gini-simpson' || type=='gs' || type == 'all'  || type == 'herfindahl-hirschman' || type == 'hh') {
    m_d <- as.data.frame(1 - rowSums(propX ^ 2, na.rm=TRUE))
    #colnames(m_d) <- c('gini.simpson')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  	if (type == 'gini-simpson' || type=='gs' || type == 'all') {
  		colnames(m_d) <- c('gini.simpson')
  	}
  	else {
  		colnames(m_d) <- c('herfindahl.hirschman')
  	}
  }
  if(type == 'simpson' || type=='s' || type == 'all' ) {
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
  	  	#print("prop_i")
  	  	#print(prop_i)
  	  	#print("prop_j")
  	  	#print(prop_j)
  	  	
  	  	if(ms_label=='disparity.sum')
  	  	{
  	  		prop_i[prop_i>0] <- 1 #binarizing proportion
  	  		prop_j[prop_j>0] <- 1 #binarizing proportion
  	  		alpha <- 1; beta<-1
  	  		entity_variety <- length(entity_data[entity_data>0])
  	  		print(paste("entity", entity, 'ent data', entity_data, 'ent variety', entity_variety))
  	  	}
  	  	#print("prop_i")
  	  	#print(prop_i)
  	  	#print("prop_j")
  	  	#print(prop_j)
  	  	p_ij <- prop_i * prop_j
				p_ij_mask <- p_ij
  			#p_ij_mask[p_ij_mask != 0] <- 1
  			p_ij_mask[upper.tri(p_ij_mask)] <- 1
  			#print("p_ij_maks")
  			#print(p_ij_mask)
  		  diag(p_ij_mask) <- 0
  		  p_ij_mask[lower.tri(p_ij_mask)] <- 0
  	  	rs_entity <- ((disX^alpha)*disX_mask) * ((p_ij^beta)*p_ij_mask) #masks ensures for proportions to use only upper triangle matrix. For distances, to use only existant distances and upper triangle matrix.
  		  m_d[entity, ms_label] <- sum(rs_entity)
  		  
  		  if(ms_label=='disparity.sum')
  		  { 
  		  	m_d[entity, 'disparity.mean'] <- sum(rs_entity)/entity_variety 
  		  }
  		  #print(entity)
  		  #print(p_ij)
  		  #print(p_ij_mask)
  	    #print((p_ij^beta)*p_ij_mask)
  	  }
  	#print(rs_entity)
  	#print((p_ij^beta)*p_ij_mask)
  	#print(propX)
  	#print(disX)
  	#print(disX_mask)
  	
      diversity <- merge(diversity,m_d, by=0, all=TRUE)
      rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  
  return(diversity)
}



#' @title Variety or Richeness
#' @description It computes the variety (number of distinct types) or simple diversity of an entity. It is also know as richeness. 
#' @param data A numeric matrix with: entities as rows, categories as columns and cells as value of abundance. It could also be a matrix with categories in rows and entities in columns, but in that case, the paramenter "entity_col" should be set to TRUE. The matrix must include proper names for rows and columns. The parameter data, can also be a dataframe with three columns in this order: entities, categories, value of abundance.
#' @param sort Indicates if results should be ordered or not. Define it to FALSE to avoid ordering.
#' @param entity_col Entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @examples 
#' dim_variety(data=pantheon)
#' dim_variety(data=pantheon, sort=FALSE)
#' @return A dataframe with values of variety for each entity.
#' @export
dim_variety <- function(data, sort=TRUE, entity_col=FALSE)
{
	vari <- diversity(data, type='v', entity_col=entity_col)
	if(sort != FALSE)
	{
		vari['category'] <- row.names(vari)
		vari <- vari[order(vari$variety, decreasing = TRUE), ]
		vari['category'] <- NULL	
	}
	
	return(vari)
}

#' @title Ubiquity of categories across entities
#' @description It computes the ubiquity or the rearnes of the categories
#' @param data A numeric matrix with: entities as rows, categories as columns and cells as value of abundance. It could also be a matrix with categories in rows and entities in columns, but in that case, the paramenter "entity_col" should be set to TRUE. The matrix must include proper names for rows and columns. The parameter data, can also be a dataframe with three columns in this order: entities, categories, value of abundance.
#' @param entity_col Entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @examples 
#' ub <- u_ubiquity(data=pantheon)
#' ub
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
#' @description It reads a file with data shaped as a matrix or as edges list. Several types of formats are allowed.
#' @param path A string representing the path to data file. If it is shaped as a matrix, the first column must include proper names of the categories. If it is shaped as edges list, it must contain three columns, which are entity, category, value. 
#' @param sep Separator field used in the file to separate columns, if it is a CSV file. Default value is comma.
#' @param type It indicates the type of data to be read. This parameter facilitate the input of diverse type of data files, such as spss or stata. Posible options are the names of the mentioned softwares. Default value is csv.
#' @param entity_col Entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @return A data frame with three columns (entity, category, value).
#' @examples 
#' #reading an edges list or panel shape, source data must include three columns
#' path <-  path_to_panel_file <- system.file("extdata", "PantheonEdges.csv", package = "diver")
#' sep <- ','
#' data <- read_data(path)
#' #reading a table
#' path <-  path_to_panel_file <- system.file("extdata", "PantheonMatrix.csv", package = "diver")
#' sep <- ','
#' data <- read_data(path)
#' #reading a table which includes the entities in the columns
#' path <-  path_to_matrix_file <- system.file("extdata", "Geese.csv", package = "diver")
#' data <- read_data(path, entity_col=TRUE)
#' @export
#' @importFrom reshape2 melt
#' @importFrom foreign read.spss read.dta
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
#' @param data A numeric matrix with: entities as rows, categories as columns and cells as value of abundance. It could also be a matrix with categories in rows and entities in columns, but in that case, the paramenter "entity_col" should be set to TRUE. The matrix must include proper names for rows and columns. The parameter data, can also be a dataframe with three columns in this order: entities, categories, value of abundance.
#' @param method A distance or disimilarity method available in "proxy" package. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). If a similarity method is invoqued (as cosine) a proper transformation to disimilarity will be retrieved. Default is Euclidean.
#' @param entity_col Entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @return A distance or disimilarity matrix
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
#' @description It takes data of abundance of categories in entities and computes the sum and average of disparities between categories PRESENT in the entity
#' @param data A numeric matrix with: entities as rows, categories as columns and cells as value of abundance. It could also be a matrix with categories in rows and entities in columns, but in that case, the paramenter "entity_col" should be set to TRUE. The matrix must include proper names for rows and columns. The parameter data, can also be a dataframe with three columns in this order: entities, categories, value of abundance.
#' @param method A distance or disimilarity method available in "proxy" package. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). If a similarity method is invoqued (as cosine) a proper transformation to disimilarity will be retrieved. Default is euclidean.
#' @param entity_col Entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @return method A data frame with resulting disparity measures of each entity in the dataset. Sum of disparities and average of disparities are computed.
#' @examples 
#' dim_disparity(pantheon)
#' @export
dim_disparity <- function(data, method='cosine', entity_col=FALSE) {
  disparity <- diversity(data=data, method=method, type='disparity')
  return(disparity)
}

#' @title Main measures of balance
#' @description A procedure to compute several measures associated to balance or eveness.
#' @param data A numeric matrix with: entities as rows, categories as columns and cells as value of abundance. It could also be a matrix with categories in rows and entities in columns, but in that case, the paramenter "entity_col" should be set to TRUE. The matrix must include proper names for rows and columns. The parameter data, can also be a dataframe with three columns in this order: entities, categories, value of abundance.
#' @param entity_col Entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @examples 
#' dim_balance(pantheon)
#' @export
dim_balance <- function(data, entity_col=FALSE )
{
	balance <- diversity(data, type='entropy', entity_col=entity_col) #first balance measure
	measures <- c( 'gini-simpson','evenness' )
	for(measure in measures)
	{
		m_b <- diversity(data, type=measure, entity_col=entity_col)
		balance <- merge(balance,m_b, by=0, all=TRUE)
		rownames(balance) <- balance$Row.names; balance$Row.names <- NULL
	}
	
	return(balance)
}

#' @title Most common measures used in Ecology to analyze biodiversity
#' @description A procedure to compute the most common measures used to analyze the biodiversity in an ecosystem, such as Berger-Parker, Entropy and Simpson with their variations.
#' @param data A numeric matrix with: entities as rows, categories as columns and cells as value of abundance. It could also be a matrix with categories in rows and entities in columns, but in that case, the paramenter "entity_col" should be set to TRUE. The matrix must include proper names for rows and columns. The parameter data, can also be a dataframe with three columns in this order: entities, categories, value of abundance.
#' @param entity_col entities are in columns. The analysis assumes that in the data matrix, the entities are in rows, but, if in the data matrix the entities are in the columns and the categories in the rows, then, the parameter entity_col, should be set to TRUE. Default is FALSE
#' @return a data frame with measures often used to measure diversity of an ecosystem
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
