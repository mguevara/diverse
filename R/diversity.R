#' @title  Main function to compute diversity measures
#' @description It receives an object with data especifying entities (entitys), categories (species) and values of presence or abundance, and calculates a number of diversity measures.
#' @param data A numeric matrix with entities as rows and categories as columns and cells as value of abundance. A dataframe with three columns: entities, categories, value of abundance.
#' @param type A mnemonic string referencing the diversity measure. List of available measures: "variety", "entropy", "gini", "simpson", "true-diversity", "herfindah-hirschman", "berger-parker", "renyi", "evenness", "rao", "rao-stirling". A list of short mnemonics for each measure: "v", "e", "g", "s", "td", 'hh', 'bp,'re', ev', 'r',and 'rs'. The default for type is "all". More information for each measure in details and examples.
#' @param dis a square matrix of distances or disimilarities between categories. It must include in the rownames the exact name used for each category in the dataset. Only the upper triangle will be used. If not matrix distance is especified, a matrix of similarities is computed by using the method defined in the parameter method. This for diversity measures that include the dimension of disparity as Rao-Stirling measure. 
#' @param method "rao-stirling" and "rao" measures, use a disparity function to measure the distance between objects. For example: "cosine", "jaccard", "euclidean". The default for method is cosine. All distance measures availables in package proxy.
#' @param entity_col entities are in columns. The analysis assumes that the entities are in rows, but entities could be listed in columns, if that is the case, this parameter should be set to TRUE. Default is FALSE
#' @param q parameter for true diversity index measure. This parameter is also used for the Renyi entropy. Default is 0.
#' @param alpha parameter for Rao-Stirling diversity measure. As default we consider alpha=1.
#' @param beta parameter for Rao-Stirling diversity measure. As default we consider beta=1.
#' @details Available diversity measures are:
#' 
#' Variables: \eqn{N} (category count), \eqn{p_i} (proportion of entity comprises category i), \eqn{d_{ij}} (disparity between \eqn{i} and \eqn{j}).
#' 
#' variety: N, category counts per entity [MacArthur 1965]
#' 
#' entropy: Shannon entropy per entity [Shannon 1948] \deqn{- \sum_i\left(p_i \log p_i\right)}
#' 
#' gini:  Gini-Simpson index per object [Gini 1912].  It is also known as the Gibbs-Martin index or the Blau index in sociology, psychology and management studies. \deqn{1 - \sum_i\left(p_i^2\right)}
#' 
#' simpson: Simpson index per object [Simpson 1949]. This measure is known as the Herfindahl-Hirschman index in economy.  \deqn{ D = \sum_i n_i(n_i-1) / N(N-1)} 
#' 
#' Also the Simpson's Index of Diversity (\eqn{1-D}) and the Reciprocal Simpson \eqn{1/D} are retrieved.
#' 
#' true-diversity: True diversity index per entity [Hill 1973]. This measure is \eqn{q} parameterized. Default for \eqn{q} is 0.  
#' 
#' berger-parker: it is equals to the maximum \eqn{p_i} value in the entity, i.e. the proportional abundance of the most abundant type. 
#'  
#' renyi: Renyi entropy per object. It is a generalization of the Shannon entropy parameterized by \eqn{q}. It corresponds to the logarithm of the true diversity. Default for \eqn{q} is 0. \deqn{\left(1-q\right)^{-1} \log\left(\sum_i p_i^q\right)}
#' 
#' evenness:  Shannon evenness per object across categories [Pielou, 1969] \deqn{-\sum_i\left(p_i \log p_i\right)/\log{N} }
#' 
#' rao-stirling: Rao-Stirling diversity per object across categories [Stirling, 2007]. As default we consider alpha=1 and beta=1.
#' As pairwise disparities (d_ij) the measure considers Jaccard, Euclidean, Cosine or others measures available in the package proxy. \deqn{\sum_{ij}{d_{ij}}^\alpha {\left(p_i p_j \right)}^\beta}
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
#' #reading csv data matrix
#' path_to_file <- system.file("extdata", "PantheonMatrix.csv", package = "diver")
#' X <- read.data(path = path_to_file)
#' diversity(data=X, type="gini")
#' diversity(data=X, type="rao-stirling", method="cosine")
#' diversity(data=X, type="all", method="jaccard")
#' 
#' #reading csv dataframe
#' path_to_file <- system.file("extdata", "PantheonEdges.csv", package = "diver")
#' X <- read.data(path = path_to_file)
#' #true diversity
#' diversity(data=X, type="td", q=1)
#' #rao stirling with differente parameters
#' diversity(data=X, type="rao-stirling", method="euclidean", alpha=0, beta=1)
#' @export
diversity <- function(data, type="all", dis=NULL, method='euclidean', entity_col=FALSE, q=0, alpha=1, beta=1){
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
  if(type == 'gini' || type=='g' || type == 'all'  || type == 'herfindahl-hirschman' || type == 'hh') {
    m_d <- as.data.frame(1 - rowSums(propX ^ 2, na.rm=TRUE))
    colnames(m_d) <- c('gini')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  	if (type == 'gini' || type=='g' || type == 'all') {
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
    colnames(m_d) <- c('rényi-entropy')
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
  		disX <- distances(X, entity_col = entity_col, method=method) #compute distances first		
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
#' @param data Data to be processed as dataframe or as matrix.
#' @param sort Indicates if results should be ordered or not. Define it to FALSE to avoid ordering.
#' @param entity_col The entities are located in rows but, if they are located in columns, then entity_col should be set to TRUE. Default is FALSE 
#' @examples 
#' variety(data=pantheon)
#' variety(data=pantheon, sort=FALSE)
#' @return a dataframe with values of variety
#' @export
variety <- function(data, sort=TRUE, entity_col=FALSE)
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

#' @title Ubiquity
#' @description It computes the ubiquity or the rearnes of the categories
#' @param data Data to be processed as dataframe or as matrix. 
#' @param entity_col The entities are located in rows but, if they are located in columns, then entity_col should be set to TRUE. Default is FALSE
#' @examples 
#' ub <- ubiquity(data=pantheon)
#' ub
#' @return a dataframe with values of frequency per category. Decreasing order
#' @export
ubiquity <- function(data, entity_col = FALSE)
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
#' @return A data frame with three columns, even when the input file is shaped as a matrix.
#' @examples 
#' #reading an edges list or panel shape, source data must include three columns
#' path <-  path_to_panel_file <- system.file("extdata", "PantheonEdges.csv", package = "diver")
#' sep <- ','
#' data <- read.data(path)
#' #reading a table
#' path <-  path_to_panel_file <- system.file("extdata", "PantheonMatrix.csv", package = "diver")
#' sep <- ','
#' data <- read.data(path)
#' #reading a table which includes the entities in the columns
#' path <-  path_to_matrix_file <- system.file("extdata", "Geese.csv", package = "diver")
#' data <- read.data(path, entity_col=TRUE)
#' @export
#' @importFrom reshape2 melt
#' @importFrom foreign read.spss read.dta
read.data <- function(path, type='csv',sep=',', entity_col=FALSE){

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
#' @param data A data frame or matrix of objects x categories
#' @param method List of available disparity methods: "cosine", "jaccard", "euclidean". The default for method is cosine.
#' @param entity_col The entities are located in rows but, if they are located in columns, then entity_col should be set to TRUE. Default is FALSE. 
#' @return A distance matrix
#' @examples 
#' Xdis <- distances(pantheon)
#' Xdis <- distances(pantheon, method="jaccard", entity_col=TRUE)
#' Xdis <- distances(pantheon, method="cosine", entity_col=TRUE)
#' @export
#' @importFrom proxy dist
distances <- function(data, method='euclidean', entity_col=FALSE){
    X <- get_data(data=data, entity_col=entity_col)
	  disX <- as.matrix(dist(t(X), method=method), diag=1) 
  	return(disX)
}

#' @title A procedure to compute the sum and average of disparities of entitys
#' @description It takes data of presence-abundance of categories in entities and computes the sum and average of disparities between categories PRESENT in the entity
#' @param data A matrix or dataframe of entities, categories and values of presence 
#' @param method a distance measure available in proxy package.
#' @param entity_col aggregation type for disparity analysis. The analysis is conducted per row but it can also be performed by column setting entity_col = TRUE. Default is FALSE. 
#' @return A data frame with disparity measures as columns for each entity of data. Sum of disparities and average of disparities are computed.
#' @examples 
#' disp <- disparity(pantheon)
#' @export
disparity <- function(data, method='cosine', entity_col=FALSE) {
  disparity <- diversity(data=data, method=method, type='disparity')
  return(disparity)
}

#' @title Main measures of balance
#' @description A procedure to compute several measures of the dimension balance of the diversity.
#' @param data A matrix of data with row and column names. Or a dataframe with three columns entity, category and value
#' @param entity_col The entities are located in rows but, if they are located in columns, then entity_col should be set to TRUE. Default is FALSE
#' @return a dataset that includes main measures of balance
#' @examples 
#' balance(pantheon)
#' @export
balance <- function(data, entity_col=FALSE )
{
	balance <- diversity(data, type='entropy', entity_col=entity_col) #first balance measure
	measures <- c( 'gini','evenness' )
	for(measure in measures)
	{
		m_b <- diversity(data, type=measure, entity_col=entity_col)
		balance <- merge(balance,m_b, by=0, all=TRUE)
		rownames(balance) <- balance$Row.names; balance$Row.names <- NULL
	}
	
	return(balance)
}

#' @title Biodiversity
#' @description A procedure to compute the most common measures used to analyze the biodiversity of a ecoentity, such as Berger-Parker, Entropy and Simpson with their variations
#' @param data A matrix of data with row and column names. Or a dataframe with three columns entity, category and value
#' @param entity_col The entities are located in rows but, if they are located in columns, then entity_col should be set to TRUE. Default is FALSE
#' @return a data frame with the main measures of biodiversity
#' @examples 
#' str(geese)
#' biodiversity(geese)
#' @export
biodiversity <- function(data, entity_col=FALSE)
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
