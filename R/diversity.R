#' @title  Diversity measures
#' @description It receives an object with data especifying entities (systems), categories (species) and values of presence or abundance, and calculates a number of diversity measures.
#' @param data A numeric matrix with entities as rows and categories as columns and cells as value of abundance. A dataframe with three columns: entities, categories, value of abundance.
#' @param type A mnemonic string referencing the diversity measure. List of available measures: "variety", "entropy", "gini", "simpson", "true", "inverse-simpson", "herfindahl–hirschman","berger-parker", "renyi", "evenness", "rao","rao-stirling". A list of short mnemonics for each measure: 'v', 'e', 'g', 's', 'td', 'is', 'hh', 'bp,'re', ev', 'r',and 'rs'. The default for type is "all". More information for each measure in details and examples. 
#' @param method "rao-stirling" and "rao" measures, use a disparity function to measure the distance between objects. For example: "cosine", "jaccard", "euclidean". The default for method is cosine. All distance measures availables in package proxy.
#' @param agg_type aggregation type for diversity analysis. The analysis is conducted per row, but it can also be conducted by column via setting agg_type = "col". Default is NULL. 
#' @param q parameter for true diversity index measure. This parameter is also used for the Rényi entropy. Default is 0.
#' @param alpha parameter for Rao-Stirling diversity measure. As default we consider alpha=1.
#' @param beta parameter for Rao-Stirling diversity measure. As default we consider beta=1.
#' @details Available diversity measures are:
#' 
#' Variables: N (category count), p_i (proportion of system comprises category i), d_ij (disparity between i and j).
#' 
#' variety: N, category counts per system [MacArthur 1965]
#' 
#' entropy: - sum(p_i log p_i), Shannon entropy per system [Shannon 1948]
#' 
#' gini: 1 - sum(p_i^2), Gini-Simpson index per object [Gini 1912].  It is also known as the Gibbs–Martin index or the Blau index in sociology, psychology and management studies.
#' 
#' simpson: sum(p_i^2), Simpson index per object [Simpson 1949]. This measure is known as the Herfindahl–Hirschman index in economy.  
#' 
#' true: (sum(p_i^q))(1-q)^-1, true diversity index per object [Hill 1973]. This measure is q parameterized. Default for q is 0.  
#' 
#' berger-parker: it is equals to the maximum p_i value in the system, i.e. the proportional abundance of the most abundant type. 
#' 
#' inverse-simpson: (sum(p_i^2))^-1, inverse simpson index per object. This measure is the true diversity at q = 2.
#' 
#' renyi: log((sum(p_i^q))(1-q)^-1), Rényi entropy per object. It is a generalization of the Shannon entropy parameterized by q. It corresponds to the logarithm of the true diversity. Default for q is 0. 
#' 
#' evenness: (-sum(p_i log p_i))/log N, Shannon evenness per object across categories [Pielou, 1969] 
#' 
#' rao-stirling: (sum((d_ij)^alpha (p_i p_j)^beta), Rao-Stirling diversity per object across categories [Stirling, 2007]. As default we consider alpha=1 and beta=1.
#' As pairwise disparities (d_ij) the measure considers Jaccard, Euclidean and Cosine. 
#'  
#' @return A data frame with diversity measures as columns for each object of data
#' @references
#' Gini, C. (1912). "Italian: Variabilità e mutabilità" 'Variability and Mutability', Memorie di metodologica statistica.
#' 
#' Hill, M. (1973). "Diversity and evenness: a unifying notation and its consequences". Ecology 54: 427–432.
#' 
#' MacArthur, R. (1965). "Patterns of Species Diversity". Biology Reviews 40: 510-533.  
#' 
#' Pielou, E. (1969). "An Introduction to Mathematical Ecology". Wiley.
#' 
#' Shannon, C. (1948). "A Mathematical Theory of Communication". Bell System Technical Journal 27 (3): 379–423.
#' 
#' Simpson, A. (1949). "Measurement of Diversity". Nature 163: 41-48.
#' 
#' Stirling, A. (2007). "A General Framework for Analysing Diversity in Science, Technology and Society". Journal of the Royal Society Interface 4: 707-719.
#' @examples
#' #reading csv data matrix
#' path_to_file <- system.file("extdata", "PantheonMatrix.csv", package = "diveR")
#' X <- read.data(path = path_to_file)
#' diversity(data=X, type="gini")
#' diversity(data=X, type="rao-stirling", method="cosine")
#' diversity(data=X, type="all", method="jaccard")
#' 
#' #reading csv dataframe
#' path_to_file <- system.file("extdata", "PantheonEdges.csv", package = "diveR")
#' X <- read.data(path = path_to_file)
#' #true diversity
#' diversity(data=X, type="td", q=1)
#' #rao stirling with differente parameters
#' diversity(data=X, type="rao-stirling", method="euclidean", alpha=0, beta=1)
#' @export
diversity <- function(data, type="all", method='euclidean', agg_type=NULL, q=0, alpha=1, beta=1){
  X <- get_data(data, agg_type)
	diversity <- data.frame(row.names=rownames(X))
	
  if (type == 'variety' || type =='v' || type== 'all') {
    m_d <- as.data.frame(rowSums(X>0, na.rm=TRUE))
    colnames(m_d) <- c('variety')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  propX <- X / rowSums(X, na.rm=TRUE)
  if(type == 'entropy' || type=='e' || type == 'all') { 
    m_d <- as.data.frame(-1 * rowSums(propX * log(propX), na.rm=TRUE))
    colnames(m_d) <- c('entropy')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'gini' || type=='g' || type == 'all') {
    m_d <- as.data.frame(1 - rowSums(propX ^ 2, na.rm=TRUE))
    colnames(m_d) <- c('gini')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'simpson' || type=='s' || type == 'all' || type == 'herfindahl-hirschman' || type == 'hh') {
    m_d <- as.data.frame(rowSums(propX ^ 2, na.rm=TRUE))
    if (type == 'simpson' || type=='s' || type == 'all') {
      colnames(m_d) <- c('simpson')
    }
    else {
      colnames(m_d) <- c('herfindahl-hirschman')
    }
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'true' || type=='td' || type == 'all') {
    p <- 1/(1-q)
    m_d <- as.data.frame((rowSums(propX ^ q, na.rm=TRUE)) ^ p)
    colnames(m_d) <- c('true-diversity')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'berger-parker' || type=='bp' || type == 'all') {
    m_d <- as.data.frame(apply(propX, 1, max))
    colnames(m_d) <- c('berger-parker')
    diversity <- merge(diversity,m_d, by=0, all=TRUE)
    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  if(type == 'inverse-simpson' || type=='is' || type == 'all') {
    q <- 2
    p <- 1/(1-q)
    m_d <- as.data.frame((rowSums(propX ^ q, na.rm=TRUE)) ^ p)
    colnames(m_d) <- c('inverse-simpson')
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
  	disX <- distances(X, agg_type = agg_type, method=method) #compute distances first	
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
  		m_d[,'rao-stirling'] <- NA; ms_label <- 'rao-stirling'	
  	}
  	if(type=='rao' || type=='r')
  	{
  		m_d[,'rao'] <- NA	 ;  ms_label <- 'rao'
  	}
  	if(type=='disparity' || type=='d')
  	{
  		m_d[,'disparity_sum'] <- NA	
  		m_d[,'disparity_mean'] <- NA	
  		ms_label <- 'disparity_sum'
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
  	  	
  	  	if(ms_label=='disparity_sum')
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
  		  
  		  if(ms_label=='disparity_sum')
  		  { 
  		  	m_d[entity, 'disparity_mean'] <- sum(rs_entity)/entity_variety 
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



#' @title Get Data
#' @description It takes data as dataframe (edges) or as matrix (table) to be exported in proper form to be used by the diversity function.
#' @param data Data to be processed as dataframe or as matrix. 
#' @param data_agg Diferent of NULL if column analysis is needed.
#' @examples 
#' X <- get_data(data=d, agg_type=NULL)
#' @export
get_data <- function(data, agg_type)
{
	if (is.data.frame(data)) {
		if(!is.null(agg_type)) {
			#diversity <- data.frame(row.names=levels(data[,2]))
			data <- droplevels(data) #delete un used levels
			X <- matrix(0, nrow=nlevels(data[,2]), ncol=nlevels(data[,1]), dimnames=list(levels(data[,2]),levels(data[,1])))
			X[cbind(data[,2], data[,1])] <- data[,3]
		}
		else {
			#diversity <- data.frame(row.names=levels(data[,1]))
			X <- matrix(0, nrow=nlevels(data[,1]), ncol=nlevels(data[,2]), dimnames=list(levels(data[,1]),levels(data[,2])))
			X[cbind(data[,1], data[,2])] <- data[,3]
		}
	}
	else {
		if (!is.null(agg_type)) {
			X <- t(data)
			#diversity <- data.frame(row.names=rownames(X))
		}
		else {
			X <- data
			#diversity <- data.frame(row.names=rownames(X))
		}  
	}
	return(X)
	
}


#' @title Variety
#' @description It computes the variety or simple diversity of a system. Number of types [!?? What do you mean here with "Number of types"? It seems to be an unfinished sentence!?]
#' @param data Data to be processed as dataframe or as matrix. 
#' @examples 
#' vari <- varity(data=d)
#' @return a dataframe with values of variety
#' @export
variety <- function(data, sort=TRUE)
{
	vari <- diversity(data, type='v')
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
#' @param data_agg Diferent of NULL if column analysis is needed.
#' @examples 
#' ub <- ubiquity(data=d)
#' @return a dataframe with values of frequency per category. Decreasing order
ubiquity <- function(data)
{
	ubiq <- diversity(data, type='v', method='euclidean' , agg_type='col')
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
#' path <-  path_to_matrix_file <- system.file("extdata", "PantheonMatrix.csv", package = "diveR")
#' sep <- ','
#' data <- read.data(path)
#' path <-  path_to_matrix_file <- system.file("extdata", "PantheonEdges.csv", package = "diveR")
#' data <- read.data(path)
read.data <- function(path, type='csv',sep=','){

	if(type=='csv')
	{
		data_temp <- read.csv(path, sep=sep, nrow=1)
		if(ncol(data_temp)>3)
		{#matrix shape
			data <- read.csv(path,sep=sep)
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
		data <- melt(data, na.rm=TRUE)
	}
	row.names(data) <- NULL
  return(data)
}

#' @title A procedure to create a disparity matrix from a data frame or a matrix
#' @description It takes a data frame or a matrix to create a disparity matrix
#' @param data A data frame or matrix of objects x categories
#' @param method List of available disparity methods: "cosine", "jaccard", "euclidean". The default for method is cosine.
#' @param agg_type aggregation type for pairwise disparity measure. The analysis is conducted per row but it can also be performed by column setting agg_type = "col". Default is row. 
#' @return A distance matrix
#' @examples 
#' Xdis <- dist_mat(data)
#' Xdis <- dist_mat(data, method="jaccard", agg_type='col')
distances <- function(data, method='cosine', agg_type=NULL){
    X <- get_data(data=data, agg_type=agg_type)
  
    if (method == 'jaccard') {
      disX <- as.matrix(dist(t(X), method="Jaccard"), diag=1)  
    }
		if (method == 'euclidean') {
      disX <- as.matrix(dist(t(X), method="euclidean"), diag=1)
    }
    if (method == "cosine") {
      disX <- as.matrix(dist(t(X), method="cosine"), diag=1)
    }

  return(disX)
}

#' @title A procedure to compute the sum and average of disparities of systems
#' @description It takes data of presence-abundance of categories in entities and computes the sum and average of disparities between categories PRESENT in the entity
#' @param data A matrix or dataframe of entities, categories and values of presence 
#' @param method a distance measure available in proxy package.
#' @param agg_type aggregation type for disparity analysis. The analysis is conducted per row but it can also be performed by column setting agg_type = "col". Default is "row". 
#' @return A data frame with disparity measures as columns for each entity of data. Sum of disparities and average of disparities are computed.
#' @examples 
#' disp <- disparity(pantheon)
disparity <- function(data, method='cosine', agg_type=NULL) {
  disparity <- diversity(data=data, method=method, type='disparity')
  return(disparity)
}

#' @title Balance
#' @description A procedure to compute several measures of the dimension balance of the diversity.
#' @param data A matrix of data with row and column names. Or a dataframe with three columns entity, category and value
#' @examples 
#' balance(data)
balance <- function(data )
{
	balance <- diversity(data, type='entropy') #first balance measure
	measures <- c( 'gini','simpson', 'berger-parker', 'inverse-simpson', 'evenness' )
	for(measure in measures)
	{
		m_b <- diversity(data, type=measure)
		balance <- merge(balance,m_b, by=0, all=TRUE)
		rownames(balance) <- balance$Row.names; balance$Row.names <- NULL
	}
	
	return(balance)
}



#' @title A procedure to plot the matrix of data as a pheatmap
#' @description It takes a matrix of data and plots a pheatmap of that matrix
#' @param data A matrix of data to be ploted
#' @param fontsize The size of the font used in the plot.
#' @param fontsize_row The font size of the labels in the rows 
#' @examples 
#' pheatmap(data)
plot.pheatmap <- function(data,title=NULL, fontsize=14, fontsize_row=7, color = c('blue','yellow'))
{
	pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE,main=paste(title),show_colnames=FALSE,legend=FALSE,fontsize = fontsize,fontsize_row=7, color=color)
	
}
