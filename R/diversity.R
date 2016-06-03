#' @title  \strong{Main} function to compute diversity measures
#' @description \strong{Main} function of the package. The diversity function computes diversity measures for a dataset with entities, categories and values.
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param type A string or a vector of strings of nemonic strings referencing to the available diversity measures. The available measures are: "variety", (Shannon) "entropy", "blau","gini-simpson", "simpson", "hill-numbers", "herfindahl-hirschman", "berger-parker", "renyi", (Pielou) "evenness", "rao", "rao-stirling". A list of short mnemonics for each measure: "v", "e", "gs", "s", "td", "hh", "bp", "re", "ev", "r", and "rs". The default for type is "all" which computes all available formulas.
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @param dis Optional square matrix of distances or dissimilarities between categories. It allows the user to provide her own matrix of dissimilarities between categories. The category names have to be both in the rows and in the columns, and these must be the exact same names used by the categories in the argument "data". Only the upper triangle will be used. If  the argument "dis" is not defined, and the user requires a measure that uses disparities (e.g. Rao), then a matrix of disparities is computed internally using the method defined by the argument 'method'. The default value is NULL.
#' @param method The "rao-stirling" and "rao"-diversity indices use a disparity function to measure the distance between objects. If the user does not provide a matrix with disparities by using the argument 'dis', then a matrix of disparities is computed using the method specified in this argument (method). Possible values for this argument are distance or dissimilarity methods available in "proxy" package as for example "Euclidean", "Kullback" or "Canberra". This argument also accepts a similarity method available in the "proxy" package, as for example: "cosine", "correlation" or "Jaccard" among others. In the latter case, a correspondent transformation to a dissimilarity measure will be retrieved. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). The default value is Euclidean distance.
#' @param q The parameter used for the hill numbers. This argument is also used for the Renyi entropy and HCDT entropy. The default value is 0.
#' @param alpha Parameter for Rao-Stirling diversity. The default value is 1.
#' @param beta Parameter for Rao-Stirling diversity. The default value is 1.
#' @param base Base of the logarithm. Used in Entropy calculations. The default value is exp(1).
#' @details  
#' Notation used in the following formulas: \eqn{N}, category count; \eqn{p_i}, proportion of entity comprises category \eqn{i}; \eqn{d_{ij}}, disparity between \eqn{i} and \eqn{j};  \eqn{q},\eqn{\alpha} and \eqn{\beta}, arguments.
#' 
#' The available diversity measures included in the package are listed above. The titles of the formulas are the possible mnemonic values that the argument "type" might take to compute that formula (i.e. diversity(data, type='variety') or diversity(data, type='v'):
#' 
#' 
#' \strong{variety, v:}
#'  Category counts per entity [MacArthur 1965] \deqn{\sum_i(p_i^0)}. 
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
#' \strong{hill-numbers, td,hn:}
#' Hill Numbers [Hill 1973]. This measure is \eqn{q} parameterized. When \eqn{q=1}, it results in the exponential of Shannon Entropy. Default for \eqn{q} is 0, this is the variety or richness. \deqn{(\sum_ip_{i}^q)^{1/(1-q)}}
#' 
#' 
#' \strong{berger-parker, bp:}
#' Berger-Parker index is equals to the maximum \eqn{p_i} value in the entity, i.e. the proportional abundance of the most abundant type. When this measure is required, the reciprocal measure is also computed.
#'  
#'  
#' \strong{renyi, re:}
#'  Renyi entropy per object. This measure is a generalization of the Shannon entropy parameterized by \eqn{q}. It corresponds to the logarithm of the hill numbers. The default value for \eqn{q} is 0. \deqn{(1-q)^{-1} \log(\sum_i p_i^q)}
#' 
#' 
#' \strong{evenness, ev:}
#'  Pielou evenness per object across categories [Pielou, 1969]. It is based in Shannon Entropy \deqn{-\sum_i(p_i \log p_i)/\log{v} }
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
#' 
#' Rafols, I., & Meyer, M. (2009). Diversity and network coherence as indicators of interdisciplinarity: case studies in bionanoscience. Scientometrics, 82(2), 263-287.
#'
#' Rafols, I. (2014). Knowledge Integration and Diffusion: Measures and Mapping of Diversity and Coherence. In Y. Ding, R. Rousseau, & D. Wolfram (Eds.), Measuring Scholarly Impact (pp. 169-190). Springer International Publishing. 
#'  
#' Chavarro, D., Tang, P., & Rafols, I. (2014). Interdisciplinarity and research on local issues: evidence from a developing country. Research Evaluation, 23(3), 195-209.

#' @examples
#' data(pantheon)
#' diversity(pantheon)
#' diversity(pantheon, type='variety')
#' diversity(geese, type='berger-parker', category_row=TRUE)
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
#' #hill numbers
#' diversity(data=X, type="td", q=1)
#' #rao stirling with differente arguments
#' diversity(data=X, type="rao-stirling", method="euclidean", alpha=0, beta=1)
#' #more than one diversity measure
#' diversity(data=X, type=c('e','ev','bp','s'))
#' @export
diversity <- function(data, type="all", category_row=FALSE, dis=NULL, method='euclidean', q=0, alpha=1, beta=1, base=exp(1)){
  X <- get_data(data, category_row)
	diversity <- data.frame(row.names=rownames(X))
	sumsX <- rowSums(X, na.rm=TRUE)
	propX <- X / rowSums(X, na.rm=TRUE)
	
	for(measure in type)
	{
		
	  if (measure == 'variety' || measure =='v' || measure== 'all') {
	    m_d <- as.data.frame(rowSums(X>0, na.rm=TRUE))
	    colnames(m_d) <- c('variety')
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
	  if(measure == 'entropy' || measure=='e' || measure == 'all') { 
	    m_d <- as.data.frame(-1 * rowSums(propX * log(propX, base=base), na.rm=TRUE))
	    colnames(m_d) <- c('entropy')
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
		
		if(measure == 'all'  || measure == 'herfindahl-hirschman' || measure == 'hh' || measure == 'hhi'){
			m_d <- as.data.frame(rowSums(propX ^ 2, na.rm=TRUE))
			colnames(m_d) <- c('HHI')
			diversity <- merge(diversity,m_d, by=0, all=TRUE)
			rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
		}
		if(measure == 'blau' || measure=='b' || measure == 'all') {
			m_d <- as.data.frame(1 - rowSums(propX ^ 2, na.rm=TRUE))
			colnames(m_d) <- c('blau.index')
			diversity <- merge(diversity,m_d, by=0, all=TRUE)
			rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
		}
	  if(measure == 'gini-simpson' || measure=='gs' || measure == 'all') {
	    m_d <- as.data.frame(1 - rowSums(propX ^ 2, na.rm=TRUE))
	  	colnames(m_d) <- c('gini.simpson')
	    #compute associated measures as concentration and inverse
	    m_d['gini.simpson.C'] <- 1-m_d$gini.simpson #concentration
	    m_d['gini.simpson.R'] <- 1/m_d$gini.simpson.C #reciprocal
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
		if(measure == 'simpson' || measure == 's' || measure == 'all')
		{
	    X_simp <- X
	  	X_simp[X_simp==0] <- NA
	  	m_d <- as.data.frame(rowSums((X_simp*(X_simp-1))/matrix(sumsX*(sumsX-1), ncol=ncol(X_simp), nrow=nrow(X_simp)), na.rm=TRUE)) 
	  	colnames(m_d) <- c('simpson.D')
	  	m_d['simpson.I'] <- 1-m_d$simpson.D #I index
	  	m_d['simpson.R'] <- 1/m_d$simpson.D #R reciprocal
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
	  if(measure == 'hn' || measure=='td' || measure == 'all') {
	  	propX_td <- propX
	  	propX_td[propX_td==0] <- NA
	  	
	  	if(q == 0) #equals to variety
	  		m_d <- as.data.frame(rowSums(X>0, na.rm=TRUE))
	  	else if(q == 1) #an approximation is computed
	    	{
	  			m_d <- -1*as.data.frame(rowSums(propX * log(propX, base=exp(1)), na.rm=TRUE))
	  			m_d <- exp(m_d) #exponential of Shannon Entropy
	  		}
	  	else
	  	{
	  		p <- 1/(1-q)
	  		m_d <- as.data.frame((rowSums(propX ^ q, na.rm=TRUE)) ^ p)
	  	}
	    colnames(m_d) <- c('hill.numbers')
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
	  if(measure == 'berger-parker' || measure=='bp' || measure == 'all') {
	    m_d <- as.data.frame(apply(propX, 1, max))
	    colnames(m_d) <- c('berger.parker.D')
	  	m_d['berger.parker.I'] <- 1/m_d$berger.parker
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
	  if(measure == 'renyi' || measure=='re' || measure == 'all') {
	  	if(q == 0) #equals to log(variety)
	  		m_d <- log(as.data.frame(rowSums(X>0, na.rm=TRUE)), base=base)
	  	else if(q==1)  #tends to Shannon Entropy
	  	{
	  		m_d <- -1*as.data.frame(rowSums(propX * log(propX, base=base), na.rm=TRUE))
	  	}
	  	else
	  	{
	  		m_d <- as.data.frame(-1* log(rowSums(propX ^ q, na.rm=TRUE) , base=base)/(q-1))
	  	}
	    
	    colnames(m_d) <- c('renyi.entropy')
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
	  if(measure == 'evenness' || measure=='ev' || measure == 'all') {
	    m_d <- as.data.frame(-1 * rowSums(propX * log(propX, base=base), na.rm=TRUE)/log(rowSums(X>0, na.rm=TRUE), base=base)) # if N == 1 -> NaN
	    colnames(m_d) <- c('evenness')
	    diversity <- merge(diversity,m_d, by=0, all=TRUE)
	    rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  }
		if(measure == 'hcdt' || measure == 'all') {
			propX_td <- propX
			propX_td[propX_td==0] <- NA
			
			if(q == 0) #equals to variety
				m_d <- as.data.frame(rowSums(X>0, na.rm=TRUE))
			else if(q == 1) #Shannon Entropy with natural log
			{
				m_d <- -1*as.data.frame(rowSums(propX * log(propX, base=exp(1)), na.rm=TRUE))
			}
			else
			{
				m_d <- as.data.frame((1 - rowSums(propX ^ q, na.rm=TRUE))/(q-1))
			}
			colnames(m_d) <- c('hcdt.entropy')
			diversity <- merge(diversity,m_d, by=0, all=TRUE)
			rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
		}
	  if(measure == 'rao-stirling' || measure=='rs' || measure == 'all' || measure=='rao' || measure=='r' || measure=='disparity' || measure=='d')
	  {
	  	if(is.null(dis)) #computing distances in case that it is needed
	  	{
	  		disX <- dis_categories(X, category_row = category_row, method=method) #compute distances first		
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
	  	if(measure == 'rao-stirling' || measure=='rs' || measure=='all')
	  	{
	  		m_d[,'rao.stirling'] <- NA; ms_label <- 'rao.stirling'	
	  		alpha_r<- alpha; beta_r <- beta
	  	}
	  	if(measure=='rao' || measure=='r')
	  	{
	  		m_d[,'rao'] <- NA	 ;  ms_label <- 'rao'
	  		alpha_r <- 1; beta_r<-1
	  	}
	  	if(measure=='disparity' || measure=='d')
	  	{
	  		m_d[,'disparity.sum'] <- NA	
	  		m_d[,'disparity.mean'] <- NA	
	  		ms_label <- 'disparity.sum'
	  		alpha_r <- 1; beta_r<-1
	  	}
	  	for(entity in row.names(propX)) #go into each entity
	  	{
	  	  	entity_data <- propX[entity,]
	  	  	prop_i <- matrix(entity_data, nrow=N,ncol=N, byrow=TRUE)
	  	  	prop_j <- matrix(entity_data, nrow=N,ncol=N, byrow=FALSE) #transpose of previous one
	  	  	  	  	
	  	  	if(ms_label=='disparity.sum')
	  	  	{
	  	  		prop_i[prop_i>0] <- 1 #binarizing proportion
	  	  		prop_j[prop_j>0] <- 1 #binarizing proportion
	  	  		entity_variety <- length(entity_data[entity_data>0])
	  	  		#print(paste("entity", entity, 'ent data', entity_data, 'ent variety', entity_variety))
	  	  	}
	  	  	
	  	  	p_ij <- prop_i * prop_j
					p_ij_mask <- p_ij
	  			p_ij_mask[upper.tri(p_ij_mask)] <- 1
	  			
	  		  diag(p_ij_mask) <- 0
	  		  p_ij_mask[lower.tri(p_ij_mask)] <- 0
	  	  	rs_entity <- ((disX^alpha_r)*disX_mask) * ((p_ij^beta_r)*p_ij_mask) #masks ensures for proportions to use only upper triangle matrix. For distances, to use only existant distances and upper triangle matrix.
	  		  m_d[entity, ms_label] <- sum(rs_entity, na.rm = TRUE)
	  		  
	  		  if(ms_label=='disparity.sum')
	  		  { 
	  		  	m_d[entity, 'disparity.mean'] <- sum(rs_entity, na.rm=TRUE)/(entity_variety * (entity_variety-1)/2)
	  		  }
	  	}#end for entities 
	  	diversity <- merge(diversity,m_d, by=0, all=TRUE)
	  	rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
	  } #end if disparity measures
  }#end loop measures
  
  return(diversity)
}


#' @title Balance or proportions
#' @description Computes the proportions or probabilities of raw values. 
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' 
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @examples 
#' balance(data=geese, category_row = TRUE)
#' @return A matrix of entities-categories with proportions.
#' @export
balance <- function(data, category_row=FALSE)
{
	X <- get_data(data, category_row)
	propX <- X / rowSums(X, na.rm=TRUE)
	return(propX)
}


#' @title Variety or Richness
#' @description Computes the variety (number of distinct types) or simple diversity of an entity. It is also known as richness. 
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param sort Indicates whether results should be ordered or not. Define it to FALSE to avoid ordering.
#' @param decreasing If argument "sort" is set to TRUE, this argument indicates descending order. The default value is TRUE. 
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @examples 
#' variety(data=pantheon)
#' variety(data=pantheon, sort=FALSE)
#' @return A dataframe with values of variety for each entity.
#' @export
variety <- function(data, sort=TRUE, decreasing=TRUE, category_row=FALSE)
{
	vari <- diversity(data, type='v', category_row=category_row)
	if(sort != FALSE)
	{
		vari['category'] <- row.names(vari)
		vari <- vari[order(vari$variety, decreasing = decreasing), ]
		vari['category'] <- NULL	
	}
	
	return(vari)
}

#' @title Pre-process the raw data
#' @description Allows to filter, binarize and/or normalize raw data.
#' Also filter and binarization is available.
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value.  
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @param norm Methods to compute normalized values. Possible values are 'p', 'proportions', 'rca', 'rca_norm' and 'ai'. RCA refers to Revealed Comparative Advantages [Balassa 1986], rca_norm normalizes the RCAs between -1 and with 1, ai refers to the Activity Index. 
#' @param filter A threshold below which values are replaced with NA.
#' @param binary A boolean value to indicate if values distinct from NA are replaced with 1.
#' @details If the three arguments 'norm', 'filter' and 'binary' are used, then the same sequential order is applied in the calculations.
#' @references Balassa, B. (1986). Comparative advantage in manufactured goods: a reappraisal. The Review of Economics and Statistics, 315-319.
#' @examples 
#' #raw values
#' values(data=pantheon)
#' values(data = scidat)
#' #proportions
#' values(data = scidat, norm='p')
#' #revealed comparative advantages
#' values(data = scidat, norm='rca')
#' values(data = scidat, norm='rca', filter=1)
#' values(data = scidat, norm='rca', filter=1, binary=TRUE)
#' @return A matrix with the raw, normalized, filtered and\/or binarized data.
#' @export
values <- function(data, category_row = FALSE, norm=NULL, filter=NULL, binary=FALSE)
{
	X <- get_data(data, category_row)
	if(!is.null(norm))
	{
		if(norm == 'p' || norm =='s' || norm == 'share' || norm == 'proportions')
		{
			X <- X / rowSums(X, na.rm=TRUE)
		}
		else if(norm == 'rca' || norm=='ai' || norm == 'rca_norm')
		{
			propX <- X / rowSums(X, na.rm=TRUE)
			total <- sum(X, na.rm=TRUE)
			sum_categ <- colSums(X, na.rm=TRUE)
			propCateg <- sum_categ / total
			propCateg_mat <- matrix(propCateg, ncol=length(colnames(propX)), nrow=length(rownames(propX)), byrow = TRUE)
			X <- propX/propCateg_mat
			if(norm == 'ai' || norm == 'rca_norm')
			{
				#X <- ((2 * (X-max(X))) / (max(X)-min(X)))+1
				X <- (X-1) /(X+1)
			}
		} #end rca normalization
	}#end normalization
	#filter
	if(!is.null(filter))
	{
		X[X <= filter] <- NA
	}
	#binary
	if(binary==TRUE)
	{
		X[!is.na(X)] <- 1
	}

	return(X)
}

#' @title Ubiquity of categories across entities
#' @description Computes the ubiquity or the rareness of the categories
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value.  
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @examples 
#' ub <- ubiquity(data=pantheon)
#' @return A dataframe with values of number of entities where the category is present. Ordered in decreasing order.
#' @export
ubiquity <- function(data, category_row = FALSE)
{
	ubiq <- diversity(data, type='v', method='euclidean' , category_row= (!category_row))
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
#' @param type It indicates the type of data to be read. This argument facilitates the input of diverse types of data files, such as spss or stata. Possible options are the names of the mentioned software. The default value is csv.
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @return A data frame with three columns (entity, category, value).
#' @examples 
#' #reading an edges list or panel shape, source data must include three columns
#' path <-   system.file("extdata", "PantheonEdges.csv", package = "diverse")
#' sep <- ","
#' data <- read_data(path)
#' #reading a table
#' path  <- system.file("extdata", "PantheonMatrix.csv", package = "diverse")
#' sep <- ","
#' data <- read_data(path)
#' #reading a table which includes the entities in the columns
#' path <- system.file("extdata", "Geese.csv", package = "diverse")
#' data <- read_data(path, category_row=TRUE)
#' @export
#' @importFrom reshape2 melt
#' @importFrom foreign read.spss read.dta
#' @importFrom utils read.csv
read_data <- function(path, type='csv',sep=',', category_row=FALSE){

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
		if(category_row==TRUE)
		{
			data <- t(data)
		}
			
		data <- melt(data, na.rm=TRUE)
	}
	row.names(data) <- NULL
  return(data)
}

#' @title A procedure to create a disparity matrix between categories.
#' @description Takes a data frame or a matrix to create a disparity matrix
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param method A distance or dissimilarity method available in "proxy" package as for example "Euclidean", "Kullback" or "Canberra". This argument also accepts a similarity method available in the "proxy" package, as for example: "cosine", "correlation" or "Jaccard" among others. In the latter case, a correspondent transformation to a dissimilarity measure will be retrieved. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). The default value is Euclidean distance.
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @return A distance or dissimilarity square matrix
#' @examples 
#' Xdis <- dis_categories(pantheon)
#' @export
#' @importFrom proxy dist
dis_categories <- function(data, method='euclidean', category_row=FALSE){
    X <- get_data(data=data, category_row=category_row)
	  disX <- as.matrix(dist(t(X), method=method, convert_distances=TRUE), diag=1) 
  	return(disX)
}


#' @title A procedure to create a disparity matrix between entities
#' @description Takes a data frame or a matrix to create a disparity matrix
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param method A distance or dissimilarity method available in "proxy" package as for example "Euclidean", "Kullback" or "Canberra". This argument also accepts a similarity method available in the "proxy" package, as for example: "cosine", "correlation" or "Jaccard" among others. In the latter case, a correspondent transformation to a dissimilarity measure will be retrieved. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). The default value is Euclidean distance.
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @return A distance or dissimilarity square matrix
#' @examples 
#' Xdis <- dis_entities(pantheon)
#' #for science dataset
#' dis_entities(scidat, method='cosine')
#' @export
#' @importFrom proxy dist
dis_entities <- function(data, method='euclidean', category_row=FALSE){
	X <- get_data(data=data, category_row=!category_row)
	disX <- as.matrix(dist(t(X), method=method, convert_distances=TRUE), diag=1) 
	return(disX)
}

#' @title A procedure to compute the sum and average of disparities
#' @description Computes the sum and the average of distances or disparities between the categories.
#' @param data A numeric matrix with entities \eqn{i} in the rows and categories \eqn{j} in the columns. Cells show the respective value (value of abundance) of entity \eqn{i} in the category \eqn{j}. It can also be a transpose of the previous matrix, that is, a matrix with categories in the rows and entities in the columns. Yet in that case, the argument "category_row" has to be set to TRUE. The matrix must include names for the rows and the columns. The argument "data", also accepts a dataframe with three columns in the following order: entity, category and value. 
#' @param method A distance or dissimilarity method available in "proxy" package as for example "Euclidean", "Kullback" or "Canberra". This argument also accepts a similarity method available in the "proxy" package, as for example: "cosine", "correlation" or "Jaccard" among others. In the latter case, a correspondent transformation to a dissimilarity measure will be retrieved. A list of available methods can be queried by using the function \code{\link[proxy]{pr_DB}}. e.g. summary(pr_DB). The default value is Euclidean distance.
#' @param category_row A flag to indicate that categories are in the rows. The analysis assumes that the categories are in the columns of the matrix. If the categories are in the rows and the entities in the columns, then the argument "category_row" has to be set to TRUE. The default value is FALSE.
#' @return A data frame with disparity measures for each entity in the dataset. Both the sum of disparities and the average of disparities are computed.
#' @examples 
#' data(pantheon)
#' disparity(data= pantheon)
#' disparity(data = pantheon, method='Canberra')
#' #For scientific publications
#' #Same disparities, since all countries authored all entities
#' disparity(scidat)
#' disparity(data= scidat, method='cosine')
#' #Creating differences by measuring Revealed Compartive Advantages 
#' disparity(values(scidat, norm='rca', filter=1))
#' #Activity Index for scientometrics
#' disparity(values(scidat, norm='ai', filter=0), method='cosine')
#' #Using binarization of values and a binary metric for dissimilarities.
#' disparity(values(scidat, norm='ai', filter=0, binary=TRUE), method='jaccard')
#' @export
disparity <- function(data, method='euclidean', category_row=FALSE) {
  disparity <- diversity(data=data, category_row=category_row, method=method, type='disparity')
  return(disparity)
}

#' @title A procedure to simulate labeled individuals for one category
#' @description Simulates a number of individuals tagged in N different categories, given a distribution such as log normal or normal.
#' @param n_categ number of categories 
#' @param size number of individuals. 
#' @param category_prefix a prefix to be used as part of the category label
#' @param type distribution name. The distribution is used to simulate how individuals are created. Use 'log-normal' for log normal distribution or 'normal' for normal distribution. Default value is 'log-normal'
#' @param mean parameter for normal or log-normal distribution. Default value is 0.
#' @param sd parameter for normal or log-normal distribution. Default value is 1.
#' @importFrom stats rnorm rlnorm
#' @return A vector of category labels. 
#' @examples 
#' sim_individuals(n_categ=50, size=10000, type='log-normal', mean=0.507, sd=1.183)
#' @export
sim_individuals <- function(n_categ, size,  category_prefix='', type = 'log-normal', mean=0, sd=1) {
		probabilities <- NA
		individuals <- ''
		if(type=='log-normal')
		{
			probabilities <- rlnorm(n_categ, meanlog = mean, sdlog = sd)
			#probabilities <- rlnorm(50, meanlog=0.507, sdlog = 1.183)
			probabilities <- probabilities/sum(probabilities)
		}
		if(type == 'normal')
		{
			probabilities <- rnorm(n_categ, mean = mean, sd = sd)
			probabilities <- probabilities/sum(probabilities)
		}
		if(is.na(probabilities[1])==FALSE)
		{
			
			while(length(unique(individuals)) != n_categ)
			{
				individuals <- sample(paste(category_prefix,(1:n_categ), sep=''), size = size, replace = TRUE, prob = probabilities )
			}
			
		}
		else
		{
			while(length(unique(individuals)) != n_categ)
			{
				individuals <- sample(paste(category_prefix,(1:n_categ), sep=''), size = size, replace = TRUE)
			}
		}
		return(individuals)
}


#' @title A procedure to simulate entities 
#' @description Simulates an entity with values of abundance for some categories. 
#' @param n_categ number of categories 
#' @param size number of individuals. Default value is 7 times n_categ.
#' @param category_prefix a prefix to be used as part of the category label
#' @param values values of abundance. This argument can be both, a distribution name or a vector of integers. The distribution is used to simulate individuals that are aggregated in frequencies or values of abundance. In the second case, an integer or a vector of integers of possible values of abundance to be used randomly. Default value is 'log-normal'
#' @param mean parameter for normal or log-normal distribution. Default value is 0.
#' @param sd parameter for normal or log-normal distribution. Default value is 1.
#' @return A data frame with two columns: category and value of abundance. 
#' @examples 
#' sim_entity(n_categ=50,  category_prefix='ctg', values=1) #equal value
#' #random numbers for values of abundance
#' sim_entity(n_categ=50,  category_prefix='ctg', values=sample(1:100, replace=TRUE)) 
#' sim_entity(n_categ=50,  category_prefix='ctg', values='log-normal') #equal value
#' @export
sim_entity <- function(n_categ, category_prefix='', values = 'log-normal', size=-1,    mean=0, sd=1) {
	data_entity = data.frame()
	if(is.numeric(values)==TRUE)
	{
		Value <- values
		if(length(values)>1)
			Value <- sample(x = values, size = n_categ, replace=TRUE)
		Category <- paste(category_prefix,(1:n_categ),sep='')
		data_entity <- data.frame(Category, Value)
	}
	else
	{
		if(size==-1)
				size_s<-n_categ*7
		else
				size_s <- size
		data_entity <- as.data.frame(table(sim_individuals(n_categ = n_categ, category_prefix = category_prefix, size = size_s, type = values, mean = mean, sd=sd)))
	}
	
	return(data_entity)
}


#' @title A procedure to simulate datasets 
#' @description Simulates a dataset with values of variety for each entity and possible values of abundance. 
#' @param n_categ a vector with number of categories for each entity. The number of entities to create is defined by the length of this vector.
#' @param size number of individuals. A number or a vector of numbers for each entity. Default value is 7 times variety. 
#' @param category_prefix a prefix to be used as part of the category label
#' @param entity_prefix a prefix to be used as part of the entity label
#' @param values values of abundance. This argument can be both, a distribution name or a vector of integers. The distribution is used to simulate individuals that are aggregated in frequencies or values of abundance. Use 'log-normal' for log normal distribution or 'normal' for normal distribution.  In the second case, an integer or a vector of integers of possible values of abundance to be used randomly. Default value is 'log-normal'
#' @param mean parameter for normal or log-normal distribution. Default value is 0.
#' @param sd parameter for normal or log-normal distribution. Default value is 1.
#' @param category_random boolean argument to determine if categories should be taken randomly (TRUE) or sequentially (FALSE). Default is FALSE
#' @return A data frame with three columns: entity, category and value of abundance. 
#' @examples 
#' sim_dataset(n_categ=50,  category_prefix='ctg', values=1) #equal value, just one entity
#' #Several entities with random values
#' n_entities <- 50
#' v_n_c <- sample(1:100, size = n_entities, replace=TRUE)
#' v_v <- sample(10:5000, size= n_entities, replace=TRUE)
#' d <- sim_dataset(n_categ = v_n_c, values= v_v, category_random = TRUE)
#' @export
sim_dataset <- function(n_categ, category_prefix='', entity_prefix='', values = 'log-normal', size=-1,mean=0, sd=1, category_random=FALSE) {
	data_set = data.frame()
	
	for(i in seq(1,length(n_categ)))
	{
		values_e <- values
		
		size_e <- size
		if(length(size)>1)
		{
			size_e <- size[i]
		}
		data_entity = sim_entity(n_categ = n_categ[i], category_prefix = category_prefix, values = values_e, size = size_e, mean = mean, sd=sd)
		Entity <- paste(entity_prefix,i,sep = "")
		data_entity<- cbind(Entity, data_entity)
		data_set <- rbind(data_set,data_entity)
	} 
	
	colnames(data_set) <- c("Entity", "Category", "Value")
	if(category_random == TRUE){
		#print(data_set$Category)
		categories <- as.vector(unique(data_set$Category))
		#print(categories)
		categ_vector <- vector()#new vector of categories
		
		for(n in n_categ)
		{
		  categ_vector <- c(categ_vector, sample(categories, n))	
		}
		#print("Random categories")
		data_set$Category <- categ_vector
	}
	
	return(data_set)
}


#' @title Transforms data to be used in Entropart package 
#' @description Transform a dataframe used in diverse to values of abundance to be used in Entropart. 
#' @param data a dataframe used in diverse with three columns, entity, category and value of abundance.
#' @return An object of type matrix of abundance to be used in entropart to create metacommunities. 
#' @examples 
#' ab <- to_entropart(sim_dataset(c(1,2))) 
#' @export
to_entropart <- function(data)
{
	Abundances <- data.frame(t(values(data)))
	Species <- rownames(Abundances)
	Communities <- colnames(Abundances)
	return(Abundances)
}