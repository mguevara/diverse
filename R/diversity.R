#' @title  Diversity measures
#' @description It takes an object x category matrix and calculates a number of diversity measures
#' @param data A numeric matrix or data frame with objects as rows and categories as columns
#' @param type A mnemonic string referencing the diversity measure. List of available measures: "variety", "entropy", "gini", "simpson", "true", "inverse-simpson", "herfindahl–hirschman", "renyi", "evenness", "rao-stirling". A list of short mnemonics for each measure: 'v', 'e', 'g', 's', 't', 'inv', 'hh', 're', ev', and 'rs'. The default for type is "all". More information for each measure in details and examples. 
#' @param method "rao-stirling" uses a disparity function between objects. List of available disparity methods: "cosine", "jaccard", "euclidean". The default for method is cosine.
#' @param agg_type aggregation type for diversity analysis. The analysis is conducted per row but it can also be performed by column setting agg_type = "col". Default is NULL. 
#' @param q parameter for true diversity index measure. This parameter is also used for the Rényi entropy. Default is 0.
#' @param alpha parameter for Rao-Stirling diversity measure. As default we consider alpha=1.
#' @param beta parameter for Rao-Stirling diversity measure. As default we consider beta=1.
#' @details Available diversity measures are (written for an object x category matrix): 
#' 
#' Variables: N (category count), p_i (proportion of system comprises category i), d_ij (disparity between i and j).
#' 
#' variety: N, category counts per object [MacArthur 1965]
#' 
#' entropy: - sum(p_i log p_i), Shannon entropy per object [Shannon 1948]
#' 
#' gini: 1 - sum(p_i^2), Gini-Simpson index per object [Gini 1912].  It is also known as the Gibbs–Martin index or the Blau index in sociology, psychology and management studies.
#' 
#' simpson: sum(p_i^2), Simpson index per object [Simpson 1949]. This measure is known as the Herfindahl–Hirschman index in economy.  
#' 
#' true: (sum(p_i^q))(1-q)^-1, true diversity index per object [Hill 1973]. This measure is q parameterized. Default for q is 0.  
#' 
#' berger-parker: it is equals to the maximum p_i value in the dataset, i.e. the proportional abundance of the most abundant type. 
#' 
#' inverse-simpson: (sum(p_i^2))^-1, inverse simpson index per object. This measure is the true diversity at q = 2.
#' 
#' renyi: log((sum(p_i^q))(1-q)^-1), Rényi entropy per object. It is a generalization of the Shannon entropy parameterized by q. It corresponds to the logarithm of the true diversity. Default for q is 0. 
#' 
#' evenness: (-sum(p_i log p_i))/log N, Shannon evenness per object across categories [Pielou 1969] 
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
#' X <- readEdges(path="~/MyDiversity/data/toy.edges", sepr=' ', we=TRUE)
#' diversity(X, type="gini")
#' diversity(X, type="rao-stirling", method="cosine")
#' diversity(X, type="all", method="jaccard")
#' 
#' data <- readCSV(path="~/MyDiversity/data/sitc_cnt_62.csv", sepr=' ', we=TRUE)
#' diversity(data, type="gini")
#' diversity(data, type="rao-stirling", method="cosine")
#' diversity(data, type="all", method="jaccard")
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
  if(type == 'true' || type=='t' || type == 'all') {
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
  if(type == 'inverse-simpson' || type=='inv' || type == 'all') {
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
  if(!is.null(method)) {
    if (method == 'jaccard') {
      disX <- as.matrix(dist(t(X), method="Jaccard"), diag=1)
      
    }
    else if (method == 'euclidean') {
      disX <- as.matrix(dist(t(X), method="euclidean"), diag=1)
    }
    else if (method == "cosine") {
      disX <- as.matrix(dist(t(X), method="cosine"), diag=1)
    }
  }
  else {
    disX <- as.matrix(dist(t(X), method="euclidean"), diag=1)
  }
  if(type == 'rao-stirling' || type=='rs' || type == 'all') {
      m_d <- as.data.frame(rowSums(propX^beta %*% disX^alpha * propX^beta))
      colnames(m_d) <- c('rao-stirling')
      diversity <- merge(diversity,m_d, by=0, all=TRUE)
      rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  }
  #if(type == 'disparity' || type=='dis' || type == 'all') {
    #  m_d <- as.data.frame(rowSums(disX))
    #  colnames(m_d) <- c('disparity')
    #  diversity <- merge(diversity,m_d, by=0, all=TRUE)
    #  rownames(diversity) <- diversity$Row.names; diversity$Row.names <- NULL
  #}
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

#' @title Ubiquity
#' @description It computes the ubiquity or the rearnes of the categories
#' @param data Data to be processed as dataframe or as matrix. 
#' @param data_agg Diferent of NULL if column analysis is needed.
#' @examples 
#' ub <- ubiquity(data=d)
#' @return a dataframe with values of frequency per category. Decreasing order
#' @export
ubiquity <- function(data)
{
	ubiq <- diversity(data, type='v', method='euclidean' , agg_type='col')
	colnames(ubiq) <- 'ubiquity'
	ubiq['category'] <- row.names(ubiq)
	ubiq <- ubiq[order(ubiq$ubiquity, decreasing = TRUE), ]
	ubiq['category'] <- NULL
	return(ubiq)
}

#' @title A procedure to read a matrix from a file
#' @description It takes a file and creates a matrix for diversity analysis
#' @param path A string representing the path to a file. Each row corresponds to a matrix row (Ex.: 1 2 0 1 1 0 represents a row a matrix with 6 cols)
#' @param cols Number of columns of the matrix
#' @param sepr Separator field used in the file to separate columns
#' @return A matrix with objects as rows and categories as rows
#' @examples 
#' path <- "~/MyDiversity/data/toy.mat"
#' cols <- 6
#' sepr <- ' '
#' X <- readMatrix(path,cols,sepr)
#' @export
readMatrix <- function(path,cols,sepr){
  X <- matrix(scan(file=path, sep=sepr), ncol=cols, byrow=TRUE)
  return(X)
}

#' @title A procedure to read a list of edges from a file
#' @description It takes a file and creates a matrix for diversity analysis
#' @param path A string representing the path to a file. Rows and columns are described by integers (Ex.: 1 2 1)
#' @param sepr Separator field used in the file to separate columns
#' @param we It indicates if the list of edges includes weights or not. Default is TRUE
#' @return A matrix with objects as rows and categories as rows
#' @examples 
#' path <- "~/MyDiversity/data/toy.edges"
#' sepr <- ' '
#' we <- TRUE
#' X <- readEdges(path,sepr,we)
readEdges <- function(path,sepr,we=TRUE){
  if (!we) col = 2 else col = 3
  L <- matrix(scan(file=path, sep=sepr), ncol=col, byrow=TRUE)
  rows <- max(L[,1])
  cols <- max(L[,2])
  X <- array(0, c(rows,cols))
  if (col == 3) 
    X[L[,1:2]] <- L[,3]
  else
    X[L[,1:2]] <- 1
  return(X)
}

#' @title A procedure to read a data frame from a csv file
#' @description It takes a file and creates a data frame for diversity analysis
#' @param path A string representing the path to a csv file. Rows and columns are described by strings (Ex.: "asind" "0010" 216000).
#' @param sepr Separator field used in the file to separate columns
#' @param we It indicates if the list of edges includes weights or not. Default is TRUE
#' @return A data frame with objects as rows and categories as cols
#' @examples 
#' path <- "~/MyDiversity/data/sitc_cnt_62.csv"
#' sepr <- ' '
#' we <- TRUE
#' data <- readCSV(path,sepr,we)
readCSV <- function(path,sepr,we=TRUE){
  data <- data.frame()
  #if (!we) col = 2 else col = 3
  col_classes <- c('factor', 'factor', NA)
  data <- rbind(data, read.csv(path, sep=sepr, colClasses=col_classes))
  return(data)
}

#' @title A procedure to write a diversity data frame to a csv file
#' @description It takes a data frame and creates a csv file for diversity analysis
#' @param frame A data frame with diversity measures
#' @param path A string representing the path to a csv file
#' @param sepr Separator field used in the file to separate columns
#' @return A csv file with objects as rows and diversity measures as cols
#' @examples 
#' writeCSV(data,"~/MyDiversity/data/data.csv",' ')
writeCSV <- function(frame, path, sepr) {
  write.table(frame, path, sep=sepr, row.names=FALSE) 
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
dist_mat <- function(data, method=NULL, agg_type=NULL){
  if (is.data.frame(data)) {
    if(is.null(agg_type)) {
      X <- matrix(0, nrow=nlevels(data[,2]), ncol=nlevels(data[,1]), dimnames=list(levels(data[,2]),levels(data[,1])))
      X[cbind(data[,2], data[,1])] <- data[,3]
    }
    else {
      X <- matrix(0, nrow=nlevels(data[,1]), ncol=nlevels(data[,2]), dimnames=list(levels(data[,1]),levels(data[,2])))
      X[cbind(data[,1], data[,2])] <- data[,3]
    }
  }
  else {
    if (!is.null(agg_type)) {
      X <- t(data)
    }
    else {
      X <- data
    }  
  }
  if(!is.null(method)) {
    if (method == 'jaccard') {
      disX <- as.matrix(dist(t(X), method="Jaccard"), diag=1)
      
    }
    else if (method == 'euclidean') {
      disX <- as.matrix(dist(t(X), method="euclidean"), diag=1)
    }
    else if (method == "cosine") {
      disX <- as.matrix(dist(t(X), method="cosine"), diag=1)
    }
  }
  else {
    disX <- as.matrix(dist(t(X), method="cosine"), diag=1)
  }
  return(disX)
}

#' @title A procedure to analyze a disparity matrix
#' @description It takes a disparity matrix to get average or sums of disparities
#' @param data A disparity matrix to be analyzed
#' @param method List of available disparity matrix analysis methods: "avg", "sum", "all". The default for method is "all".
#' @param agg_type aggregation type for disparity analysis. The analysis is conducted per row but it can also be performed by column setting agg_type = "col". Default is "row". 
#' @return A data frame with disparity measures as columns for each object of data
#' @examples 
#' Xdis <- dist_mat(data, method="euclidean")
#' disp <- disparity(Xdis)
#' disp <- disparity(Xdis, method="avg", agg_type='col')
disparity <- function(X, method="all", agg_type="row") {
  if (agg_type == "col") {
    X <- t(X)
    disparity <- data.frame(row.names=rownames(X))
  }
  else {
    disparity <- data.frame(row.names=rownames(X))
  }
  if (method == 'sum' || method == 'v' || method == 'all') {
    m_d <- as.data.frame(rowSums(X, na.rm=TRUE))
    colnames(m_d) <- c('sum')
    disparity <- merge(disparity, m_d, by=0, all=TRUE)
    rownames(disparity) <- disparity$Row.names; disparity$Row.names <- NULL
  }
  if (method == 'avg' || method =='a' || method == 'all') {
    m_d <- as.data.frame(rowMeans(X, na.rm=TRUE))
    colnames(m_d) <- c('avg')
    disparity <- merge(disparity, m_d, by=0, all=TRUE)
    rownames(disparity) <- disparity$Row.names; disparity$Row.names <- NULL
  }
  return(disparity)
}