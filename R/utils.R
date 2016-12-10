# Get Data
# It takes data as dataframe (edges) or as matrix (table) to be exported in proper form to be used by the diversity function.
# data Data to be processed as dataframe or as matrix. 
# category_row TRUE if column analysis is needed.
get_data <- function(data, category_row=FALSE)
{
	if (is.data.frame(data)) {
		entities <- as.vector(unique(data[,1]))
		categories <- as.vector(unique(data[,2]))
		n_ent <- length(entities)
		n_categ <- length(categories)

		data['id_ctg'] <- match(data[,2], categories) #getting numeric index for
		data['id_ent'] <- match(data[,1], entities) #getting numeric index for
		
		if(category_row==TRUE) {
			X <- matrix(0, nrow=n_categ, ncol=n_ent, dimnames=list(categories,entities))
			X[cbind(data[,4], data[,5])] <- data[,3]
		}
		else {
			X <- matrix(0, nrow=n_ent, ncol=n_categ, dimnames=list(entities,categories))
			X[cbind(data[,5], data[,4])] <- data[,3]
		}
	}
	else {
		if (category_row==TRUE) {
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