####
## THIS LIBRARY IS A COLLECTION OF VARIOUS USEFUL TOOLS I WROTE
#


###
# This object can count occurrences of different items
#
new.counter <- function() {

	a <- list()

	###
	# Add x occurrences to item id
	#
	add <- function(id,x) {
		if (!is.character(id)) {
			stop("Illegal argument:",id)
		}
		if (is.null(a[[id]])) {
			a[[id]] <<- x
		} else {
			a[[id]] <<- a[[id]] + x
		}
	}

	# increase counter for item id by 1
	inc <- function(id) add(id,1)

	# get the counter state for id
	get <- function(id) a[[id]]

	# list counts for all ids
	ls <- function() a

	# export counter state as a string
	export <- function() {
		paste(lapply(names(a), function(id) paste(id,"=",a[[id]],sep="") ), collapse=",")
	}

	# import counter state from string
	import.add <- function(strs) { 
		lapply(strsplit(strs,","), function(eqs) {
			lapply(strsplit(eqs,"="), function(vals) {
				add(vals[[1]],as.numeric(vals[[2]]))
			})
		})
		invisible()
	}

	structure(
		list(
			inc = inc,
			add = add,
			get = get,
			ls = ls,
			export = export,
			import.add = import.add
		),
		class="yogicounter"
	)
}


#Function for returning the i'th ranked item in a list
ith.rank <- function(values, i) sort(values,decreasing=TRUE)[i]

###
# Matthew's correlation coefficient (MCC)
# 
mcc <- function(t, scores, truth) {

	# exclude <- is.na(scores) | is.na(truth)
	# scores <- scores[-exclude]
	# truth <- truth[-exclude]

	.truth <- truth == 1
	.calls <- scores >= t

	tp <- sum(.truth & .calls)
	tn <- sum(!.truth & !.calls)
	fp <- sum(.calls & !.truth)
	fn <- sum(.truth & !.calls)

	# mcc <- (tp * tn - fp * fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
	#this formula prevents integer overflow errors
	mcc <- exp( log(tp*tn - fp*fn) - ( log(tp+fp) + log(tp+fn) + log(tn+fp) + log(tn+fn) )/2 )
	prec <- tp/(tp+fp)
	recall <- tp/(tp+fn)
	c(mcc=mcc,prec=prec,recall=recall)
}



###
# assemble list of matrices into 3d matrix
#
matrix.3d <- function(matrices) {
	nr <- nrow(matrices[[1]])
	nc <- ncol(matrices[[1]])
	nm <- length(matrices)
	mat <- rep(0,nr*nc*nm)
	dim(mat) <- c(nr,nc,nm)
	rownames(mat) <- rownames(matrices[[1]])
	colnames(mat) <- colnames(matrices[[1]])
	for(i in 1:nm) {
		mat[,,i] <- matrices[[i]]
	}
	mat
}



###
# This function forms all possible subsets of a given set
#
combo <- function(l) {
	do.call(c,lapply(1:length(l),function(n){
		tab <- combn(l,n)
		lapply(1:ncol(tab),function(i)tab[,i])
	}))
}

