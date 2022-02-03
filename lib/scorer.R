###
# This class calculates BFGY2H s and s' scores
#
new.scorer <- function(barcodes.file="res/barcodes.csv") {

	# barcodes.file <- "res/barcodes.csv"
	barcodes <- read.csv(barcodes.file,header=FALSE,stringsAsFactors=FALSE)
	colnames(barcodes) <- c("ad.db","type","id","symbol","version","bc.name","upseq","dnseq")

	if (any(duplicated(barcodes$bc.name))) {
		stop("Error in barcode descriptor: barcode names (bc.name) must be unique!")
	}

	.matrix.template <- matrix(
		0,
		nrow=sum(barcodes$ad.db=="AD"),
		ncol=sum(barcodes$ad.db=="DB"),
		dimnames=list(
			rownames=with(barcodes,bc.name[ad.db=="AD"]),
			colnames=with(barcodes,bc.name[ad.db=="DB"])
		)
	)

	get.matrix.template <- function() {
		.matrix.template
	}


	#####
	# normalize by +HIS
	#
	calc.s <- function(cs, cp, alpha=1) {
		.cs <- cs + alpha
		.cp <- cp + alpha

		.fs <- .cs / sum(.cs)

		.f.ad <- apply(.cp,1,sum) / sum(.cp)
		.f.db <- apply(.cp,2,sum) / sum(.cp)

		.f.est <- get.matrix.template()
		for (i in 1:nrow(.cp)) {
			for (j in 1:ncol(.cp)) {
				.f.est[i,j] <- .f.ad[i] * .f.db[j]
			}
		}
		return( .fs / .f.est )
	}

	###
	# calculate s'
	#
	calc.s.prime <- function(.s, rho=0.75) {
		apply(.s, 2, function(col) {
			med <- median(col)
			beta <- quantile(col[col > med] - med, rho)
			(col-med)/beta
		})
	}


	structure(list(
		calc.s = calc.s,
		calc.s.prime = calc.s.prime,
		get.matrix.template = get.matrix.template
		# collapse = collapse
	), class="scorer")

}
