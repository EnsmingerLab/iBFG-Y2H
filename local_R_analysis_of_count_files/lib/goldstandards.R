library("hash",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)
source("lib/synergizer.R")

###
# This class contains a number of methods for dealing with Gold standard data
#
new.goldstandards <- function() {

	#Get HI3
	load.HI3 <- function(hi3.file) {
		hi3 <- read.delim(hi3.file,stringsAsFactors=FALSE)
		hi.genes <- unique(c(hi3[,1],hi3[,2]))

		#Translate HI3 to gene names
		syn <- new.synergizer()
		transl <- syn$translate(
			authority="ensembl",
			species="Homo sapiens",
			domain="entrezgene",
			range="hgnc_symbol",
			ids=hi.genes
		)
		transl.h <- hash(keys=names(transl), values=transl)
		do.call(rbind,apply(hi3,1,function(ia) {
			list(transl.h[[ as.character(ia[1]) ]], transl.h[[ as.character(ia[2]) ]] )
		}))
	}

	load.HI.II.14 <- function(hi.file) {
		hi <- read.delim(hi.file, stringsAsFactors=FALSE)
		hi[,c(2,4)]
	}


	#Get PRS 
	load.PRS <- function(prs.file) {
		read.delim(prs.file,stringsAsFactors=FALSE,header=FALSE)
	}

	#Adapt Gold standard to test matrix
	#linear: gold standard has linear test space (didn't test all-by-all)
	adapt <- function(pairs, genewise.matrix, linear=FALSE,sep=NULL) {

		#create hash structure to quickly find all partners
		hi3.hash <- hash()
		for (i in 1:nrow(pairs)) {
			a <- unlist(pairs[i,1])
			b <- unlist(pairs[i,2])
			if ( !is.null(a) && !is.null(b) &&
			     length(a) > 0 && length(b) > 0 ) {
				for (.a in a) {
					if (!has.key(.a, hi3.hash)) {
						hi3.hash[[.a]] <- b
					} else {
						hi3.hash[[.a]] <- union(hi3.hash[[.a]], b)
					}
				}
				for (.b in b) {
					if (!has.key(.b, hi3.hash)) {
						hi3.hash[[.b]] <- a
					} else {
						hi3.hash[[.b]] <- union(hi3.hash[[.b]], a)
					}
				}
			}
		}

		#index barcode ids -> gene
		symbols <- if (!is.null(sep)) sapply(strsplit(barcodes$symbol,sep),function(x)x[2]) else barcodes$symbol
		bc2gene <- hash(keys=barcodes$id, values=symbols)

		#create gold standard for test space from hi3
		gold.matrix <- genewise.matrix
		for (i in 1:nrow(gold.matrix)) {
			id.i <- rownames(gold.matrix)[i]
			gene.i <- bc2gene[[id.i]]
			for (j in 1:ncol(gold.matrix)) {
				id.j <- colnames(gold.matrix)[j]
				gene.j <- bc2gene[[id.j]]
				if (has.key(gene.i, hi3.hash) && has.key(gene.j,hi3.hash)) {
					if (gene.j %in% hi3.hash[[gene.i]]) {
						gold.matrix[i,j] <- 1
					} else {
						gold.matrix[i,j] <- ifelse(linear,NA,0)
					}
				} else {
					gold.matrix[i,j] <- NA
				}
			}
		}
		gold.matrix
	}

	structure(
		list(
			adapt=adapt,
			load.HI3=load.HI3,
			load.PRS=load.PRS,
			load.HI.II.14=load.HI.II.14
		),
		class="GoldStandards"
	)

}