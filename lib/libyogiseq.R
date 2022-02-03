
#Function to turn string into character array
to.char.array <- function (str) sapply(1:nchar(str), function(i) substr(str,i,i))
char.at <- function(str,i) substr(str,i,i)


new.sequence <- function(sequence, qual=NULL, id=NULL) {

	.seq <- sequence
	.qual <- qual
	.id <- id

	toString <- function() {
		.seq
	}
	getQuality <- function(is) {
		if (!is.null(.qual)) {
			.qual[is]
		} else {
			# warning("This sequence has no quality track.")
			NULL
		}
	}
	getID <- function() {
		if (!is.null(.id)) {
			.id
		} else {
			# warning("This sequence has no ID field.")
			NULL
		}
	}

	structure(list(
		toString=toString,
		getQuality=getQuality,
		getID=getID
	),class="yogiseq")
}
print.yogiseq <- function(s) print(paste("<YogiSeq:",s$getID(),">"))
summary.yogiseq <- function(s) c(id=s$getID(),sequence=s$toString(),phred=paste(s$getQuality(),collapse=","))
length.yogiseq <- function(s) nchar(s$toString())

reverseComplement <- function(seq) {		
	trans <- c(A='T',C='G',G='C',T='A',N='N',R='Y',Y='R',S='S',W='W',K='M',M='K')
	if (any(class(seq) == "yogiseq")) {
		revSeq <- paste(rev(sapply(to.char.array(seq$toString()), function(nc) trans[nc])),collapse="")
		revQual <- rev(seq$getQuality())
		new.sequence(revSeq,qual=revQual,id=seq$getID())
	} else {
		paste(rev(sapply(to.char.array(seq), function(nc) trans[nc])),collapse="")
	}
}

subseq <- function(s,from,to) {
	if (!any(class(s) == "yogiseq")) stop("First argument must be a YogiSeq object")
	new.sequence(
		substr(s$toString(),from,to), 
		if (!is.null(s$getQuality())) s$getQuality(from:to) else NULL,
		s$getID()
	)
}


writeFASTA <- function(con,seqs) {
	for (i in 1:length(seqs)) {
		s <- seqs[[i]]
		if (class(s) == "yogiseq") {
			writeLines(c(
				paste(">",s$getID(),sep=""),
				s$toString()
			),con)
		} else if (class(s) == "character") {
			writeLines(c(
				paste(">",names(seqs)[i],sep=""),
				s
			),con)
		} else {
			warning("Skipping unsupported data type",class(s))
		}
	}
}

readFASTA <- function(con) {
	out <- list()
	id <- NULL
	seq <- NULL
	i <- 0
	while(length(line <- readLines(con, n=1)) > 0) {
		if (substr(line,1,1)==">") {
			#if old sequence exists, add it to the output
			if (!is.null(id)) {
				out[[length(out)+1]] <- new.sequence(seq,id=id)
				cat(paste("\r Read",i <- i+1,"sequences    "))
			}
			#new sequence
			id <- substr(line,2,nchar(line))
			seq <- ""
		} else {
			seq <- paste(seq,line,sep="")
		}
	}
	#add last sequence to output
	if (!is.null(id)) {
		out[[length(out)+1]] <- new.sequence(seq,id=id)
		cat(paste("\r Read",i <- i+1,"sequences    \n"))
	}
	out
}


writeFASTQ <- function(con, seqs) {

	#function for decoding phred quality scores
	qualScale <- to.char.array("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
	qual2string <- function(qual) paste(qualScale[qual-32],collapse="")

	writeLines(unlist(lapply(seqs, function(s) {
		c(
			paste("@",s$getID(),sep=""),
			s$toString(),
			"+",
			qual2string(s$getQuality())
		)
	})),con)

}


#creates a new fastq parser object
new.fastq.parser <- function(con) {

	.con <- con

	#function for decoding phred quality scores
	qualScale <- to.char.array("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
	string2phred <- function(string) {
		out <- sapply(to.char.array(string), function(x) which(qualScale == x))
		names(out) <- NULL
		out+32
	}

	#function for parsing the next n entries from the open fastq file (or less if less than n remain)
	parse.next <- function(n=10,ignore.quality=FALSE) {

		contents <- list()
		i <- 0

		while ((i <- i+1) <= n && length(lines <- readLines(.con, n=4)) > 0) {

			if (length(lines) < 4 || substr(lines[1],1,1) != "@" || substr(lines[3],1,1) != "+") {
				stop("Corrupt read:\n",paste(lines,collapse="\n"))
			}

			id <- strsplit(substr(lines[1],2,nchar(lines[1])), " ", fixed=TRUE)[[1]][1]
			sequence <- lines[2]

			quality <- if (ignore.quality) NULL else string2phred(lines[4])

			contents[[length(contents)+1]] <- new.sequence(sequence,id=id,qual=quality)

		}

		contents
	}

	structure(list(parse.next=parse.next),class="yogi.fastq.parser")
}



#alignment algorithm requires bitwise operations
library(bitops)

##
# Needleman-Wunsch global alignment algorithm
#
new.alignment <- function(s1, s2) {

	if (any(class(s1)=="yogiseq")) {
		c1 <- c("$",to.char.array(s1$toString()))
	} else {
		c1 <- c("$",to.char.array(s1))
	}
	if (any(class(s2)=="yogiseq")) {
		c2 <- c("$",to.char.array(s2$toString()))
	} else {
		c2 <- c("$",to.char.array(s2))
	}

	#init score matrix
	mat <- matrix(nrow=length(c1), ncol=length(c2))
	mat[1,] <- 1:length(c2) - (c1[1] == c2[1])
	mat[,1] <- 1:length(c1) - (c1[1] == c2[1])

	#init trace matrix
	trace <- matrix(0, nrow=length(c1), ncol=length(c2))
	trace[1,] <- 4
	trace[,1] <- 2
	trace[1,1] <- 0

	#compute alignment matrix
	for (i in 2:length(c1)) {
		for (j in 2:length(c2)) {
			options <- c(
				rep = mat[i-1,j-1] + (c1[i] != c2[j]),
				del = mat[i-1,j] + 1,
				ins = mat[i,j-1] + 1
			)
			mat[i,j] <- min(options)

			tr.bitmasks <- 2^(which(options == min(options))-1)
			for (mask in tr.bitmasks) {
				trace[i,j] <- bitOr(trace[i,j],mask)
			}
		}
	}

	getMatrix <- function() {
		mat
	}

	getDistance <- function() {
		mat[length(c1),length(c2)]
	}

	.mutations <- NULL
	.mapping <- NULL
	run.trace <- function() {

		rep <- 1
		del <- 2
		ins <- 4

		muts <- list()
		map <- list()

		i <- length(c1)
		j <- length(c2)

		while (i > 1 || j > 1) {
			if (bitAnd(trace[i,j], rep) > 0) {
				if (c1[i] != c2[j]) {
					muts[[length(muts)+1]] <- c(c1[i], i-1, j-1, c2[j])
				}
				map[[length(map)+1]] <- c(i-1, j-1)
				i <- i-1
				j <- j-1
			} else if (bitAnd(trace[i,j], del)) {
				muts[[length(muts)+1]] <- c(c1[i], i-1, j-1, "-")
				map[[length(map)+1]] <- c(i-1, NA)
				i <- i-1
			} else if (bitAnd(trace[i,j], ins)) {
				muts[[length(muts)+1]] <- c("-", i-1, j-1, c2[j])
				map[[length(map)+1]] <- c(NA, j-1)
				j <- j-1
			} else {
				stop("uninitialized trace at ",i,j)
			}
		}
		# if (c1[1] != c2[1]) {
		# 	muts[[length(muts)+1]] <- c(c1[1],i,c2[1])
		# }
		.mapping <<- do.call(rbind,rev(map))
		.mutations <<- do.call(rbind,rev(muts))
	}

	getMutations <- function() {
		if (is.null(.mutations)) run.trace()
		.mutations
	}

	getMappings <- function() {
		if (is.null(.mapping)) run.trace()
		.mapping
	}

	printAlignment <- function() {

		if (is.null(.mapping)) run.trace()

		chars <- do.call(cbind,lapply(1:nrow(.mapping), function(k) {
			i <- .mapping[k,1]
			j <- .mapping[k,2]
			char1 <- if (is.na(c1[i+1])) '-' else c1[i+1]
			char2 <- if (is.na(c2[j+1])) '-' else c2[j+1]
			matchChar <- if (is.na(c1[i+1]) || is.na(c2[j+1])) " " 
				else if (c1[i+1] == c2[j+1]) "|" else "."
			c(char1,matchChar,char2)
		}))

		cat("\nLevenstein distance:",getDistance(),"\n")

		for (wrap in 0:(ncol(chars)/70)) {
			startcol <- wrap*70 + 1
			endcol <- if (startcol+69 > ncol(chars)) ncol(chars) else startcol+69
			cat("\n",paste(apply(chars[,startcol:endcol],1,paste,collapse=""),collapse="\n"),"\n",sep="")

		}

	}

	structure(list(
		getMatrix=getMatrix,
		getDistance=getDistance,
		getMutations=getMutations,
		getMappings=getMappings,
		printAlignment=printAlignment
	),class="yogialign")

}


##
# Creates a new translator object for translating Nucleotide strings to Amino acid strings.
#
init.translator <- function(ctable.file="codontable.txt") {

	##
	# Creates a new codon table object
	#
	init.codon.table <- function(con) {

		nc2single <- list()
		nc2triple <- list()

		while (length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
			cols <- strsplit(line,"\t")[[1]]
			aa3 <- cols[1]
			aa1 <- cols[2]
			codons <- strsplit(cols[3], "|", fixed=TRUE)[[1]]
			for (codon in codons) {
				nc2single[codon] <- aa1
				nc2triple[codon] <- aa3
			}
		}
		
		# Return single-letter code of aminoacid that is encoded by the given codon
		getSingleForCodon <- function(codon) {
			nc2single[[codon]]
		}

		structure(list(
			getSingleForCodon=getSingleForCodon
		),class="codonTable")
	}

	tryCatch({

		con <- file(ctable.file, open="r")
		codons <- init.codon.table(con)

	},
	error = function(ex) {
		cat(ex)
	},
	finally = {
		if (isOpen(con)) {
			close(con)
		}
	})

	# translates a given nucleotide sequence
	translate <- function(nucl) {

		aaQual <- NULL
		if (any(class(nucl)=="yogiseq")) {
			ncSeq <- nucl$toString()
			if (!is.null(nucl$getQuality())) {
				aaQual <- sapply(seq(1,length(nucl),3), function(i) min(nucl$getQuality(i:(i+2))))
			}
		} else {
			ncSeq <- nucl
		}

		if (nchar(ncSeq) == 0) stop("translate: empty string! ",ncSeq)

		aa <- paste(sapply(
			seq(1,nchar(ncSeq),3),
			function(i) {
				a <- codons$getSingleForCodon(substr(ncSeq,i,i+2))
				if(is.null(a)) "" else a
			}
		), collapse="")

		aaseq <- to.char.array(aa)

		if (any(aaseq == "*")) {
			cutoff <- min(which(aaseq == "*"))
			aa <- paste(aaseq[1:cutoff], collapse="")
		}

		if (!is.null(aaQual)) {
			attr(aa,"quality") <- aaQual
		}
		aa
	}

	list(translate=translate)	
}



#input: list of lists of mutatiion descriptor strings 
# (e.g 'A20G', 'truncation', 'nonsense' or 'silent')
#output: matrix of mutations
mutlist2matrix <- function(mutations, num.aa) {

	#initialize the change matrix
	change.matrix <- matrix(0,nrow=21,ncol=num.aa,
		dimnames=list(
			c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'),
			1:num.aa
		)
	)
	
	#remove any mutations that come with trucations and nonsense
	muts <- unlist(mutations[sapply(mutations, {
		function(x) length(x) > 0 && !any(x == "truncation" | x == "nonsense")
	})])

	# Mark mutations in the matrix
	for (sac in muts) {

		if (sac == "silent") next

		pos <- as.numeric(substr(sac,2,nchar(sac)-1))
		aa <- substr(sac,nchar(sac),nchar(sac))

		change.matrix[aa,pos] <- change.matrix[aa,pos] + 1
	}

	change.matrix
}

# Plots the mutation coverage for a given change matrix
plotMutCoverage <- function(change.matrix, sequence, translator, all=FALSE, main="") {

	.sequence <- ifelse (any(class(sequence) == "yogiseq"), sequence$toString(), sequence)
	#translate sequence to protein
	protein <- trans$translate(sequence)

	#a function that returns the i'th codon from the template
	codon.at <- function(dna, i) substr(dna,3*i-2, 3*i)

	#init reachability matrix
	reach.matrix <- matrix(NA,nrow=21,ncol=ncol(change.matrix),
		dimnames=list(
			c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'),
			1:nchar(protein)
		)
	)
	# fill reachability matrix according to sequence
	for (i in 1:ncol(change.matrix)) {
		codon <- codon.at(.sequence,i)
		for (pos in 1:3) {
			for (nc in c('A','C','G','T')) {
				mut.codon <- codon
				substr(mut.codon,pos,pos) <- nc
				aa <- translator$translate(mut.codon)
				reach.matrix[aa,i] <- 0
			}
		}
		reach.matrix[char.at(protein,i),i] <- -1
	}


	# compute coverage
	if (!all) {
		coverage <- (apply(change.matrix,2,function(x) sum(na.omit(x) > 0)) 
			/ (apply(reach.matrix,2,function(x) sum(!is.na(x)))-1))
	} else {
		coverage <- apply(change.matrix,2,function(x) sum(na.omit(x) > 0)) / 20
	}

	#define drawing layot, set drawing color to gray, adjust margins
	layout(matrix(c(1,2),ncol=1), heights=c(1,3))
	op <- par(fg="gray",mar=c(0,4.1,4.1,2.1))	

	# draw a bar plot for coverage
	ylim <- if (max(coverage) > 1) c(0,max(coverage)) else c(0,1)
	barplot(coverage,
		main=main, 
		xlab="Position",
		ylab="Coverage", 
		ylim=ylim,
		border=NA,
		names.arg=NA,
		col="darkolivegreen3"
	)

	# Compute a color gradient to represent the mutation counts
	maxVal <- max(apply(change.matrix,1,function(x) max(na.omit(x))))
	colors <- colorRampPalette(c("white", "orange"))(5)

	### Draw the diagram
	# use horizontal axis labels
	op <- c(op,par(las=1))
	par(mar=c(5.1,4.1,0,2.1))
	# create an empty plot
	plot(0,
		type='n',
		axes=FALSE,
		xlim=c(0,ncol(change.matrix)), 
		ylim=c(0,21),
		xlab="Position",
		ylab="Amino acid"
	)
	# iterate over each matrix entry and draw the contents on the plot
	for (x in 1:ncol(change.matrix)) {
		for (y in 1:21) {
			if (change.matrix[y,x] > 0) {
				#observed mutations are drawn in a color shade corresponding to their count
				col <- colors[ceiling(5*change.matrix[y,x]/maxVal)+1]
				rect(x-1,22-y,x,21-y,col=col, lty="blank")
			}
		}
	}
	for (x in 1:ncol(change.matrix)) {
		for (y in 1:21) {
			if (!is.na(reach.matrix[y,x])) {
				if (reach.matrix[y,x] == -1) {
					#original amino acids are marked in gray
					rect(x-1,22-y,x,21-y,col="gray")
				} else if (!all) {
					#reachable aminoacids are marked with dotted outline
					rect(x-1,22-y,x,21-y, lty="dashed",lwd=2)
				}
			}
		}
	}
	# draw axes
	axis(1, at=c(1,seq(5,ncol(change.matrix),5))-.5, labels=c(1,seq(5,ncol(change.matrix),5)))
	axis(2, at=(1:21)-.5, labels=rev(rownames(change.matrix)) )

	par(op)
}


##
# Mutagenic PCR simulation function
# sequence = original DNA sequence, should start with a start codon and end with a stop codon
# cycles = number of PCR cycles to simulate. Should be > 1, but too large numbers will affect runtime and memory usage exponentially
# init.amount = Initial amount of template molecules to use in the simulation
# etr = Enzyme-to-Template ratio. Defaults to 1/1
# mutation.rate = Mutations per bp introduced by the enzyme per replication process.
#
pcr.sim <- function(sequence, translator, cycles=10, init.amount=100, etr=1, mut.rate=1/2000) {

	.sequence <- ifelse(any(class(sequence) == "yogiseq"), sequence$toString(), sequence) 

	enzyme.amount <- round(init.amount * etr)

	#calculate sampling bias based on Mutazyme II bias and sequence bias
	pol.bias <- c(A=0.2675, C=0.2325, G=0.2325, T=0.2675)
	seq.bias <- table(to.char.array(.sequence)) / nchar(.sequence)
	bias <- pol.bias * seq.bias / sum(pol.bias * seq.bias)
	cbias <- c(bias[1],sapply(2:4, function(i) sum(bias[1:i])))
	names(cbias) <- names(bias)

	#make index of nucleotide positions
	nuc.positions <- sapply(c('A','C','G','T'), function(nuc) which(to.char.array(.sequence) == nuc))

	#mutation transition matrix based on Mutazyme II
	mut <- cbind(
		rbind(
			A=c(A=0,   C=.047,G=.175,T=.285),
			C=c(A=.141,C=0,   G=.041,T=.255),
			G=c(A=.255,C=.041,G=0,   T=.141),
			T=c(A=.285,C=.175,G=.047,T=0   )
		) * .5,
		DEL=rep(.048,4)/4,
		INS=rep(.008,4)/4
	) * 4
	mut <- mut / apply(mut,1,sum)
	cmut <- cbind(mut[,1],sapply(2:ncol(mut), function(i) apply(mut[,1:i],1,sum)))
	dimnames(cmut) <- dimnames(mut)

	#seed molecule pool with templates
	pool <- list()
	for (i in 1:init.amount) pool[[i]] <- list()

	#perform PCR cycles
	for (c in 1:cycles) {

		num.reactions <- min(length(pool),enzyme.amount)
		templates <- sample(pool, num.reactions)
		num.muts <- rpois(num.reactions, nchar(.sequence) * mut.rate)

		new.mutations <- sapply(num.muts, function(num.mut) {

			if (num.mut == 0) {
				return(list())
			}

			# use bias table to figure out how many of each nucleotide to pick for mutating
			to.sample <- table(sapply(1:num.mut, function(i) {
				names(which.min(which(runif(1,0,1) < cbias)))
			}))

			#pick positions to mutate
			to.mutate <- sapply(names(to.sample), function(nuc) {
				sample(nuc.positions[[nuc]], to.sample[nuc])
			})

			#implement mutations
			unlist(sapply(names(to.mutate), function(nuc) {
				sapply(to.mutate[[nuc]], function(pos) {
					#sample mutation
					to.nuc <- names(which.min(which(runif(1,0,1) < cmut[nuc,])))

					if (to.nuc == "DEL" || to.nuc == "INS") {
						return("nonsense")
					} else {

						codon.number <- floor((pos-1) / 3) + 1
						codon.start <- 3*codon.number - 2
						from.codon <- substr(.sequence,codon.start,codon.start+2)
						
						change.pos <- pos - codon.start + 1
						to.codon <- from.codon
						substr(to.codon,change.pos,change.pos) <- to.nuc

						from.aa <- translator$translate(from.codon)
						to.aa <- translator$translate(to.codon)

						if (from.aa == to.aa) {
							return("silent")
						} else if (to.aa == "*") {
							return("truncation")
						} else {
							return(paste(from.aa,codon.number,to.aa,sep=""))
						}
					}
				})
			}))

		})
		names(new.mutations) <- NULL
	
		#add mutagenized copies to pool
		pool[(length(pool)+1):(length(pool)+num.reactions)] <- sapply(1:num.reactions, function(i) {
			c(templates[[i]], new.mutations[[i]])
		})

	}

	#return pool without original templates
	pool[-(1:init.amount)]
}


default.error <- function(ex) {
	print(ex)
	traceback(ex)
	stop()
}
protect <- function(filename, fun, mode="r",error=default.error) {
	tryCatch({
		con <- file(filename, open=mode)
		fun(con)
	},
	error = error,
	finally = {
		if (exists("con") && isOpen(con)) {
			close(con)
		}
	})
}


# processFile <- function(file,f) {
	# tryCatch({
	# 	con <- file(file, open="r")
	# 	f(con)
	# },
	# error = function(ex) {
	# 	traceback(ex)
	# },
	# finally = {
	# 	if (exists("con") && isOpen(con)) {
	# 		close(con)
	# 	}
	# })
# }

# test1 <- NULL
# processFile("test1.fastq",function(con) {
# 	test1 <<- parseFASTQ(con)
# })

# test2 <- NULL
# processFile("test2.fastq",function(con) {
# 	test2 <<- parseFASTQ(con)
# })







