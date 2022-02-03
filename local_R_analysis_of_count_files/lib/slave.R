#!/software/R/bin/Rscript

##### 
# This script as a worker instance that uses BLAST to identify barcodes.
#

library("Biostrings",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)
library("hash",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)

source("lib/liblogging.R")   #Logger
source("lib/libyogitools.R") #Various tools
source("lib/cliargs.R")      #Command-line argument processing


#get the command arguments
# R1 reads chunk file. Example:
#/home/nyachie/proj/cresample2/Data/proj-test_Jochen/fragmented_fasta/PhiX_S1_L001_R1_001_1037.fna
r1.file <- getArg("r1",required=TRUE)
# R2 reads chunk file.
r2.file <- getArg("r2",required=TRUE)
# output file.
out.file <- getArg("out",required=TRUE)
# log file
log.file <- getArg("log",required=TRUE)

#Start a logger instance for the given log file.
logger <- new.logger(log.file)

#Create a counter to track alignment problems
qual.report <- new.counter()

#Barcode sequences. Taken from:
#/home/nyachie/proj/cresample2/Data/db/CNT_Ver1/fasta/Barcode_CNT_Ver1_08162012.fna
barcode.fa   <- getArg("barcodeSeq",default="res/barcodes.fasta")

#Primer sequences. Taken from:
#/home/nyachie/proj/cresample2/Data/db/CNT_Ver1/fasta/const-seq.fna
primer.fa    <- getArg("vectorSeq",default="res/vector.fasta")

# Muxtag sequence map. Taken from:
#/home/nyachie/proj/cresample2/Data/bar2num.txt
bar2num.file <- getArg("muxIDs",default="res/bar2num.txt")

#Column names for blast output file
blast.out.fmt <- c(
	"read","template","bit","4","5","6",
	"str.on.read","end.on.read","str","end","evalue","12"
)

# Load barcode sequences
logger$info("Reading fasta files")
seqs <- append(
	readDNAStringSet(barcode.fa,format="fasta"),
	readDNAStringSet(primer.fa,format="fasta")
)
# Load Read information.
reads1 <- readDNAStringSet(r1.file,format="fasta")
reads2 <- readDNAStringSet(r2.file,format="fasta")

#Record the number of input reads
qual.report$add("reads",length(reads1)+length(reads2))

###
# Check for coherence of read files!
#
if (!all(names(reads1) == names(reads2))) {
	logger$fatal("R1 and R2 files do not match each other!")
	stop()
}

###
# This file contains the multiplexing tags (muxtags) and their ID numbers
#
bar2num <- read.csv(bar2num.file,header=FALSE,stringsAsFactors=FALSE)
#format muxtag IDs so they always have 2 digits
bar2num[,1] <- sprintf("%02d",bar2num[,1])
colnames(bar2num) <- c("id","seq")

# get.seqs <- function(set, ids) {
# 	lapply(ids, function(id) set[[which(names(set)==id)]])
# }

###
# This function calls a BLAST instance to map the barcodes,
# applies some filter criteria on the results and extracts 
# Multiplex tag sequences from the reads
#
process.reads <- function(read.file, direction) {

	#run BLAST
	logger$info(paste("Running BLAST on",direction))

	con <- pipe(paste(
		"/home/rothlab/jweile/bin/blastn",
		"-task blastn-short",
		"-strand plus",
		"-db",barcode.fa,
		"-outfmt 10",
		"-evalue 1e-8",
		"-query",read.file
	),open="r")
	#read BLAST output from stdout stream
	bl.out <- read.csv(con,header=FALSE,stringsAsFactors=FALSE)
	#bl.out <- read.csv("blast/out.BC_blast/PhiX_S1_L001_R1_001_1037.blast",header=FALSE,stringsAsFactors=FALSE)
	close(con)

	colnames(bl.out) <- blast.out.fmt

	#report loss of reads
	qual.report$add("unaligned",length(reads1)-length(unique(bl.out$read)))

	logger$info("Filtering")

	#retrieve lengths of original templates for each hit.
	.sl <- width(seqs)
	names(.sl) <- names(seqs)
	bl.out$tmpl.length <- sapply(bl.out$template, function(x) .sl[[x]])

	#apply Nozomu's filter criteria
	valid <- bl.out$str.on.read < bl.out$end.on.read
	valid <- valid & bl.out$bit >= 90
	valid <- valid & (bl.out$str <= 20)
	valid <- valid & (bl.out$end >= bl.out$tmpl.length-20)

	#report filter loss
	qual.report$add("invalid",sum(!valid))

	bl.out <- bl.out[valid,]

	#for each read, only keep the hit with the highest bitscore
	bl.out <- bl.out[order(bl.out$read, bl.out$bit, decreasing=TRUE),]
	bl.out <- bl.out[!duplicated(bl.out$read),]

	bl.out <- data.frame(
		read = bl.out$read,
		direction = direction,
		template = sapply(bl.out$template, function(tmpl) {
			if (regexpr("TAG$",tmpl) > 0) "BARCODE" else tmpl
		}),
		name = bl.out$template,
		str = bl.out$str.on.read - bl.out$str + 1,
		end = bl.out$end.on.read + bl.out$tmpl.length - bl.out$end,
		evalue = bl.out$evalue,
		bitscore = bl.out$bit,
		stringsAsFactors=FALSE
	)

	#report non-barcode hits
	qual.report$add("non-barcode",sum(bl.out$template != "BARCODE"))

	#drop non-barcode hits
	bl.out <- bl.out[bl.out$template == "BARCODE",]

	#extract multiplexing tag
	rseqs <- if (direction=="R1") reads1 else reads2
	rseq.index <- hash(keys=names(rseqs),values=1:length(rseqs))
	bl.out$mux.tag <- apply(bl.out,1, function(row) {
		rseq <- rseqs[[ rseq.index[[ row[["read"]] ]] ]]
		start <- as.numeric(row["str"])-as.numeric(row['str']) + 1
		end <- as.numeric(row["str"])-as.numeric(row['str']) + 9
		if (start > 0 && start <= 20 && end < nchar(rseq)) {
			as.character(subseq(rseq,start=start,end=end))
		} else {
			NA
		}
	})

	return(bl.out)
}

###
# This function takes muxtag sequences from the reads and identifies
# them using pairwise alignments
call.mux.tags <- function(sequences, bar2num) {

	logger$info("Analyzing multiplexing tag")
	
	#alignment score matrix
	score.matrix <- nucleotideSubstitutionMatrix(match = 0, mismatch = -1, baseOnly = FALSE)

	#to speed things up we only check unique sequences and then extrapolate the results
	uq.seqs <- na.omit(unique(sequences))
	hits <- lapply(uq.seqs,function(qry) {
		#Try our luck with string matching, for speed's sake
		hit <- which(bar2num$seq == qry)
		#Only if that fails we use alignments
		if (length(hit) == 0) {
			scores <- sapply(bar2num$seq, function(ref) {
				pairwiseAlignment(pattern=qry,subject=ref, 
					substitutionMatrix=score.matrix, gapOpening=-1, gapExtension=-1,
					scoreOnly=TRUE
				)
			})
			#Even the best alignment score has to be better than -5 to be accepted.
			max.score <- max(scores)
			if (max.score < -5 || sum(scores==max.score) > 1) {
				hit <- NA
			} else {
				hit <- which(scores == max.score)
			}
		}

		if (is.na(hit)) NA else bar2num$id[[hit]]
	})
	#Use a hashmap to translate the unique result back to all results.
	seq2hit <- hash(keys=uq.seqs,values=hits)
	sapply(sequences, function(sequence) if(is.na(sequence)) NA else seq2hit[[sequence]])
}

###
# This function counts the combinations of barcodes
#
count.combinations <- function(data) {

	logger$info("Counting combinations")

	#Read IDs can be duplicated in the result table
	reads <- unique(data$read)
	qual.report$add("pairs",length(reads1))
	qual.report$add("eliminated.pairs",length(reads1)-length(reads))

	# A function to standardize barcode IDs
	processName <- function(name) {
		#remove leading "c" if exists
		name <- sub("^c","",name)
		split <- strsplit(name,"-")[[1]]
		c(id=split[1],type=split[2])
	}

	# New hashmap to count the combinations.
	counts <- hash()

	for (read in reads) {

		#For each readpair id find index of R1 and R2 in the table
		idx1 <- which(data$read == read & data$direction == "R1")
		idx2 <- which(data$read == read & data$direction == "R2")

		#detect and discard missing partners
		if (length(idx1) == 0 || length(idx2) == 0) {
			qual.report$inc("incomplete.pairs")
			next
			# return(rep(NA,4))
		}

		#re-format barcode ids
		name1 <- processName(data[idx1,"name"])
		name2 <- processName(data[idx2,"name"])

		#check for correct fusion and identify fusion type.
		if (name1["type"] == name2["type"]) {
			fusion.type <- if (name1["type"] == "UPTAG") "UpUp" else "DnDn"
		} else {
			#report invalid fusion
			qual.report$inc("bad.fusion")
			next
		}

		#detect and discard if read pair is missing muxtags
		if (any(is.na(data[c(idx1,idx2),"mux.call"]))) {
			#log missing mux tags
			qual.report$inc("missing.muxtag")
			next
		}

		#build a tag that summarizes both multiplex tags
		mux.combo <- paste("P", data[idx1,"mux.call"], "-P", data[idx2,"mux.call"], sep="")

		#add the count to the hashmap
		key <- paste(
			mux.combo,
			name1[["id"]],
			name2[["id"]],
			fusion.type,
			sep="_"
		)
		counts[[key]] <- if (is.null(counts[[key]])) 1 else counts[[key]]+1

	}

	#turn the hashmap into a table and return it
	out <- data.frame(do.call(rbind,strsplit(keys(counts),"_")),values(counts))
	colnames(out) <- c("mux","bc1","bc2","fusion","count")
	return(out)
	
}

####
# Use the functions above to do the actual work:
#

#process both read files and merge the two resulting tables
data <- rbind(
	process.reads(r1.file,"R1"),
	process.reads(r2.file,"R2")
)
#identify multiplex tags and add them to the result table.
data$mux.call <- call.mux.tags(data$mux.tag,bar2num)

#do the counting
counts <- count.combinations(data)

#write the results to file
logger$info("Writing output file")
write.table(counts,out.file,sep="\t",quote=FALSE,row.names=FALSE)

#write the quality report to the log file.
logger$info(paste("Quality report:",qual.report$export()))

###
# Delete input files when we're done, to save HD space.
#
logger$info("Cleaning up sequence fragments")
file.remove(r1.file)
file.remove(r2.file)
