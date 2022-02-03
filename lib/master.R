#!/software/R/bin/Rscript

library("hash",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)

source("lib/sge2.R")         #SUN Grid Engine
source("lib/libyogiseq.R")   #FastQ parser
source("lib/scorer.R")       #Scoring methods
source("lib/cliargs.R")      #Command-line argument processing
source("lib/liblogging.R")   #Writing log files
source("lib/libyogitools.R") #Some handy tools

###
# Get command line arguments
#

# Files containing the R1 and R2 reads. Example file:
# /home/nyachie/proj/cresample2/Data/proj-test_Jochen/miseq_data/PhiX_S1_L001_R1_001.fastq
r1.files <- getArg("r1",required=TRUE)
r2.files <- getArg("r2",required=TRUE)
if (length(strsplit(r1.files,",")[[1]]) != length(strsplit(r2.files,","))) {
	stop("There needs to be an equal number of r1 and r2 files!")
}
rfile.table <- cbind(r1=strsplit(r1.files,",")[[1]],r2=strsplit(r2.files,",")[[1]])


# The session tag is used to name the ouptut directory together with a timestamp.
session.tag <- getArg("session",default="bfgy2h")

# The chuck size determines how many reads are processed by each slave script.
chunk.size <- as.numeric(getArg("chunksize",default=20000))

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

# This option determines how many jobs will be submitted to the cluster at maximum.
# When more jobs are available, they get buffered internally until enough room exists in the queue.
max.queue <- as.numeric(getArg("maxQueue",default=60))

# This file contains information about the barcodes. Example:
# /home/nyachie/proj/cresample2/Data/db/CNT_Ver1/csv/MainDB.CNT_Ver1_08162012.csv
barcodes.descriptor.file <- getArg("barcodeDescriptor",default="res/barcodes.csv")
# We only really need it to passi it on to the constructor of the scorer
scorer <- new.scorer(barcodes.descriptor.file)

#Barcode sequences. Taken from:
#/home/nyachie/proj/cresample2/Data/db/CNT_Ver1/fasta/Barcode_CNT_Ver1_08162012.fna
barcode.seq   <- getArg("barcodeSeq",default="res/barcodes.fasta")

#Primer sequences. Taken from:
#/home/nyachie/proj/cresample2/Data/db/CNT_Ver1/fasta/const-seq.fna
vector.seq    <- getArg("vectorSeq",default="res/vector.fasta")

# Muxtag sequence map. Taken from:
#/home/nyachie/proj/cresample2/Data/bar2num.txt
mux.id.file <- getArg("muxIDs",default="res/bar2num.txt")

# This file identifies the different multiplex tags. Example: 
# /home/nyachie/proj/cresample2/Data/db/CNT_Ver1/csv/Multiplexing_test-Jochen.csv
muxtags.file <- getArg("multiplexDescriptor",default="res/muxtags.csv")


###
# Read the multiplex file and assign proper headers. There are two different versions.
#
muxtags <- read.csv(muxtags.file,header=FALSE,stringsAsFactors=FALSE)
if (ncol(muxtags) == 5) {
	colnames(muxtags) <- c("tag","category","selection","version","comment")
} else if (ncol(muxtags) == 7) {
	colnames(muxtags) <- c("tag","category","selection","version","comment","date","project")
} else {
	stop("Unsupported format for muxtags.csv")
}

###
# Create output directory
#
timestamp <- format(Sys.time(),format='%Y-%m-%d_%H-%M-%S')
out.dir <- paste(session.tag,"_",timestamp,"/", sep="")
dir.create(out.dir, mode="0755")

# Set up log file
logger <- new.logger(paste(out.dir,"master.log",sep=""))



###
# PHASE 1: CREATE CHUNKS OF READS ON-THE-FLY AND START SLAVE JOBS ON THEM
#

###
# This function uses an already open FASTQ parser and an index 
# to create a new FASTA chunk file. Returns a FASTA file and a 
# flag saying whether no more chunks can be produced.
#
# parser = FastQ parser instance
# i = chunk number
# direction = R1 or R2
#
make.chunk <- function(parser, i, direction) {

	outfile <- paste(out.dir,direction,"-",i,".fasta",sep="")
	outcon <- file(outfile,open="w")

	curr.chunk.size <- 0
	while (length(seqs <- parser$parse.next(ignore.quality=TRUE)) > 0 
			&& curr.chunk.size < chunk.size) {
		writeFASTA(outcon,seqs)
		curr.chunk.size <- curr.chunk.size+length(seqs)
	}
	close(outcon)

	list(file=outfile, last=(curr.chunk.size < chunk.size))
}

###
# Function to find out whether file is GZip file
#
is.gz <- function(f) {
	substr(f,nchar(f)-2,nchar(f))==".gz" || 
		regexpr("gzip compressed data", system(paste("file",f),intern=TRUE) ) > 0
}

###
# work through read file pairs
#
counts.files <- list()
# Construct SunGridEngine object with designated maximum queue size.
sge <- new.sge(max.queue.length=max.queue, logger=logger, debug=debug.mode)

invisible(apply(rfile.table, 1, function(rfiles) {

	r1.file <- rfiles[["r1"]]
	r2.file <- rfiles[["r2"]]

	###
	# Open connections to sequencing results files
	# And open FastQ parsers on them.
	#
	con.r1 <- file(r1.file,open="r")
	if (is.gz(r1.file)) {
		con.r1 <- gzcon(con.r1)
	}
	con.r2 <- file(r2.file,open="r")
	if (is.gz(r2.file)) {
		con.r2 <- gzcon(con.r2)
	}
	parser.r1 <- new.fastq.parser(con.r1)
	parser.r2 <- new.fastq.parser(con.r2)


	done <- FALSE
	i <- 0
	while (!done) {
		i <- i+1
		# Make chunks for R1 and R2
		r1.chunk <- make.chunk(parser.r1,i,"R1")
		r2.chunk <- make.chunk(parser.r2,i,"R2")

		#Designate an output file for them
		counts.file <- paste(out.dir,"counts_",i,".tsv",sep="")
		counts.files[[length(counts.files)+1]] <<- counts.file

		#Create a job id
		job.id <- paste(session.tag,timestamp,i,sep="_")
		#Designate a log file.
		slave.log <- paste(out.dir,"slave_",i,".log",sep="")

		#Submit Slave job to SunGridEngine
		sge$enqueue(
			id=job.id,
			command="/software/R/bin/Rscript",
			arguments=list(
				"lib/slave.R",
				paste("r1=",r1.chunk$file,sep=""),
				paste("r2=",r2.chunk$file,sep=""),
				paste("out=",counts.file,sep=""),
				paste("log=",slave.log,sep=""),
				paste("barcodeSeq=",barcode.seq,sep=""),
				paste("vectorSeq=",vector.seq,sep=""),
				paste("muxIDs=",mux.id.file,sep="")
			)
		)

		#We're done if we run out of reads to process
		done <- r1.chunk$last || r2.chunk$last
	}
	#Close file connections
	close(con.r1)
	close(con.r2)

}))

#Wait for the remaining jobs to finish
sge$wait(verbose=TRUE)



####
# PHASE 2: CONSOLIDATE JOB RESULTS
#

# Construct a metadata table of all the different count matrices
matrices <- cbind(
	rbind(muxtags,muxtags),
	updn=c(rep("UpUp",nrow(muxtags)),rep("DnDn",nrow(muxtags))),
	stringsAsFactors=FALSE
)
# Construct a list of matrices (AD barcodes X DB barcodes). Each list entry 
# corresponds to a row in the "matrices" metadata table.
counts <- lapply(1:nrow(matrices),function(dummy) scorer$get.matrix.template())

# Create hashes to quickly find the matrix indices corresponding to the barcodes.
# This is just to speed up the 
ad.idx <- hash(keys=rownames(counts[[1]]),values=1:nrow(counts[[1]]))
db.idx <- hash(keys=colnames(counts[[1]]),values=1:ncol(counts[[1]]))

# A shortcut function to check whether a value exists and is valid
is.valid <- function(x) {
	!is.null(x) && !is.na(x) && (length(x) == 1)
}

# Construct a counter instance to count instances of alignment problems
qual.report <- new.counter()

####
# join count tables together
#
logger$info("Consolidating results")
k <- 0
pb <- txtProgressBar(max=length(counts.files),style=3)
#for each data file
for (counts.file in counts.files) {

	# logger$info("Evaluating",counts.file,"\n")

	#read in the data table from the file
	sub.counts <- read.delim(counts.file,stringsAsFactors=FALSE)

	#iterate over rows
	for (i in 1:nrow(sub.counts)) {

		#if the multiplexing tag is valid
		if (sub.counts$mux[[i]] %in% muxtags$tag) {

			#identify the index in the tag table corresponding to the mux tag
			# mux.idx <- which(muxtags$tag == sub.counts[i,"mux"])
			#get the barcode ids
			ad.bc <- sub.counts$bc1[[i]]
			db.bc <- sub.counts$bc2[[i]]

			#identify the matrix corresponding to the mux tag and fusion type
			mat.idx <- with(matrices, which(
				tag == sub.counts$mux[[i]] &
				updn == sub.counts$fusion[[i]]
			))
			
			#if matrix was found, add the counts to it
			if (is.valid(mat.idx) && is.valid(ad.bc) && is.valid(db.bc) 
				&& has.key(ad.bc,ad.idx) && has.key(db.bc,db.idx) ) {
				ad.i <- ad.idx[[ad.bc]]
				db.j <- db.idx[[db.bc]]
				counts[[mat.idx]][ad.i,db.j] <- counts[[mat.idx]][ad.i,db.j] + sub.counts$count[[i]]
			} else {
				#log invalid barcode combination
				qual.report$inc("incompatible.barcodes")
			}
		} else {
			#log wrong muxtag
			qual.report$inc("incompatible.muxtags")
		}
	}
	setTxtProgressBar(pb,k<-k+1)
}

logger$info(paste("Quality report:",qual.report$export()))

###
# write count matrices to file
# TODO: This is unsafe. Information should not be conveyed via file names.
#       Should export list of matrices and descriptor table instead.
#
logger$info("Writing results to file")
sapply(1:nrow(matrices), function(i) {
	if (ncol(matrices) > 6) {
		tag <- paste(matrices[i,c("category","selection","version","updn","project")],collapse="_")
	} else {
		tag <- paste(matrices[i,c("category","selection","version","updn")],collapse="_")
	}
	outfile <- paste(out.dir,tag,"_counts.tsv",sep="")
	write.table(counts[[i]], outfile, sep="\t", quote=FALSE)
	outfile
})


###
# Delete input files when we're done, to save HD space.
#
logger$info("Cleaning up count fragments")
invisible(lapply(counts.files, function(cfile) {
	file.remove(cfile)
}))


logger$info("Done!")
