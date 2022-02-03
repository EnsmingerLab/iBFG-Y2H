#!/software/R/bin/Rscript

###
# Import libraries
#
library("parallel",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)
library("hash",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)
source("lib/cliargs.R")
source("lib/scorer.R")
source("lib/libyogiplot.R")
source("lib/cytoscape.R")
source("lib/liblogging.R")
source("lib/libyogitools.R")
source("lib/goldstandards.R")

###
# Get command line parameters
#

# HI3 interactome file to be used as a Gold Standard
hi3.file <- getArg("HI3","res/HI3.tsv")

# Positive reference set (PRS) of protein interactions (also Gold Standard)
prs.file <- getArg("PRS","res/prs.tsv")

# HI-II-14 Gold standard file
hi.ii.file <- getArg("HIII14","res/HI-II-14.tsv")

# The barcode descriptor file.
barcodes.file <- getArg("barcodeDescriptor","res/barcodes.csv")

# working directory
dir <- getArg("dir","./")
if (regexpr("/$",dir) < 0) dir <- paste(dir,"/",sep="")

# Log file
logger <- new.logger(paste(dir,"bfgy2h_analysis.log",sep=""))

# rho-parameter for s' score. Previously called beta_th
rho <- as.numeric(getArg("rho",0.75))
# alpha-parameter for s score. 
alpha <- as.numeric(getArg("alpha",1))
# id of permissive condition
permissive.id <- getArg("permissiveId","+His")

#separator for genenames in the barcodes descriptor file
genesep <- getArg("geneSep",NULL)

# select sub-project to process
subproject <- getArg("project",NULL)

#Number of CPU cores to use
n.cores <- as.numeric(getArg("cores",4))

#Gene space
gene.space.file <- getArg("space",NULL)
gene.space <- if (is.null(gene.space.file)) {
	NULL
} else {
	read.delim(gene.space.file,stringsAsFactors=FALSE)
}

###
# load count tables from file system
# TODO: This is just way too error-prone in case people use spaces or underscores
# in their project names. Should use rdata binaries instead.
#
logger$info("Loading data")
file.names <- paste(dir,list.files(dir,pattern="_counts.tsv"),sep="")
if (length(file.names)==0) {
	logger$fatal("No files found!")
	stop("No files found!")
}
matrix.descriptors <- as.data.frame(do.call(rbind,strsplit(gsub(".+/","",file.names),"_|\\.")),stringsAsFactors=FALSE)
if (ncol(matrix.descriptors) == 6) {
	colnames(matrix.descriptors) <- c("category","selection","version","updn","matrix.type","format")
} else if (ncol(matrix.descriptors) == 7) {
	colnames(matrix.descriptors) <- c("category","selection","version","updn","project","matrix.type","format")
} else {
	logger$fatal("Count files seem to be named incorrectly.")
	stop("Count files seem to be named incorrectly.")
}
matrices <- lapply(lapply(file.names, read.delim, stringsAsFactors=FALSE), as.matrix)

# The list of matrices containing the different counts
counts <- matrices[which(matrix.descriptors$matrix.type=="counts")]
# A table of metadata for the matrices. Each row in the table corresponds to one matrix
counts.descriptor <- matrix.descriptors[matrix.descriptors$matrix.type=="counts",]

if (ncol(matrix.descriptors) == 7) {
	if (!is.null(subproject)) {
		if (any(matrix.descriptors$project==subproject)) {
			counts <- matrices[which(matrix.descriptors$project==subproject)]
			counts.descriptor <- matrix.descriptors[which(matrix.descriptors$project==subproject),]
		} else {
			logger$fatal("Project ID does not match with any project data!")
			stop("Project ID does not match with any project data!")
		}
	} else if (length(unique(matrix.descriptors$project)) > 1) {
		logger$fatal("Multiple projects found, but no project id specified!")
		stop("Multiple projects found, but no project id specified!")
	}
}



###
# STEP 1: Draw UpTag-DnTag correlations as dot plots
#
outfile <- paste(dir,"count_corr.png",sep="")
if (!file.exists(outfile)) {
	logger$info("Drawing correlation plots")
	ups <- which(counts.descriptor$updn == "UpUp")
	dns <- sapply(ups,with(counts.descriptor, function(up) which(
		category==category[up] & 
		selection==selection[up] & 
		version==version[up] & 
		updn=="DnDn"
	)))
	grid <- best.grid(length(ups))
	png(outfile,350*grid[2],350*grid[1])
	op <- par(mfrow=grid)
	apply(cbind(up=ups,dn=dns),1,function(ud){
		x <- as.vector(counts[[ ud[["up"]] ]])
		y <- as.vector(counts[[ ud[["dn"]] ]])
		r <- format(cor(x,y),digits=2)
		label <- paste(counts.descriptor[ud[["up"]],c("category","selection","version")],collapse=" ")
		plot(
			log(x+1), 
			log(y+1),
			xlab="log(UpTag counts)",
			ylab="log(DnTag counts)",
			main=paste(label,"R =",r),
			pch=20,
			cex=.5,
			col="royalblue"
		)
	})
	par(op)
	invisible(dev.off())
}


###
# STEP 2: CALCULATE SCORES
#

#Load scorer tool
scorer <- new.scorer(barcodes.file)

###
# set up descriptor table for score matrices with index pointers to 
# the correct selective and permissive count matrices that serve as inputs.
#
selective <- which(counts.descriptor$selection != permissive.id)
permissive <- which(counts.descriptor$selection == permissive.id)
score.descriptor <- as.data.frame(do.call(rbind,lapply(selective, function(selective.idx) {
	permissive.idx <- with(counts.descriptor,which(
		selection == permissive.id
		& category == category[selective.idx]
		& version == version[selective.idx]
		& updn == updn[selective.idx]
	))
	c(
		selective=selective.idx, #pointer to -his/+3at matrix
		permissive=permissive.idx, #pointer to corresponding +his matrix
		counts.descriptor[selective.idx,c("category","selection","version","updn")]
	)
})),stringsAsFactors=FALSE)

###
#calculate scores
#
logger$info("Calculating scores")
s <- mclapply(1:nrow(score.descriptor), function(i) {
	scorer$calc.s(
		counts[[ score.descriptor[[i,"selective"]] ]], 
		counts[[ score.descriptor[[i,"permissive"]] ]],
		alpha = alpha
	)
},mc.cores=n.cores)
s.prime <- lapply(s, scorer$calc.s.prime, rho=rho)



###
# STEP 3: COLLAPSE SCORE MATRICES TO GENE-WISE REPRESENTATION
#

###
# Load barcode table
#
barcodes <- read.csv(barcodes.file,header=FALSE,stringsAsFactors=FALSE)
colnames(barcodes) <- c("ad.db","type","id","symbol","version","bc.name","upseq","dnseq")
#delete sequences. they just slow us down
barcodes$upseq <- NULL
barcodes$dnseq <- NULL

real.symbol <- if (is.null(genesep)) {
	barcodes$symbol
} else {
	sapply(strsplit(barcodes$symbol,genesep),function(x) if (length(x)==2) x[[2]] else "")
}

#list of AD genes from barcode table 
ad.genes <- if (is.null(gene.space)) {
	with(barcodes,unique(id[ad.db=="AD"]))
} else {
	with(barcodes,unique(id[ad.db=="AD" & (real.symbol %in% gene.space$symbol)]))
}
#list of DB genes from barcode table
db.genes <- if (is.null(gene.space)) {
	with(barcodes,unique(id[ad.db=="DB"]))
} else {
	with(barcodes,unique(id[ad.db=="DB" & (real.symbol %in% gene.space$symbol)]))
}

###
# Construct a descriptor table for gene-wise matrices 
# Listing all possible barcode combinations for each gene pair
# 
genewise.descriptor <- with(score.descriptor,expand.grid(
	category=unique(category),
	selection=unique(selection),
	version=unique(version),
	updn=unique(updn),
	ad.bc.v=unique(barcodes$version),#AD barcode version
	db.bc.v=unique(barcodes$version),#DB barcode version
	stringsAsFactors=FALSE
))

###
# Template for a gene-wise matrix
#
genewise.template <- matrix(
	0, nrow=length(ad.genes), ncol=length(db.genes),
	dimnames=list(
		rownames=ad.genes, colnames=db.genes
	)
)

###
# Re-format the score matrices into a gene-wise scoring matrices
# Using the above descriptor table
#
logger$info("Reformatting matrix for gene-wise scores")
hickups <- 0
genewise.scores <- mclapply(1:nrow(genewise.descriptor), function(i) {

	mat.in <- s.prime[[which(
		score.descriptor$category==genewise.descriptor$category[[i]] & 
		score.descriptor$selection==genewise.descriptor$selection[[i]] &
		score.descriptor$version==genewise.descriptor$version[[i]] &
		score.descriptor$updn==genewise.descriptor$updn[[i]]
	)]]

	ad.bc.v <- genewise.descriptor$ad.bc.v[[i]]
	db.bc.v <- genewise.descriptor$db.bc.v[[i]]

	ad.idx <- with(barcodes,which(ad.db=="AD" & version == ad.bc.v))
	ad.hash <- if (length(ad.idx) > 0) hash(keys=barcodes$id[ad.idx],values=barcodes$bc.name[ad.idx]) else hash()
	db.idx <- with(barcodes,which(ad.db=="DB" & version == db.bc.v))
	db.hash <- if (length(db.idx) > 0) hash(keys=barcodes$id[db.idx],values=barcodes$bc.name[db.idx]) else hash()
	

	mat.out <- genewise.template
	for (ad.gene in ad.genes) {
		ad.bc <- ad.hash[[ad.gene]]
		for (db.gene in db.genes) {
			db.bc <- db.hash[[db.gene]]
			if (length(ad.bc) == 1 && length(db.bc) == 1) {
				mat.out[ad.gene,db.gene] <- mat.in[ad.bc,db.bc]
			} else {
				mat.out[ad.gene,db.gene] <- NA
				hickups <<- hickups +1
			}
		}
	}
	mat.out
}
,mc.cores=n.cores
)

if (hickups > 0) {
	logger$warn(paste(hickups,"unassigned values!"))
}

#Turn list of matrices into 3D-matrix
scores.3d <- matrix.3d(genewise.scores)


####
# STEP 4: USE GOLD STANDARDS TO PICK BEST IPS-SCORES
#

###
# Compute a list of IPS matrices across parameter space
#
logger$info("Calculationg all possible data consolidations")
#all categories (-AA, ALL)
cats <- unique(unlist(score.descriptor$category))
#all selections (-His, +3AT)
sels <- unique(unlist(score.descriptor$selection))
#the parameter space to be searched:
parameter.space <- expand.grid(#all possible combinations of the below
	categories=combo(cats),#all possible subsets of categories
	selections=combo(sels),# all possible subsets of selections
	consolidator=append(lapply(1:8,function(r) list("rank",r)),"average"), #ranks 1-8 + average
	stringsAsFactors=FALSE
)
ips.scores <- mclapply(1:nrow(parameter.space), function(i) {

	cats <- parameter.space[i,"categories"][[1]]
	sels <- parameter.space[i,"selections"][[1]]
	consolidator <- parameter.space[i,"consolidator"][[1]]

	#Find submatrix for the selected criteria
	submatrix <- with(genewise.descriptor, 
		scores.3d[,,(category %in% cats) & (selection %in% sels)]
	)
	
	switch(consolidator[[1]],
		rank={
			apply(submatrix,c(1,2),ith.rank, consolidator[[2]])
		},
		average={
			apply(submatrix,c(1,2), mean, na.rm=TRUE)
		}
	)

},mc.cores=n.cores)


###
# Obtain gold standards
#
gs <- new.goldstandards()
logger$info("Loading and translating HI3")
hi3 <- gs$load.HI3(hi3.file)
logger$info("Loading PRS")
prs <- gs$load.PRS(prs.file)
logger$info("Loading HI-II-14")
hi.ii <- gs$load.HI.II.14(hi.ii.file)

logger$info("Adapting Gold Standards to Search Space")
hi3.matrix <- gs$adapt(hi3,genewise.template,sep=genesep)
prs.matrix <- gs$adapt(prs,genewise.template,linear=FALSE,sep=genesep)
hi.ii.matrix <- gs$adapt(hi.ii,genewise.template,sep=genesep)
unionGS.matrix <- union.matrix <- apply(
	matrix.3d(list(hi3.matrix,prs.matrix)),
	c(1,2),
	function(z) {
		if(is.na(z[[1]])) {
			z[[2]] 
		}
		else if (is.na(z[[2]])) {
			z[[1]] 
		}
		else{ 
			max(z)
		}
	}
)

###
# Evaluate the maximum Matthew Correlation Coefficient (MCC) for all the results from 
# above against three gold standards: HI3, PRS and union of HI3 and PRS
#
logger$info("Evaluating consolidation options against gold standards")
mccs <- lapply(list(hi3=hi3.matrix,prs=prs.matrix,union=unionGS.matrix,hi.ii=hi.ii.matrix), function(gold.matrix) {

	do.call(rbind,mclapply(1:nrow(parameter.space), function(i) {
		ips <- ips.scores[[i]]

		pairs <- list()
		for (i in 1:nrow(ips)) {
			for (j in 1:ncol(ips)) {
				if (!is.na(gold.matrix[i,j])) {
					pairs[[length(pairs)+1]] <- c(score=ips[i,j],gold=gold.matrix[i,j])
				}
			}
		}
		pairs <- do.call(rbind,pairs)

		if (any(is.na(ips)) || min(ips)==max(ips)) {
			return(c(t=NA,mcc=NA))
		}

		ts <- seq(min(ips),max(ips),length.out=1000)
		mccs <- t(sapply(ts, function(thr) mcc(thr, pairs[,"score"], pairs[,"gold"]) ))
		idx <- which.max(mccs[,"mcc"])
		c(best.t=ts[idx],mccs[idx,])

	},mc.cores=n.cores))

})

#Pick the parameters with the best MCC in the Union set.
best.param <- which.max(mccs$hi.ii[,"mcc"])
#And pick the corresponding IPS matrix.
best.ips.matrix <- ips.scores[[best.param]]


###
# STEP 5: OUTPUT 
#

###
# Save data to hard drive
# 
logger$info("Exporting matrices")
save(scores.3d, file=paste(dir,"s-prime.rdata",sep=""))
save(best.ips.matrix, file=paste(dir,"ips.rdata",sep=""))

###
#Draw MCCs plots
#
logger$info("Drawing MCC plots")
#all categories (-AA, ALL)
cats <- unique(unlist(score.descriptor$category))
#all selections (-His, +3AT)
sels <- unique(unlist(score.descriptor$selection))
cat.sel <- expand.grid(#all possible combinations of the below
	categories=combo(cats),#all combinations of categories
	selections=combo(sels),#all combinations of selections
	stringsAsFactors=FALSE
)
mfrow <- best.grid(nrow(cat.sel))
plot.cols <- c(colorRampPalette(c("steelblue1","royalblue4"))(length(unique(parameter.space[,3]))-1),"black")
png(paste(dir,"MCCs.png",sep=""),400*mfrow[2],300*mfrow[1])
op <- par(mfrow=mfrow,las=3,cex=1.2)
invisible(lapply(1:nrow(cat.sel), function(i) {
	cat <- cat.sel[i,"categories"][[1]]
	sel <- cat.sel[i,"selections"][[1]]
	cat.match <- sapply(parameter.space$categories,function(x)all(x==cat))
	sel.match <- sapply(parameter.space$selections,function(x)all(x==sel))
	idxs <- which(cat.match & sel.match)
	names <- sapply(parameter.space[idxs,"consolidator"],paste,collapse="=")
	title <- paste(paste(cat,collapse=", "),paste(sel,collapse=", "), sep="; ")
	data <- do.call(c,lapply(mccs, function(.mccs) .mccs[idxs,"mcc"]))
	n <- length(mccs)
	xs <- barplot(
		data,
		names.arg=rep(names,n),
		main=title,
		ylim=c(0,.6),
		ylab="max(MCC)",
		col=plot.cols,
		border=NA
	)
	grid(NA,NULL)
	a <- xs[1]
	b <- xs[length(xs)]
	mtext(toupper(names(mccs)),side=1,line=0,las=0,at=a+(1+2*(0:(n-1)))*(b-a)/(2*n))
	if (best.param %in% idxs) {
		x <- xs[3*length(names)+which(idxs==best.param),]
		cat(x,"\n")
		points(x,.59,pch="*",cex=2)
	}
}))
par(op)
invisible(dev.off())


###
# Draw a sorted barchart for the distribution of the top-ranking interactions
#
logger$info("Drawing IPS barplot")
# Make a table of interactions
ips.table <- do.call(rbind,lapply(1:nrow(best.ips.matrix),function(i) {
	do.call(rbind,lapply(1:ncol(best.ips.matrix),function(j) {
		list(AD=ad.genes[[i]],DB=db.genes[[j]],IPS=best.ips.matrix[i,j])
	}))
}))
#format as proper frame
ips.table <- data.frame(
	AD=unlist(ips.table[,"AD"]),
	DB=unlist(ips.table[,"DB"]),
	IPS=unlist(ips.table[,"IPS"]),
	stringsAsFactors=FALSE
)
#sort by score
ips.table <- ips.table[order(ips.table$IPS, decreasing=TRUE),]
best.t <- mccs$union[best.param,"best.t"]
best.cutoff <- max(which(ips.table$IPS > best.t))

bc2gene <- hash(keys=barcodes$id, values=barcodes$symbol)

png(paste(dir,"ips_barplot.png",sep=""),1000,500)
layout(matrix(1:3,ncol=1),heights=c(1,.25,.5))

op <- par(mar=c(0,4,4,1)+.1,las=3,cex=.7)
n.ints <- best.cutoff*2
xs <- barplot(
	ips.table$IPS[1:n.ints],
	log="y",
	ylab="IS score",
	col="steelblue3",
	space=0
	# names.arg=paste(db.sym,"-",ad.sym)
)
x <- xs[best.cutoff,]
y <- max(ips.table$IPS) * .9
lines(x=c(x,x),y=c(ips.table$IPS[best.cutoff],y))
arrows(x0=x,x1=xs[best.cutoff-10,],y0=y,length=.1)
text(xs[best.cutoff-5,],y,paste("MCC-optimal (",best.cutoff," pairs)",sep=""),pos=1)
par(op)

adapt.gs <- function(gs.matrix) {
	sapply(1:n.ints, function(i) {
		val <- gs.matrix[ips.table$AD[i],ips.table$DB[i]]
		if (is.na(val)) {
			"gray"
		} else if (val==0) {
			"white"
		} else {
			"steelblue3"
		}
	})
}

op <- par(mar=c(0,4,0,1)+.1,cex=.7,las=2)
plot(0,type="n",axes=FALSE,xlab="",ylab="",ylim=c(0,3),xlim=c(0,n.ints))
axis(2,at=1:3-.5,labels=c("Union","HI3","PRS"),tick=FALSE)
rect(1:n.ints-1,2,1:n.ints,3,col=adapt.gs(prs.matrix))
rect(1:n.ints-1,1,1:n.ints,2,col=adapt.gs(hi3.matrix))
rect(1:n.ints-1,0,1:n.ints,1,col=adapt.gs(unionGS.matrix))
par(op)

op <- par(mar=c(0,4,0,1)+.1,cex=.7)
plot(0,type="n",axes=FALSE,xlab="",ylab="",ylim=c(0,2),xlim=c(0,n.ints))
ad.sym <- sapply(ips.table$AD[1:n.ints],function(x)bc2gene[[x]])
db.sym <- sapply(ips.table$DB[1:n.ints],function(x)bc2gene[[x]])
text((1:n.ints), 2, paste(db.sym,"-",ad.sym),srt=90,pos=2)
# text((1:n.ints)-1.3, 1.1, ad.sym,srt=90,pos=4)
# text((1:n.ints), 0.9, db.sym,srt=90,pos=2)
# axis(2,at=c(.5,1.5),labels=c("DB-X","AD-Y"),tick=FALSE)
par(op)
invisible(dev.off())


###
# Draw heatmap of best IPS score matrix
#
logger$info("Drawing IPS heatmap")
yp <- new.yogi.plotter()

png(paste(dir,"ips_heatmap.png",sep=""),2400,2000)

ad.names <- sapply(rownames(best.ips.matrix),function(name) with(barcodes,unique(symbol[id==name])))
db.names <- sapply(colnames(best.ips.matrix),function(name) with(barcodes,unique(symbol[id==name])))
ad.prs.idx <- max(which(sapply(rownames(best.ips.matrix),function(name) with(barcodes,unique(type[id==name])))=="PRS"))
db.prs.idx <- max(which(sapply(colnames(best.ips.matrix),function(name) with(barcodes,unique(type[id==name])))=="PRS"))

yp$heatmap(
	best.ips.matrix,
	row.names=ad.names,
	col.names=db.names,
	cex=.7,
	log=FALSE,
	row.separator=ad.prs.idx,
	col.separator=db.prs.idx
)

invisible(dev.off())


###
# Use IPS score matrix to construct netwoork for cytoscape
#
logger$info("Exporting data to Cytoscape format")
cy <- new.cytoscape()
cy$from.matrix(
	best.ips.matrix,
	mccs$union[best.param,"best.t"],
	out.file=paste(dir,"cytoscape.xgmml",sep=""),
	row.names=ad.names,col.names=db.names
)





###
# Quality analysis
#
logger$info("Generating report on read-pair fates")
outfile <- paste(dir,"quality.png",sep="")

if (!file.exists(outfile)) {

	con <- pipe(paste("grep Quality ",dir,"slave_*.log",sep=""), open="r")
	lines <- readLines(con)
	close(con)

	qual <- new.counter()
	qual$import.add(sub(".*Quality report: ","",lines))
	raw.nums <- qual$ls()

	con <- pipe(paste("grep Quality ",dir,"master.log",sep=""), open="r")
	lines <- readLines(con)
	close(con)

	qual$import.add(sub(".*Quality report: ","",lines))
	raw.nums <- qual$ls()

	if (length(raw.nums) > 0) {
		readwise <- with(raw.nums, list(unaligned=unaligned/(unaligned+invalid), invalid=invalid/(unaligned+invalid)))
		error.percent <- with(raw.nums, c(
				unaligned=(eliminated.pairs+incomplete.pairs)*readwise$unaligned/pairs,
				invalid=(eliminated.pairs+incomplete.pairs)*readwise$invalid/pairs,
				#bad.fusion=bad.fusion/pairs,
				missing.muxtag=missing.muxtag/pairs,
				incompatible.muxtags=incompatible.muxtags/pairs,
				incompatible.barcodes=incompatible.barcodes/pairs
			)
		)
		error.percent[["normal"]] <- 1-sum(error.percent)

		png(outfile,800,700)
		labels <- paste(names(error.percent)," ",sapply(100*error.percent,format,digits=2),"%",sep="")
		op <- par(mar=c(0,3,0,3))
		pie(error.percent,labels=labels,clockwise=TRUE,col=grey.colors(length(error.percent)),cex=.7)
		par(op)
		invisible(dev.off())
	}
}


###
# HTML output
#
interpolate <- function(string, values) {
	for (name in names(values)) {
		string <- gsub(paste("%",name,sep=""),values[[name]],string)
	}
	string
}

con <- file("html/resultpage.html",open="r")
template <- paste(readLines(con),collapse="\n")
close(con)

logger$info("Generating report page")
page <- interpolate(template, c(
	id=sub("/","",dir)#,
	# alpha=alpha,
	# rho=rho,
	# selection=select.code,
	# consolidator.options=paste(
	# 	"<option value=\"rank\" ",ifelse(consolidator=="rank","selected",""),">rank</option>\n",
	# 	"<option value=\"average\" ",ifelse(consolidator=="average","selected",""),">average</option>\n",
	# 	sep=""
	# ),
	# rank=c.rank
))

con <- file(paste(dir,"result.html",sep=""),open="w")
writeLines(page,con)
close(con)


logger$info("Done!")
