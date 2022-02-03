#!/usr/bin/Rscript

new.yogi.plotter <- function() {

	###
	# this function draws a heatmap. 
	#
	heatmap <- function(data,cex=.5,color.stops=c("white","black"),
		row.names=rownames(data),col.names=colnames(data),log=FALSE,
		col.separator=NULL, row.separator=NULL) {
		#define the number of colors to use
		ncolors <- 10

		#create color palette
		colors <- colorRampPalette(color.stops)(ncolors)

		#SET UP PLOT LAYOUT
		#set up layout with main plot taking 90% of the space and legend 10%
		layout(matrix(c(1,2),nrow=1),widths=c(9,1))

		#MAIN PLOT
		if (log) {
			data <- log10(data)
		}
		#set margins for the main plot, and set axis labels to perpendicular mode
		op <- par(mar=c(2,6,6,0),las=2)
		#draw heatmap with correct orientation
		image(t(data)[,nrow(data):1], col=colors, axes=FALSE)
		#draw separators
		if (!is.null(row.separator)) {
			b <- 1/(nrow(data)-1)
			abline(h=(nrow(data)-row.separator)*b-b/2,lty="dashed",col="gray")
		}
		if (!is.null(col.separator)) {
			b <- 1/(ncol(data)-1)
			abline(v=col.separator*b-b/2,lty="dashed",col="gray")
		}
		#add axis labels
		axis(3,at=seq(0,1,length.out=ncol(data)),labels=col.names,cex.axis=cex)
		axis(2,at=seq(0,1,length.out=nrow(data)),labels=rev(row.names),cex.axis=cex)
		#restore graphics parameters
		par(op)

		#LEGEND

		#set margins and label orientation for legend
		op <- par(las=1,mar=c(2,1,4,4))
		#create a new plot that will become the legend
		plot(0,type="n",xlim=c(0,1),ylim=c(0,ncolors),axes=FALSE,ylab="",xlab="")
		#draw colored rectangles along the y axis of the new plot
		rect(0,0:(ncolors-1),1,1:ncolors,col=colors,border=NA)
		#add an axis label that associates the rectangles with the values, round them to 2 decimal digits
		labels <- signif(seq(min(data),max(data),length.out=ncolors),2)
		if (log) {
			labels <- parse(text=sprintf("10^%s",labels))
		}
		axis(4,at=(1:ncolors)-.5,labels=labels,cex.axis=cex)
		#restore graphics parameters
		par(op)

	}

	structure(list(
		heatmap=heatmap
	),class="yogiplot")

}


###
# This function computes the best 2D layout for a given number
# of images to be drawn together.
# n = number of images
# returns: a vector of two integers representing rows and columns
# of the best layout.
#
best.grid <- function(n,l=1) {
	r <- ceiling(sqrt(n))
	factors <- do.call(rbind,lapply((r-l):(r+l),function(rows) {
		do.call(rbind,lapply(rows:(r+l), function(cols) {
			val <- rows*cols - n
			c(rows=rows, cols=cols, fit=ifelse(val<0, Inf, val) )
		}))
	}))
	factors[which.min(factors[,"fit"]),c("rows","cols")]
}
