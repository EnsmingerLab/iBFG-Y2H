library("XML")

###
# This object can be used to crate cytoscape files from matrices
new.cytoscape <- function() {

	to.edge.list <- function(mat,t, 
			row.names=rownames(mat),col.names=colnames(mat)) {

		edges <- list()
		for (i in 1:nrow(mat)) {
			for (j in 1:ncol(mat)) {
				if (mat[i,j] > t) {
					edges[[length(edges)+1]] <- c(
						n1=rownames(mat)[i],
						n1.name=row.names[[i]],
						n2=colnames(mat)[j],
						n2.name=col.names[[j]],
						weight=mat[i,j]
					)
				}
			}
		}
		edges <- do.call(rbind,edges)
		w <- as.numeric(edges[,"weight"])
		width <- 1 + (w-min(w)) * 39 / (max(w)-min(w))
		cbind(edges,width=width)
	}

	grid.layout <- function(edges) {

		nodes <- unique(c(edges[,"n1"],edges[,"n2"]))
		lookup <- rbind(edges[,c("n1","n1.name")], edges[,c("n2","n2.name")])
		names <- sapply(nodes, function(id) lookup[lookup[,1]==id,2][[1]])

		d <- 100
		r <- floor(sqrt(length(nodes)))
		i <- 0:(length(nodes)-1)
		x <- d * floor(i / r)
		y <- d * (i %% r)

		data.frame(id=i,canonical=nodes,label=names,x=x,y=y)
	}

	from.matrix <- function(mat, t, out.file="network.xgmml", 
			row.names=rownames(mat),col.names=colnames(mat), 
			graph.name="Network") {

		edges <- to.edge.list(mat,t, row.names=row.names, col.names=col.names)
		nodes <- grid.layout(edges)
		canon2id <- data.frame(id=nodes[,"id"],row.names=nodes[,"canonical"])

		ns <- "http://www.cs.rpi.edu/XGMML"
		nss <- c(
			ns,
			dc="http://purl.org/dc/elements/1.1/", 
			xlink="http://www.w3.org/1999/xlink", 
			rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#", 
			cy="http://www.cytoscape.org"
		)
		doc <- newXMLDoc(namespaces=nss)

		z <- xmlTree("graph", attrs=c(label=graph.name, directed="1"), namespaces=nss, doc=doc)
		z$addNode("att", attrs=c(name="documentVersion", value="1.1"))
		z$addNode("att", attrs=c(name="networkMetadata"), close=FALSE)
			z$addNode("RDF",namespace="rdf", close=FALSE)
				z$addNode("Description",namespace="rdf", attrs=c(`rdf:about`="http://www.cytoscape.org/"),close=FALSE)
					z$addNode("type","Protein-Protein Interaction",namespace="dc")
					z$addNode("description","N/A",namespace="dc")
					z$addNode("identifier","N/A",namespace="dc")
					z$addNode("date",format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),namespace="dc")
					z$addNode("title",graph.name,namespace="dc")
					z$addNode("source","http://www.cytoscape.org/",namespace="dc")
					z$addNode("format","Cytoscape-XGMML",namespace="dc")
				z$closeTag()
			z$closeTag()
		z$closeTag()
		z$addNode("att", attrs=c(type="string", name="backgroundColor", value="#ffffff"))
		z$addNode("att", attrs=c(type="real", name="GRAPH_VIEW_ZOOM", value="1.0"))
		z$addNode("att", attrs=c(type="real", name="GRAPH_VIEW_CENTER_X", value="0.0"))
		z$addNode("att", attrs=c(type="real", name="GRAPH_VIEW_CENTER_Y", value="0.0"))
		z$addNode("att", attrs=c(type="boolean", name="NODE_SIZE_LOCKED", value="true"))

		apply(nodes, 1, function(node) {
			z$addNode("node", attrs=c(label=node[["canonical"]], id=node[["id"]]), close=FALSE)
				z$addNode("att", attrs=c(type="string",name="NODE_TYPE",value="DefaultNode"))
				z$addNode("att", attrs=c(type="string",name="canonicalName",value=node[["canonical"]]))
				z$addNode("att", attrs=c(type="string",name="label",value=node[["label"]]))
				z$addNode("graphics", attrs=c(type="ELLIPSE", 
					h="40.0", w="40.0", 
					x=node[["x"]],y=node[["y"]],
					fill="#cccccc", width="1", outline="#666666"
				))
			z$closeTag()
			node
		})

		apply(edges,1, function(edge) {
			label <- paste(edge["n1"],"(DirectedEdge)",edge["n2"])
			z$addNode("edge", attrs=c(label=label[[1]], source=canon2id[edge["n1"],1], target=canon2id[edge["n2"],1]), close=FALSE)
				z$addNode("att", attrs=c(type="string",name="canonicalName",value=label))
				z$addNode("att", attrs=c(type="real",name="score",value=edge[["weight"]]))
				z$addNode("att", attrs=c(type="string",name="interaction",value="DirectedEdge"))
				z$addNode("graphics", attrs=c(width=edge[["width"]],fill="#000000"))
			z$closeTag()
			edge
		})

		saveXML(doc,file=out.file)

	}

	structure(list(
		from.matrix = from.matrix
	),class="cytoscape")
}

