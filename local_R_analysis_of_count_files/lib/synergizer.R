#/usr/bin/Rscript

library("RJSONIO")
library("RCurl")

new.synergizer <- function() {

	host <- "http://llama.mshri.on.ca/cgi/synergizer/serv?"

	post <- function(query) {

		.query <- toJSON(query)

		reply <- postForm(host,
	        .opts = list(postfields=.query,
		          httpheader=c(
		          	'Content-Type'='application/json', 
		            'Accept'='application/json'
		          ),
		          ssl.verifypeer = FALSE
		     )
		)
		if (substr(reply,1,1)!= "{") {
			stop(reply)
		} else {
			result <- fromJSON(reply)
			if (is.null(result$error)) {
				result$result
			} else {
				stop(result$error)
			}
		}
	}

	version <- function() {
		post(list(
			method="version",
			params=list(),
			id=0
		))
	}

	translate <- function(authority,species,domain,range,ids) {
		result <- post(list(
			method="translate",
			params=list(list(
				authority=authority,
				species=species,
				domain=domain,
				range=range,
				ids=ids
			)),
			id=0
		))
		out <- lapply(result,function(x) unlist(x[-1]))
		names(out) <- lapply(result,function(x) x[[1]])
		out
	}

	availableAuthorities <- function() {
		post(list(
			method="available_authorities",
			params=list(),
			id=0
		))
	}

	availableSpecies <- function(authority) {
		post(list(
			method="available_species",
			params=list(
				authority
			),
			id=0
		))
	}

	availableDomains <- function(authority,species) {
		post(list(
			method="available_domains",
			params=list(
				authority,
				species
			),
			id=0
		))
	}

	availableRanges <- function(authority,species,domain) {
		post(list(
			method="available_ranges",
			params=list(
				authority,
				species,
				domain
			),
			id=0
		))
	}

	structure(list(
		version=version,
		availableAuthorities=availableAuthorities,
		availableSpecies=availableSpecies,
		availableDomains=availableDomains,
		availableRanges=availableRanges,
		translate=translate
	),class="synergizer")

}


test.synergizer <- function() {

	syn <- new.synergizer()

	cat(syn$version(),"\n")

	syn$translate(
		"ncbi","Homo sapiens","uniprot","entrezgene",
		list(
			"P00338","P00367","P00374","P00390","P00403","P00439"
		)
	) 

}

