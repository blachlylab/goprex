package main

import "regexp"

const (
	ENST = iota
	ENSG
	UCSC
	REFSEQ
	SYMBOL
	UNKNOWN
)

var idDescriptions = map[int]string{
	ENST:    "ensembl! transcript",
	ENSG:    "ensembl! gene",
	UCSC:    "UCSC transcript id",
	REFSEQ:  "NCBI Refseq id",
	SYMBOL:  "Official gene symbol",
	UNKNOWN: "Unknown identifier type",
}

// TO DO: figure out what column to pull in case of UCSC or REFSEQ identifiers
var idGff3Names = map[int]string{
	ENST:   "transcript_id",
	ENSG:   "gene_id",
	UCSC:   "TBD",
	REFSEQ: "TBD",
	SYMBOL: "gene_name",
}

// decode identifier (a string like "ENST00000123" or "PLCG2"
// and return an int representing the type of identifier
func decodeId(identifier string) int {
	matched := false

	// BUG(james) need to fix regexes to allow ENSTNNN.X or ENSGNNN.Y
	matched, _ = regexp.MatchString("^ENST[0-9]{11}", identifier)
	if matched {
		return ENST
	}

	matched, _ = regexp.MatchString("^ENSG[0-9]{11}", identifier)
	if matched {
		return ENSG
	}

	matched, _ = regexp.MatchString("^uc[0-9]{3}[a-z]{3}\\.", identifier)
	if matched {
		return UCSC
	}

	matched, _ = regexp.MatchString("^[NX][GM]_", identifier)
	if matched {
		return REFSEQ
	}

	matched, _ = regexp.MatchString("[A-Z0-9][A-Za-z0-9]{1,}", identifier)
	if matched {
		return SYMBOL
	}

	warn("I could not understand your gene id: " + identifier)
	return UNKNOWN

}
