package main

type Region struct {
	chrom  string
	start  int
	end    int
	strand byte
}

//methods bound to the Region struct

// is the Region totally uninitialized?
func (r Region) isEmpty() bool {
	if r.start == 0 && r.end == 0 && r.chrom == "" && r.strand == 0 {
		return true
	}
	return false
}

// is the Region "greater than" the passed Region?
// BUG(karl) I am not sure this is complete?
func (r Region) greaterThan(inR Region) bool {
	//if r's Region is bigger than inR's Region, return true
	if r.chrom == inR.chrom && r.strand == inR.strand {
		if r.start < inR.start {
			return true
		} else if r.end > inR.end {
			return true
		}
	} else if inR.isEmpty() {
		// if the checked Region is empty, r has to be bigger
		return true
	}
	return false
}

// take the union of r and inR to create the largest possible Region
func (r Region) expandTo(inR Region) Region {

	// take the union of r and inR to create the largest possible Region
	if r.chrom == inR.chrom && r.strand == inR.strand {
		if r.start > inR.start {
			r.start = inR.start
		} else if r.end < inR.end {
			r.end = inR.end
		}
	} else if r.start == 0 && r.end == 0 {
		r = inR
	}
	return r
}

func appendIfNew(refList []Region, addition Region) []Region {
	for _, v := range refList {
		if v == addition {
			return refList
		}
	}
	return append(refList, addition)
}

func expandIfNew(runningRegion, addition Region) Region {
	if addition.greaterThan(runningRegion) {
		return runningRegion.expandTo(addition)
	} else {
		return runningRegion
	}
}
