package main

import (
	// native imports
	"bufio"
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"strconv"
	"strings"

	// external imports
	"github.com/blachlylab/gff3"
	"github.com/mitchellh/go-homedir"
)

// GetAppris will parse the principal isoform number out of tag field
// Cannot use the gff3 'get' methods since this can be part of a comma-separated list
func GetAppris(myRecord *gff3.Record) int {
	keyword := "appris_principal_"
	tag, ok := myRecord.AttributesField["tag"]
	if !ok || !strings.Contains(tag, keyword) {
		return 99
	}
	i := strings.Index(tag, keyword) + len(keyword)
	out, err := strconv.Atoi(tag[i : i+1])
	if err != nil {
		log.Fatal("cannot convert " + tag[i:i+1] + " to int")
	}
	return out
}

// GetTrump will determine which of two GFF3 records is more likely the "principal" isoform
// Checks principal tag, transcript support level, and level
// Based on info from here: http://www.gencodegenes.org/faq.html
// TODO: this might be faster returning a bool rather than a record
//       but then the order of records input will be important
func GetTrump(rec1, rec2 *gff3.Record) *gff3.Record {
	// is either of them empty?
	if !rec1.Complete {
		return rec2
	} else if !rec2.Complete {
		return rec1
	}
	// is one of them a lower principal isoform?
	var principal [2]int
	principal[0] = GetAppris(rec1)
	principal[1] = GetAppris(rec2)
	if principal[0] > principal[1] {
		return rec2
	} else if principal[1] > principal[0] {
		return rec1
	}
	// does one of them have lower transcript support level?
	var support [2]int
	// missing key means TSL = 0; there is no such thing in the specs
	// http://useast.ensembl.org/Help/Glossary?id=492
	support[0], _ = strconv.Atoi(rec1.AttributesField["transcript_support_level"])
	support[1], _ = strconv.Atoi(rec2.AttributesField["transcript_support_level"])
	if support[0] > support[1] {
		return rec2
	} else if support[1] > support[0] {
		return rec1
	}
	// is one of them lower level?
	var level [2]int
	level[0], _ = strconv.Atoi(rec1.AttributesField["level"])
	level[1], _ = strconv.Atoi(rec2.AttributesField["level"])
	if level[0] > level[1] {
		return rec2
	} else if level[1] > level[0] {
		return rec1
	}

	// are they effectively equivalent?
	if rec1.StartField == rec2.StartField && rec1.EndField == rec2.EndField {
		// if they are completely identical, including start and stop, arbitrarily pick one
		return rec1
	}

	// couldnt determine a winner
	panic("could not determine a winner")
}

// readGff3File reads a Gff3 file,
// and returns a map of symbols (gene names) to gff3.Records
func readGff3File(filename string, passGenes map[string]string) (map[string]*gff3.Record, error) {
	// function to read a(n optionally gzipped) gff3 file
	// load gff3 reader
	var myReader *gff3.Reader
	fi, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer fi.Close()
	if strings.HasSuffix(filename, ".gz") {
		// gzipped file needs to be decompressed
		fz, err := gzip.NewReader(fi)
		if err != nil {
			return nil, err
		}
		defer fz.Close()
		myReader = gff3.NewReader(fz)
	} else {
		// no gzip; pass reader to gff3 directly
		myReader = gff3.NewReader(fi)
	}
	// reader is ready; begin parsing
	out := make(map[string]*gff3.Record)
	for featName, _ := range passGenes {
		out[featName] = &gff3.Record{}
	}
	var res *gff3.Record
	// name this loop so it can be continued from the nested loop
recordLoop:
	for myRecord, err := myReader.Read(); err != io.EOF; myRecord, err = myReader.Read() {
		if err != nil {
			log.Fatal(err)
		}
		if !myRecord.FilterByField("type", "transcript").Complete {
			// if this isn't a start codon, skip immediately
			continue recordLoop
		}
		// is this record associated with any of the input features?
		for featName, featKind := range passGenes {
			// warning: this filter modifies the record in memory
			// 			must reset the Complete flag before checking the next gene
			myRecord.FilterByAttribute(featKind, featName)
			if myRecord.Complete {
				// this record matches this feature
				// identify the winning isoform encountered so far
				res = GetTrump(out[featName], myRecord)
				if res.Complete {
					out[featName] = res
				} else if out[featName].Complete && !res.Complete {

					// this gene has two principal isoforms
					warn("cannot process " + featName + " having multiple principal isoforms!")
				}

				// we know where this record belongs; continue on to the next gff3 row
				continue recordLoop
			} else {
				//reset the complete flag and check the next gene
				myRecord.Complete = true
			}
		}
	}
	return out, nil
}

// expandRegion takes a gff3 Record, r; and a number of nucleotides
// upstream and downstream by which to expand the Region def'n
// it returns an expanded Region struct
func expandRegion(r *gff3.Record, up int, down int) Region {
	bedStart := 0
	bedEnd := 0
	if r.StrandField == '+' {
		bedStart = r.StartField - up
		bedEnd = r.StartField + down
	} else if r.StrandField == '-' {
		bedStart = r.EndField - down
		bedEnd = r.EndField + up
	} else {
		fmt.Println(r)
		warn("no strand found!")
	}
	out := Region{chrom: r.SeqidField, start: bedStart, end: bedEnd, strand: r.StrandField}
	return out
}

// doBedStuff builds a temporary BED file containing the Region of interest
// and executes bedtools' getfasta command
func doBedStuff(r Region, fastaIn string, fastaOut string, name string) {
	log.Println("doBedStuff() name: " + name)
	log.Println("doBedStuff() r.start: " + strconv.Itoa(r.start))
	tempDir := os.TempDir()
	tempFile, err := ioutil.TempFile(tempDir, "prex_")
	if err != nil {
		abort(err)
	}
	defer os.Remove(tempFile.Name())

	bedName := name + ";" + r.chrom + ":" + strconv.Itoa(r.start) + "-" + strconv.Itoa(r.end) + "(" + string(r.strand) + ")"
	bedString := strings.Join([]string{r.chrom, strconv.Itoa(r.start), strconv.Itoa(r.end), bedName, ".", string(r.strand)}, "\t")
	err = ioutil.WriteFile(tempFile.Name(), []byte(bedString+"\n"), 600)
	if err != nil {
		abort(err)
	}
	// if the region length is 0, bedtools will skip the feature and write an empty fasta
	_, err = exec.Command("bedtools", "getfasta", "-name", "-s", "-fi", fastaIn, "-bed", tempFile.Name(), "-fo", fastaOut).Output()
	if err != nil {
		abort(err)
	}
	//wg.Done()
}

func loadConfig() map[string]string {
	var config struct {
		Gff3  string
		Fasta string
	}
	file, err := os.Open("./prex.json")
	if err != nil {
		abort(err)
	}
	jsonParser := json.NewDecoder(file)
	if err = jsonParser.Decode(&config); err != nil {
		abort(err)
	}
	config.Fasta, _ = homedir.Expand(config.Fasta)
	config.Gff3, _ = homedir.Expand(config.Gff3)
	return map[string]string{"fasta": config.Fasta, "gff3": config.Gff3}
}

//var wg sync.WaitGroup

func main() {
	flagGff3 := flag.String("gff3", "", "gtf annotation")
	flagFasta := flag.String("fasta", "", "fasta sequence file")
	flagUp := flag.Int("up", 0, "upstream distance")
	flagDown := flag.Int("down", 0, "downstream distance")
	flag.Parse()
	inGenes := flag.Args()

	if len(inGenes) < 1 {
		warn("No arguments found! Pass some feature names!")
		os.Exit(1)
	} else if len(inGenes) == 1 {
		// is this a file?
		fi, err := os.Open(inGenes[0])
		if err == nil {
			// this appears to be a file
			defer fi.Close()
			var geneFileGenes []string
			scanner := bufio.NewScanner(fi)
			for scanner.Scan() {
				line := scanner.Text()
				if strings.TrimSpace(line) != "" {
					geneFileGenes = append(geneFileGenes, line)
				}
				inGenes = geneFileGenes
			}
		}
		// otherwise, assume this is a gene identifier
	}

	if *flagUp == 0 && *flagDown == 0 {
		warn("Must define upstream and/or downstream; default to start codon")
		*flagDown = 3
	}

	config := loadConfig()
	fasta := config["fasta"]
	if *flagFasta != "" {
		// if command line flag is provided, take it
		fasta = *flagFasta
	}
	if _, err := os.Stat(fasta); err != nil {
		abort(err)
	}
	gff3 := config["gff3"]
	if *flagGff3 != "" {
		// if command line flag is provided, take it
		gff3 = *flagGff3
	}

	passGenes := make(map[string]string)
	for _, v := range inGenes {
		passGenes[v] = idGff3Names[decodeId(v)]
	}
	info("reading " + gff3)
	f, err := readGff3File(gff3, passGenes)
	if err != nil {
		abort(err)
	}
	info("probing genes ...")
	//inGenes := []string{"GATA2", "DNMT3A","RUNX1","ASXL1","MADE_UP_GENE"}
	for featName, _ := range passGenes {
		info(featName + " → " + idDescriptions[decodeId(featName)])
		if f[featName].Complete {
			info(featName + " found")
			outFasta := strings.Join([]string{featName, "fa"}, ".")
			doBedStuff(expandRegion(f[featName], *flagDown, *flagUp), fasta, outFasta, featName)
			info("\tdone!")
			fmt.Println()
		} else {
			warn("no gene found for " + featName)
		}
	}
}

// Atoi returns also an err parameter
// Need this mustAtoi form to use inline
func mustAtoi(s string) int {
	i, err := strconv.ParseInt(s, 0, 0)
	if err != nil {
		warn("mustAtoi() " + s)
		abort(err)
	}
	return int(i)
}

func info(message string) {
	log.Println("[ok] " + message)
}

func warn(message string) {
	log.Println("[* ] " + message)
}

func abort(err error) {
	log.Fatalln("[!!] " + err.Error())
}
