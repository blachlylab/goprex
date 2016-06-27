package main

// native imports
import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"os/exec"
	"strconv"
	"strings"
	//"sync"
	"flag"
	"log"

	// external imports
	"github.com/blachlylab/gff3"
	"github.com/mitchellh/go-homedir"
)

// GetAppris will parse the principal isoform number out of tag field
// Cannot use the gff3 'get' methods since this can be part of a comma-separated list
func GetAppris(myRecord *gff3.Record) int {
	keyword := "appris_principal_"
	tag, ok := myRecord.AttributesField["tag"]
	if !ok {
		return -1
	}
	if !strings.Contains(tag, keyword) {
		return -1
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
	// couldnt determine a winner
	return &gff3.Record{}
}

// readGzFile reads a specified gzipped Gff3 file,
// (uncompressed Gff3 files are not supported at this time)
// and returns a map of symbols (gene names) to Region definitions
func readGzFile(filename string, passGenes map[string]string) (map[string]Region, error) {
	// function to read (gzipped agnostic) files
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
	out := make(map[string]*gff3.Record)
	for featName, _ := range passGenes {
		out[featName] = &gff3.Record{}
	}
	var res *gff3.Record
	for myRecord, err := myReader.Read(); err != io.EOF; myRecord, err = myReader.Read() {
		if err != nil {
			log.Fatal(err)
		}
		// record the TSS associated with the principal isoform of every gene of interest
		for featName, featKind := range passGenes {
			if myRecord.FilterByField("type", "start_codon").FilterByAttribute(featKind, featName).Complete {
				// this record matches this feature
				res = GetTrump(out[featName], myRecord)
				if res.Complete {
					out[featName] = res
				} else {
					fmt.Println("failed ot identify principal")
				}

			}
		}
		// // if strings.Contains(myRecord.AttributesField["tag"], "appris_principal") {
		// // 	// a principal isoform
		// // 	// out[myRecord.AttributesField[field]] = Region{chrom: myRecord.SeqidField, start: myRecord.StartField, end: myRecord.EndField, strand: myRecord.StrandField}

		// }

	}
	fmt.Println(out["RUNX1"])
	return make(map[string]Region), nil
}

// validateID takes a lookup table , f (from readGzFile); and an identifier
// (i.e. gene name), v. It returns true if that id was found in the
// annotation, or false if not
func validateID(f map[string]Region, v string) bool {
	if _, ok := f[v]; ok {
		return true
	} else {
		warn("invalid or unknown identifier: " + v)
		return false
	}
	/*
		    resCount := f[v]
			if resCount == 0 {
				warn("nothing found for " + v)
				return false
			} else if resCount > 1 {
				warn("too many primary isoforms for " + v)
				//fmt.Println(f[v])
				return false
			}
			return true
	*/
}

// expandRegion takes a Region, r; and a number of nucleotides
// upstream and downstream by which to expand the Region def'n
// it returns the expanded Region
func expandRegion(r Region, up int, down int) Region {
	bedStart := 0
	bedEnd := 0
	if r.strand == '+' {
		bedStart = r.start - up
		bedEnd = r.start + down
	} else if r.strand == '-' {
		bedStart = r.end - down
		bedEnd = r.end + up
	} else {
		fmt.Println(r)
		warn("no strand found!")
	}
	out := Region{chrom: r.chrom, start: bedStart, end: bedEnd, strand: r.strand}
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
		warn("Must define upstream and/or downstream")
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
	f, err := readGzFile(gff3, passGenes)
	if err != nil {
		abort(err)
	}
	info("probing genes ...")
	//inGenes := []string{"GATA2", "DNMT3A","RUNX1","ASXL1","MADE_UP_GENE"}
	panic("oh no")
	for featName, _ := range passGenes {
		info(featName + " → " + idDescriptions[decodeId(featName)])
		if true { //validateID(f, featName)
			//wg.Add(1)
			info(featName + " found")
			outFasta := strings.Join([]string{featName, "fa"}, ".")
			doBedStuff(expandRegion(f[featName], *flagDown, *flagUp), fasta, outFasta, featName)
			info("\tdone!")
			fmt.Println()
		} else {
			fmt.Println()
		}
	}
	//wg.Wait()
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
