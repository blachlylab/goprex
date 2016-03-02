package main
import (
    "fmt"
    "bufio"
    "strings"
    "strconv"
    "compress/gzip"
    "os"
    "os/exec"
    "io/ioutil"
    "encoding/json"
    "sync"
    "flag"
)

type feat struct {
    reg region
    feature string 
    strand string 
    geneID string
    geneName string
    transcriptID string 
    appris int 
} 

type region struct {
    chrom string 
    start int 
    end int 
    strand string
}

func getGene(line string) feat {
    spl     := strings.Split(line,"\t")
    chrom   := spl[0]
    feature := spl[2]
    start,_ := strconv.Atoi(spl[3])
    end,_   := strconv.Atoi(spl[4])
    strand  := spl[6]
    geneID  := ""
    transcriptID := ""
    geneName := ""
    appris := 0
    tags := strings.Split(spl[8],";")
    for _,v := range(tags) {
        tag := strings.Split(v, "=")
        if tag[0] == "gene_id" {
            geneID = tag[1]
        } else if tag[0] == "transcript_id" {
            transcriptID = tag[1]
        } else if tag[0] == "gene_name" {
            geneName = tag[1]
        } else if strings.Contains(tag[1], "appris_principal") {
            for _,field := range(strings.Split(tag[1],",")) {
                if strings.Contains(field, "appris_principal") {
                    appris,_ = strconv.Atoi(strings.Split(field, "_")[2])
                }
            }
             
        }
    }
    thisFeat := feat{
        reg: region{chrom: chrom, start: start, end: end, strand:strand},
        feature: feature,
        geneName: geneName,
        geneID: geneID,
        transcriptID: transcriptID,
        appris: appris,
    }
    return thisFeat 
}

func appendIfNew(refList []region, addition region ) []region{
    for _,v := range(refList) {
        if v == addition {
            return refList 
        }
    }
    return append(refList, addition)
}

func readGzFile(filename string) (map[string][]region, error) {
    // function to read gzipped files
    fi, err := os.Open(filename)
    if err != nil {
        return nil, err
    }
    defer fi.Close()
    fz, err := gzip.NewReader(fi)
    // TODO: detect file format (i.e. .gz) and dynamically open or gzip.open as necessary 
    if err != nil {
        return nil, err
    }
    defer fz.Close()
    scanner := bufio.NewScanner(fz)
    out := make(map[string][]region )
    for scanner.Scan() {
        line := scanner.Text()
        if !strings.HasPrefix(line, "#") {
            thisFeat := getGene(line)
            if (thisFeat.appris > 0) && (thisFeat.feature == "start_codon") {
                // add redundant copies of this feature to the map with 
                // gene_id, transcript_id, and gene_name keys  
                // TODO: only add the features that were detected? 
                //       make a channel for each feature type? 
		//	 maybe use pointers to avoid duplicating data in memory?
                out[thisFeat.geneID] = appendIfNew(out[thisFeat.geneID], thisFeat.reg )
                out[thisFeat.transcriptID] = appendIfNew(out[thisFeat.transcriptID], thisFeat.reg )
                out[thisFeat.geneName] = appendIfNew(out[thisFeat.geneName], thisFeat.reg )
            }
        }
    }  
    return out, nil
}

func validateID(f map[string][]region, v string) bool {
    resCount := len(f[v])
    if resCount == 0 {
        warn("nothing found for " + v)
        return false
    } else if resCount > 1 {
        warn("too many primary isoforms for " + v)
        //fmt.Println(f[v])
        return false
    }
    return true
}

func doGff3Stuff(r region, up int, down int) region {
    bedStart := 0
    bedEnd := 0
    if r.strand == "+" {
        bedStart = r.start - up
        bedEnd   = r.start + down 
    } else if r.strand == "-" {
        bedStart = r.start - down 
        bedEnd   = r.start + up
    } else {
        warn("no strand found!")
    }
    out := region{chrom:r.chrom, start:bedStart, end:bedEnd, strand:r.strand}
    return out
}

func doBedStuff(r region, fastaIn string, fastaOut string, name string) {
     tempDir := os.TempDir()
     tempFile, err := ioutil.TempFile(tempDir, "prex_")
     if err != nil {
     	abort(err)
     }
     defer os.Remove(tempFile.Name())
     bedName := name + ";" + r.chrom + ":" + strconv.Itoa(r.start) + "-" + strconv.Itoa(r.end) + "(" + r.strand + ")"
     bedString := strings.Join([]string{r.chrom, strconv.Itoa(r.start), strconv.Itoa(r.end), bedName, ".", r.strand}, "\t")
     err = ioutil.WriteFile(tempFile.Name(), []byte(bedString + "\n"), 600)
     if err != nil {
     	abort(err)
     }
     _, err = exec.Command("bedtools", "getfasta", "-name", "-s", "-fi", fastaIn, "-bed", tempFile.Name(), "-fo", fastaOut).Output()
     if err != nil {
     	abort(err)
     }
     wg.Done()
}

func loadConfig() map[string]string {
     var config struct {
     	 Gff3 string
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
     return map[string]string{"fasta": config.Fasta, "gff3":config.Gff3}
}

var wg sync.WaitGroup

func main() {
    flagGff3  := flag.String("gff3", "", "gtf annotation")
    flagFasta := flag.String("fasta", "", "fasta sequence file")     
    flagUp    := flag.Int("up", 0, "upstream distance")
    flagDown  := flag.Int("down", 0, "downstream distance")
    flag.Parse()
    inGenes := flag.Args()

    if len(inGenes) < 1 {
       warn("No arguments found! Pass some feature names!")    
       os.Exit(1)
    }
    if *flagUp == 0 *flagDown == 0 {
        warn("Must define upstream and/or downstream")
    }

    config := loadConfig()
    fasta  := config["fasta"]
    if *flagFasta != "" {
       // if command line flag is provided, take it
       fasta = *flagFasta
    }
    if _, err := os.Stat(fasta); err != nil {
       abort(err)
    }
    gff3   := config["gff3"]
    if *flagGff3 != "" {
       // if command line flag is provided, take it
       gff3 = *flagGff3 
    }
    info("reading " + gff3)
    f, err := readGzFile(gff3)
    if err != nil {
        abort(err)
    }
    info("probing genes ...")
    //inGenes := []string{"GATA2", "DNMT3A","RUNX1","ASXL1","MADE_UP_GENE"} 
    
    for _,v := range(inGenes) {
       if validateID(f,v) {
       	   wg.Add(1)
           info(v)
	   outFasta := strings.Join([]string{v, "fa"},".")
           go doBedStuff(doGff3Stuff(f[v][0], *flagDown, *flagUp), fasta, outFasta, v)
           info("\tdone!")
	   fmt.Println()
       } else {
       	   fmt.Println()
       }
    }
    wg.Wait()
}

func info(message string) {
    fmt.Println("[ok] " + message)
}
    
func warn(message string) {
    fmt.Println("[* ] " + message)
}

func abort(message error) {
    fmt.Println("[!!] " + message.Error())
    os.Exit(1)
}
