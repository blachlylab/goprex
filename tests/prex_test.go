package tests

import (
	"os"
	"reflect"
	"testing"

	prex "../../goprex"
	"github.com/blachlylab/gff3"
)

var myRecord = gff3.Record{
	Complete:    true,
	SeqidField:  "chr1",
	SourceField: "HAVANA",
	TypeField:   "CDS",
	StartField:  3499796,
	EndField:    3499924,
	ScoreField:  0,
	StrandField: '-',
	PhaseField:  2,
	AttributesField: map[string]string{
		"tag":                      "basic,appris_principal_1,CCDS",
		"gene_id":                  "ENSG00000162591.15",
		"transcript_id":            "ENST00000356575.8",
		"gene_type":                "protein_coding",
		"transcript_type":          "protein_coding",
		"exon_number":              "22",
		"havana_gene":              "OTTHUMG00000000611.7",
		"gene_name":                "MEGF6",
		"transcript_status":        "KNOWN",
		"exon_id":                  "ENSE00001206792.1",
		"protein_id":               "ENSP00000348982.4",
		"transcript_support_level": "1",
		"havana_transcript":        "OTTHUMT00000354866.1",
		"ID":                       "CDS:ENST00000356575.8:22",
		"Parent":                   "ENST00000356575.8",
		"gene_status":              "KNOWN",
		"transcript_name":          "MEGF6-007",
		"level":                    "2",
		"ccdsid":                   "CCDS41237.1",
	},
}

var newRecord = gff3.Record{
	Complete:    true,
	SeqidField:  "chr1",
	SourceField: "HAVANA",
	TypeField:   "CDS",
	StartField:  3499796,
	EndField:    3499924,
	ScoreField:  0,
	StrandField: '-',
	PhaseField:  2,
	AttributesField: map[string]string{
		"tag":                      "basic,appris_principal_3,CCDS",
		"gene_id":                  "ENSG00000162591.15",
		"transcript_id":            "ENST00000356575.8",
		"gene_type":                "protein_coding",
		"transcript_type":          "protein_coding",
		"exon_number":              "22",
		"havana_gene":              "OTTHUMG00000000611.7",
		"gene_name":                "MEGF6",
		"transcript_status":        "KNOWN",
		"exon_id":                  "ENSE00001206792.1",
		"protein_id":               "ENSP00000348982.4",
		"transcript_support_level": "1",
		"havana_transcript":        "OTTHUMT00000354866.1",
		"ID":                       "CDS:ENST00000356575.8:22",
		"Parent":                   "ENST00000356575.8",
		"gene_status":              "KNOWN",
		"transcript_name":          "MEGF6-007",
		"level":                    "2",
		"ccdsid":                   "CCDS41237.1",
	},
}

func TestPrincipalIsoformIdentification(t *testing.T) {
	if prex.GetAppris(&myRecord) != 1 {
		t.Errorf("error parsing principal isoform tag")
	}
}

func TestPrincipalTrumpRule(t *testing.T) {
	res := prex.GetTrump(&newRecord, &myRecord)
	if !reflect.DeepEqual(res, &myRecord) || !res.Complete {
		t.Errorf("GetTrump failed to identify the principal isoform")
	}
	newRecord.AttributesField["tag"] = myRecord.AttributesField["tag"]
	newRecord.AttributesField["transcript_support_level"] = "3"
	res = prex.GetTrump(&newRecord, &myRecord)
	if !reflect.DeepEqual(res, &myRecord) || !res.Complete {
		t.Errorf("GetTrump failed to identify the principal isoform")
	}
	newRecord.AttributesField["transcript_support_level"] = myRecord.AttributesField["transcript_support_level"]
	newRecord.AttributesField["level"] = "1"
	res = prex.GetTrump(&newRecord, &myRecord)
	if !reflect.DeepEqual(res, &newRecord) || !res.Complete {
		t.Errorf("GetTrump failed to identify the principal isoform")
	}
	newRecord.AttributesField["level"] = myRecord.AttributesField["level"]
	res = prex.GetTrump(&newRecord, &myRecord)
	if !reflect.DeepEqual(res, &gff3.Record{}) || res.Complete {
		t.Errorf("GetTrump incorrectly found a principal isoform between two, identical structs")
	}
}

func TestMain(m *testing.M) {
	os.Exit(m.Run())
}
