## Gonetics

Gonetics is a bioinformatics library providing basic data structures for handling genetic data. The documentation is available [here](https://godoc.org/github.com/pbenner/gonetics).

### GRanges

Create a GRanges object with three ranges on the first chromosome:

```go
  seqnames := []string{"chr1", "chr1", "chr1"}
  from     := []int{100000266, 100000271, 100000383}
  to       := []int{100000291, 100000296, 100000408}
  strand   := []byte{'+', '+', '-'}

  granges  := NewGRanges(seqnames, from, to, strand)
  fmt.Println(granges)
```
	  seqnames                 ranges strand
	1     chr1 [100000266, 100000291)      +
	2     chr1 [100000271, 100000296)      +
	3     chr1 [100000383, 100000408)      -

Add some meta data to the GRanges object:

```go
  granges.AddMeta("data", []float64{1.0, 2.0, 3.0})
```
	  seqnames                 ranges strand |          data
	1     chr1 [100000266, 100000291)      + |      1.000000
	2     chr1 [100000271, 100000296)      + |      2.000000
	3     chr1 [100000383, 100000408)      - |      3.000000

Find overlaps of two GRanges objects:

```go
  rSubjects := NewGRanges(
    []string{"chr4", "chr4", "chr4", "chr4"},
    []int{100, 200, 300, 400},
    []int{150, 250, 350, 450},
    []byte{})
  rQuery := NewGRanges(
    []string{"chr1", "chr4", "chr4", "chr4", "chr4", "chr4"},
    []int{100, 110, 190, 340, 390, 450},
    []int{150, 120, 220, 360, 400, 500},
    []byte{})

  queryHits, subjectHits := FindOverlaps(rQuery, rSubjects)
```
	  queryHits: [1 2 3 4 5]
	subjectHits: [0 1 2 3 3]

### Genes

Download gene list from UCSC and export it to file:

```go
  genes := ImportGenesFromUCSC("hg19", "ensGene")
  genes.WriteTable("hg19.knownGene.txt", true, false)
  fmt.Println(genes)
```

	                 names seqnames          transcripts                  cds strand
	     1 ENST00000456328     chr1 [   11868,    14409) [   14409,    14409)      +
	     2 ENST00000515242     chr1 [   11871,    14412) [   14412,    14412)      +
	     3 ENST00000518655     chr1 [   11873,    14409) [   14409,    14409)      +
	     4 ENST00000450305     chr1 [   12009,    13670) [   13670,    13670)      +
	     5 ENST00000423562     chr1 [   14362,    29370) [   29370,    29370)      -
	                   ...      ...                  ...                  ...       
	204936 ENST00000420810     chrY [28695571, 28695890) [28695890, 28695890)      +
	204937 ENST00000456738     chrY [28732788, 28737748) [28737748, 28737748)      -
	204938 ENST00000435945     chrY [28740997, 28780799) [28780799, 28780799)      -
	204939 ENST00000435741     chrY [28772666, 28773306) [28773306, 28773306)      -
	204940 ENST00000431853     chrY [59001390, 59001635) [59001635, 59001635)      +

Import expression data from a GTF file:

```go
  genes.ReadGTF("genesExpr_test.gtf.gz", "transcript_id", "FPKM", false)
```

	                 names seqnames          transcripts                  cds strand |          expr
	     1 ENST00000456328     chr1 [   11868,    14409) [   14409,    14409)      + |      0.073685
	     2 ENST00000515242     chr1 [   11871,    14412) [   14412,    14412)      + |      0.000000
	     3 ENST00000518655     chr1 [   11873,    14409) [   14409,    14409)      + |      0.000000
	     4 ENST00000450305     chr1 [   12009,    13670) [   13670,    13670)      + |      0.000000
	     5 ENST00000423562     chr1 [   14362,    29370) [   29370,    29370)      - |     10.413931
	                   ...      ...                  ...                  ...        |           ...
	204936 ENST00000420810     chrY [28695571, 28695890) [28695890, 28695890)      + |      0.000000
	204937 ENST00000456738     chrY [28732788, 28737748) [28737748, 28737748)      - |      0.000000
	204938 ENST00000435945     chrY [28740997, 28780799) [28780799, 28780799)      - |      0.000000
	204939 ENST00000435741     chrY [28772666, 28773306) [28773306, 28773306)      - |      0.000000
	204940 ENST00000431853     chrY [59001390, 59001635) [59001635, 59001635)      + |      0.000000

### Peaks

Import peaks from a MACS xls file:

```go
  peaks := ReadXlsPeaks("peaks_test.xls")
```
	   seqnames             ranges strand |  abs_summit     pileup -log10(pvalue) fold_enrichment -log10(qvalue)
	 1       2L [   5757,    6001)      * |        5865  33.000000      19.809300        6.851880      17.722200
	 2       2L [  47233,   47441)      * |       47354  36.000000      19.648200        6.263150      17.566200
	 3       2L [  66379,   67591)      * |       66957 252.000000     350.151250       50.986050     346.525450
	 4       2L [  72305,   72838)      * |       72525 170.000000     208.558240       34.460930     205.734390
	 5       2L [  72999,   73218)      * |       73130  25.000000      12.711700        5.239670      10.700880
	        ...                ...        |         ...        ...            ...             ...            ...
	12       2R [3646319, 3646794)      * |     3646442  37.000000      23.176910        7.455710      21.063850
	13       2R [3666770, 3668041)      * |     3667119 215.000000     279.229060       41.551060     276.108090
	14       2R [3668231, 3668441)      * |     3668363  22.000000       9.943110        4.476950       7.976070
	15       2R [3670063, 3670393)      * |     3670180  38.000000      19.474590        5.901360      17.393440
	16       2R [3670470, 3670927)      * |     3670719 227.000000     305.243350       45.831180     301.974760

### Track

Import ChIP-seq reads from bed files and create a track with the normalized signal:

```go
  fmt.Fprintf(os.Stderr, "Parsing reads (treatment) ...\n")
  treatment1 := GRanges{}
  treatment1.ReadBed6("SRR094207.bed")
  treatment2 := GRanges{}
  treatment2.ReadBed6("SRR094208.bed")
  fmt.Fprintf(os.Stderr, "Parsing reads (control)   ...\n")
  control1   := GRanges{}
  control1.ReadBed6("SRR094215.bed")
  control2   := GRanges{}
  control2.ReadBed6("SRR094216.bed")

  genome  := Genome{}
  genome.ReadFile("Data/hg19.genome")
  d       := 200 // d=200 (see *_peaks.xls)
  binsize := 100 // binsize of the track
  pcounts := 1   // pseudocounts
  track := NormalizedTrack("H3K4me3",
    []GRanges{treatment1, treatment2}, []GRanges{control1, control2},
    genome, d, binsize, pcounts, pcounts, false)
```

Export track to wig or bigWig:
```go
  track.WriteWiggle("track.wig", "track description", true)
  track.WriteBigWig("track.bw",  "track description", true)
```
