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

	   seqnames             ranges strand |  abs_summit -log10(pvalue) fold_enrichment
	 1       2L [   5757,    6000)      * |        5865      19.809300        6.851880
	 2       2L [  47233,   47440)      * |       47354      19.648200        6.263150
	 3       2L [  66379,   67590)      * |       66957     350.151250       50.986050
	 4       2L [  72305,   72837)      * |       72525     208.558240       34.460930
	 5       2L [  72999,   73217)      * |       73130      12.711700        5.239670
	        ...                ...        |         ...            ...             ...
	12       2R [3646319, 3646793)      * |     3646442      23.176910        7.455710
	13       2R [3666770, 3668040)      * |     3667119     279.229060       41.551060
	14       2R [3668231, 3668440)      * |     3668363       9.943110        4.476950
	15       2R [3670063, 3670392)      * |     3670180      19.474590        5.901360
	16       2R [3670470, 3670926)      * |     3670719     305.243350       45.831180
