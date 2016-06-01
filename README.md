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
