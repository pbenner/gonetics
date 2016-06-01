## Gonetics

Documentation is available [here](https://godoc.org/github.com/pbenner/gonetics).

### GRanges

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

