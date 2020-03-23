$ cat test1.bed 
chr1    9720000 9745600 S1      0       .
chr1    9745600 9747200 S2      0       .
chr1    9747200 9747400 S2      0       .
chr1    9747400 9749000 S1      0       .
$ cat test2.bed 
chr1    9720000 9745400 S1      0       .
chr1    9745400 9747200 S2      0       .
chr1    9747200 9747400 S2      0       .
chr1    9747400 9749000 S1      0       .
$ cat test3.bed 
chr1    9720000 9747200 S2      0       .
chr1    9747200 9749000 S1      0       .
$ ./segmentationDifferential S1 test{1,2,3}.bed
seqnames    from      to  name labels
    chr1 9720000 9745600    S1    0,1
    chr1 9747200 9749000    S1  0,1,2
