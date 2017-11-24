
## Apply a function to a set of tracks

### Example: Combine tracks for differential analysis

Define mapping `f.go`:

```go
package main

import "math"

func F(seqname string, position int, values ...float64) float64 {

  x1 := values[0]
    x2 := values[1]

  return math.Abs(x1-x2)*(1.0+x1)/(1.0+x2)

}
```

Compile with `go build -buildmode=plugin f.go`. Execute with `bigWigMap -v --bin-size=200 f.so result.bw input1.bw input2.bw`
