/* Copyright (C) 2016 Philipp Benner
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package gonetics

/* -------------------------------------------------------------------------- */

import "sort"

/* -------------------------------------------------------------------------- */

// Container for genes. Tx contains the transcription start and end
// positions. Cds specifies the coding regions.
type Genes struct {
  Names      []string
  Seqnames   []string
  Tx         []Range
  Cds        []Range
  Strand     []byte
  index      map[string]int
  Meta
}

/* constructors
 * -------------------------------------------------------------------------- */

func newGenes(names, seqnames []string, tx, cds []Range, strand []byte) Genes {
  n := len(names)
  if len(tx)       != n || len(cds)    != n  ||
    (len(seqnames) != n || len(strand) != n) {
    panic("NewGenes(): invalid arguments!")
  }
  index := map[string]int{}
  for i := 0; i < n; i++ {
    // check if strand is valid
    if strand[i] != '+' && strand[i] != '-' {
      panic("NewGenes(): Invalid strand!")
    }
    index[names[i]] = i
  }
  return Genes{names, seqnames, tx, cds, strand, index, Meta{}}
}

func NewGenes(names, seqnames []string, txFrom, txTo, cdsFrom, cdsTo []int, strand []byte) Genes {
  n := len(names)
  if len(txFrom)   != n || len(txTo)   != n  ||
    (len(cdsFrom)  != n || len(cdsTo)  != n) ||
    (len(seqnames) != n || len(strand) != n) {
    panic("NewGenes(): invalid arguments!")
  }
  tx  := make([]Range, n)
  cds := make([]Range, n)
  for i := 0; i < n; i++ {
    // create tx and cds ranges
    tx [i] = NewRange( txFrom[i],  txTo[i])
    cds[i] = NewRange(cdsFrom[i], cdsTo[i])
  }
  return newGenes(names, seqnames, tx, cds, strand)
}

func (g *Genes) Clone() Genes {
  n := g.Length()
  names    := make([]string, n)
  seqnames := make([]string, n)
  cds      := make([]Range, n)
  tx       := make([]Range, n)
  strand   := make([]byte, n)
  copy(names,    g.Names)
  copy(seqnames, g.Seqnames)
  copy(cds,      g.Cds)
  copy(tx,       g.Tx)
  copy(strand,   g.Strand)
  result := newGenes(names, seqnames, cds, tx, strand)
  result.Meta = g.Meta.Clone()
  return result
}

/* -------------------------------------------------------------------------- */

// Number of genes in this object.
func (g Genes) Length() int {
  return len(g.Names)
}

// Returns the index of a gene.
func (g Genes) FindGene(name string) (int, bool) {
  i, ok := g.index[name]
  return i, ok
}

func (g1 Genes) Append(g2 Genes) Genes {
  names    := append(g1.Names,    g2.Names...)
  seqnames := append(g1.Seqnames, g2.Seqnames...)
  cds      := append(g1.Cds,      g2.Cds...)
  tx       := append(g1.Tx,       g2.Tx...)
  strand   := append(g1.Strand,   g2.Strand...)

  result := newGenes(names, seqnames, cds, tx, strand)
  result.Meta = g1.Meta.Append(g2.Meta)

  return result
}

func (g Genes) Remove(indices []int) Genes {
  if len(indices) == 0 {
    return g.Clone()
  }
  indices = removeDuplicatesInt(indices)
  sort.Ints(indices)

  n := g.Length()
  m := n - len(indices)
  // convert indices to subset indices
  idx := make([]int, m)
  for i, j, k := 0, 0, 0; i < g.Length(); i++ {
    for k < len(indices)-1 && i > indices[k] {
      k++
    }
    if i != indices[k] {
      idx[j] = i
      j++
    }
  }
  result := g.Subset(idx)
  result.Meta = g.Meta.Subset(idx)

  return result
}

// Returns a subset of the genes given by indices. This function may
// also be used to alter the order of genes.
func (g Genes) Subset(indices []int) Genes {
  n := len(indices)
  names    := make([]string, n)
  seqnames := make([]string, n)
  txFrom   := make([]int, n)
  txTo     := make([]int, n)
  cdsFrom  := make([]int, n)
  cdsTo    := make([]int, n)
  strand   := make([]byte, n)

  for i := 0; i < n; i++ {
    names   [i] = g.Names   [indices[i]]
    seqnames[i] = g.Seqnames[indices[i]]
    txFrom  [i] = g.Tx      [indices[i]].From
    txTo    [i] = g.Tx      [indices[i]].To
    cdsFrom [i] = g.Cds     [indices[i]].From
    cdsTo   [i] = g.Cds     [indices[i]].To
    strand  [i] = g.Strand  [indices[i]]
  }
  result := NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand)
  result.Meta = g.Meta.Subset(indices)

  return result
}

// Returns the set of genes starting from row ifrom till ito (row ito
// is not included).
func (g Genes) Slice(ifrom, ito int) Genes {
  n := ito-ifrom
  names    := make([]string, n)
  seqnames := make([]string, n)
  txFrom   := make([]int, n)
  txTo     := make([]int, n)
  cdsFrom  := make([]int, n)
  cdsTo    := make([]int, n)
  strand   := make([]byte, n)

  for i := ifrom; i < ito; i++ {
    names   [i-ifrom] = g.Names   [i]
    seqnames[i-ifrom] = g.Seqnames[i]
    txFrom  [i-ifrom] = g.Tx      [i].From
    txTo    [i-ifrom] = g.Tx      [i].To
    cdsFrom [i-ifrom] = g.Cds     [i].From
    cdsTo   [i-ifrom] = g.Cds     [i].To
    strand  [i-ifrom] = g.Strand  [i]
  }
  result := NewGenes(names, seqnames, txFrom, txTo, cdsFrom, cdsTo, strand)
  result.Meta = g.Meta.Slice(ifrom, ito)

  return result
}

// Sort rows using data from a meta column. The parameter name gives
// the name of the columns. If reverse is true, rows are sorted in
// descending order.
func (g Genes) Sort(name string, reverse bool) (Genes, error) {
  j, err := g.sortedIndices(name, reverse)
  if err != nil {
    return Genes{}, err
  }
  return g.Subset(j), nil
}

/* convert to GRanges
 * -------------------------------------------------------------------------- */

// Export Genes object to GRanges. The transcription region is exported
// as the standard ranges field. Ranges of coding sequences and gene
// names are exported as meta data.
func (g *Genes) GRanges() *GRanges {
  n := g.Length()
  names   := make([]string, n)
  txFrom  := make([]int, n)
  txTo    := make([]int, n)
  cdsFrom := make([]int, n)
  cdsTo   := make([]int, n)
  for i := 0; i < n; i++ {
    names[i]   = g.Names[i]
    txFrom[i]  = g.Tx[i].From
    txTo[i]    = g.Tx[i].To
    cdsFrom[i] = g.Cds[i].From
    cdsTo[i]   = g.Cds[i].To
  }
  result := NewGRanges(g.Seqnames, txFrom, txTo, g.Strand)
  result.Meta = g.Meta.Clone()
  result.AddMeta("names",   names)
  result.AddMeta("cdsFrom", cdsFrom)
  result.AddMeta("cdsTo",   cdsTo)
  return &result
}

/* convert to string
 * -------------------------------------------------------------------------- */

func (genes Genes) String() string {
  return genes.PrettyPrint(10)
}
