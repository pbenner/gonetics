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

//import "fmt"

/* -------------------------------------------------------------------------- */

// Container for genes. Tx contains the transcription start and end
// positions. Cds specifies the coding regions.
type Genes struct {
  GRanges
  // pointer to some meta columns
  Names []string
  Cds   []Range
  index map[string]int
}

/* constructors
 * -------------------------------------------------------------------------- */

func newGenes(granges GRanges) Genes {
  names := granges.GetMetaStr("names")
  if len(names) == 0 {
    panic("NewGenes(): names column is missing")
  }
  cds := granges.GetMetaRange("cds")
  if len(cds) == 0 {
    panic("NewGenes(): cds column is missing")
  }
  index := map[string]int{}
  for i := 0; i < granges.Length(); i++ {
    // check if strand is valid
    if granges.Strand[i] != '+' && granges.Strand[i] != '-' {
      panic("NewGenes(): Invalid strand!")
    }
    index[names[i]] = i
  }
  return Genes{granges, names, cds, index}
}

func NewGenes(names, seqnames []string, txFrom, txTo, cdsFrom, cdsTo []int, strand []byte) Genes {
  granges := NewGRanges(seqnames, txFrom, txTo, strand)
  // construct cds ranges
  cds := make([]Range, granges.Length())
  for i := 0; i < granges.Length(); i++ {
    cds[i] = NewRange(cdsFrom[i], cdsTo[i])
  }
  granges.AddMeta("names", names)
  granges.AddMeta("cds",   cds)
  return newGenes(granges)
}

func (g *Genes) Clone() Genes {
  return newGenes(g.GRanges.Clone())
}

/* -------------------------------------------------------------------------- */

func (obj Genes) Remove(indices []int) Genes {
  r := obj.GRanges.Remove(indices)
  return newGenes(r)
}

func (obj Genes) RemoveOverlapsWith(subject Genes) Genes {
  r := obj.GRanges.RemoveOverlapsWith(subject.GRanges)
  return newGenes(r)
}

func (obj Genes) KeepOverlapsWith(subject Genes) Genes {
  r := obj.GRanges.KeepOverlapsWith(subject.GRanges)
  return newGenes(r)
}

func (obj Genes) Subset(indices []int) Genes {
  r := obj.GRanges.Subset(indices)
  return newGenes(r)
}

func (obj Genes) Slice(ifrom, ito int) Genes {
  r := obj.GRanges.Slice(ifrom, ito)
  return newGenes(r)
}

func (obj Genes) Intersection(subject Genes) Genes {
  r := obj.GRanges.Intersection(subject.GRanges)
  return newGenes(r)
}

func (obj Genes) Sort(name string, reverse bool) (Genes, error) {
  if r, err := obj.GRanges.Sort(name, reverse); err != nil {
    return Genes{}, err
  } else {
    return newGenes(r), nil
  }
}

func (obj Genes) FilterGenome(genome Genome) Genes {
  r := obj.GRanges.FilterGenome(genome)
  return newGenes(r)
}

/* -------------------------------------------------------------------------- */

// Returns the index of a gene.
func (g Genes) FindGene(name string) (int, bool) {
  i, ok := g.index[name]
  return i, ok
}

/* convert to string
 * -------------------------------------------------------------------------- */

func (genes Genes) String() string {
  return genes.PrintPretty(10)
}
