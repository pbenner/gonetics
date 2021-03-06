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
import "math/rand"

/* -------------------------------------------------------------------------- */

type GenomeRng struct {
  Weights []float64
  Genome  Genome
  MaxLen  int
}

func NewGenomeRng(genome Genome) GenomeRng {
  maxLen := 0
  for i := 0; i < genome.Length(); i++ {
    if maxLen < genome.Lengths[i] {
      maxLen = genome.Lengths[i]
    }
  }
  // each chromosome is weighted by its length
  weights := make([]float64, genome.Length())
  sum     := 0.0
  for i := 0; i < genome.Length(); i++ {
    weights[i] = float64(genome.Lengths[i])
    sum       += weights[i]
  }
  // compute cumulative probabilities
  weights[0] /= sum
  for i := 1; i < genome.Length(); i++ {
    weights[i] = weights[i-1] + weights[i]/sum
  }
  return GenomeRng{weights, genome, maxLen}
}

func (rng GenomeRng) Draw(wsize int) (int, int) {
  if wsize > rng.MaxLen {
    return -1, 0
  }
  k := 0
  // draw a chromosome
  t := 0.0
  for {
    p := rand.Float64()
    for i := 0; i < len(rng.Weights); i++ {
      if t <= p && p < rng.Weights[i] {
        k = i; break
      }
      t = rng.Weights[i]
    }
    if rng.Genome.Lengths[k] - wsize >= 0 {
      break
    }
  }
  i := rand.Intn(rng.Genome.Lengths[k] - wsize + 1)

  return k, i
}

/* -------------------------------------------------------------------------- */

// Generate a GRanges object where each genomic range is generated at random.
// Each chromosome is weighted by its length.
func RandomGRanges(n, wsize int, genome Genome, useStrand bool) GRanges {
  rng := NewGenomeRng(genome)
  // allocate data for the granges object
  seqnames := make([]string, n)
  from     := make([]int, n)
  to       := make([]int, n)
  strand   := []byte{}
  if useStrand {
    strand = make([]byte, n)
  }
  for i := 0; i < n; i++ {
    j, position := rng.Draw(wsize)
    if j == -1 {
      return NewEmptyGRanges(0)
    }
    seqnames[i] = genome.Seqnames[j]
    from    [i] = position
    to      [i] = position + wsize
    if useStrand {
      k := rand.Intn(2)
      strand[i] = []byte{'+', '-'}[k]
    }
  }
  return NewGRanges(seqnames, from, to, strand)
}

/* -------------------------------------------------------------------------- */

func (r GRanges) RandomPermutation() GRanges {
  idx := rand.Perm(r.Length())
  return r.Subset(idx)
}
