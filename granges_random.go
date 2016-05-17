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

package bioinf

/* -------------------------------------------------------------------------- */

import "math/rand"

/* -------------------------------------------------------------------------- */

// Generate a GRanges object where each genomic range is generated at random.
// Genomic bounds are not checked. (This is required by the GRanges
// ImportTrack() method!)
func RandomGRanges(n, wsize int, genome Genome, useStrand bool) GRanges {

  seqnames := make([]string, n)
  from     := make([]int, n)
  to       := make([]int, n)
  strand   := []byte{}
  if useStrand {
    strand = make([]byte, n)
  }
  n_seq := len(genome.Seqnames)

  for i := 0; i < n; i++ {
    j := rand.Intn(n_seq)
    if genome.Lengths[j]-wsize < 0 {
      panic("window size is too large")
    }
    position := rand.Intn(genome.Lengths[j]-wsize)
    seqnames[i] = genome.Seqnames[j]
    from[i]     = position
    to[i]       = position + wsize
    if useStrand {
      k := rand.Intn(2)
      strand[i] = []byte{'+', '-'}[k]
    }
  }
  return NewGRanges(seqnames, from, to, strand)
}
