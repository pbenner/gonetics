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

import "fmt"

/* -------------------------------------------------------------------------- */

func (granges *GRanges) ReadBam(filename string, args... interface{}) error {
  var reader *BamReader

  if r, err := NewBamReader(filename, args...); err != nil {
    return err
  } else {
    reader = r
  }
  seqnames := []string{}
  from     := []int{}
  to       := []int{}
  strand   := []byte{}
  sequence := []string{}
  mapq     := []int{}
  cigar    := []string{}
  flag     := []int{}

  for block := range reader.ReadBlocks() {
    if block.Error != nil {
      return block.Error
    }
    if block.RefID == -1 {
      continue
    }
    if block.Flag.Unmapped() {
      // read is unmapped, skip...
      continue
    }
    if block.RefID < 0 || int(block.RefID) > len(reader.Genome.Seqnames) {
      return fmt.Errorf("bam file `%s' contains invalid RefID `%d'", filename, block.RefID)
    }
    seqnames = append(seqnames, reader.Genome.Seqnames[block.RefID])
    from     = append(from,     int(block.Position))
    to       = append(to,       int(block.Position + block.LSeq))
    // check if read is unmapped
    if block.Flag.Revcomp() {
      strand = append(strand,   '-')
    } else {
      strand = append(strand,   '+')
    }
    mapq     = append(mapq,     int(block.MapQ))
    flag     = append(flag,     int(block.Flag))
    sequence = append(sequence, block.Seq.String())
    cigar    = append(cigar,    block.Cigar.String())
  }
  *granges = NewGRanges(seqnames, from, to, strand)
  granges.AddMeta("sequence", sequence)
  granges.AddMeta("flag", flag)
  granges.AddMeta("mapq", mapq)
  granges.AddMeta("cigar", cigar)
  
  return nil
}
