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

  for block := range reader.ReadSingleEnd() {
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
    if block.Flag.ReverseStrand() {
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

func (granges *GRanges) ReadBamPairedEnd(filename string, args... interface{}) error {
  var reader *BamReader

  if r, err := NewBamReader(filename, args...); err != nil {
    return err
  } else {
    reader = r
  }
  seqnames  := []string{}
  from      := []int{}
  to        := []int{}
  strand    := []byte{}
  // meta data
  sequence1 := []string{}
  sequence2 := []string{}
  mapq1     := []int{}
  mapq2     := []int{}
  cigar1    := []string{}
  cigar2    := []string{}
  flag1     := []int{}
  flag2     := []int{}

  for r := range reader.ReadPairedEnd() {
    if r.Error != nil {
      return r.Error
    }
    block1 := &r.Block1
    block2 := &r.Block2
    if block1.RefID == -1 || block1.RefID != block2.RefID {
      continue
    }
    if block1.RefID < 0 || int(block1.RefID) > len(reader.Genome.Seqnames) {
      return fmt.Errorf("bam file `%s' contains invalid RefID `%d'", filename, block1.RefID)
    }
    seqnames  = append(seqnames,  reader.Genome.Seqnames[block1.RefID])
    from      = append(from,      int(block1.Position))
    to        = append(to,        int(block2.Position + block2.LSeq))
    strand    = append(strand,    '*')
    // parse meta data
    mapq1     = append(mapq1,     int(block1.MapQ))
    mapq2     = append(mapq2,     int(block2.MapQ))
    flag1     = append(flag1,     int(block1.Flag))
    flag2     = append(flag2,     int(block2.Flag))
    sequence1 = append(sequence1, block1.Seq.String())
    sequence2 = append(sequence2, block2.Seq.String())
    cigar1    = append(cigar1,    block1.Cigar.String())
    cigar2    = append(cigar2,    block2.Cigar.String())
  }
  *granges = NewGRanges(seqnames, from, to, strand)
  granges.AddMeta("sequence1", sequence1)
  granges.AddMeta("sequence2", sequence2)
  granges.AddMeta("flag1", flag1)
  granges.AddMeta("flag2", flag2)
  granges.AddMeta("mapq1", mapq1)
  granges.AddMeta("mapq2", mapq2)
  granges.AddMeta("cigar1", cigar1)
  granges.AddMeta("cigar2", cigar2)
  
  return nil
}
