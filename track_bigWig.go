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
import "encoding/binary"
import "math"
import "strings"

/* -------------------------------------------------------------------------- */

func (track *Track) parseBlock(buffer []byte, genome Genome) error {
  header := BigWigDataHeader{}
  header.ReadBuffer(buffer)
  // crop header from buffer
  buffer = buffer[24:]

  if idx := int(header.ChromId); idx < 0 || idx > genome.Length() {
    return fmt.Errorf("invalid chromosome id")
  }
  seq := track.Data[genome.Seqnames[int(header.ChromId)]]

  switch header.Type {
  default:
    return fmt.Errorf("unsupported block type")
  case 3:
    if len(seq) != 4*len(buffer) {
      return fmt.Errorf("data block has invalid length")
    }
    for i := 0; i < len(buffer); i += 4 {
      value1  := binary.LittleEndian.Uint32(buffer[i:i+4])
      value2  := math.Float32frombits(value1)
      seq[i/4] = float64(value2)
    }
  }

  return nil
}

func (track *Track) parseBWIndex(bwf *BigWigFile, vertex *RTreeVertex, genome Genome) error {

  if vertex.IsLeaf != 0 {
    for i := 0; i < int(vertex.NChildren); i++ {
      if block, err := vertex.GetBlock(bwf.Fptr, bwf.Header, i); err != nil {
        return err
      } else {
        track.parseBlock(block, genome)
      }
    }
  } else {
    for i := 0; i < int(vertex.NChildren); i++ {
      if err := track.parseBWIndex(bwf, &vertex.Children[i], genome); err != nil {
        return err
      }
    }
  }
  return nil
}

func (track *Track) ReadBigWig(filename, description string, binsize int) error {

  bwf := new(BigWigFile)
  if err := bwf.Open(filename); err != nil {
    return err
  }
  defer bwf.Close()

  seqnames := make([]string, len(bwf.ChromData.Keys))
  lengths  := make([]int,    len(bwf.ChromData.Keys))

  for i := 0; i < len(bwf.ChromData.Keys); i++ {
    if len(bwf.ChromData.Values[i]) != 8 {
      return fmt.Errorf("reading `%s' failed: invalid chromosome list", filename)
    }
    idx := int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][0:4]))
    if idx >= len(bwf.ChromData.Keys) {
      return fmt.Errorf("reading `%s' failed: invalid chromosome index", filename)
    }
    seqnames[idx] = strings.TrimRight(string(bwf.ChromData.Keys[i]), "\x00")
    lengths [idx] = int(binary.LittleEndian.Uint32(bwf.ChromData.Values[i][4:8]))
  }
  genome := NewGenome(seqnames, lengths)

  *track = AllocTrack(description, genome, binsize)

  if err := track.parseBWIndex(bwf, &bwf.Index.Root, genome); err != nil {
    return fmt.Errorf("reading `%s' failed: %v", filename, err)
  }
  return nil
}
