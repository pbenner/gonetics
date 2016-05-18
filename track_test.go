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

//import   "fmt"
import   "testing"

/* -------------------------------------------------------------------------- */

func TestTrack1(t *testing.T) {

  genome := ReadGenome("Data/hg19.genome")
  track  := NewTrack("test", genome, 10)

  track.Set("chrX", 10, 13.0)
  track.Add("chrX", 10, 10.0)

  if v, _ := track.At("chrX", 18); v != 23 {
    t.Error("TestTrack1 failed!")
  }

}

func TestTrack2(t *testing.T) {

  genome := ReadGenome("Data/hg19.genome")
  track  := NewTrack("test", genome, 10)

  seqnames := []string{"chr1", "chr1"}
  from     := []int { 6, 17}
  to       := []int {23, 28}
  strand   := []byte{'+', '+'}
  reads    := NewGRanges(seqnames, from, to, strand)

  track.AddReads(reads, 0)

  if v, _ := track.At("chr1",  0); v != 0.5 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 10); v != 1.4 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 20); v != 1.1 {
    t.Error("TestTrack2 failed!")
  }

}
