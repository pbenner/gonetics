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
import   "math"
import   "testing"

/* -------------------------------------------------------------------------- */

func TestTrack1(t *testing.T) {

  genome := ReadGenome("Data/hg19.genome")
  track  := NewTrack("test", genome, 100)

  track.Set("chrX", 101, 13.0)
  track.Add("chrX", 101, 10.0)

  if v, _ := track.At("chrX", 180); v != 23 {
    t.Error("TestTrack1 failed!")
  }

}

func TestTrack2(t *testing.T) {

  genome := ReadGenome("Data/hg19.genome")
  track  := NewTrack("test", genome, 100)

  // bin                                         bin                                         bin                                         bin
  // 00:   010 020 030 040 050 060 070 080 090   01:   110 120 130 140 150 160 170 180 190   02:   210 220 230 240 250 260 270 280 290   03:
  // |||+---+---+---+---+---+---+---+---+---+---+|||+---+---+---+---+---+---+---+---+---+---+|||+---+---+---+---+---+---+---+---+---+---+|||
  //      (nucleotides #99 and #100) 99-100+1 = 2                                         100                                          31                     
  // r1:                                       |-------------------------------------------------------------|
  //                                                                             200-174+1=27                                          86
  // r2:                                                                          |-----------------------------------------------|
  //                                           33
  // r3:|------------|
  seqnames := []string{"chr1", "chr1", "chr1"}
  from     := []int { 99, 174, 01}
  to       := []int {231, 286, 33}
  strand   := []byte{'+', '+', '+'}
  reads    := NewGRanges(seqnames, from, to, strand)

  track.AddReads(reads, 0)
  
  if v, _ := track.At("chr1",   1); math.Abs(v - 0.35) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 100); math.Abs(v - 0.35) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 101); math.Abs(v - 1.27) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 200); math.Abs(v - 1.27) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 201); math.Abs(v - 1.17) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 300); math.Abs(v - 1.17) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }

}
