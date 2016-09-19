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

import   "fmt"
import   "math"
import   "testing"

/* -------------------------------------------------------------------------- */

func TestTrack1(t *testing.T) {

  genome := Genome{}
  genome.ReadFile("Data/hg19.genome")
  track  := AllocTrack("test", genome, 100)

  track.Set("chrX", 100, 13.0)
  track.Add("chrX", 100, 10.0)

  if v, _ := track.At("chrX", 180); v != 23 {
    t.Error("TestTrack1 failed!")
  }

}

func TestTrack2(t *testing.T) {

  genome := Genome{}
  genome.ReadFile("Data/hg19.genome")
  track  := AllocTrack("test", genome, 100)

  // bin                                         bin                                         bin                                         bin
  // 00:   010 020 030 040 050 060 070 080 090   01:   110 120 130 140 150 160 170 180 190   02:   210 220 230 240 250 260 270 280 290   03:
  // |||+---+---+---+---+---+---+---+---+---+---+|||+---+---+---+---+---+---+---+---+---+---+|||+---+---+---+---+---+---+---+---+---+---+|||
  // r1:  (nucleotides #98 and #199)     100-98=2                                         100                                          31
  //                                           |-------------------------------------------------------------| [98, 231)
  // r2:                                                                           200-173=27                                          86
  //                                                                              |-----------------------------------------------| [173, 286)
  // r3:                                       33
  //    |------------| [0,33)
  seqnames := []string{"chr1", "chr1", "chr1"}
  from     := []int { 98, 173, 00}
  to       := []int {231, 286, 33}
  strand   := []byte{'+', '+', '+'}
  reads    := NewGRanges(seqnames, from, to, strand)

  track.AddReads(reads, 0)
  
  if v, _ := track.At("chr1",   0); math.Abs(v - 0.35) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1",  99); math.Abs(v - 0.35) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 100); math.Abs(v - 1.27) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 199); math.Abs(v - 1.27) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 200); math.Abs(v - 1.17) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v, _ := track.At("chr1", 299); math.Abs(v - 1.17) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }

}

func TestTrack3(t *testing.T) {

  filename := "track_test.1.wig"

  genome := NewGenome([]string{"test1", "test2"}, []int{100, 200})
  track1 := AllocTrack("Test Track", genome, 10)
  track1.Data["test1"] = []float64{0.0,0.0,0.0,0.0,4.5,5.6,0.0,7.8,8.9,0.0}
  track1.Data["test2"] = []float64{0.1,1.2,2.3,3.4,4.5,5.6,math.NaN(),math.NaN(),8.9,9.0,0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.0}

  track1.WriteWiggle(filename, "test description")

  track2 := AllocTrack("", genome, 10)
  track2.Map(func(x float64) float64 { return math.NaN() })
  if err := track2.ReadWiggle(filename); err != nil {
    t.Error(err)
  }
  if track1.Name != track2.Name {
    t.Error("TestTrack3 failed")
  }
  for name, seq1 := range track1.Data {
    seq2 := track2.Data[name]
    if seq2 == nil {
      t.Error("TestTrack3 failed")
    }
    fmt.Println("seq1:",seq1)
    fmt.Println("seq2:",seq2)
    for i := 0; i < len(seq1); i++ {
      if math.IsNaN(seq1[i]) && math.IsNaN(seq2[i]) {
        continue
      }
      if seq1[i] != seq2[i] {
        t.Error("TestTrack3 failed")
      }
    }
  }
}

func TestTrack4(t *testing.T) {

  track := NewTrack("",
    []string{"chr1"},
    [][]float64{{10, 10, 10, 20, 20, 10, 5, 5, 2}},
    100)
  result := []float64{10.0, 10.0, 10.0, 20.0, 20.0, 15.0, 20.0/3.0, 5.5, 5.5}

  track.Smoothen(20, []int{1,2,3,4})

  seq := track.Data["chr1"]

  for i := 0; i < len(seq); i++ {
    if math.Abs(seq[i] - result[i]) > 1e-8 {
      t.Error("TestTrack4 failed")
    }
  }
}

func TestTrack5(t *testing.T) {

  track := NewTrack("",
    []string{"chr1"},
    [][]float64{{10, 2, 5, 2, 2, 1, 5, 5, 2}},
    100)
  result := []float64{4.2, 4.2, 4.2, 2.4, 3.0, 3.0, 3.0, 3.0, 3.0}

  track.Smoothen(20, []int{1,2,3,4,5})

  seq := track.Data["chr1"]

  for i := 0; i < len(seq); i++ {
    if math.Abs(seq[i] - result[i]) > 1e-8 {
      t.Error("TestTrack5 failed")
    }
  }
}

func TestTrack6(t *testing.T) {

  track := NewTrack("",
    []string{"chr1"},
    [][]float64{{1, 2, 5, 2, 2, 1, 1, 1, 2}},
    100)
  result := []float64{17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0}

  track.Smoothen(20, []int{1,2,3,4,5,6,7,8,9,10})

  seq := track.Data["chr1"]

  for i := 0; i < len(seq); i++ {
    if math.Abs(seq[i] - result[i]) > 1e-8 {
      t.Error("TestTrack6 failed")
    }
  }
}

func TestTrack7(t *testing.T) {

  track := Track{}

  if err := track.ReadBigWig("track_test.1.bw", ""); err != nil {
    t.Error(err)
  }
  if err := track.WriteBigWig("track_test.1.bw.tmp", ""); err != nil {
    t.Error(err)
  }
  // reset track and read file again
  track = Track{}
  if err := track.ReadBigWig("track_test.1.bw.tmp", ""); err != nil {
    t.Error(err)
  }
}

func TestTrack8(t *testing.T) {

  track1 := Track{}
  track2 := Track{}

  if err := track1.ReadBigWig("track_test.3.bw", ""); err != nil {
    t.Error(err)
  }
  if err := track1.WriteBigWig("track_test.3.bw.tmp", ""); err != nil {
    t.Error(err)
  }
  if err := track2.ReadBigWig("track_test.3.bw.tmp", ""); err != nil {
    t.Error(err)
  }
}

func TestTrack9(t *testing.T) {

  filename1 := "track_test.4.bw"
  filename2 := "track_test.5.bw"

  genome := NewGenome([]string{"test1", "test2"}, []int{100, 200})
  track1 := AllocTrack("Test Track", genome, 10)
  track1.Data["test1"] = []float64{0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.0}
  track1.Data["test2"] = []float64{math.NaN(),1.2,2.3,3.4,4.5,5.6,math.NaN(),math.NaN(),8.9,9.0,0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,math.NaN()}

  track1.WriteBigWig(filename1, "test description", true)
  track1.WriteBigWig(filename2, "test description", false)

  track2 := AllocTrack("", genome, 10)
  track2.Map(func(x float64) float64 { return math.NaN() })
  track3 := AllocTrack("", genome, 10)
  track3.Map(func(x float64) float64 { return math.NaN() })
  err1 := track2.ReadBigWig(filename1, "")
  err2 := track3.ReadBigWig(filename2, "")
  if err1 != nil {
    t.Error(err1)
  }
  if err2 != nil {
    t.Error(err2)
  }
  for name, seq1 := range track1.Data {
    seq2 := track2.Data[name]
    seq3 := track3.Data[name]
    if seq2 == nil {
      t.Error("TestTrack3 failed")
    }
    if seq3 == nil {
      t.Error("TestTrack3 failed")
    }
    for i := 0; i < len(seq1); i++ {
      if math.IsNaN(seq1[i]) && math.IsNaN(seq2[i]) && math.IsNaN(seq3[i]) {
        continue
      }
      if math.Abs(seq1[i] - seq2[i]) > 1e-4 {
        t.Errorf("TestTrack3 failed for sequence `%s' at position `%d'", name, i)
      }
      if math.Abs(seq1[i] - seq3[i]) > 1e-4 {
        t.Errorf("TestTrack3 failed for sequence `%s' at position `%d'", name, i)
      }
    }
  }
}
