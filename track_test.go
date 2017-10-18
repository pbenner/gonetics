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

func TestTrack2(t *testing.T) {

  genome := Genome{}
  genome.Import("Data/hg19.genome")
  track  := AllocSimpleTrack("test", genome, 100)

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

  GenericMutableTrack{track}.AddReads(reads, 0, true)

  seq, _ := track.GetSequence("chr1")
  
  if v := seq.At(  0); math.Abs(v - 0.35) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v := seq.At( 99); math.Abs(v - 0.35) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v := seq.At(100); math.Abs(v - 1.27) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v := seq.At(199); math.Abs(v - 1.27) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v := seq.At(200); math.Abs(v - 1.17) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }
  if v := seq.At(299); math.Abs(v - 1.17) > 1e-8 {
    t.Error("TestTrack2 failed!")
  }

}

func TestTrack3(t *testing.T) {

  filename := "track_test.1.wig"

  genome := NewGenome([]string{"test1", "test2"}, []int{100, 200})
  track1 := AllocSimpleTrack("Test Track", genome, 10)
  track1.Data["test1"] = []float64{0.1,1.2,2.3,3.4,4.5,5.6,math.NaN(),math.NaN(),8.9,9.0}
  track1.Data["test2"] = []float64{0.1,1.2,2.3,3.4,4.5,5.6,math.NaN(),math.NaN(),8.9,9.0,0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.0}

  track1.WriteWiggle(filename, "test description")

  track2 := AllocSimpleTrack("", genome, 10)
  // initialize all values to NaN
  GenericMutableTrack{track2}.Map(track2, func(name string, i int, x float64) float64 { return math.NaN() })

  if err := track2.ImportWiggle(filename); err != nil {
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
    for i := 0; i < len(seq1); i++ {
      if math.IsNaN(seq1[i]) != math.IsNaN(seq2[i]) {
        t.Errorf("test failed for sequence `%s' at position `%d'", name, i)
      }
      if !math.IsNaN(seq1[i]) && seq1[i] != seq2[i] {
        t.Errorf("test failed for sequence `%s' at position `%d'", name, i)
      }
    }
  }
}

func TestTrack4(t *testing.T) {

  track, _ := NewSimpleTrack("",
    []string{"chr1"},
    [][]float64{{10, 10, 10, 20, 20, 10, 5, 5, 2}},
    100)
  result := []float64{10.0, 10.0, 10.0, 20.0, 20.0, 15.0, 20.0/3.0, 5.5, 5.5}

  GenericMutableTrack{track}.Smoothen(20, []int{1,2,3,4})

  seq := track.Data["chr1"]

  for i := 0; i < len(seq); i++ {
    if math.Abs(seq[i] - result[i]) > 1e-8 {
      t.Error("TestTrack4 failed")
    }
  }
}

func TestTrack5(t *testing.T) {

  track, _ := NewSimpleTrack("",
    []string{"chr1"},
    [][]float64{{10, 2, 5, 2, 2, 1, 5, 5, 2}},
    100)
  result := []float64{4.2, 4.2, 4.2, 2.4, 3.0, 3.0, 3.0, 3.0, 3.0}

  GenericMutableTrack{track}.Smoothen(20, []int{1,2,3,4,5})

  seq := track.Data["chr1"]

  for i := 0; i < len(seq); i++ {
    if math.Abs(seq[i] - result[i]) > 1e-8 {
      t.Error("TestTrack5 failed")
    }
  }
}

func TestTrack6(t *testing.T) {

  track, _ := NewSimpleTrack("",
    []string{"chr1"},
    [][]float64{{1, 2, 5, 2, 2, 1, 1, 1, 2}},
    100)
  result := []float64{17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0, 17.0/9.0}

  GenericMutableTrack{track}.Smoothen(20, []int{1,2,3,4,5,6,7,8,9,10})

  seq := track.Data["chr1"]

  for i := 0; i < len(seq); i++ {
    if math.Abs(seq[i] - result[i]) > 1e-8 {
      t.Error("TestTrack6 failed")
    }
  }
}

func TestTrack7(t *testing.T) {

  track := SimpleTrack{}

  genome := Genome{}
  genome.Import("track_test.1.genome")

  if err := track.ReadBigWig("track_test.1.bw", "", BinMean, 10, 0, math.NaN()); err != nil {
    t.Error(err)
  }
  if err := track.WriteBigWig("track_test.1.bw.tmp", genome); err != nil {
    t.Error(err)
  }
  // reset track and read file again
  track = SimpleTrack{}
  if err := track.ReadBigWig("track_test.1.bw.tmp", "", BinMean, 10, 0, math.NaN()); err != nil {
    t.Error(err)
  }
}

func TestTrack8(t *testing.T) {

  track1 := SimpleTrack{}
  track2 := SimpleTrack{}

  genome := Genome{}
  genome.Import("track_test.3.genome")

  if err := track1.ReadBigWig("track_test.3.bw", "", BinMean, 10, 0, math.NaN()); err != nil {
    t.Error(err)
  }
  if err := track1.WriteBigWig("track_test.3.bw.tmp", genome); err != nil {
    t.Error(err)
  }
  if err := track2.ReadBigWig("track_test.3.bw.tmp", "", BinMean, 10, 0, math.NaN()); err != nil {
    t.Error(err)
  }
}

func TestTrack9(t *testing.T) {

  filename := "track_test.4.bw"

  genome := NewGenome([]string{"test1", "test2"}, []int{100, 200})
  track1 := AllocSimpleTrack("Test Track", genome, 10)
  track1.Data["test1"] = []float64{0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.0}
  track1.Data["test2"] = []float64{math.NaN(),1.2,2.3,3.4,4.5,5.6,math.NaN(),math.NaN(),8.9,9.0,0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,math.NaN()}

  track1.WriteBigWig(filename, genome)

  track2 := AllocSimpleTrack("", genome, 10)

  err1 := track2.ReadBigWig(filename, "", BinMean, 10, 0, math.NaN())
  if err1 != nil {
    t.Error(err1)
  }
  for name, seq1 := range track1.Data {
    seq2 := track2.Data[name]
    if seq2 == nil {
      t.Error("test failed")
    }
    for i := 0; i < len(seq1); i++ {
      if math.IsNaN(seq1[i]) != math.IsNaN(seq2[i]) {
        t.Errorf("test failed for sequence `%s' at position `%d'", name, i)
      }
      if !math.IsNaN(seq1[i]) && math.Abs(seq1[i] - seq2[i]) > 1e-4 {
        t.Errorf("test failed for sequence `%s' at position `%d'", name, i)
      }
    }
  }
}

func TestTrack10(t *testing.T) {

  filename := "track_test.1.wig"

  genome := NewGenome([]string{"test1", "test2"}, []int{100, 200})
  track1 := AllocSimpleTrack("Test Track", genome, 10)
  track1.Data["test1"] = []float64{0.0,0.0,0.0,0.0,4.5,5.6,0.0,7.8,8.9,0.0}
  track1.Data["test2"] = []float64{0.1,1.2,2.3,3.4,4.5,5.6,0,0,8.9,9.0,0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.0}

  track1.WriteWiggle(filename, "test description")

  track2 := AllocSimpleTrack("", genome, 10)
  if err := track2.ImportWiggle(filename); err != nil {
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
    for i := 0; i < len(seq1); i++ {
      if math.IsNaN(seq1[i]) && math.IsNaN(seq2[i]) {
        continue
      }
      if seq1[i] != seq2[i] {
        t.Error("TestTrack10 failed")
      }
    }
  }
}

func TestTrack11(t *testing.T) {

  filename := "track_test.4.bw"

  binSize    := 10
  binOverlap := 1

  genome := NewGenome([]string{"test1", "test2"}, []int{100, 200})
  track1 := AllocSimpleTrack("Test Track", genome, binSize)
  track1.Data["test1"] = []float64{0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.0}
  track1.Data["test2"] = []float64{math.NaN(),1.2,2.3,3.4,4.5,5.6,math.NaN(),math.NaN(),8.9,9.0,0.1,1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,math.NaN()}

  track1.WriteBigWig(filename, genome)

  track2 := AllocSimpleTrack("", genome, binSize)

  err1 := track2.ReadBigWig(filename, "", BinMax, binSize, binOverlap, math.NaN())
  if err1 != nil {
    t.Error(err1)
  }
  for name, seq1 := range track1.Data {
    seq2 := track2.Data[name]
    if seq2 == nil {
      t.Error("test failed")
    }
    for i := 0; i < len(seq1); i++ {
      value := math.NaN()
      // find maximum value
      for j := i-binOverlap; j <= i+binOverlap; j++ {
        if j < 0 || j >= len(seq1) {
          continue
        }
        if math.IsNaN(value) || value < seq1[j] {
          value = seq1[j]
        }
      }
      if math.IsNaN(value) != math.IsNaN(seq2[i]) {
        t.Errorf("test failed for sequence `%s' at position `%d'", name, i)
      }
      if !math.IsNaN(value) && math.Abs(value - seq2[i]) > 1e-4 {
        t.Errorf("test failed for sequence `%s' at position `%d'", name, i)
      }
    }
  }
}
