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

func TestGRanges1(t *testing.T) {

  granges := NewEmptyGRanges(0)
  granges.ReadBed6("granges_test.bed")
  granges.AddMeta("TestMeta",
    []int{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20})

  //fmt.Println(granges)

  if granges.Length() != 20 {
    t.Error("TestGRanges1 failed!")
  }
  if granges.MetaLength() != 1 {
    t.Error("TestGRanges1 failed!")
  }

}

func TestGRanges2(t *testing.T) {

  granges := NewEmptyGRanges(0)
  granges.ReadBed6("granges_test.bed")
  granges.AddMeta("TestMeta1",
    [][]int{
      { 1, 2},{ 2, 3},{ 3, 4},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 8},{ 8, 9},{ 9,10},{10,11},
      {11,12},{12,13},{13,14},{14,15},{15,16},{16,17},{17,18},{18,19},{19,20},{20,21}})
  granges.AddMeta("TestMeta2",
    []int{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20})

  //fmt.Println(granges)

  if granges.Length() != 20 {
    t.Error("TestGRanges2 failed!")
  }
  if granges.MetaLength() != 2 {
    t.Error("TestGRanges2 failed!")
  }
}

func TestGRanges3(t *testing.T) {
  seqnames := []string{"chr1", "chr1", "chr1"}
  from     := []int{100000266, 100000271, 100000383}
  to       := []int{100000291, 100000296, 100000408}
  strand   := []byte{'+', '+', '-'}

  granges  := NewGRanges(seqnames, from, to, strand)

  if granges.Length() != 3 {
    t.Error("TestGRanges3 failed!")
  }
}

func TestGRangesRandom(t *testing.T) {
  genome  := Genome{}
  genome.ReadFile("Data/hg19.genome")
  granges := RandomGRanges(1000, 10000, genome, true)

  if granges.Length() != 1000 {
    t.Error("TestGRangesRandom failed!")
  }
}
