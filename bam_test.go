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

func TestBam1(t *testing.T) {

  granges := GRanges{}
  if err := granges.ImportBamSingleEnd("bam_test.1.bam"); err != nil {
    t.Error(err)
  }
  if granges.Length() != 12 {
    t.Error("TestBam1 failed")
  }
  if granges.GetMetaStr("cigar")[3] != "6M14N1I5M" {
    t.Error("TestBam1 failed")
  }
}

func TestBam2(t *testing.T) {

  granges := GRanges{}
  if err := granges.ImportBamPairedEnd("bam_test.2.bam"); err != nil {
    t.Error(err)
  }
  if granges.Length() != 2335 {
    t.Error("TestBam2 failed")
  }
}
