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

func TestPromoters1(t *testing.T) {

  genes, err1 := ReadUCSCGenes("Data/hg19.knownGene.txt.gz")
  if err1 != nil {
    t.Error("TestGenes1 failed!")
  }

  promoters, err2 := Promoters(genes, 500, 500)
  if err2 != nil {
    t.Error("TestGenes1 failed!")
  }
  if promoters.Length() != 51384 {
    t.Error("TestGenes1 failed!")
  }
}
