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

package bioinf

/* -------------------------------------------------------------------------- */

//import   "fmt"
import   "testing"

/* -------------------------------------------------------------------------- */

func TestGenes1(t *testing.T) {

  genes := ReadUCSCGenes("../Data/hg19.knownGene.txt.gz")
  //genes := ReadUCSCGenes("../Data/hg19.ensGene.txt.gz")

  if genes.Length() != 82960 {
    t.Error("TestGenes1 failed!")
  }
}

// func TestGenes2(t *testing.T) {

//   genes := ImportGenesFromUCSC("hg19", "ensGene")
//   fmt.Println(genes)

// }
