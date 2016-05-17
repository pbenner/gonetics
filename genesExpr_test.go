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
import   "math"
import   "testing"

/* -------------------------------------------------------------------------- */

func TestGenesExpr1(t *testing.T) {

  genes := ReadUCSCGenes("../Data/hg19.ensGene.txt.gz")
  genes  = ReadGTF("genesExpr_test.gtf.gz", "transcript_id", "FPKM", genes, false)

  if len(genes.GetMetaFloat("expr")) != 204940 {
    t.Error("TestGenesExpr1 failed!")
  }

  if math.Abs(genes.GetMetaFloat("expr")[4] - 10.413931) > 1e-4 {
    t.Error("TestGenesExpr1 failed!")
  }
}
