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
//import   "math"
import   "testing"

/* -------------------------------------------------------------------------- */

func TestGenesExpr1(t *testing.T) {

//  genes := ImportGenesFromUCSC("hg19", "ensGene")
//  genes.WriteTable("Data/hg19.ensGene.txt.gz", false, true)

  genes, err1 := ReadUCSCGenes("Data/hg19.ensGene.txt.gz")
  if err1 != nil {
    t.Error("TestGenesExpr1 failed!")
  }

  err2 := genes.ReadGTFExpr("genesExpr_test.gtf.gz", "transcript_id", "FPKM")
  if err2 != nil {
    t.Error("TestGenesExpr1 failed!")
  }

  if len(genes.GetMetaFloat("expr")) != 204940 {
    t.Error("TestGenesExpr1 failed!")
  }
}
