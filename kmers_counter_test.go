/* Copyright (C) 2019 Philipp Benner
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

//import "fmt"
import "testing"

/* -------------------------------------------------------------------------- */

func TestKmersCounter1(test *testing.T) {
  r := []string{
    "acgt|tgca|tgca|acgt",
    "acnt|tgna|tnca|angt",
    "agcg|tcgc|gcga|cgct",
    "agng|tcnc|gnga|cnct",
    "ancg|tngc|gcna|cgnt",
    "anng|tnnc|gnna|cnnt",
    "annt|tnna|tnna|annt",
    "cagc|gtcg|cgac|gctg",
    "canc|gtng|cnac|gntg",
    "cgcg|gcgc|gcgc|cgcg",
    "cgtc|gcag|ctgc|gacg",
    "cgnc|gcng|cngc|gncg",
    "cgng|gcnc|gngc|cncg",
    "ctnc|gang|cntc|gnag",
    "cnnc|gnng|cnnc|gnng",
    "cnng|gnnc|gnnc|cnng",
    "acgtc|tgcag|ctgca|gacgt",
    "acgnc|tgcng|cngca|gncgt",
    "acntc|tgnag|ctnca|gangt",
    "acnnc|tgnng|cnnca|gnngt",
    "agcgc|tcgcg|cgcga|gcgct",
    "agcnc|tcgng|cncga|gngct",
    "agngc|tcncg|cgnga|gcnct",
    "agnnc|tcnng|cnnga|gnnct",
    "ancgc|tngcg|cgcna|gcgnt",
    "ancnc|tngng|cncna|gngnt",
    "angtc|tncag|ctgna|gacnt",
    "angnc|tncng|cngna|gncnt",
    "anngc|tnncg|cgnna|gcnnt",
    "anntc|tnnag|ctnna|gannt",
    "annnc|tnnng|cnnna|gnnnt",
    "cagcg|gtcgc|gcgac|cgctg",
    "cagng|gtcnc|gngac|cnctg",
    "cancg|gtngc|gcnac|cgntg",
    "canng|gtnnc|gnnac|cnntg",
    "cgacg|gctgc|gcagc|cgtcg",
    "cgang|gctnc|gnagc|cntcg",
    "cgcng|gcgnc|gncgc|cngcg",
    "cgtng|gcanc|gntgc|cnacg",
    "cgncg|gcngc|gcngc|cgncg",
    "cgnng|gcnnc|gnngc|cnncg",
    "cnang|gntnc|gnanc|cntng",
    "cncng|gngnc|gncnc|cngng",
    "cnnng|gnnnc|gnnnc|cnnng" }
  s := []int{
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
    1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 2, 1, 1, 2 }

  kmersCounter, _ := NewKmersCounter(4, 5, true, true, true, nil, GappedNucleotideAlphabet{})
  counts := kmersCounter.CountKmers([]byte("acgtcgcg"))
  for i, it := 0, counts.Iterate(); it.Ok(); it.Next() {
    if it.GetName() != r[i] {
      test.Error("test failed")
    }
    if it.GetCount() != s[i] {
      test.Error("test failed")
    }
    i++
  }
}
