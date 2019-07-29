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

func TestKmerGraph1(test *testing.T) {
  kmers, _ := NewKmerCatalogue(2, 6, false, false, true, nil, GappedNucleotideAlphabet{})
  kmers.GetKmerClass("at")
  kmers.GetKmerClass("tc")
  kmers.GetKmerClass("gctc")
  kmers.GetKmerClass("gcta")
  kmers.GetKmerClass("annnc")
  kmers.GetKmerClass("atnnc")
  kmers.GetKmerClass("anntc")
  kmers.GetKmerClass("anctc")
  kmers.GetKmerClass("angtc")
  kmers.GetKmerClass("agctc")
  kmers.GetKmerClass("aagntc")
  kmers.GetKmerClass("aagctc")
  kmers.GetKmerClass("agctca")

  graph := NewKmerGraphFromCatalogue(kmers)

  r := make(map[string][]string)
  // rank 1
  r["tc"] = []string{
    "anntc|gannt" }
  r["at"] = []string{
    "atnnc|gnnat" }
  r["annnc"] = []string{
    "atnnc|gnnat",
    "anntc|gannt" }
  // rank 2
  r["atnnc"] = []string{
    "at|at",
    "annnc|gnnnt" }
  r["anntc"] = []string{
    "ga|tc",
    "anctc|gagnt",
    "angtc|gacnt",
    "annnc|gnnnt" }
  // rank 3
  r["gctc"] = []string{
    "agctc|gagct" }
  r["angtc"] = []string{
    "anntc|gannt" }
  r["anctc"] = []string{
    "agctc|gagct",
    "anntc|gannt" }
  // rank 4
  r["agctc"] = []string{
    "gagc|gctc",
    "anctc|gagnt",
    "aagctc|gagctt",
    "agctca|tgagct" }
  r["aagntc"] = []string{
    "aagctc|gagctt" }
  // rank 5
  r["aagctc"] = []string{
    "agctc|gagct",
    "aagntc|ganctt" }
  r["agctca"] = []string{
    "agctc|gagct" }

  for query, result := range r {
    if s := graph.RelatedKmers(query); len(s) != len(result) {
      test.Errorf("test failed for %v -> %v", query, s)
    } else {
      for i, kmer := range s {
        if result[i] != kmer.String() {
          test.Errorf("test failed for %v -> %v", query, result[i])
        }
      }
    }
  }
}
