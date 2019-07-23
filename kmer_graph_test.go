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
  kmers, _ := NewKmerCatalogue(4, 6, false, false, true, nil, GappedNucleotideAlphabet{})
  kmers.GetKmerClass("gctc")
  kmers.GetKmerClass("gcta")
  kmers.GetKmerClass("anntc")
  kmers.GetKmerClass("anctc")
  kmers.GetKmerClass("agctc")
  kmers.GetKmerClass("aagctc")
  kmers.GetKmerClass("agctca")

  graph := NewKmerGraphFromCatalogue(kmers)

  r1 := []string{
    "gagc|gctc",
    "anctc|gagnt",
    "anntc|gannt",
    "aagctc|gagctt",
    "agctca|tgagct" }
  r2 := []string{
    "agctc|gagct",
    "anntc|gannt" }

  for i, kmer := range graph.RelatedKmers("agctc") {
    if r1[i] != kmer.String() {
      test.Error("test failed")
    }
  }
  for i, kmer := range graph.RelatedKmers("anctc") {
    if r2[i] != kmer.String() {
      test.Error("test failed")
    }
  }
  if r := graph.RelatedKmers("agatc"); len(r) != 0 {
    test.Error("test failed")
  }
}
