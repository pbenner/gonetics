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

func TestKmersIterator1(t *testing.T) {
  c := "cgta"
  r := []string{
    "acgtaa",
    "acgtac",
    "acgtag",
    "acgtat",
    "ccgtaa",
    "ccgtac",
    "ccgtag",
    "ccgtat",
    "gcgtaa",
    "gcgtac",
    "gcgtag",
    "gcgtat",
    "tcgtaa",
    "tcgtac",
    "tcgtag",
    "tcgtat" }

  i := 0
  for it := NewKmersCylinderIterator(6, 2, NucleotideAlphabet{}, 1, c); it.Ok(); it.Next() {
    if r[i] != it.Get() {
      t.Error("test failed")
    }
    i++
  }
}

func TestKmersIterator2(test *testing.T) {
  r := []string{
    "ctaaa",
    "ctaca",
    "ctaga",
    "ctata",
    "ctcaa",
    "ctcca",
    "ctcga",
    "ctcta",
    "ctgaa",
    "ctgca",
    "ctgga",
    "ctgta",
    "cttaa",
    "cttca",
    "cttga",
    "cttta" }

  it := NewKmersInstantiationIterator(GappedNucleotideAlphabet{}, "ctnna", false)
  for i := 0; it.Ok(); it.Next() {
    if it.Get() != r[i] {
      test.Error("test failed")
    }
    i++
  }
}

func TestKmersIterator3(test *testing.T) {
  r := []string{
    "ctnna",
    "ctnaa",
    "ctnca",
    "ctnga",
    "ctnta",
    "ctana",
    "ctaaa",
    "ctaca",
    "ctaga",
    "ctata",
    "ctcna",
    "ctcaa",
    "ctcca",
    "ctcga",
    "ctcta",
    "ctgna",
    "ctgaa",
    "ctgca",
    "ctgga",
    "ctgta",
    "cttna",
    "cttaa",
    "cttca",
    "cttga",
    "cttta" }

  it := NewKmersInstantiationIterator(GappedNucleotideAlphabet{}, "ctnna", true)
  for i := 0; it.Ok(); it.Next() {
    if it.Get() != r[i] {
      test.Error("test failed")
    }
    i++
  }
}
