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

func TestKmersCounter1(t *testing.T) {
  c := []byte("cgta")
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
    if r[i] != string(it.Get()) {
      t.Error("test failed")
    }
    i++
  }
}

