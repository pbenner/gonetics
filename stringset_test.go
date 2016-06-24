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

func TestStringSet1(t *testing.T) {

  ss  := EmptyStringSet()
  err := ss.ReadFasta("stringset_test.fa")

  if err != nil {
    t.Error("TestStringSet1 failed")
  }
  if len(ss["chr1"]) != 2011 {
    t.Error("TestStringSet1 failed")
  }
  if len(ss["chr2"]) != 1582 {
    t.Error("TestStringSet1 failed")
  }
}
