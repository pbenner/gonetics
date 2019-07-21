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

func TestKmerCatalogue1(test *testing.T) {
  kmers, _ := NewKmerCatalogue(4, 5, false, false, true, nil, GappedNucleotideAlphabet{})
  j := kmers.GetId  ("anntc")
  i := kmers.GetId  ("anntc")
  s := kmers.GetName("anntc")
  t := kmers.GetName("gannt")
  if i != j {
    test.Error("test failed")
  }
  if s != t || s != "anntc|gannt"{
    test.Error("test failed")
  }
}

func TestKmerCatalogue2(test *testing.T) {
  kmers, _ := NewKmerCatalogue(4, 5, false, false, true, nil, GappedNucleotideAlphabet{})
  j := kmers.GetId  ("gannt")
  i := kmers.GetId  ("anntc")
  s := kmers.GetName("anntc")
  t := kmers.GetName("gannt")
  if i != j {
    test.Error("test failed")
  }
  if s != t || s != "anntc|gannt"{
    test.Error("test failed")
  }
}

func TestKmerCatalogue3(test *testing.T) {
  kmers, _ := NewKmerCatalogue(4, 5, true, true, true, nil, GappedNucleotideAlphabet{})
  j := kmers.GetId  ("tnnag")
  i := kmers.GetId  ("anntc")
  s := kmers.GetName("anntc")
  t := kmers.GetName("gannt")
  if i != j {
    test.Error("test failed")
  }
  if s != t || s != "anntc|tnnag|ctnna|gannt"{
    test.Error("test failed")
  }
}

func TestKmerCatalogue4(test *testing.T) {
  kmers, _ := NewKmerCatalogue(4, 5, true, true, true, nil, GappedNucleotideAlphabet{})
  j := kmers.GetId  ("ctnna")
  i := kmers.GetId  ("anntc")
  s := kmers.GetName("anntc")
  t := kmers.GetName("gannt")
  if i != j {
    test.Error("test failed")
  }
  if s != t || s != "anntc|tnnag|ctnna|gannt"{
    test.Error("test failed")
  }
}
