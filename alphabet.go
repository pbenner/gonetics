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

import "fmt"

/* -------------------------------------------------------------------------- */

type Alphabet interface {
  Bases (i byte) ([]byte, error)
  Code  (i byte) (  byte, error)
  Decode(i byte) (  byte, error)
  Length()       int
}

type ComplementableAlphabet interface {
  Bases (i byte) ([]byte, error)
  Code           (i byte) (byte, error)
  Decode         (i byte) (byte, error)
  Length         ()       int
  Complement     (i byte) (byte, error)
  ComplementCoded(i byte) (byte, error)
}

/* -------------------------------------------------------------------------- */

type NucleotideAlphabet struct {
}

func (NucleotideAlphabet) Bases(i byte) ([]byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return []byte{'a'}, nil
  case 'C': fallthrough
  case 'c': return []byte{'c'}, nil
  case 'G': fallthrough
  case 'g': return []byte{'g'}, nil
  case 'T': fallthrough
  case 't': return []byte{'t'}, nil
  default:  return nil, fmt.Errorf("Code(): `%c' is not part of the alphabet", i)
  }
}

func (NucleotideAlphabet) Code(i byte) (byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return 0, nil
  case 'C': fallthrough
  case 'c': return 1, nil
  case 'G': fallthrough
  case 'g': return 2, nil
  case 'T': fallthrough
  case 't': return 3, nil
  default:  return 0xFF, fmt.Errorf("Code(): `%c' is not part of the alphabet", i)
  }
}

func (NucleotideAlphabet) Decode(i byte) (byte, error) {
  switch i {
  case 0:  return 'a', nil
  case 1:  return 'c', nil
  case 2:  return 'g', nil
  case 3:  return 't', nil
  default: return 0xFF, fmt.Errorf("Code(): `%d' is not a code of the alphabet", int(i))
  }
}

func (NucleotideAlphabet) Length() int {
  return 4
}

func (NucleotideAlphabet) ComplementCoded(i byte) (byte, error) {
  switch i {
  case 0:  return 3, nil
  case 1:  return 2, nil
  case 2:  return 1, nil
  case 3:  return 0, nil
  default: return 0xFF, fmt.Errorf("ComplementCoded(): `%d' is not a code of the alphabet", int(i))
  }
}

func (NucleotideAlphabet) Complement(i byte) (byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return 't', nil
  case 'C': fallthrough
  case 'c': return 'g', nil
  case 'G': fallthrough
  case 'g': return 'c', nil
  case 'T': fallthrough
  case 't': return 'a', nil
  default:  return 0xFF, fmt.Errorf("Complement(): `%c' is not part of the alphabet", i)
  }
}

/* -------------------------------------------------------------------------- */

type GappedNucleotideAlphabet struct {
}

func (GappedNucleotideAlphabet) Bases(i byte) ([]byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return []byte{'a'}, nil
  case 'C': fallthrough
  case 'c': return []byte{'c'}, nil
  case 'G': fallthrough
  case 'g': return []byte{'g'}, nil
  case 'T': fallthrough
  case 't': return []byte{'t'}, nil
  case 'N': fallthrough
  case 'n': return []byte{'a', 'c', 'g', 't'}, nil
  default:  return nil, fmt.Errorf("Code(): `%c' is not part of the alphabet", i)
  }
}

func (GappedNucleotideAlphabet) Code(i byte) (byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return 0, nil
  case 'C': fallthrough
  case 'c': return 1, nil
  case 'G': fallthrough
  case 'g': return 2, nil
  case 'T': fallthrough
  case 't': return 3, nil
  case 'N': fallthrough
  case 'n': return 4, nil
  default:  return 0xFF, fmt.Errorf("Code(): `%c' is not part of the alphabet", i)
  }
}

func (GappedNucleotideAlphabet) Decode(i byte) (byte, error) {
  switch i {
  case 0:  return 'a', nil
  case 1:  return 'c', nil
  case 2:  return 'g', nil
  case 3:  return 't', nil
  case 4:  return 'n', nil
  default: return 0xFF, fmt.Errorf("Code(): `%d' is not a code of the alphabet", int(i))
  }
}

func (GappedNucleotideAlphabet) Length() int {
  return 5
}

func (GappedNucleotideAlphabet) ComplementCoded(i byte) (byte, error) {
  switch i {
  case 0:  return 3, nil
  case 1:  return 2, nil
  case 2:  return 1, nil
  case 3:  return 0, nil
  case 4:  return 4, nil
  default: return 0xFF, fmt.Errorf("ComplementCoded(): `%d' is not a code of the alphabet", int(i))
  }
}

func (GappedNucleotideAlphabet) Complement(i byte) (byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return 't', nil
  case 'C': fallthrough
  case 'c': return 'g', nil
  case 'G': fallthrough
  case 'g': return 'c', nil
  case 'T': fallthrough
  case 't': return 'a', nil
  case 'N': fallthrough
  case 'n': return 'n', nil
  default:  return 0xFF, fmt.Errorf("Complement(): `%c' is not part of the alphabet", i)
  }
}
