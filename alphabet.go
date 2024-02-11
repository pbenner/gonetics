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
  Bases            (i byte) ([]byte, error)
  Matching         (i byte) ([]byte, error)
  Code             (i byte) (  byte, error)
  Decode           (i byte) (  byte, error)
  IsAmbiguous      (i byte) (  bool, error)
  IsWildcard       (i byte) (  bool, error)
  Length           ()       int
  LengthUnambiguous()       int
  String           ()       string
}

type ComplementableAlphabet interface {
  Bases            (i byte) ([]byte, error)
  Matching         (i byte) ([]byte, error)
  Code             (i byte) (byte, error)
  Decode           (i byte) (byte, error)
  Complement       (i byte) (byte, error)
  ComplementCoded  (i byte) (byte, error)
  IsAmbiguous      (i byte) (bool, error)
  IsWildcard       (i byte) (bool, error)
  Length           ()       int
  LengthUnambiguous()       int
  String           ()       string
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
  default:  return nil, fmt.Errorf("Bases(): `%c' is not part of the alphabet", i)
  }
}

func (NucleotideAlphabet) Matching(i byte) ([]byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return []byte{'a'}, nil
  case 'C': fallthrough
  case 'c': return []byte{'c'}, nil
  case 'G': fallthrough
  case 'g': return []byte{'g'}, nil
  case 'T': fallthrough
  case 't': return []byte{'t'}, nil
  default:  return nil, fmt.Errorf("Matching(): `%c' is not a non-ambiguous letter of the alphabet", i)
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
  default: return 0xFF, fmt.Errorf("Decode(): `%d' is not a code of the alphabet", int(i))
  }
}

func (NucleotideAlphabet) IsAmbiguous(i byte) (bool, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return false, nil
  case 'C': fallthrough
  case 'c': return false, nil
  case 'G': fallthrough
  case 'g': return false, nil
  case 'T': fallthrough
  case 't': return false, nil
  default:  return false, fmt.Errorf("IsAmbiguous(): `%c' is not part of the alphabet", i)
  }
}

func (NucleotideAlphabet) IsWildcard(i byte) (bool, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return false, nil
  case 'C': fallthrough
  case 'c': return false, nil
  case 'G': fallthrough
  case 'g': return false, nil
  case 'T': fallthrough
  case 't': return false, nil
  default:  return false, fmt.Errorf("IsWildcard(): `%c' is not part of the alphabet", i)
  }
}

func (NucleotideAlphabet) Length() int {
  return 4
}

func (NucleotideAlphabet) LengthUnambiguous() int {
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

func (NucleotideAlphabet) String() string {
  return "nucleotide alphabet"
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
  default:  return nil, fmt.Errorf("Bases(): `%c' is not part of the alphabet", i)
  }
}

func (GappedNucleotideAlphabet) Matching(i byte) ([]byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return []byte{'a', 'n'}, nil
  case 'C': fallthrough
  case 'c': return []byte{'c', 'n'}, nil
  case 'G': fallthrough
  case 'g': return []byte{'g', 'n'}, nil
  case 'T': fallthrough
  case 't': return []byte{'t', 'n'}, nil
  default:  return nil, fmt.Errorf("Matching(): `%c' is not a non-ambiguous letter of the alphabet", i)
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
  default: return 0xFF, fmt.Errorf("Decode(): `%d' is not a code of the alphabet", int(i))
  }
}

func (GappedNucleotideAlphabet) Length() int {
  return 5
}

func (GappedNucleotideAlphabet) LengthUnambiguous() int {
  return 4
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

func (GappedNucleotideAlphabet) IsAmbiguous(i byte) (bool, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return false, nil
  case 'C': fallthrough
  case 'c': return false, nil
  case 'G': fallthrough
  case 'g': return false, nil
  case 'T': fallthrough
  case 't': return false, nil
  case 'N': fallthrough
  case 'n': return true, nil
  default:  return false, fmt.Errorf("IsAmbiguous(): `%c' is not part of the alphabet", i)
  }
}

func (GappedNucleotideAlphabet) IsWildcard(i byte) (bool, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return false, nil
  case 'C': fallthrough
  case 'c': return false, nil
  case 'G': fallthrough
  case 'g': return false, nil
  case 'T': fallthrough
  case 't': return false, nil
  case 'N': fallthrough
  case 'n': return true, nil
  default:  return false, fmt.Errorf("IsAmbiguous(): `%c' is not part of the alphabet", i)
  }
}

func (GappedNucleotideAlphabet) String() string {
  return "gapped nucleotide alphabet"
}

/* -------------------------------------------------------------------------- */

type AmbiguousNucleotideAlphabet struct {
}

func (AmbiguousNucleotideAlphabet) Bases(i byte) ([]byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return []byte{'a'}, nil
  case 'C': fallthrough
  case 'c': return []byte{'c'}, nil
  case 'G': fallthrough
  case 'g': return []byte{'g'}, nil
  case 'T': fallthrough
  case 't': return []byte{'t'}, nil
  case 'W': fallthrough
  case 'w': return []byte{'a', 't'}, nil
  case 'S': fallthrough
  case 's': return []byte{'c', 'g'}, nil
  case 'M': fallthrough
  case 'm': return []byte{'a', 'c'}, nil
  case 'K': fallthrough
  case 'k': return []byte{'g', 't'}, nil
  case 'R': fallthrough
  case 'r': return []byte{'a', 'g'}, nil
  case 'Y': fallthrough
  case 'y': return []byte{'c', 't'}, nil
  case 'B': fallthrough
  case 'b': return []byte{'c', 'g', 't'}, nil
  case 'D': fallthrough
  case 'd': return []byte{'a', 'g', 't'}, nil
  case 'H': fallthrough
  case 'h': return []byte{'a', 'c', 't'}, nil
  case 'V': fallthrough
  case 'v': return []byte{'a', 'c', 'g'}, nil
  case 'N': fallthrough
  case 'n': return []byte{'a', 'c', 'g', 't'}, nil
  default:  return nil, fmt.Errorf("Bases(): `%c' is not part of the alphabet", i)
  }
}

func (AmbiguousNucleotideAlphabet) Matching(i byte) ([]byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return []byte{'a', 'w', 'm', 'r', 'd', 'h', 'v', 'n'}, nil
  case 'C': fallthrough
  case 'c': return []byte{'c', 's', 'm', 'y', 'b', 'h', 'v', 'n'}, nil
  case 'G': fallthrough
  case 'g': return []byte{'g', 's', 'k', 'r', 'b', 'd', 'v', 'n'}, nil
  case 'T': fallthrough
  case 't': return []byte{'t', 'w', 'k', 'y', 'b', 'd', 'h', 'n'}, nil
  default:  return nil, fmt.Errorf("Matching(): `%c' is not a non-ambiguous letter of the alphabet", i)
  }
}

func (AmbiguousNucleotideAlphabet) Code(i byte) (byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return 0, nil
  case 'C': fallthrough
  case 'c': return 1, nil
  case 'G': fallthrough
  case 'g': return 2, nil
  case 'T': fallthrough
  case 't': return 3, nil
  case 'W': fallthrough
  case 'w': return 4, nil
  case 'S': fallthrough
  case 's': return 5, nil
  case 'M': fallthrough
  case 'm': return 6, nil
  case 'K': fallthrough
  case 'k': return 7, nil
  case 'R': fallthrough
  case 'r': return 8, nil
  case 'Y': fallthrough
  case 'y': return 9, nil
  case 'B': fallthrough
  case 'b': return 10, nil
  case 'D': fallthrough
  case 'd': return 11, nil
  case 'H': fallthrough
  case 'h': return 12, nil
  case 'V': fallthrough
  case 'v': return 13, nil
  case 'N': fallthrough
  case 'n': return 14, nil
  default:  return 0xFF, fmt.Errorf("Code(): `%c' is not part of the alphabet", i)
  }
}

func (AmbiguousNucleotideAlphabet) Decode(i byte) (byte, error) {
  switch i {
  case  0:  return 'a', nil
  case  1:  return 'c', nil
  case  2:  return 'g', nil
  case  3:  return 't', nil
  case  4:  return 'w', nil
  case  5:  return 's', nil
  case  6:  return 'm', nil
  case  7:  return 'k', nil
  case  8:  return 'r', nil
  case  9:  return 'y', nil
  case 10:  return 'b', nil
  case 11:  return 'd', nil
  case 12:  return 'h', nil
  case 13:  return 'v', nil
  case 14:  return 'n', nil
  default: return 0xFF, fmt.Errorf("Decode(): `%d' is not a code of the alphabet", int(i))
  }
}

func (AmbiguousNucleotideAlphabet) Length() int {
  return 15
}

func (AmbiguousNucleotideAlphabet) LengthUnambiguous() int {
  return 4
}

func (AmbiguousNucleotideAlphabet) ComplementCoded(i byte) (byte, error) {
  switch i {
  case  0:  return  3, nil
  case  1:  return  2, nil
  case  2:  return  1, nil
  case  3:  return  0, nil
  case  4:  return  4, nil
  case  5:  return  5, nil
  case  6:  return  7, nil
  case  7:  return  6, nil
  case  8:  return  9, nil
  case  9:  return  8, nil
  case 10:  return 13, nil
  case 11:  return 12, nil
  case 12:  return 11, nil
  case 13:  return 10, nil
  case 14:  return 14, nil
  default: return 0xFF, fmt.Errorf("ComplementCoded(): `%d' is not a code of the alphabet", int(i))
  }
}

func (AmbiguousNucleotideAlphabet) Complement(i byte) (byte, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return 't', nil
  case 'C': fallthrough
  case 'c': return 'g', nil
  case 'G': fallthrough
  case 'g': return 'c', nil
  case 'T': fallthrough
  case 't': return 'a', nil
  case 'W': fallthrough
  case 'w': return 'w', nil
  case 'S': fallthrough
  case 's': return 's', nil
  case 'M': fallthrough
  case 'm': return 'k', nil
  case 'K': fallthrough
  case 'k': return 'm', nil
  case 'R': fallthrough
  case 'r': return 'y', nil
  case 'Y': fallthrough
  case 'y': return 'r', nil
  case 'B': fallthrough
  case 'b': return 'v', nil
  case 'D': fallthrough
  case 'd': return 'h', nil
  case 'H': fallthrough
  case 'h': return 'd', nil
  case 'V': fallthrough
  case 'v': return 'b', nil
  case 'N': fallthrough
  case 'n': return 'n', nil
  default:  return 0xFF, fmt.Errorf("Complement(): `%c' is not part of the alphabet", i)
  }
}

func (AmbiguousNucleotideAlphabet) IsAmbiguous(i byte) (bool, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return false, nil
  case 'C': fallthrough
  case 'c': return false, nil
  case 'G': fallthrough
  case 'g': return false, nil
  case 'T': fallthrough
  case 't': return false, nil
  case 'W': fallthrough
  case 'w': return true, nil
  case 'S': fallthrough
  case 's': return true, nil
  case 'M': fallthrough
  case 'm': return true, nil
  case 'K': fallthrough
  case 'k': return true, nil
  case 'R': fallthrough
  case 'r': return true, nil
  case 'Y': fallthrough
  case 'y': return true, nil
  case 'B': fallthrough
  case 'b': return true, nil
  case 'D': fallthrough
  case 'd': return true, nil
  case 'H': fallthrough
  case 'h': return true, nil
  case 'V': fallthrough
  case 'v': return true, nil
  case 'N': fallthrough
  case 'n': return true, nil
  default:  return false, fmt.Errorf("IsAmbiguous(): `%c' is not part of the alphabet", i)
  }
}

func (AmbiguousNucleotideAlphabet) IsWildcard(i byte) (bool, error) {
  switch i {
  case 'A': fallthrough
  case 'a': return false, nil
  case 'C': fallthrough
  case 'c': return false, nil
  case 'G': fallthrough
  case 'g': return false, nil
  case 'T': fallthrough
  case 't': return false, nil
  case 'W': fallthrough
  case 'w': return false, nil
  case 'S': fallthrough
  case 's': return false, nil
  case 'M': fallthrough
  case 'm': return false, nil
  case 'K': fallthrough
  case 'k': return false, nil
  case 'R': fallthrough
  case 'r': return false, nil
  case 'Y': fallthrough
  case 'y': return false, nil
  case 'B': fallthrough
  case 'b': return false, nil
  case 'D': fallthrough
  case 'd': return false, nil
  case 'H': fallthrough
  case 'h': return false, nil
  case 'V': fallthrough
  case 'v': return false, nil
  case 'N': fallthrough
  case 'n': return true, nil
  default:  return false, fmt.Errorf("IsWildcard(): `%c' is not part of the alphabet", i)
  }
}

func (AmbiguousNucleotideAlphabet) String() string {
  return "ambiguous nucleotide alphabet"
}
