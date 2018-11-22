/* Copyright (C) 2018 Philipp Benner
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

// Gardiner-Garden and Frommer, J. Mol. Biol. (1987) 196 (2), 261-282:
//   observed * length / (number of C * number of G)
func ObservedOverExpectedCpG(sequence []byte) float64 {
  n_c   := 0
  n_g   := 0
  n_cpg := 0
  for j := 0; j < len(sequence); j++ {
    if sequence[j] == 'c' || sequence[j] == 'C' {
      n_c += 1
    } else
    if sequence[j] == 'g' || sequence[j] == 'G' {
      n_g += 1
    }
  }
  for j := 0; j < len(sequence)-1; j++ {
    if (sequence[j] == 'c' || sequence[j] == 'C') && (sequence[j+1] == 'g' || sequence[j+1] == 'G') {
      n_cpg += 1
    }
  }
  if n_cpg != 0 {
    return float64(n_cpg*len(sequence))/float64(n_c*n_g)
  } else {
    return 0.0
  }
}
