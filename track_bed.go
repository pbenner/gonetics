/* Copyright (C) 2016, 2017 Philipp Benner
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
import "io"

/* -------------------------------------------------------------------------- */

func (track GenericTrack) WriteBed(w io.Writer) error {
  r, err := track.GRanges(); if err != nil {
    return err
  }
  // write to file
  return r.WriteBed6(w)
}

func (track GenericTrack) ExportBed(filename string, compress bool) error {
  r, err := track.GRanges(); if err != nil {
    return err
  }
  // write to file
  return r.ExportBed6(filename, compress)
}
