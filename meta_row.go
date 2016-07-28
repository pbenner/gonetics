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

//import "fmt"
import "errors"

/* -------------------------------------------------------------------------- */

type MetaRow struct {
  meta Meta
  idx  int
}

/* -------------------------------------------------------------------------- */

func NewMetaRow(m Meta, i int) MetaRow {
  return MetaRow{m, i}
}

/* -------------------------------------------------------------------------- */

func (m MetaRow) GetMeta(name string) interface{} {
  tmp := m.meta.MetaData[m.idx]
  if tmp == nil {
    return nil
  }
  switch v := tmp.(type) {
  case [][]string:  return v[m.idx]
  case   []string:  return v[m.idx]
  case [][]float64: return v[m.idx]
  case   []float64: return v[m.idx]
  case [][]int:     return v[m.idx]
  case   []int:     return v[m.idx]
  default: panic("Row(): invalid type!")
  }
}

func (m MetaRow) GetMetaStr(name string) (string, error) {
  r := m.GetMeta(name)
  if r != nil {
    return r.(string), nil
  }
  return "", errors.New("invalid meta name")
}

func (m MetaRow) GetMetaFloat(name string) (float64, error) {
  r := m.GetMeta(name)
  if r != nil {
    return r.(float64), nil
  }
  return 0.0, errors.New("invalid meta name")
}

func (m MetaRow) GetMetaInt(name string) (int, error) {
  r := m.GetMeta(name)
  if r != nil {
    return r.(int), nil
  }
  return 0, errors.New("invalid meta name")
}
