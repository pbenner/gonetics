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

import "bufio"
import "fmt"
import "os"

/* -------------------------------------------------------------------------- */

// Export GRanges as GTF file. Required GTF fields should be provided
// as meta columns named sources, features, scores, and frames. All other
// meta columns are exported as optional fields.
func (granges GRanges) WriteGTF(filename string) {
  f, err := os.Create(filename); check(err)
  defer f.Close()

  w := bufio.NewWriter(f)
  defer w.Flush()

  sources  := granges.GetMetaStr("sources")
  features := granges.GetMetaStr("features")
  scores   := granges.GetMetaInt("scores")
  frames   := granges.GetMetaInt("frames")

  for i := 0; i < granges.Length(); i++ {
    fmt.Fprintf(w, "%s", granges.Seqnames[i])
    if len(sources) > 0 {
      fmt.Fprintf(w, "\t%s", sources[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(features) > 0 {
      fmt.Fprintf(w, "\t%s", features[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].From)
    fmt.Fprintf(w, "\t%d", granges.Ranges[i].To)
    if len(scores) > 0 {
      fmt.Fprintf(w, "\t%d", scores[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }
    if len(granges.Strand) > 0 && granges.Strand[i] != '*' {
      fmt.Fprintf(w, "\t%c", granges.Strand[i])
    } else {
      fmt.Fprintf(w, "\t%c", '.')
    }
    if len(frames) > 0 {
      fmt.Fprintf(w, "\t%d", frames[i])
    } else {
      fmt.Fprintf(w, "\t%s", ".")
    }

    if granges.MetaLength() != 0 {
      printedTab := false
      for k := 0; k < granges.MetaLength(); k++ {
        if granges.MetaName[k] == "sources" {
          continue
        }
        if granges.MetaName[k] == "features" {
          continue
        }
        if granges.MetaName[k] == "scores" {
          continue
        }
        if granges.MetaName[k] == "frames" {
          continue
        }
        if printedTab {
          w.WriteString(" ")
        } else {
          w.WriteString("\t")
          printedTab = true
        }
        // print name of the meta data
        fmt.Fprintf(w, "%s ", granges.MetaName[k])
        // print data
        switch v := granges.MetaData[k].(type) {
        case []string : fmt.Fprintf(w, "\"%s\"", v[i])
        case []float64: fmt.Fprintf(w, "\"%f\"", v[i])
        case []int    : fmt.Fprintf(w, "\"%d\"", v[i])
        case [][]string:
          w.WriteString("\"")
          for j := 0; j < len(v[i]); j++ {
            if j != 0 {
              w.WriteString(" ")
            }
            fmt.Fprintf(w, "%s", v[i][j])
          }
          w.WriteString("\"")
        case [][]float64:
          w.WriteString("\"")
          for j := 0; j < len(v[i]); j++ {
            if j != 0 {
              w.WriteString(" ")
            }
            fmt.Fprintf(w, "%f", v[i][j])
          }
          w.WriteString("\"")
        case [][]int:
          w.WriteString("\"")
          for j := 0; j < len(v[i]); j++ {
            if j != 0 {
              w.WriteString(" ")
            }
            fmt.Fprintf(w, "%d", v[i][j])
          }
          w.WriteString("\"")
        }
        w.WriteString(";")
      }
    }
    w.WriteString("\n")
  }
}
