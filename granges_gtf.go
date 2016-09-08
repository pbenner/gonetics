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
import "bufio"
import "compress/gzip"
import "os"
import "strconv"
import "strings"
import "unicode"

/* i/o
 * -------------------------------------------------------------------------- */

type gtfOptional map[string]interface{}
type gtfTypeMap  map[string]string

func readGTFParseOptional(fields []string, gtfOpt gtfOptional, typeMap gtfTypeMap, length int) (gtfOptional, error) {
  if len(fields) % 2 == 1 {
    return nil, fmt.Errorf("ReadGTF(): invalid file format!")
  }
  // loop through list
  for i := 0; i < len(fields); i += 2 {
    name     := fields[i]
    valueStr := fields[i+1]
    if _, ok := typeMap[name]; ok {
      // get the data vector
      switch v := gtfOpt[name].(type) {
      case []int:
        if value, err := strconv.ParseInt(valueStr, 10, 64); err != nil {
          return nil, err
        } else {
          v = append(v, int(value))
        }
      case []float64:
        if value, err := strconv.ParseFloat(valueStr, 64); err != nil {
          return nil, err
        } else {
          v = append(v, value)
        }
      case []string:
        v = append(v, valueStr)
      }
    }
  }
  // check that all optional fields are available
  for name, values := range gtfOpt {
    thisLength := 0
    switch v := values.(type) {
    case []int    : thisLength = len(v)
    case []float64: thisLength = len(v)
    case []string : thisLength = len(v)
    }
    if thisLength == length-1 {
      return nil, fmt.Errorf("optional field `%s' is missing at line `%d'", name, length+1)
    }
  }
  return gtfOpt, nil
}

func readGTFParseLine(line string) []string {
  // if quoted
  q := false
  f := func(r rune) bool {
    if r == '"' {
      q = !q
    }
    // A quote is treated as a white space so that it is removed from the
    // line. Otherwise a white space is removed only if q (quote) is false.
    return r == '"' || ((unicode.IsSpace(r) || r == ';') && q == false)
  }
  return strings.FieldsFunc(line, f)
}

// Parse expression data from a GTF file (gene transfer format). The data
// is added as a meta column named "expr" to the gene list. Parameters:
//  geneIdName: Name of the optional field containing the gene id
//  exprIdName: Name of the optional field containing the expression data
//  genes: List of query genes
func (granges *GRanges) ReadGTF(filename string, optNames, optTypes []string, verbose bool) error {
  var scanner *bufio.Scanner
  // open file
  f, err := os.Open(filename)
  check(err)
  defer f.Close()
  // check if file is gzipped
  if isGzip(filename) {
    g, err := gzip.NewReader(f)
    check(err)
    defer g.Close()
    scanner = bufio.NewScanner(g)
  } else {
    scanner = bufio.NewScanner(f)
  }
  if len(optNames) != len(optTypes) {
    return fmt.Errorf("ReadGTF(): invalid arguments")
  }
  // construct type map
  seqname  := []string{}
  source   := []string{}
  feature  := []string{}
  start    := []int{}
  end      := []int{}
  score    := []float64{}
  strand   := []byte{}
  frame    := []int{}
  gtfOpt   := make(gtfOptional)
  typeMap  := make(gtfTypeMap)
  for i := 0; i < len(optNames); i++ {
    typeMap[optNames[i]] = optTypes[i]
  }
  for name, typeStr := range typeMap {
    switch typeStr {
    case "[]int"    : gtfOpt[name] = []int{}
    case "[]float64": gtfOpt[name] = []float64{}
    case "[]string" : gtfOpt[name] = []string{}
    default:
      return fmt.Errorf("ReadGTF(): invalid type `%s' for optional field `%s'", typeStr, name)
    }
  }

  for i := 0; scanner.Scan(); i++ {
    if err := scanner.Err(); err != nil {
      return err
    }
    fields := readGTFParseLine(scanner.Text())
    if len(fields) == 0 {
      continue
    }
    if len(fields) < 8 {
      panic("File must have at least eight columns!")
    }
    seqname = append(seqname, fields[0])
    source  = append(source,  fields[1])
    feature = append(feature, fields[2])
    if v, err := strconv.ParseInt(fields[3], 10, 64); err != nil {
      return err
    } else {
      start = append(start, int(v))
    }
    if v, err := strconv.ParseInt(fields[4], 10, 64); err != nil {
      return err
    } else {
      end = append(end, int(v))
    }
    if fields[5] == "." {
      score = append(score, 0)
    } else {
      if v, err := strconv.ParseFloat(fields[5], 64); err != nil {
        return err
      } else {
        score = append(score, v)
      }
    }
    strand = append(strand, fields[6][0])
    if fields[7] == "." {
      frame = append(frame, -1)
    } else {
      if v, err := strconv.ParseInt(fields[7], 10, 64); err != nil {
        return err
      } else {
        frame = append(frame, int(v))
      }
    }
    // parse optional fields
    if tmp, err := readGTFParseOptional(fields[8:len(fields)], gtfOpt, typeMap, i); err != nil {
      return err
    } else {
      gtfOpt = tmp
    }
  }
  for name, values := range gtfOpt {
    granges.AddMeta(name, values)
  }
  return nil
}

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
