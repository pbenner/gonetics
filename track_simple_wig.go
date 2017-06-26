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
import "compress/gzip"
import "errors"
import "fmt"
import "io"
import "math"
import "os"
import "strconv"
import "strings"

/* -------------------------------------------------------------------------- */

func (track SimpleTrack) writeWiggle_fixedStep(w io.Writer, seqname string, sequence []float64) error {
  for i, gap := 0, true; i < len(sequence); i++ {
    if !math.IsNaN(sequence[i]) {
      if gap {
        if _, err := fmt.Fprintf(w, "fixedStep chrom=%s start=%d span=%d step=%d\n", seqname, i*track.BinSize+1, track.BinSize, track.BinSize); err != nil {
          return err
        }
        gap = false
      }
      if _, err := fmt.Fprintf(w, "%f\n", sequence[i]); err != nil {
        return err
      }
    } else {
      gap = true
    }
  }
  return nil
}

func (track SimpleTrack) writeWiggle_variableStep(w io.Writer, seqname string, sequence []float64) error {
  if _, err := fmt.Fprintf(w, "variableStep chrom=%s span=%d\n", seqname, track.BinSize); err != nil {
    return err
  }
  for i := 0; i < len(sequence); i++ {
    if !math.IsNaN(sequence[i]) {
      if _, err := fmt.Fprintf(w, "%d %f\n", i*track.BinSize+1, sequence[i]); err != nil {
        return err
      }
    }
  }
  return nil
}

// Export the track to wiggle format. For sparse tracks (more than half
// of the values are zero), variable step formatting is used.
func (track SimpleTrack) WriteWiggle(filename, description string) error {
  f, err := os.Create(filename)
  if err != nil {
    return err
  }
  defer f.Close()

  w := bufio.NewWriter(f)
  defer w.Flush()

  if _, err := fmt.Fprintf(w, "track type=wiggle_0 name=\"%s\" description=\"%s\"\n", track.Name, description); err != nil {
    return err
  }

  for seqname, sequence := range track.Data {
    n := 0
    // check if sequence is dense or sparse
    for i := 0; i < len(sequence); i++ {
      if math.IsNaN(sequence[i]) || sequence[i] == 0 {
        n++
      }
    }
    if n >= len(sequence)/2 {
      // sparse data track
      if err := track.writeWiggle_variableStep(w, seqname, sequence); err != nil {
        return err
      }
    } else {
      // dense data track
      if err := track.writeWiggle_fixedStep(w, seqname, sequence); err != nil {
        return err
      }
    }
  }
  return nil
}

func readWiggle_header(scanner *bufio.Scanner, result *SimpleTrack) error {
  fields := fieldsQuoted(scanner.Text())

  for i := 1; i < len(fields); i++ {
    headerFields := strings.FieldsFunc(fields[i], func(r rune) bool { return r == '=' })
    if len(headerFields) != 2 {
      return errors.New("invalid declaration line")
    }
    switch headerFields[0] {
    case "name":
      result.Name = removeQuotes(headerFields[1])
    case "type":
      if removeQuotes(headerFields[1]) != "wiggle_0" {
        return errors.New("unsupported wiggle format")
      }
    }
  }
  return nil
}

func readWiggle_fixedStep(scanner *bufio.Scanner, result *SimpleTrack) error {
  fields   := fieldsQuoted(scanner.Text())
  seqname  := ""
  sequence := []float64{}
  position := 0
  // parse header
  for i := 1; i < len(fields); i++ {
    headerFields := strings.FieldsFunc(fields[i], func(r rune) bool { return r == '=' })
    if len(headerFields) != 2 {
      return errors.New("invalid declaration line")
    }
    switch headerFields[0] {
    case "chrom": seqname = removeQuotes(headerFields[1])
    case "start":
      t, err := strconv.ParseInt(headerFields[1], 10, 64)
      if err != nil {
        return err
      }
      position = result.Index(int(t)-1)
    case "step":
      t, err := strconv.ParseInt(headerFields[1], 10, 64)
      if err != nil {
        return err
      }
      if result.BinSize != int(t) {
        return errors.New("step sizes does not match the binSize of the track")
      }
    }
  }
  if seqname == "" {
    return errors.New("declaration line is missing the chromosome name")
  }
  // parse data
  if position < 0 {
    return errors.New("declaration line defines invalid start position")
  }
  sequence, ok := result.Data[seqname]
  // if the sequence is not available in the track, continue parsing
  // the file
  for scanner.Scan() {
    fields = strings.Fields(scanner.Text())
    if len(fields) != 1 {
      break
    }
    t, err := strconv.ParseFloat(fields[0], 64)
    if err != nil {
      return err
    }
    if ok && position < len(sequence) {
      sequence[position] = t
      position++
    }
  }
  if ok {
    result.Data[seqname] = sequence
  }
  return nil
}

func readWiggle_variableStep(scanner *bufio.Scanner, result *SimpleTrack) error {
  fields   := fieldsQuoted(scanner.Text())
  seqname  := ""
  // parse header
  for i := 1; i < len(fields); i++ {
    headerFields := strings.FieldsFunc(fields[i], func(r rune) bool { return r == '=' })
    if len(headerFields) != 2 {
      return errors.New("invalid declaration line")
    }
    switch headerFields[0] {
    case "chrom": seqname = removeQuotes(headerFields[1])
    case "span":
      t, err := strconv.ParseInt(headerFields[1], 10, 64)
      if err != nil {
        return err
      }
      if result.BinSize != int(t) {
        return errors.New("span does not match the binSize of the track")
      }
    }
  }
  if seqname == "" {
    return errors.New("declaration line is missing the chromosome name")
  }
  // parse data
  sequence, ok := result.Data[seqname]
  for scanner.Scan() {
    fields = strings.Fields(scanner.Text())
    if len(fields) != 2 {
      break
    }
    t1, err := strconv.ParseInt(fields[0], 10, 64)
    if err != nil {
      return err
    }
    t2, err := strconv.ParseFloat(fields[1], 64)
    if err != nil {
      return err
    }
    if t1 <= 0 {
      return errors.New("invalid chromosomal position")
    }
    position := result.Index(int(t1)-1)
    if ok && position < len(sequence) {
      sequence[position] = t2
    }
  }
  if ok {
    result.Data[seqname] = sequence
  }
  return nil
}

// Import data from wiggle files.
func (track *SimpleTrack) ReadWiggle(reader io.Reader) error {

  header  := false
  fields  := []string{}
  scanner := bufio.NewScanner(reader)

  if !scanner.Scan() {
    return nil
  }
  for {
    fields = strings.Fields(scanner.Text())
    if len(fields) == 0 {
      break
    }
    // header?
    if fields[0] == "track" {
      if header == false {
        header = true
        err := readWiggle_header(scanner, track)
        if err != nil {
          return err
        }
        if !scanner.Scan() {
          return nil
        }
      } else {
        return errors.New("file contains more than one track definition line")
      }
    } else if fields[0] == "browser" {
      // skip any browser options
      if !scanner.Scan() {
        return nil
      }
    } else if fields[0] == "fixedStep" {
      err := readWiggle_fixedStep(scanner, track)
      if err != nil {
        return err
      }
    } else if fields[0] == "variableStep" {
      err := readWiggle_variableStep(scanner, track)
      if err != nil {
        return err
      }
    } else {
      return errors.New("unknown sequence type (i.e. not fixedStep or variableStep)")
    }
  }
  return nil
}

func (track *SimpleTrack) ImportWiggle(filename string) error {
  var r io.Reader
  // open file
  f, err := os.Open(filename)
  if err != nil {
    return err
  }
  defer f.Close()
  // check if file is gzipped
  if isGzip(filename) {
    g, err := gzip.NewReader(f)
    if err != nil {
      return err
    }
    defer g.Close()
    r = g
  } else {
    r = f
  }
  return track.ReadWiggle(r)
}
