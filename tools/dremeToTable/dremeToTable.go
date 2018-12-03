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

package main

/* -------------------------------------------------------------------------- */

import   "fmt"
import   "encoding/xml"
import   "io/ioutil"
import   "log"
import   "math"
import   "strconv"
import   "strings"
import   "os"
import   "unicode"

import   "github.com/pborman/getopt"

import . "github.com/pbenner/gonetics"

/* -------------------------------------------------------------------------- */

type Config struct {
  Alpha    float64
  Format   string
  Verbose  int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* ------------------------------------------------------------------------- */

type Dreme struct {
  XMLName xml.Name   `xml:"dreme"`
  Model       Model  `xml:"model"`
  Motifs    []Motif  `xml:"motifs>motif"`
}

type Model struct {
  XMLName xml.Name        `xml:"model"`
  Alphabet   AlphabetType `xml:"alphabet"`
  Background Background   `xml:"background"`
}

type AlphabetType struct {
  XMLName xml.Name   `xml:"alphabet"`
  Name        string `xml:"name,attr"`
}

type Background struct {
  XMLName xml.Name    `xml:"background"`
  Attrs   []xml.Attr  `xml:",any,attr"`
}

type Motif struct {
  XMLName xml.Name `xml:"motif"`
  Pos       []Pos  `xml:"pos"`
}

type Pos struct {
  XMLName   xml.Name    `xml:"pos"`
  Attrs   []xml.Attr    `xml:",any,attr"`
}


/* ------------------------------------------------------------------------- */

func (obj Model) GetAlphabet() (Alphabet, error) {
  switch obj.Alphabet.Name {
  case "DNA":
    return NucleotideAlphabet{}, nil
  default:
    return nil, fmt.Errorf("invalid alphabet `%s'", obj.Alphabet.Name)
  }
}

/* ------------------------------------------------------------------------- */

func (obj Pos) GetValue(letter byte) (float64, error) {
  for k := 0; k < len(obj.Attrs); k++ {
    if len(obj.Attrs[k].Name.Local) != 1 || unicode.ToLower(rune(letter)) != unicode.ToLower(rune(obj.Attrs[k].Name.Local[0])) {
      continue
    }
    v, err := strconv.ParseFloat(obj.Attrs[k].Value, 64); if err != nil {
      return 0.0, err
    }
    return v, nil
  }
  return 0.0, fmt.Errorf("invalid letter `%c'", letter)
}

func (obj Background) GetValue(letter byte) (float64, error) {
  for k := 0; k < len(obj.Attrs); k++ {
    if len(obj.Attrs[k].Name.Local) != 1 || unicode.ToLower(rune(letter)) != unicode.ToLower(rune(obj.Attrs[k].Name.Local[0])) {
      continue
    }
    v, err := strconv.ParseFloat(obj.Attrs[k].Value, 64); if err != nil {
      return 0.0, err
    }
    return v, nil
  }
  return 0.0, fmt.Errorf("invalid background letter `%c'", letter)
}

/* ------------------------------------------------------------------------- */

func (obj Motif) AsPPM(model Model, alpha float64) (TFMatrix, error) {
  alphabet, err := model.GetAlphabet(); if err != nil {
    return TFMatrix{}, err
  }
  // allocate values matrix
  values := make([][]float64, alphabet.Length())
  for i := 0; i < alphabet.Length(); i++ {
    values[i] = make([]float64, len(obj.Pos))
  }
  for i := 0; i < alphabet.Length(); i++ {
    c, err := alphabet.Decode(byte(i))
    if err != nil {
      panic("internal error")
    }
    for j, pos := range obj.Pos {
      x, err := pos.GetValue(c)
      if err != nil {
        return TFMatrix{}, err
      }
      values[i][j] = (x + alpha)/(1.0 + float64(alphabet.Length())*alpha)
    }
  }
  return TFMatrix{values}, nil
}

func (obj Motif) AsPWM(model Model, alpha float64) (TFMatrix, error) {
  alphabet, err := model.GetAlphabet(); if err != nil {
    return TFMatrix{}, err
  }
  // allocate values matrix
  values := make([][]float64, alphabet.Length())
  for i := 0; i < alphabet.Length(); i++ {
    values[i] = make([]float64, len(obj.Pos))
  }
  for i := 0; i < alphabet.Length(); i++ {
    c, err := alphabet.Decode(byte(i))
    if err != nil {
      panic("internal error")
    }
    for j, pos := range obj.Pos {
      x, err := pos.GetValue(c)
      if err != nil {
        return TFMatrix{}, err
      }
      y, err := model.Background.GetValue(c)
      if err != nil {
        return TFMatrix{}, err
      }
      // normalize x
      x = (x + alpha)/(1.0 + float64(alphabet.Length())*alpha)
      values[i][j] = math.Log(x/y)
    }
  }
  return TFMatrix{values}, nil
}

/* ------------------------------------------------------------------------- */

func dremeToTable(config Config, filename string, basename string) {
  xmlFile, err := os.Open(filename); if err != nil {
		log.Fatal(err)
	}
  byteValue, err := ioutil.ReadAll(xmlFile); if err != nil {
    log.Fatal(err)
  }
  xmlFile.Close()

  dreme    := Dreme{}
  tfmatrix := TFMatrix{}

  if err := xml.Unmarshal(byteValue, &dreme); err != nil {
    panic(err)
  }
  for i := 0; i < len(dreme.Motifs); i++ {
    PrintStderr(config, 1, "Parsing motif %d...\n", i+1)
    
    switch config.Format {
    case "pwm":
      if pwm, err := dreme.Motifs[i].AsPWM(dreme.Model, config.Alpha); err != nil {
        log.Fatal(err)
      } else {
        tfmatrix = pwm
      }
      case "ppm":
      if ppm, err := dreme.Motifs[i].AsPPM(dreme.Model, config.Alpha); err != nil {
        panic(err)
      } else {
        tfmatrix = ppm
      }
    default:
      panic("internal error")
    }
    if file, err := os.Create(fmt.Sprintf("%s-%04d.table", basename, i)); err != nil {
      log.Fatal(err)
    } else {
      if err := tfmatrix.WriteMatrix(file); err != nil {
        log.Fatal(err)
      }
    }
  }
}

/* ------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optHelp     := options.   BoolLong("help",              'h',        "print help")
  optVerbose  := options.CounterLong("verbose",           'v',        "be verbose")
  optAlpha    := options. StringLong("pseudo-probability", 0 , "0.0", "pseudo probability mass added to the PPM")
  optFormat   := options. StringLong("format",             0 , "pwm", "output format [PWM (default), PPM]")

  options.SetParameters("<INPUT.bw> <BASENAME>")
  options.Parse(os.Args)

  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if len(options.Args()) != 2 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  } 
  if v, err := strconv.ParseFloat(*optAlpha, 64); err != nil {
    log.Fatal(err)
  } else {
    config.Alpha = v
  }
  config.Format  = strings.ToLower(*optFormat)
  config.Verbose = *optVerbose

  switch config.Format {
  case "pwm":
  case "ppm":
  default:
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filename := options.Args()[0]
  basename := options.Args()[1]

  dremeToTable(config, filename, basename)
}