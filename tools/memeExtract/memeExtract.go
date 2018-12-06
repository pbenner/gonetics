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
  Alpha        float64
  FilterEValue float64
   InputFormat string
  OutputFormat string
  OutputType   string
  Verbose      int
}

/* -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* ------------------------------------------------------------------------- */

type Motif interface {
  AsPWM(alphabet Alphabet, background []float64, alpha float64) (TFMatrix, error)
  AsPPM(alphabet Alphabet, background []float64, alpha float64) (TFMatrix, error)
}

/* meme xml format
 * ------------------------------------------------------------------------- */

type Meme struct {
  XMLName xml.Name           `xml:"MEME"`
  Alphabet    MemeAlphabet   `xml:"training_set>alphabet"`
  Model       MemeModel      `xml:"model"`
  Motifs    []MemeMotif      `xml:"motifs>motif"`
}

type MemeAlphabet struct {
  XMLName xml.Name   `xml:"alphabet"`
  Name        string `xml:"name,attr"`
}

type MemeAlphabetMatrix struct {
  XMLName xml.Name              `xml:"alphabet_matrix"`
  Arrays    []MemeAlphabetArray `xml:"alphabet_array"`
}

type MemeAlphabetArray struct {
  XMLName xml.Name              `xml:"alphabet_array"`
  Values    []MemeAlphabetValue `xml:"value"`
}

type MemeAlphabetValue struct {
  XMLName xml.Name    `xml:"value"`
  Letter      string  `xml:"letter_id,attr"`
  Value       string  `xml:",innerxml"`
}

type MemeModel struct {
  XMLName    xml.Name              `xml:"model"`
  Background     MemeAlphabetArray `xml:"background_frequencies>alphabet_array"`
}

type MemeMotif struct {
  XMLName       xml.Name               `xml:"motif"`
  PValue            string             `xml:"p_value,attr"`
  EValue            string             `xml:"e_value,attr"`
  Scores            MemeAlphabetMatrix `xml:"scores>alphabet_matrix"`
  Probabilities     MemeAlphabetMatrix `xml:"probabilities>alphabet_matrix"`
}

/* ------------------------------------------------------------------------- */

func (obj Meme) GetAlphabet() (Alphabet, error) {
  switch strings.ToUpper(obj.Alphabet.Name) {
  case "DNA":
    return NucleotideAlphabet{}, nil
  default:
    return nil, fmt.Errorf("invalid alphabet `%s'", obj.Alphabet.Name)
  }
}

func (obj MemeAlphabetArray) GetValues(alphabet Alphabet) ([]float64, error) {
  r := make([]float64, alphabet.Length())
  for _, value := range obj.Values {
    if len(value.Letter) != 1 {
      return nil, fmt.Errorf("background has invalid letter")
    }
    i, err := alphabet.Code(byte(value.Letter[0]))
    if err != nil {
      return nil, fmt.Errorf("background has invalid letter")
    }
    v, err := strconv.ParseFloat(value.Value, 64)
    if err != nil {
      return nil, fmt.Errorf("background has invalid value: %v", err)
    }
    r[i] = v
  }
  return r, nil
}

func (obj MemeAlphabetMatrix) GetValues(alphabet Alphabet) ([][]float64, error) {
  r := make([][]float64, alphabet.Length())
  for k := 0; k < alphabet.Length(); k++ {
    r[k] = make([]float64, len(obj.Arrays))
  }
  for j, array := range obj.Arrays {
    for _, value := range array.Values {
      if len(value.Letter) != 1 {
        return nil, fmt.Errorf("background has invalid letter")
      }
      i, err := alphabet.Code(byte(value.Letter[0]))
      if err != nil {
        return nil, fmt.Errorf("background has invalid letter")
      }
      v, err := strconv.ParseFloat(value.Value, 64)
      if err != nil {
        return nil, fmt.Errorf("background has invalid value: %v", err)
      }
      r[i][j] = v
    }
  }
  return r, nil
}

func (obj MemeMotif) AsPPM(alphabet Alphabet, background []float64, alpha float64) (TFMatrix, error) {
  v, err := obj.Probabilities.GetValues(alphabet)
  if err != nil {
    return TFMatrix{}, err
  }
  if alpha != 0.0 {
    for i := 0; i < len(v); i++ {
      for j := 0; j < len(v[i]); j++ {
        v[i][j] = (v[i][j] + alpha)/(1.0 + float64(alphabet.Length())*alpha)
      }
    }
  }
  return TFMatrix{v}, nil
}

func (obj MemeMotif) AsPWM(alphabet Alphabet, background []float64, alpha float64) (TFMatrix, error) {
  v, err := obj.Probabilities.GetValues(alphabet)
  if err != nil {
    return TFMatrix{}, err
  }
  for i := 0; i < len(v); i++ {
    for j := 0; j < len(v[i]); j++ {
      // normalize x
      if alpha < 0.0 {
        // add pseudoprobability mass as computed by meme
        v[i][j] = v[i][j] - background[i]*alpha
      } else {
        v[i][j] = (v[i][j] + alpha)/(1.0 + float64(alphabet.Length())*alpha)
      }
      v[i][j] = math.Log2(v[i][j]/background[i])
    }
  }
  return TFMatrix{v}, nil
}

/* dreme xml format
 * ------------------------------------------------------------------------- */

type Dreme struct {
  XMLName xml.Name        `xml:"dreme"`
  Model       DremeModel  `xml:"model"`
  Motifs    []DremeMotif  `xml:"motifs>motif"`
}

type DremeModel struct {
  XMLName xml.Name           `xml:"model"`
  Alphabet   DremeAlphabet   `xml:"alphabet"`
  Background DremeBackground `xml:"background"`
}

type DremeAlphabet struct {
  XMLName xml.Name   `xml:"alphabet"`
  Name        string `xml:"name,attr"`
}

type DremeBackground struct {
  XMLName xml.Name    `xml:"background"`
  Attrs   []xml.Attr  `xml:",any,attr"`
}

type DremeMotif struct {
  XMLName xml.Name     `xml:"motif"`
  Pos       []DremePos `xml:"pos"`
}

type DremePos struct {
  XMLName   xml.Name    `xml:"pos"`
  Attrs   []xml.Attr    `xml:",any,attr"`
}


/* ------------------------------------------------------------------------- */

func (obj DremeModel) GetAlphabet() (Alphabet, error) {
  switch strings.ToUpper(obj.Alphabet.Name) {
  case "DNA":
    return NucleotideAlphabet{}, nil
  default:
    return nil, fmt.Errorf("invalid alphabet `%s'", obj.Alphabet.Name)
  }
}


func (obj DremePos) GetValue(letter byte) (float64, error) {
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

func (obj DremeBackground) GetValues(alphabet Alphabet) ([]float64, error) {
  r := make([]float64, alphabet.Length())
  for k := 0; k < len(obj.Attrs); k++ {
    if len(obj.Attrs[k].Name.Local) != 1 {
      return nil, fmt.Errorf("background has invalid letter")
    }
    i, err := alphabet.Code(byte(obj.Attrs[k].Name.Local[0]))
    if err != nil {
      return nil, fmt.Errorf("background has invalid letter")
    }
    v, err := strconv.ParseFloat(obj.Attrs[k].Value, 64)
    if err != nil {
      return nil, fmt.Errorf("background has invalid value: %v", err)
    }
    r[i] = v
  }
  return r, nil
}

func (obj DremeMotif) AsPPM(alphabet Alphabet, background []float64, alpha float64) (TFMatrix, error) {
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

func (obj DremeMotif) AsPWM(alphabet Alphabet, background []float64, alpha float64) (TFMatrix, error) {
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
      y      := background[i]
      if err != nil {
        return TFMatrix{}, err
      }
      // normalize x
      if alpha < 0.0 {
        // add pseudoprobability mass as computed by meme
        x = x - y*alpha
      } else {
        x = (x + alpha)/(1.0 + float64(alphabet.Length())*alpha)
      }
      values[i][j] = math.Log2(x/y)
    }
  }
  return TFMatrix{values}, nil
}

/* ------------------------------------------------------------------------- */

func writeTFMatrix(config Config, tfmatrix TFMatrix, filename string) {
  f, err := os.Create(filename)
  if err != nil {
    log.Fatal(err)
  }
  defer f.Close()

  switch config.OutputFormat {
  case "table":
    if err := tfmatrix.WriteMatrix(f); err != nil {
      log.Fatal(err)
    }
  case "jaspar":
    if err := tfmatrix.WriteJaspar(f); err != nil {
      log.Fatal(err)
    }
  default:
    panic("internal error")
  }
}

/* ------------------------------------------------------------------------- */

func getTFMatrix(config Config, motif Motif, background []float64, alphabet Alphabet) TFMatrix {
  switch config.OutputType {
  case "pwm":
    if pwm, err := motif.AsPWM(alphabet, background, config.Alpha); err != nil {
      log.Fatal(err)
    } else {
      return pwm
    }
  case "ppm":
    if ppm, err := motif.AsPPM(alphabet, background, config.Alpha); err != nil {
      log.Fatal(err)
    } else {
      return ppm
    }
  default:
    panic("internal error")
  }
  return TFMatrix{}
}

/* ------------------------------------------------------------------------- */

func memeExtract(config Config, filename string, basename string) {
  xmlFile, err := os.Open(filename); if err != nil {
		log.Fatal(err)
	}
  byteValue, err := ioutil.ReadAll(xmlFile); if err != nil {
    log.Fatal(err)
  }
  xmlFile.Close()

  tfmatrices := []TFMatrix{}

  if config.InputFormat == "meme" {
    meme := Meme{}

    if err := xml.Unmarshal(byteValue, &meme); err != nil {
      log.Fatalf("parsing file `%s' failed: %v", filename, err)
    }
    alphabet, err := meme.GetAlphabet(); if err != nil {
      log.Fatal(err)
    }
    background, err := meme.Model.Background.GetValues(alphabet); if err != nil {
      log.Fatal(err)
    }
    tfmatrices = make([]TFMatrix, len(meme.Motifs))

    for i := 0; i < len(meme.Motifs); i++ {
      v, err := strconv.ParseFloat(meme.Motifs[i].EValue, 64); if err != nil {
        log.Fatalf("Parsing motif `%d' failed: %v", i, err)
      }
      if v > config.FilterEValue {
        PrintStderr(config, 1, "Skipping motif %d...\n", i+1)
      } else {
        PrintStderr(config, 1, "Parsing motif %d...\n", i+1)
        tfmatrices[i] = getTFMatrix(config, meme.Motifs[i], background, alphabet)
      }
    }
  }
  if config.InputFormat == "dreme" {
    dreme := Dreme{}

    if err := xml.Unmarshal(byteValue, &dreme); err != nil {
      log.Fatalf("parsing file `%s' failed: %v", filename, err)
    }
    alphabet, err := dreme.Model.GetAlphabet(); if err != nil {
      log.Fatal(err)
    }
    background, err := dreme.Model.Background.GetValues(alphabet); if err != nil {
      log.Fatal(err)
    }
    tfmatrices := make([]TFMatrix, len(dreme.Motifs))

    for i := 0; i < len(dreme.Motifs); i++ {
      PrintStderr(config, 1, "Parsing motif %d...\n", i+1)
      tfmatrices[i] = getTFMatrix(config, dreme.Motifs[i], background, alphabet)
    }
  }
  for i := 0; i < len(tfmatrices); i++ {
    if tfmatrices[i].Values == nil {
      continue
    }
    filename := fmt.Sprintf("%s-%04d.table", basename, i)
    writeTFMatrix(config, tfmatrices[i], filename)
  }
}

/* ------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  optHelp         := options.   BoolLong("help",              'h',             "print help")
  optVerbose      := options.CounterLong("verbose",           'v',             "be verbose")
  optAlpha        := options. StringLong("pseudo-probability", 0 , "-0.00001", "pseudo probability mass added to the PPM (default: 10^-5 times background probability)")
  optFilterEValue := options. StringLong("filter-e-value",     0 , "0.05",     "filter motifs by their e-value")
  optInputFormat  := options. StringLong( "input-format",      0 , "meme",     " input format [meme  (default), dreme]")
  optOutputFormat := options. StringLong("output-format",      0 , "table",    "output format [table (default), jaspar]")
  optOutputType   := options. StringLong("output-type",        0 , "pwm",      "output type   [PWM   (default), PPM]")

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
  if v, err := strconv.ParseFloat(*optFilterEValue, 64); err != nil {
    log.Fatal(err)
  } else {
    if v < 0.0 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    config.FilterEValue = v
  }
  config.InputFormat  = strings.ToLower(*optInputFormat)
  config.OutputFormat = strings.ToLower(*optOutputFormat)
  config.OutputType   = strings.ToLower(*optOutputType)
  config.Verbose      = *optVerbose

  switch config.InputFormat {
  case "meme":
  case "dreme":
  default:
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  switch config.OutputFormat {
  case "table":
  case "jaspar":
  default:
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  switch config.OutputType {
  case "pwm":
  case "ppm":
  default:
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  filename := options.Args()[0]
  basename := options.Args()[1]

  memeExtract(config, filename, basename)
}
