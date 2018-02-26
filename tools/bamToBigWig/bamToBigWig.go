/* Copyright (C) 2016-2017 Philipp Benner
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
import   "bufio"
import   "log"
import   "math"
import   "path/filepath"
import   "strconv"
import   "strings"
import   "os"

import   "github.com/pborman/getopt"
import . "github.com/pbenner/gonetics"

import   "gonum.org/v1/plot"
import   "gonum.org/v1/plot/plotter"
import   "gonum.org/v1/plot/plotutil"
import   "gonum.org/v1/plot/vg"

/* -------------------------------------------------------------------------- */

type Config struct {
  Verbose                int
  BinSummaryStatistics   string
  BWZoomLevels         []int
  BinningMethod          string
  BinSize                int
  BinOverlap             int
  NormalizeTrack         string
  ShiftReads          [2]int
  PairedAsSingleEnd      bool
  LogScale               bool
  Pseudocounts        [2]float64
  FraglenRange        [2]int
  FilterMapQ             int
  FilterReadLengths   [2]int
  FilterDuplicates       bool
  FilterStrand           byte
  FilterPairedEnd        bool
  FilterSingleEnd        bool
  SmoothenControl        bool
  SmoothenSizes        []int
  SmoothenMin            float64
  SaveFraglen            bool
  SaveCrossCorr          bool
  SaveCrossCorrPlot      bool
}

func DefaultConfig() Config {
  config := Config{}
  // set default values
  config.BinSummaryStatistics = "mean"
  config.BWZoomLevels         = nil   // zoom levels are determined automatically
  config.BinningMethod        = "simple"
  config.BinSize              = 10
  config.BinOverlap           = 0
  config.FraglenRange         = [2]int{-1, -1}
  config.FilterReadLengths    = [2]int{0,0}
  config.FilterMapQ           = 0
  config.FilterDuplicates     = false
  config.FilterStrand         = '*'
  config.FilterPairedEnd      = false
  config.FilterSingleEnd      = false
  config.LogScale             = false
  config.Pseudocounts         = [2]float64{0.0, 0.0}
  config.SmoothenControl      = false
  config.SmoothenSizes        = []int{}
  config.SmoothenMin          = 20.0
  config.SaveFraglen          = false
  config.SaveCrossCorr        = false
  config.SaveCrossCorrPlot    = false
  return config
}

/* i/o
 * -------------------------------------------------------------------------- */

func PrintStderr(config Config, level int, format string, args ...interface{}) {
  if config.Verbose >= level {
    fmt.Fprintf(os.Stderr, format, args...)
  }
}

/* utility
 * -------------------------------------------------------------------------- */

func parseFilename(filename string) (string, int) {
  if tmp := strings.Split(filename, ":"); len(tmp) == 2 {
    t, err := strconv.ParseInt(tmp[1], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    return tmp[0], int(t)
  } else
  if len(tmp) >= 2 {
    log.Fatal("invalid input file description `%s'", filename)
  }
  return filename, 0
}

/* read filters
 * -------------------------------------------------------------------------- */

// treat all paired end reads as single end reads, this allows
// to extend/crop paired end reads when adding them to the track
// with AddReads()
func filterPairedAsSingleEnd(config Config, chanIn ReadChannel) ReadChannel {
  if config.PairedAsSingleEnd == false {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    for r := range chanIn {
      r.PairedEnd = false; chanOut <- r
    }
    close(chanOut)
  }()
  return chanOut
}

func filterPairedEnd(config Config, chanIn ReadChannel) ReadChannel {
  if config.FilterPairedEnd == false {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if r.PairedEnd {
        chanOut <- r; m++
      }
      n++
    }
    PrintStderr(config, 1, "Filtered out %d unpaired reads (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    close(chanOut)
  }()
  return chanOut
}

func filterSingleEnd(config Config, veto bool, chanIn ReadChannel) ReadChannel {
  if config.FilterSingleEnd == false && !veto {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if !r.PairedEnd {
        chanOut <- r; m++
      }
      n++
    }
    PrintStderr(config, 1, "Filtered out %d paired reads (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    close(chanOut)
  }()
  return chanOut
}

func filterDuplicates(config Config, chanIn ReadChannel) ReadChannel {
  if config.FilterDuplicates == false {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if !r.Duplicate {
        chanOut <- r; m++
      }
      n++
    }
    PrintStderr(config, 1, "Filtered out %d duplicates (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    close(chanOut)
  }()
  return chanOut
}

func filterStrand(config Config, chanIn ReadChannel) ReadChannel {
  if config.FilterStrand == '*' {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if r.Strand == config.FilterStrand {
        chanOut <- r; m++
      }
      n++
    }
    PrintStderr(config, 1, "Filtered out %d reads not on strand %v (%.2f%%)\n", n-m, config.FilterStrand, 100.0*float64(n-m)/float64(n))
    close(chanOut)
  }()
  return chanOut
}

func filterMapQ(config Config, chanIn ReadChannel) ReadChannel {
  if config.FilterMapQ <= 0 {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      if r.MapQ >= config.FilterMapQ {
        chanOut <- r; m++
      }
      n++
    }
    PrintStderr(config, 1, "Filtered out %d reads with mapping quality lower than %d (%.2f%%)\n", n-m, config.FilterMapQ, 100.0*float64(n-m)/float64(n))
    close(chanOut)
  }()
  return chanOut
}

func filterReadLength(config Config, chanIn ReadChannel) ReadChannel {
  if config.FilterReadLengths[0] == 0 && config.FilterReadLengths[1] == 0 {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    n := 0
    m := 0
    for r := range chanIn {
      len := r.Range.To - r.Range.From
      if len >= config.FilterReadLengths[0] &&
        (len <= config.FilterReadLengths[1] || config.FilterReadLengths[1] == 0) {
        chanOut <- r; m++
      }
      n++
    }
    PrintStderr(config, 1, "Filtered out %d reads with non-admissible length (%.2f%%)\n", n-m, 100.0*float64(n-m)/float64(n))
    close(chanOut)
  }()
  return chanOut
}

func shiftReads(config Config, chanIn ReadChannel) ReadChannel {
  if config.ShiftReads[0] == 0 && config.ShiftReads[1] == 0 {
    return chanIn
  }
  chanOut := make(chan Read)
  go func() {
    for r := range chanIn {
      if r.Strand == '+' {
        r.Range.From += config.ShiftReads[0]
        r.Range.To   += config.ShiftReads[0]
      } else
      if r.Strand == '-' {
        r.Range.From += config.ShiftReads[1]
        r.Range.To   += config.ShiftReads[1]
      }
      if r.Range.From < 0 {
        r.Range.To   -= r.Range.From
        r.Range.From  = 0
      }
      chanOut <- r
    }
    PrintStderr(config, 1, "Shifted reads (forward strand: %d, reverse strand: %d)\n",
      config.ShiftReads[0], config.ShiftReads[1])
    close(chanOut)
  }()
  return chanOut
}

/* fragment length estimation
 * -------------------------------------------------------------------------- */

func saveFraglen(config Config, filename string, fraglen int) {
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.txt", basename)

  f, err := os.Create(filename)
  if err != nil {
    log.Fatalf("opening `%s' failed: %v", filename, err)
  }
  defer f.Close()

  fmt.Fprintf(f, "%d\n", fraglen)

  PrintStderr(config, 1, "Wrote fragment length estimate to `%s'\n", filename)
}

func saveCrossCorr(config Config, filename string, x []int, y []float64) {
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.table", basename)

  f, err := os.Create(filename)
  if err != nil {
    log.Fatalf("opening `%s' failed: %v", filename, err)
  }
  defer f.Close()

  for i := 0; i < len(x); i++ {
    fmt.Fprintf(f, "%d %f\n", x[i], y[i])
  }
  PrintStderr(config, 1, "Wrote crosscorrelation to `%s'\n", filename)
}

func saveCrossCorrPlot(config Config, filename string, fraglen int, x []int, y []float64) {
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.pdf", basename)

  // draw cross-correlation
  xy := make(plotter.XYs, len(x))
  for i := 0; i < len(x); i++ {
    xy[i].X = float64(x[i])+1
    xy[i].Y = y[i]
  }
  p, err := plot.New()
  if err != nil {
    log.Fatal(err)
  }
  p.Title.Text = ""
  p.X.Label.Text = "shift"
  p.Y.Label.Text = "cross-correlation"

  err = plotutil.AddLines(p, xy)
  if err != nil {
    log.Fatal(err)
  }

  if fraglen != -1 {
    // determine cross-correlation maximum
    max_value := 0.0
    for i := 0; i < len(x); i++ {
      if y[i] > max_value {
        max_value = y[i]
      }
    }
    // draw vertical line at fraglen estimate
    fr := make(plotter.XYs, 2)
    fr[0].X = float64(fraglen)
    fr[0].Y = 0.0
    fr[1].X = float64(fraglen)
    fr[1].Y = max_value

    err = plotutil.AddLines(p, fr)
    if err != nil {
      log.Fatal(err)
    }
  }
  if err := p.Save(8*vg.Inch, 4*vg.Inch, filename); err != nil {
    log.Fatal(err)
  }
  PrintStderr(config, 1, "Wrote cross-correlation plot to `%s'\n", filename)
}

func importFraglen(config Config, filename string, genome Genome) int {
  // try reading the fragment length from file
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.txt", basename)
  if f, err := os.Open(filename); err != nil {
    return -1
  } else {
    defer f.Close()
    PrintStderr(config, 1, "Reading fragment length from `%s'... ", filename)
    scanner := bufio.NewScanner(f)
    if scanner.Scan() {
      if fraglen, err := strconv.ParseInt(scanner.Text(), 10, 64); err == nil {
        PrintStderr(config, 1, "done\n")
        return int(fraglen)
      }
    }
    PrintStderr(config, 1, "failed\n")
    return -1
  }
}

func estimateFraglen(config Config, filename string, genome Genome) int {
  var reads ReadChannel

  PrintStderr(config, 1, "Reading tags from `%s'...\n", filename)
  if bam, err := OpenBamFile(filename, BamReaderOptions{}); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    defer bam.Close()
    reads = bam.ReadSimple(false)
  }

  // first round of filtering
  reads = filterSingleEnd(config, true, reads)
  reads = filterReadLength(config, reads)
  reads = filterDuplicates(config, reads)
  reads = filterMapQ(config, reads)

  // estimate fragment length
  PrintStderr(config, 1, "Estimating mean fragment length... ")
  if fraglen, x, y, err := EstimateFragmentLength(reads, genome, 2000, config.BinSize, config.FraglenRange); err != nil {
    PrintStderr(config, 1, "failed\n")
    if x != nil && y != nil && config.SaveCrossCorr {
      saveCrossCorr(config, filename, x, y)
    }
    if x != nil && y != nil && config.SaveCrossCorrPlot {
      saveCrossCorrPlot(config, filename, -1, x, y)
    }
    log.Fatalf("estimating read length failed: %v", err)
    return 0
  } else {
    PrintStderr(config, 1, "done\n")
    PrintStderr(config, 1, "Estimated mean fragment length: %d\n", fraglen)

    if config.SaveFraglen {
      saveFraglen(config, filename, fraglen)
    }
    if config.SaveCrossCorr {
      saveCrossCorr(config, filename, x, y)
    }
    if config.SaveCrossCorrPlot {
      saveCrossCorrPlot(config, filename, fraglen, x, y)
    }
    return fraglen
  }
}

/* -------------------------------------------------------------------------- */

func bamToBigWig(config Config, filenameTrack string, filenamesTreatment, filenamesControl []string, fraglenTreatment, fraglenControl []int, genome Genome) {

  // treatment data
  track1 := AllocSimpleTrack("treatment", genome, config.BinSize)

  // number of reads
  n_treatment := 0
  n_control   := 0

  for i, filename := range filenamesTreatment {
    fraglen := fraglenTreatment[i]

    var treatment ReadChannel
    PrintStderr(config, 1, "Reading treatment tags from `%s'...\n", filename)
    if bam, err := OpenBamFile(filename, BamReaderOptions{}); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      defer bam.Close()
      treatment = bam.ReadSimple(!config.PairedAsSingleEnd)
    }

    // first round of filtering
    treatment = filterPairedAsSingleEnd(config, treatment)
    treatment = filterPairedEnd(config, treatment)
    treatment = filterSingleEnd(config, false, treatment)
    treatment = filterReadLength(config, treatment)
    treatment = filterDuplicates(config, treatment)
    treatment = filterMapQ(config, treatment)
    // second round of filtering
    treatment = filterStrand(config, treatment)
    treatment = shiftReads(config, treatment)

    switch config.BinningMethod {
    case "simple":
      n_treatment += GenericMutableTrack{track1}.AddReads(treatment, fraglen, false)
    case "overlap":
      n_treatment += GenericMutableTrack{track1}.AddReads(treatment, fraglen, true)
    default:
      log.Fatal("invalid binning method `%s'", config.BinningMethod)
    }
  }
  if config.NormalizeTrack == "rpm" {
    PrintStderr(config, 1, "Normalizing treatment track (rpm)... ")
    c := float64(1000000)/float64(n_treatment)
    GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 {
      return c*x
    })
    // adapt pseudocounts!
    config.Pseudocounts[0] *= c
    PrintStderr(config, 1, "done\n")
  }

  if len(filenamesControl) > 0 {
    // control data
    track2 := AllocSimpleTrack("control", genome, config.BinSize)

    for i, filename := range filenamesControl {
      fraglen  := fraglenControl[i]

      var control ReadChannel
      PrintStderr(config, 1, "Reading treatment tags from `%s'...\n", filename)
      if bam, err := OpenBamFile(filename, BamReaderOptions{}); err != nil {
        PrintStderr(config, 1, "failed\n")
        log.Fatal(err)
      } else {
        defer bam.Close()
        control = bam.ReadSimple(true)
      }

      // first round of filtering
      control = filterPairedAsSingleEnd(config, control)
      control = filterPairedEnd(config, control)
      control = filterSingleEnd(config, false, control)
      control = filterReadLength(config, control)
      control = filterDuplicates(config, control)
      control = filterMapQ(config, control)
      // second round of filtering
      control = filterStrand(config, control)
      control = shiftReads(config, control)

      switch config.BinningMethod {
      case "simple":
        n_control += GenericMutableTrack{track2}.AddReads(control, fraglen, false)
      case "overlap":
        n_control += GenericMutableTrack{track2}.AddReads(control, fraglen, true)
      default:
        log.Fatal("invalid binning method `%s'", config.BinningMethod)
      }
    }
    if config.NormalizeTrack == "rpm" {
      PrintStderr(config, 1, "Normalizing control track (rpm)... ")
      c := float64(1000000)/float64(n_control)
      GenericMutableTrack{track2}.Map(track2, func(name string, i int, x float64) float64 {
        return c*x
      })
      // adapt pseudocounts!
      config.Pseudocounts[1] *= c
      PrintStderr(config, 1, "done\n")
    }
    if config.SmoothenControl {
      GenericMutableTrack{track2}.Smoothen(config.SmoothenMin, config.SmoothenSizes)
    }
    PrintStderr(config, 1, "Combining treatment and control tracks... ")
    if err := (GenericMutableTrack{track1}).Normalize(track1, track2, config.Pseudocounts[0], config.Pseudocounts[1], config.LogScale); err != nil {
      PrintStderr(config, 1, "failed\n")
      log.Fatalf("normalizing track failed: %v", err)
    }
    PrintStderr(config, 1, "done\n")
  } else {
    // no control data
    if config.Pseudocounts[0] != 0.0 {
      PrintStderr(config, 1, "Adding pseudocount `%f'... ", config.Pseudocounts[0])
      GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 { return x+config.Pseudocounts[0] })
      PrintStderr(config, 1, "done\n")
    }
    if config.LogScale {
      PrintStderr(config, 1, "Log-transforming data... ")
      GenericMutableTrack{track1}.Map(track1, func(name string, i int, x float64) float64 { return math.Log(x) })
      PrintStderr(config, 1, "done\n")
    }
  }
  PrintStderr(config, 1, "Writing track `%s'... ", filenameTrack)
  parameters := DefaultBigWigParameters()
  parameters.ReductionLevels = config.BWZoomLevels
  if err := (GenericTrack{track1}).ExportBigWig(filenameTrack, parameters); err != nil {
    PrintStderr(config, 1, "failed\n")
    log.Fatal(err)
  } else {
    PrintStderr(config, 1, "done\n")
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config := DefaultConfig()

  options := getopt.New()

  // bigWig options
  optBWZoomLevels      := options. StringLong("bigwig-zoom-levels",         0 , "", "comma separated list of BigWig zoom levels")
  // read options
  optShiftReads        := options. StringLong("shift-reads",                0 , "", "shift reads on the positive strand by `x' bps and those on the negative strand by `y' bps [format: x,y]")
  optPairedAsSingleEnd := options.   BoolLong("paired-as-single-end",       0 ,     "treat paired as single end reads")
  // options for filterering reads
  optFilterStrand      := options. StringLong("filter-strand",              0 , "", "use reads on either the forward `+' or reverse `-' strand")
  optReadLength        := options. StringLong("filter-read-lengths",        0 , "", "feasible range of read-lengths [format: min:max]")
  optFilterMapQ        := options.    IntLong("filter-mapq",                0 ,  0, "filter reads for minimum mapping quality (default: 0)")
  optFilterDuplicates  := options.   BoolLong("filter-duplicates",          0 ,     "remove reads marked as duplicates")
  optFilterPairedEnd   := options.   BoolLong("filter-paired-end",          0 ,     "remove all single end reads")
  optFilterSingleEnd   := options.   BoolLong("filter-single-end",          0 ,     "remove all paired end reads")
  // track options
  optBinningMethod     := options. StringLong("binning-method",             0 , "", "binning method (i.e. simple [default] or overlap)")
  optBinSize           := options.    IntLong("bin-size",                   0 ,  0, "track bin size [default: 10]")
  optNormalizeTrack    := options. StringLong("normalize-track",            0 , "", "normalize track with the specified method (i.e. rpm)")
  optPseudocounts      := options. StringLong("pseudocounts",               0 , "", "pseudocounts added to treatment and control signal (default: `0,0')")
  optSmoothenControl   := options.   BoolLong("smoothen-control",           0 ,     "smoothen control with an adaptive window method")
  optSmoothenSizes     := options. StringLong("smoothen-window-sizes",      0 , "", "feasible window sizes for the smoothening method [format: s1,s2,...]")
  optSmoothenMin       := options. StringLong("smoothen-min-counts",        0 , "", "minimum number of counts for the smoothening method")
  optLogScale          := options.   BoolLong("log-scale",                  0 ,     "log-transform data")
  // options for estimating and setting fragment lengths
  optFraglen           := options.    IntLong("fragment-length",            0 ,  0, "fragment length for all input files (reads are extended to the given length)")
  optFraglenRange      := options. StringLong("fragment-length-range",      0 , "", "feasible range of fragment lengths (format from:to)")
  optEstimateFraglen   := options.   BoolLong("estimate-fragment-length",   0 ,     "use crosscorrelation to estimate the fragment length")
  optSaveFraglen       := options.   BoolLong("save-fraglen",               0 ,     "save estimated fragment length in a file named <BAM_BASENAME>.fraglen.txt")
  optSaveCrossCorr     := options.   BoolLong("save-crosscorrelation",      0 ,     "save crosscorrelation between forward and reverse strands in a file named <BAM_BASENAME>.fraglen.table")
  optSaveCrossCorrPlot := options.   BoolLong("save-crosscorrelation-plot", 0 ,     "save crosscorrelation plot in a file names <BAM_BASENAME>.fraglen.pdf")
  // generic options
  optVerbose           := options.CounterLong("verbose",                   'v',     "verbose level [-v or -vv]")
  optHelp              := options.   BoolLong("help",                      'h',     "print help")

  options.SetParameters("<TREATMENT1.bam[:FRAGLEN],TREATMENT2.bam[:FRAGLEN],...> [<CONTROL1.bam[:FRAGLEN],CONTROL2.bam[:FRAGLEN],...>] <RESULT.bw>")
  options.Parse(os.Args)

  // parse options
  //////////////////////////////////////////////////////////////////////////////
  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if *optVerbose != 0 {
    config.Verbose = *optVerbose
  }
  if len(options.Args()) != 2 && len(options.Args()) != 3 {
    options.PrintUsage(os.Stderr)
    os.Exit(1)
  }
  if *optBinSize != 0 {
    if *optBinSize < 1 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    } else {
      config.BinSize = *optBinSize
    }
  }
  if *optBinningMethod != "" {
    config.BinningMethod = *optBinningMethod
  }
  if *optPseudocounts != "" {
    tmp := strings.Split(*optPseudocounts, ",")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseFloat(tmp[0], 64)
    if err != nil {
      log.Fatal(err)
    }
    t2, err := strconv.ParseFloat(tmp[1], 64)
    if err != nil {
      log.Fatal(err)
    }
    config.Pseudocounts[0] = t1
    config.Pseudocounts[1] = t2
  }
  if *optReadLength != "" {
    tmp := strings.Split(*optReadLength, ":")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseInt(tmp[0], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    t2, err := strconv.ParseInt(tmp[1], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    if t1 > t2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    config.FilterReadLengths[0] = int(t1)
    config.FilterReadLengths[1] = int(t2)
  }
  if *optFilterMapQ < 0 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
  } else {
    config.FilterMapQ = *optFilterMapQ
  }
  if *optFilterStrand != "" {
    switch *optFilterStrand {
    case "+": config.FilterStrand = '+'
    case "-": config.FilterStrand = '-'
    default:
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
  }
  if *optShiftReads != "" {
    tmp := strings.Split(*optShiftReads, ",")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseInt(tmp[0], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    t2, err := strconv.ParseInt(tmp[1], 10, 64)
    if err != nil {
      log.Fatal(err)
    }
    config.ShiftReads[0] = int(t1)
    config.ShiftReads[1] = int(t2)
  }
  if *optSmoothenControl {
    config.SmoothenControl = true
  }
  if *optSmoothenSizes != "" {
    config.SmoothenSizes = []int{}
    tmp := strings.Split(*optSmoothenSizes, ",")
    for i := 0; i < len(tmp); i++ {
      t, err := strconv.ParseInt(tmp[i], 10, 64)
      if err != nil {
        log.Fatal(err)
      }
      config.SmoothenSizes = append(config.SmoothenSizes, int(t))
    }
  }
  if *optSmoothenMin != "" {
    t, err := strconv.ParseFloat(*optSmoothenMin, 64)
    if err != nil {
      log.Fatal(err)
    }
    config.SmoothenMin = t
  }
  if *optBWZoomLevels != "" {
    tmp := strings.Split(*optBWZoomLevels, ",")
    config.BWZoomLevels = []int{}
    for i := 0; i < len(tmp); i++ {
      if t, err := strconv.ParseInt(tmp[i], 10, 64); err != nil {
        log.Fatal(err)
      } else {
        config.BWZoomLevels = append(config.BWZoomLevels, int(t))
      }
    }
  }
  if *optNormalizeTrack != "" {
    switch *optNormalizeTrack {
    case "rpm":
    default:
      log.Fatal("invalid normalization method `%s'", *optNormalizeTrack)
    }
    config.NormalizeTrack = *optNormalizeTrack
  }
  if *optFraglenRange != "" {
    tmp := strings.Split(*optFraglenRange, ":")
    if len(tmp) != 2 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
    }
    t1, err := strconv.ParseInt  (tmp[0], 10, 64)
    if err != nil {
      log.Fatalf("parsing fragment length range failed: %v", err)
    }
    t2, err := strconv.ParseInt  (tmp[1], 10, 64)
    if err != nil {
      log.Fatalf("parsing fragment length range failed: %v", err)
    }
    config.FraglenRange[0] = int(t1)
    config.FraglenRange[1] = int(t2)
  }
  if *optFilterPairedEnd && *optEstimateFraglen {
    log.Fatal("cannot estimate fragment length for paired end reads")
  }
  if *optFilterPairedEnd && *optFilterSingleEnd {
    log.Fatal("cannot filter for paired and single end reads")
  }
  config.LogScale          = *optLogScale
  config.PairedAsSingleEnd = *optPairedAsSingleEnd
  config.FilterDuplicates  = *optFilterDuplicates
  config.FilterPairedEnd   = *optFilterPairedEnd
  config.FilterSingleEnd   = *optFilterSingleEnd
  config.SaveFraglen       = *optSaveFraglen
  config.SaveCrossCorr     = *optSaveCrossCorr
  config.SaveCrossCorrPlot = *optSaveCrossCorrPlot

  // parse arguments
  //////////////////////////////////////////////////////////////////////////////
  filenamesTreatment := strings.Split(options.Args()[0], ",")
  filenamesControl   := []string{}
  filenameTrack      := ""
  if len(options.Args()) == 3 {
    filenamesControl = strings.Split(options.Args()[1], ",")
    filenameTrack = options.Args()[2]
  } else {
    filenameTrack = options.Args()[1]
  }

  fraglenTreatment := make([]int, len(filenamesTreatment))
  fraglenControl   := make([]int, len(filenamesControl))
  // split filename:fraglen
  for i, filename := range filenamesTreatment {
    filenamesTreatment[i], fraglenTreatment[i] = parseFilename(filename)
  }
  for i, filename := range filenamesControl {
    filenamesControl[i], fraglenControl[i] = parseFilename(filename)
  }
  if *optFraglen != 0 {
    for i := 0; i < len(fraglenTreatment); i++ {
      fraglenTreatment[i] = *optFraglen
    }
    for i := 0; i < len(fraglenControl); i++ {
      fraglenControl[i] = *optFraglen
    }
  }

  // read genome
  //////////////////////////////////////////////////////////////////////////////
  var genome Genome

  for _, filename := range append(filenamesTreatment, filenamesControl...) {
    g, err := BamImportGenome(filename); if err != nil {
      log.Fatal(err)
    }
    if genome.Length() == 0 {
      genome = g
    } else {
      if !genome.Equals(g) {
        log.Fatal("bam genomes are not equal")
      }
    }
  }

  // fragment length estimation
  //////////////////////////////////////////////////////////////////////////////
  if *optEstimateFraglen {
    for i, filename := range filenamesTreatment {
      if fraglenTreatment[i] != 0 {
        continue
      }
      if fraglen := importFraglen(config, filename, genome); fraglen != -1 {
        fraglenTreatment[i] = fraglen
      } else {
        fraglenTreatment[i] = estimateFraglen(config, filename, genome)
      }
    }
    for i, filename := range filenamesControl {
      if fraglenControl[i] != 0 {
        continue
      }
      if fraglen := importFraglen(config, filename, genome); fraglen != -1 {
        fraglenControl[i] = fraglen
      } else {
        fraglenControl[i] = estimateFraglen(config, filename, genome)
      }
    }
  }

  // bam -> bigWig
  //////////////////////////////////////////////////////////////////////////////
  bamToBigWig(config, filenameTrack, filenamesTreatment, filenamesControl, fraglenTreatment, fraglenControl, genome)
}
