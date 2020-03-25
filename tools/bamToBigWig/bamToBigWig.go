/* Copyright (C) 2016-2018 Philipp Benner
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
  BWZoomLevels     []int
  SaveFraglen        bool
  SaveCrossCorr      bool
  SaveCrossCorrPlot  bool
  Verbose            int
}

/* i/o
 * -------------------------------------------------------------------------- */

func printStderr(config Config, level int, format string, args ...interface{}) {
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
    log.Fatalf("invalid input file description `%s'", filename)
  }
  return filename, -1
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

  printStderr(config, 1, "Wrote fragment length estimate to `%s'\n", filename)
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
  printStderr(config, 1, "Wrote crosscorrelation to `%s'\n", filename)
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
  printStderr(config, 1, "Wrote cross-correlation plot to `%s'\n", filename)
}

func importFraglen(config Config, filename string) int {
  // try reading the fragment length from file
  basename := strings.TrimRight(filename, filepath.Ext(filename))
  filename  = fmt.Sprintf("%s.fraglen.txt", basename)
  if f, err := os.Open(filename); err != nil {
    return -1
  } else {
    defer f.Close()
    printStderr(config, 1, "Reading fragment length from `%s'... ", filename)
    scanner := bufio.NewScanner(f)
    if scanner.Scan() {
      if fraglen, err := strconv.ParseInt(scanner.Text(), 10, 64); err == nil {
        printStderr(config, 1, "done\n")
        return int(fraglen)
      }
    }
    printStderr(config, 1, "failed\n")
    return -1
  }
}

/* -------------------------------------------------------------------------- */

func main() {

  config  := Config{}
  options := getopt.New()

  // bigWig options
  optBWZoomLevels      := options. StringLong("bigwig-zoom-levels",         0 , "", "comma separated list of BigWig zoom levels")
  // read options
  optShiftReads        := options. StringLong("shift-reads",                0 , "", "shift reads on the positive strand by `x' bps and those on the negative strand by `y' bps [format: x,y]")
  optPairedAsSingleEnd := options.   BoolLong("paired-as-single-end",       0 ,     "treat paired as single end reads")
  optPairedEndStrand   := options.   BoolLong("paired-end-strand-specific", 0 ,     "strand specific paired-end sequencing")
  // options for filterering reads
  optFilterStrand      := options. StringLong("filter-strand",              0 , "", "use reads on either the forward `+' or reverse `-' strand")
  optReadLength        := options. StringLong("filter-read-lengths",        0 , "", "feasible range of read-lengths [format: min:max]")
  optFilterMapQ        := options.    IntLong("filter-mapq",                0 ,  0, "filter reads for minimum mapping quality [default: 0]")
  optFilterDuplicates  := options.   BoolLong("filter-duplicates",          0 ,     "remove reads marked as duplicates")
  optFilterPairedEnd   := options.   BoolLong("filter-paired-end",          0 ,     "remove all single end reads")
  optFilterSingleEnd   := options.   BoolLong("filter-single-end",          0 ,     "remove all paired end reads")
  optFilterChroms      := options. StringLong("filter-chromosomes",         0 , "", "remove all reads on the given chromosomes [comma separated list]")
  // track options
  optBinningMethod     := options. StringLong("binning-method",             0 , "", "binning method [`default' (increment the value of each bin by one " +
                                                                                    "that overlaps a read), `overlap' (increment the value of each bin that " +
                                                                                    "overlaps the read by the number of overlapping nucleotides), or `mean overlap' " +
                                                                                    "(increment the value of each bin that overlaps a read by the mean number of overlapping nucleotides]")
  optBinSize           := options.    IntLong("bin-size",                   0 ,  0, "track bin size [default: 10]")
  optNormalizeTrack    := options. StringLong("normalize-track",            0 , "", "normalize track with the specified method [i.e. rpkm (reads per kilobase " +
                                                                                    "per million mapped reads, i.e. {bin read count}/({total number of reads in millions}*{bin size})), " +
                                                                                    "or cpm (counts per million mapped reads, i.e. {bin read count}/{total number of reads in millions})]")
  optPseudocounts      := options. StringLong("pseudocounts",               0 , "", "pseudocounts added to treatment and control signal [default: `0,0']")
  optSmoothenControl   := options.   BoolLong("smoothen-control",           0 ,     "smoothen control with an adaptive window method")
  optSmoothenSizes     := options. StringLong("smoothen-window-sizes",      0 , "", "feasible window sizes for the smoothening method [format: s1,s2,...]")
  optSmoothenMin       := options. StringLong("smoothen-min-counts",        0 , "", "minimum number of counts for the smoothening method")
  optLogScale          := options.   BoolLong("log-scale",                  0 ,     "log-transform data")
  // options for estimating and setting fragment lengths
  optFraglen           := options.    IntLong("fragment-length",            0 ,  0, "fragment length for all input files [reads are extended to the given length]")
  optFraglenRange      := options. StringLong("fragment-length-range",      0 , "", "feasible range of fragment lengths [format from:to]")
  optFraglenBinSize    := options.    IntLong("fragment-length-bin-size",   0 ,  0, "bin size used when estimating the fragment length [default: 10]")
  optEstimateFraglen   := options.   BoolLong("estimate-fragment-length",   0 ,     "use crosscorrelation to estimate the fragment length")
  optSaveFraglen       := options.   BoolLong("save-fraglen",               0 ,     "save estimated fragment length in a file named <BAM_BASENAME>.fraglen.txt")
  optSaveCrossCorr     := options.   BoolLong("save-crosscorrelation",      0 ,     "save crosscorrelation between forward and reverse strands in a file named <BAM_BASENAME>.fraglen.table")
  optSaveCrossCorrPlot := options.   BoolLong("save-crosscorrelation-plot", 0 ,     "save crosscorrelation plot in a file names <BAM_BASENAME>.fraglen.pdf")
  // generic options
  optVerbose           := options.CounterLong("verbose",                   'v',     "verbose level [-v or -vv]")
  optHelp              := options.   BoolLong("help",                      'h',     "print help")

  options.SetParameters("<TREATMENT1.bam[:FRAGLEN],TREATMENT2.bam[:FRAGLEN],...> [<CONTROL1.bam[:FRAGLEN],CONTROL2.bam[:FRAGLEN],...>] <RESULT.bw>")
  options.Parse(os.Args)

  optionsList := []interface{}{}

  // parse options
  //////////////////////////////////////////////////////////////////////////////
  if *optHelp {
    options.PrintUsage(os.Stdout)
    os.Exit(0)
  }
  if *optVerbose != 0 {
    config.Verbose = *optVerbose
    optionsList = append(optionsList, OptionLogger{log.New(os.Stderr, "", 0)})
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
      optionsList = append(optionsList, OptionBinSize{*optBinSize})
    }
  }
  if *optBinningMethod != "" {
    switch *optBinningMethod {
    case "simple":
    case "default":
    case "overlap":
    case "mean overlap":
    default:
      log.Fatalf("invalid binning method `%s'", *optBinningMethod)
    }
    optionsList = append(optionsList, OptionBinningMethod{*optBinningMethod})
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
    optionsList = append(optionsList, OptionPseudocounts{[2]float64{t1, t2}})
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
    optionsList = append(optionsList, OptionFilterReadLengths{[2]int{int(t1), int(t2)}})
  }
  if *optFilterMapQ < 0 {
      options.PrintUsage(os.Stderr)
      os.Exit(1)
  } else {
    optionsList = append(optionsList, OptionFilterMapQ{*optFilterMapQ})
  }
  if *optFilterStrand != "" {
    switch *optFilterStrand {
    case "+":
      optionsList = append(optionsList, OptionFilterStrand{'+'})
    case "-":
      optionsList = append(optionsList, OptionFilterStrand{'-'})
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
    optionsList = append(optionsList, OptionShiftReads{[2]int{int(t1), int(t2)}})
  }
  if *optSmoothenControl {
    optionsList = append(optionsList, OptionSmoothenControl{true})
  }
  if *optSmoothenSizes != "" {
    smoothenSizes := []int{}
    tmp := strings.Split(*optSmoothenSizes, ",")
    for i := 0; i < len(tmp); i++ {
      t, err := strconv.ParseInt(tmp[i], 10, 64)
      if err != nil {
        log.Fatal(err)
      }
      smoothenSizes = append(smoothenSizes, int(t))
    }
    optionsList = append(optionsList, OptionSmoothenSizes{smoothenSizes})
  }
  if *optSmoothenMin != "" {
    t, err := strconv.ParseFloat(*optSmoothenMin, 64)
    if err != nil {
      log.Fatal(err)
    }
    optionsList = append(optionsList, OptionSmoothenMin{t})
  }
  if *optBWZoomLevels != "" {
    tmp := strings.Split(*optBWZoomLevels, ",")
    bwZoomLevels := []int{}
    for i := 0; i < len(tmp); i++ {
      if t, err := strconv.ParseInt(tmp[i], 10, 64); err != nil {
        log.Fatal(err)
      } else {
        bwZoomLevels = append(bwZoomLevels, int(t))
      }
    }
    config.BWZoomLevels = bwZoomLevels
  }
  if *optNormalizeTrack != "" {
    switch strings.ToLower(*optNormalizeTrack) {
    case "rpkm":
    case "cpm":
    default:
      log.Fatalf("invalid normalization method `%s'", *optNormalizeTrack)
    }
    optionsList = append(optionsList, OptionNormalizeTrack{strings.ToLower(*optNormalizeTrack)})
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
    optionsList = append(optionsList, OptionFraglenRange{[2]int{int(t1), int(t2)}})
  }
  if *optFraglenBinSize > 0 {
    optionsList = append(optionsList, OptionFraglenBinSize{*optFraglenBinSize})
  }
  if *optFilterPairedEnd && *optEstimateFraglen {
    log.Fatal("cannot estimate fragment length for paired end reads")
  }
  if *optFilterPairedEnd && *optFilterSingleEnd {
    log.Fatal("cannot filter for paired and single end reads")
  }
  if *optFilterChroms != "" {
    optionsList = append(optionsList, OptionFilterChroms{strings.Split(*optFilterChroms, ",")})
  }
  optionsList = append(optionsList, OptionEstimateFraglen{*optEstimateFraglen})
  optionsList = append(optionsList, OptionLogScale{*optLogScale})
  optionsList = append(optionsList, OptionPairedAsSingleEnd{*optPairedAsSingleEnd})
  optionsList = append(optionsList, OptionPairedEndStrandSpecific{*optPairedEndStrand})
  optionsList = append(optionsList, OptionFilterDuplicates{*optFilterDuplicates})
  optionsList = append(optionsList, OptionFilterPairedEnd{*optFilterPairedEnd})
  optionsList = append(optionsList, OptionFilterSingleEnd{*optFilterSingleEnd})
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

  // import fragment length
  //////////////////////////////////////////////////////////////////////////////
  for i, filename := range filenamesTreatment {
    if fraglenTreatment[i] == 0 {
      fraglenTreatment[i] = importFraglen(config, filename)
    }
  }
  for i, filename := range filenamesControl {
    if fraglenControl[i] == 0 {
      fraglenControl[i] = importFraglen(config, filename)
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  result, fraglenTreatmentEstimate, fraglenControlEstimate, err := BamCoverage(filenameTrack, filenamesTreatment, filenamesControl, fraglenTreatment, fraglenControl, optionsList...)

  // save fraglen estimates
  //////////////////////////////////////////////////////////////////////////////
  if *optEstimateFraglen {
    for i, estimate := range fraglenTreatmentEstimate {
      filename := filenamesTreatment[i]
      if config.SaveFraglen && estimate.Error == nil {
        saveFraglen(config, filename, estimate.Fraglen)
      }
      if config.SaveCrossCorr && estimate.X != nil && estimate.Y != nil {
        saveCrossCorr(config, filename, estimate.X, estimate.Y)
      }
      if config.SaveCrossCorrPlot && estimate.X != nil && estimate.Y != nil {
        saveCrossCorrPlot(config, filename, estimate.Fraglen, estimate.X, estimate.Y)
      }
    }
    for i, estimate := range fraglenControlEstimate {
      filename := filenamesControl[i]
      if config.SaveFraglen && estimate.Error == nil {
        saveFraglen(config, filename, estimate.Fraglen)
      }
      if config.SaveCrossCorr && estimate.X != nil && estimate.Y != nil {
        saveCrossCorr(config, filename, estimate.X, estimate.Y)
      }
      if config.SaveCrossCorrPlot && estimate.X != nil && estimate.Y != nil {
        saveCrossCorrPlot(config, filename, estimate.Fraglen, estimate.X, estimate.Y)
      }
    }
  }

  // process result
  //////////////////////////////////////////////////////////////////////////////
  if err != nil {
    log.Fatal(err)
  } else {
    printStderr(config, 1, "Writing track `%s'... ", filenameTrack)
    parameters := DefaultBigWigParameters()
    parameters.ReductionLevels = config.BWZoomLevels
    if err := (GenericTrack{result}).ExportBigWig(filenameTrack, parameters); err != nil {
      printStderr(config, 1, "failed\n")
      log.Fatal(err)
    } else {
      printStderr(config, 1, "done\n")
    }
  }
}
