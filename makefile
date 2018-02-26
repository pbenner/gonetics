
SUBDIRS = \
	. \
	tools/bamCheckBin \
	tools/bamGenome \
	tools/bamView \
	tools/bamToBigWig \
	tools/bigWigEditChromNames \
	tools/bigWigExtract \
	tools/bigWigExtractChroms \
	tools/bigWigGenome \
	tools/bigWigMap \
	tools/bigWigPositive \
	tools/bigWigQuantileNormalize \
	tools/bigWigQuery \
	tools/bigWigQuerySequence \
	tools/bigWigStatistics \
	tools/chromHmmTablesToBigWig \
	tools/drawGenomicRegions \
	tools/gtfToBed \
	tools/pwmScanSequences

all:

build:
	@for i in $(SUBDIRS); do \
		echo "Building $$i"; (cd $$i && go build); \
	done

install:
	@for i in $(SUBDIRS); do \
		echo "Installing $$i"; (cd $$i && go install); \
	done

test:
	@for i in $(SUBDIRS); do \
		echo "Testing $$i"; (cd $$i && go test); \
	done
