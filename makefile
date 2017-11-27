
SUBDIRS = \
	. \
	tools/bamToBigWig \
	tools/bigWigEditChromNames \
	tools/bigWigExtract \
	tools/bigWigExtractChroms \
	tools/bigWigGenome \
	tools/bigWigMap \
	tools/bigWigPositive \
	tools/bigWigQuery \
	tools/bigWigQuerySequence \
	tools/bigWigStatistics \
	tools/chromHmmTablesToBigWig \
	tools/drawGenomicRegions \
	tools/gtfToBed

all:

test:
	@for i in $(SUBDIRS); do \
		echo "Testing $$i"; (cd $$i && go test); \
	done
