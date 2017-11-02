
SUBDIRS = \
	. \
	tools/bigWigExtractChroms \
	tools/bigWigQuerySequence \
	tools/chromHmmTablesToBigWig \
	tools/bigWigGenome \
	tools/bigWigEditChromNames \
	tools/gtfToBed \
	tools/bigWigQuery \
	tools/bamToBigWig

all:

test:
	@for i in $(SUBDIRS); do \
		echo "Testing $$i"; (cd $$i && go test); \
	done
