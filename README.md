# hash10x
A toolset for efficient analysis of 10X Genomics linked read data sets, in particular for de novo assembly.

The premise of hash10x is that there is a large amount of information present in an unmapped 10X data set that is not used by supernova, or easily accesible via it (for more background information on 10X, supernova etc. see below).  The core idea is to construct a set of minhash kmer hashes from all the reads, and to build data structures with sorted hashes for each barcode and sorted barcodes for each hash.  From this we aim to do three main things:
1. cluster the hashes within each barcode so that each cluster comes from one molecule,
2. provide estimates of sequence depth (and hence genome size), molecule depth (and hence fraction of sequence coverage per molecule), molecule length distribution, etc.
3. for each hash identify whether the corresponding kmer is an *error* not present in the source genome, *heterozygous* only present in one copy in the diploid genome, *homozygous* present in two copies in the source genome that are homologous, or *multicopy* present at more than one distinct location in the source genome.
Downstream, I would like to use the heterozygous and homozygous kmers to classify candidate long read overlaps prior to long read assembly.  Also there is potential for building new direct assembly maps from the hashes.

hash10x is a work in progress. To get and build the current version:
```
git clone https://github.com/richarddurbin/hash10x.git
cd hash10x
make
```
This should make two executables, hash10x (the main program) and fq2b, which is used to turn the initial fastq data set into a compact binary record-based format that can be sorted on key by [bsort](https://github.com/pelotoncycle/bsort) for import into hash10x.

*Note: For real, vertebrate genome-sized data sets it is good to compile a multi-threaded version using OMP. This parallelizes the cluster command. To do this you need a -DOMP compiler option.  I haven't set this all up in the Makefile yet.*

To run:
```
fq2b -10x goodcodes -o L1.fqb DATA_S1_L001_R1_001.fastq.gz DATA_S1_L001_R2_001.fastq.gz
fq2b -10x goodcodes -o L2.fqb DATA_S1_L002_R1_001.fastq.gz DATA_S1_L002_R2_001.fastq.gz
etc.
cat L*.fqb > DATA.fqb           // make binary file packing 4 bases per byte, 2*(40 bytes sequence + 20 bytes 1-bit qual)
bsort -k 4 -r 120 DATA.fqb      // sort 120 byte binary records based on first 4 bytes (16 bases) as key
hash10x --readFQB DATA.fqb --writeHash DATA.hash    // construct minhashes and save data structures in binary .hash format
hash10x --readHash DATA.hash --hashStats DATA.hashStats -o DATA.codeStats --codeStats    // some histograms and summary stats
hash10x --readHash DATA.hash --hashDepthRange 30 100 --cluster 1 0 --clusterSplit --writeHash DATA.30-100.hash
// cluster hashes with depth between 30 and 100 into candidate molecules, map these to new "barcodes" and save
```
Entering just `hash10x` without arguments lists possible commands.  Any sequence of commands with their arguments can be chained within a single invocation of hash10x.

The `goodcodes` file used above should be a clean list of barcodes, one per line, either from a white list provided by 10X or, as in the yeast example below, obtained by identifying all the initial 16mers present at least 100 times.

## Example yeast simulation

To test hash10x I use a simulated data set made with [lrsim](https://github.com/aquaskyline/LRSIM) from two yeast genomes.  Where a diploid reference genome is known as in this case there are commands to build a "crib" which allows evaluation of the accuracy of clustering and potentially other operations.  Here are relevant commands:
```
wget ftp://ftp.sanger.ac.uk/pub/rd/S288c.fa.gz
wget ftp://ftp.sanger.ac.uk/pub/rd/SK1.fa.gz
gunzip *.fa.gz
perl ../LRSIM/simulateLinkedReads.pl -g S288c.fa,SK1.fa -p ./yeast -x 5 -f 50 -t 10 -m 10 -o
gzip -dc yeast_S1_L001_R1_001.fastq.gz yeast_S1_L002_R1_001.fastq.gz | perl -ne 'if (++$i%4 == 2) { print substr($_,0,16) ; print "\n" ;}' | sort | uniq -c | awk '($1 >= 100) { print $2 }' > goodcodes
fq2b -10x goodcodes -o L1.fqb yeast_S1_L001_R1_001.fastq.gz yeast_S1_L001_R2_001.fastq.gz
fq2b -10x goodcodes -o L2.fqb yeast_S1_L002_R1_001.fastq.gz yeast_S1_L002_R2_001.fastq.gz
cat L1.fqb L2.fqb > yeast.fqb
bsort -k 4 -r 120 yeast.fqb
hash10x -B 24 --readFQB yeast.fqb --writeHash yeast.hash -o yeast.hashStats --hashStats
hash10x -B 24 --readHash yeast.hash --hashDepthRange 30 100 --cluster 1 10 --cribBuild S288c.fa SK1.fa --clusterReport 1 10 // clusters codes 1..10
hash10x -B 24 --readHash yeast.hash --hashDepthRange 30 100 --cluster 1 0 --writeHash yeast.30-100.hash -o yeast.cluster30-100.codeStats --codeStats
```
This generates 2.5M read pairs (151+151bp) over ~10k barcodes, each with 10-20 ~50kb molecules, for an ~60x average depth data set.  Most hash10x operations for this case take seconds, with the clustering taking around 5 minutes single threaded on a MacBook Pro 2.5 GHz Intel Core i7.

The `-B 24` option uses a smaller hash table size (24 bits = 16M) for minhashes from this small (12Mb) genome - we end up storing ~34M hashes in 2.3M bins.  Default is 28 bits for 256M hash bins, appropriate for a 500Mb to 1Gb genome.  For larger genomes, e.g. mammalian, you could explore going up to `-B 30`.  Maximum is `-B 34`.

I could clean up the generation of goodcodes so as to work from a binary file, which would speed things up.  e.g. aiming for the following
```
fq2b -o L1.fqb yeast_S1_L001_R1_001.fastq.gz yeast_S1_L001_R2_001.fastq.gz
fq2b -o L2.fqb yeast_S1_L002_R1_001.fastq.gz yeast_S1_L002_R2_001.fastq.gz
cat L1.fqb L2.fqb > yeast-raw.fqb
fq2b -10xThresh 100 -o yeast.fqb yeast.fqb  // find barcodes present at least 100 times and fix to them, printing out stats
```

## Background

10X Genomics linked read data consist of sets of barcoded read pairs, where the barcode assigns the reads to one or more "GEM" droplets, each of which contained some number of long DNA molecules. Reads with the same barcode typically come from tens to hundreds of molecules, each of which might be several tens of kilobases long (up to several hundred kb), but with relatively low read coverage per molecule, e.g. 10%.  There are usually of the order of 1 million distinct barcodes, and a good data set for assembly purposes should be at least 50x deep. 10X Genomics as a company provides extensive tools (lariat, longranger, loupe etc.) for analysing data sets against a reference, but for de novo sequencing there is only an assembler supernova, which, while it can give impressive results for some genomes, is quite heavy to run.

In detail, the first 16bp of the first read is the barcode, with the next 7bp being spacer.
