# post-assembly-intro
Tutorial for common steps post denovo genome assembly


***Disclaimer***
To follow the demo and make the most of it, it helps if you have some basic skills with running software tools and manipulating files using the Unix shell command line. It assumes you have Docker installed on your computer (tested with Docker version 18.09.7, build 2d0083d; on Ubuntu 18.04).

## Introduction
Once you have assembled a genome you usually want to assess the quality based on certain metrics, such as contiguity, completeness, perhaps even evaluate if your assembly is contaminated with sequences from 'non-target' organisms.

In the following I'll demonstrate a few common steps. The data required for these comes with the repository - it's in the `data/` directory. Ideally you would start by cloning/downloading this repository, e.g. using `git`, like so:

```bash
(user@host)-$ git clone https://github.com/chrishah/post-assembly-intro.git
```

Which will get you a directory that you can move into and follow along with the exercise within it.

```bash
(user@host)-$ cd post-assembly-intro
```
Checkout what we have in the `data/` directory:

```bash
(user@host)-$ ls data
```
From the file names/extensions you might already have a good idea about what most of these files may be (that's the intention). I'll get to the details in little while.

## Post Assembly

__1.) Assess assembly contiguity__

At this stage you have probably heard about the most common assembly stats for describing contiguity - N50, etc.

[Quast](http://quast.sourceforge.net/quast) is a very convenient tool for calculating these stats.

We have a toy genome assembly in your `data/` directory - let's run quast on it.

```bash
docker run --rm \
-v $(pwd):/in -w /in \
reslp/quast:5.0.2 \
quast.py data/genome_assembly.fasta
```

This produces a directory `quast_results` containing the stats. You can have a quick look at those in a simple text file.
```bash
cat quast_results/latest/report.txt
```
Or check out the html version of the report (`quast_results/latest/report.html`) which you can open in your web browser.

Quast gives a whole lot of other useful things but for these I refer you to the quast [homepage](http://quast.sourceforge.net/quast). 

__2.) Assess assembly completeness__

The idea here is that there exists a set of conserved genes that are expected to be present in the genomes of most organisms. Since they are conserved they should be reasonably easy to identify using some form of similarity search to reference sequences for these genes. The percentage of these conserved genes identified in the genome you have assembled is then used as a proxy for the overall completeness of the gene space in the genome.

There are currently two tools available:
 - [BUSCO](https://busco.ezlab.org/)
 - [CEGMA](http://korflab.ucdavis.edu/datasets/cegma/)

The latter is unfortunately not maintained any more (some [history](http://www.acgt.me/blog/2015/5/18/goodbye-cegma-hello-busco)), but it can still be used if you can get it installed on your system. I have made a docker container for it, since I am planning to keep using it - see [here](https://hub.docker.com/r/chrishah/cegma).
Anyway, BUSCO is the 'new kid' (well, not so new any more) and works also very well. 

Now, let's run BUSCO (will take about 15 minutes):
```bash
(user@host)-$ docker run \
              --rm -v $(pwd):/in -w /in \
              ezlabgva/busco:v5.2.1_cv1 \
              busco -i data/genome_assembly.fasta \
              -o busco -m genome -l eukaryota \
              -c 1 \
              --augustus --augustus_species schistosoma
```
If you ran BUSCO as above it will have created one directory called `busco` (because you said `-o busco` above). 

A detailed exlanation of the parameters and the BUSCO output you can also find as part of a different [session](https://github.com/chrishah/phylogenomics-intro).

Usually, the most relevant files are:
 - `busco/run_eukaryota_odb10/short_summary.txt`
 - `busco/run_eukaryota_odb10/full_table.tsv`

and fasta files in the directory:
 - `busco/run_eukaryota_odb10/busco_sequences/`

CEGMA, once installed, or containerized, is simple to run (it has a lot of options that you can explore in your own time) - takes a while though:

***ATTENTION***
> The next step (`cegma`) will run for about an hour, so if you are in a rush, you can also skip this and look at the output that we have deposited with the repo (see below).

```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/cegma:2.5 \
cegma --threads 1 -g data/genome_assembly.fasta
```
While it is running we can skip to the next part and talk about mapping reads to genomes.

Once it's done the thing you want to be looking at is the CEGMA report in `output.completeness_report`. A gff file with the genes cegma has predicted can be found at `output.cegma.gff`.

Note that if for some reason you want to skip running CEGMA we have an example output deposited for you as part of this repo at: `data/outputs/cegma/output.completeness_report`.

Let's have a look at the completeness report:
```bash
(user@host)-$ cat output.completeness_report
```

__3.) Mapping reads to genomes__

Read mapping is covered by many online tutorials, so I'll just show you how it could be done with a tool called [bowtie-2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

First index your genome file.
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
reslp/bowtie2:2.3.5 \
bowtie2-build data/genome_assembly.fasta my_genome.index -q
```
Check out which new files have been created.
```bash
(user@host)-$ ls -hlrt
```

Then, map the reads to the indexed genome.
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
reslp/bowtie2:2.3.5 \
bowtie2 -1 data/reads.1.fastq.gz -2 data/reads.2.fastq.gz --threads 2 -q --phred33 --fr -x my_genome.index -S my_mapped_reads.sam
```
This has produced a file called `my_mapped_reads.sam`. This is simple text file formatted in [sam](https://samtools.github.io/hts-specs/SAMv1.pdf) format, that contains information on where in the genome certain reads mapped (if at all).

You can look into the file, if you dare.. - just the first 1000 lines.
```bash
(user@host)-$ head -n 1000 my_mapped_reads.sam
```
Since SAM is just a text file and for large amounts of data these files may get very big the developers have established a binary (non-human readable) version of SAM, which is called BAM.  

The next step will be to convert the SAM to a BAM file.
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
reslp/samtools:1.9 \
samtools view -bS my_mapped_reads.sam -o my_mapped_reads.bam -@ 2
```

Check out the size of the newest file `my_mapped_reads.bam` as compared to the original SAM file. Note that we have not lost any information - we have just compressed the data.
```bash
(user@host)-$ ls -hlrt
```

Two more steps that are usually being done are sorting and indexing the bam file - this is the convention, and what most downstream tools expect.
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
reslp/samtools:1.9 \
samtools sort -o my_mapped_reads.sorted.bam my_mapped_reads.bam -@ 2

(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
reslp/samtools:1.9 \
samtools index my_mapped_reads.sorted.bam -@ 2

```

A common step that I want to at least mention is the removal of duplicates. [Picard](https://broadinstitute.github.io/picard/) offers a good option there.

```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
broadinstitute/picard:2.20.6 \
java -jar /usr/picard/picard.jar MarkDuplicates \
INPUT=my_mapped_reads.sorted.bam OUTPUT=my_mapped_reads.sorted.duprmvd.bam METRICS_FILE=my_mapped_reads.sorted.duprmvd.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
```

I encourage you to inspect the assembly and reads mapping to it visually. A possible tools is the Integrative Genomics Viewer [igv](https://software.broadinstitute.org/software/igv/). You can install it at some point, but for now we're going to use through a webapp - go [here](https://igv.org/app/).

First we need to index our reference genome.
Index the genome for viewing.
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
reslp/samtools:1.9 \
samtools faidx data/genome_assembly.fasta
```

If you followed until here you're going to need the following files downloaded locally:
 - `data/genome_assembly.fasta` (came with the repo)
 - `data/genome_assembly.fasta.fai` (index prodcued in the last step above)
 - `my_mapped_reads.sorted.bam` (sorted bam file from above)
 - `my_mapped_reads.sorted.bam.bai` (index of sorted bam file from above)

In the web [app](https://igv.org/app/), go *Genome->Local File->*, then make sure to select both files, so `genome_assembly.fasta` together with `genome_assembly.fasta.fai`.

Now, to get the reads visible, go to *Tracks->Local File* and select `my_mapped_reads.sorted.bam` and `my_mapped_reads.sorted.bam.bai` again at the same time before clicking *Open*.

***ATTENTION***
> Note that not all contigs will have reads mapping to it, because this is a reduced set of reads. You can select for example conig `NODE_288_length_24172_cov_71.622494` in IGV for exploration.


Enjoy!

__4.) Blobtools__

A nice tool for assessing contamination in your genome assembly is [blobtools](https://blobtools.readme.io/docs). It summarizes aggregate properties of the assembly (GC content, coverage of each contig/scaffold) which sometimes reveals interesting patterns in assemblies.

The basic things you need are
 - an assembly - fasta file
 - information about coverage

Blobtools can extract coverage information from bam files (gladly we made one above). It also can parse the needed information directly from the fasta headers, if you have used certain assemblers, e.g. Platanus or SPAdes. Our assembly was done with SPADes, so we can try that.

Blobtools needs to be run in three steps - do consult the manual on the Blobtools [webpage](https://blobtools.readme.io/docs) to get more info on what the individual steps are doing.
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools create -i data/genome_assembly.fasta -y spades -o blobtools_spades

(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools view -i blobtools_spades.blobDB.json

(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools plot -i blobtools_spades.blobDB.json
```

The file you want to look at first of all is: `blobtools_spades.blobDB.json.bestsum.phylum.p8.span.100.blobplot.spades.png`, but there is lots more to explore on your own.

Now, let's assume you hadn't used SPAdes as your assembler, you can still use blobtools, but in this case you need to give the coverage information in a different way, e.g. a bam file.

```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools create -i data/genome_assembly.fasta -b my_mapped_reads.sorted.bam -o blobtools_bam

(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools view -i blobtools_bam.blobDB.json

(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools plot -i blobtools_bam.blobDB.json

```
Checkout `blobtools_bam.blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png` and `blobtools_bam.blobDB.json.bestsum.phylum.p8.span.100.blobplot.read_cov.bam0.png` - there shouldn't be much difference to the first result.

Finally, one of the nicest features of blobtools is that the visualizations can be taxonomically annotated - see [here](https://blobtools.readme.io/docs/taxonomic-annotation). What you'll need is a so-called 'hits files'. This is essentially a text files obtained via comparing the assembly against a reference database using `blast` or other tools - see [here](https://blobtools.readme.io/docs/taxonomy-file).

So, you could download the entirety of NCBI's nt (nucleotide) database. BLAST you assembly against it and use the info you get to annotate you blobs.

***Attention***
> Do not do this as part of the course (if you are in one right now), unless you are specifically asked. The next steps will involve downloading (currently) some 100GB worth of data and subsequently a BLAST search that might take several days, if not parallelized in a smart way.
 
So, for completeness sake, you could download the entire `nt` database and decompress it.
```bash
(user@host)-$ mkdir db
(user@host)-$ cd db
(user@host)-$ wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz"
(user@host)-$ for a in nt.*.tar.gz; do tar xzf $a; done
(user@host)-$ cd ..
```

Then, you would use BLAST to compare your genome against the database - the blobtools people suggest a certain way to use blast for this. You could change the params however you see fit of course, but the main thing is that you get the output in the right format.
```bash
(user@host)-$ ASSEMBLY=data/genome_assembly.fasta
(user@host)-$ DB=db/nt
(user@host)-$ blastn \
 -query $ASSEMBLY \
 -db $DB \
 -outfmt "6 qseqid staxids bitscore std" \
 -max_target_seqs 1 \
 -max_hsps 1 \
 -evalue 1e-25 \
 -num_threads 10 \
 -out blastn.fmt6.out.txt
```

***continue here, as part of the course***


An example file comes with the repository - check it out.
```bash
(user@host)-$ cat data/blastn.fmt6.out.txt
```


Now, try to give this additional info to blobtools. Last thing you need is a so called *taxdump*, this is a textfile that contains the taxonomic information for all entries in the Genbank databases. You can download it and unpack.
```bash
(user@host)-$ wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
(user@host)-$ tar xzfv taxdump.tar.gz
```
From this you get a file `nodes.dmp` and a `names.dmp` file - these you need in the next step.

```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools create -i data/genome_assembly.fasta -b my_mapped_reads.sorted.bam \
--nodes nodes.dmp --names names.dmp --hitsfile data/blastn.fmt6.out.txt -o blobtools_tax

(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools view -i blobtools_tax.blobDB.json

(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools plot -i blobtools_tax.blobDB.json

```

What you want to look at initially is:
 - `blobtools_tax.blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png`
 - `blobtools_tax.blobDB.json.bestsum.phylum.p8.span.100.blobplot.read_cov.bam0.png`

Nice, no?

Now what to do with this info?

You could for example take all reads that contributed to contigs that were classified as 'Chordata', and reassemble them, if that happens to be your target.

First get all contig/scaffold ids that were classified as Chordata.
```bash
(user@host)-$ grep "Chordata" blobtools_tax.blobDB.table.txt | cut -f 1 > Chordata.list.txt
```

Then, use another tool from the blobtools suite (see [here](https://blobtools.readme.io/docs/bamfilter)) to extract the relevant reads from the original bam file.
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/blobtools:v1.1.1 \
blobtools bamfilter -b my_mapped_reads.sorted.bam -i Chordata.list.txt --read_format fq --noninterleaved -o reads_Chordata
```

Thanks for joining us today!

If you have any questions, comments, feedback (good OR bad), let me know!

__What do you think?__

# Contact
Christoph Hahn - <christoph.hahn@uni-graz.at>



