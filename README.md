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
Anyway, BUSCO is the 'new kid' (well not so new any more) and works also very well. We'll run it as part of a different session, for now let's stick to CEGMA

CEGMA, once installed, or containerized, is simple to run:
```bash
(user@host)-$ docker run --rm \
-v $(pwd):/in -w /in \
chrishah/cegma:2.5 \
cegma
```

So, let's run BUSCO on our assembly.
