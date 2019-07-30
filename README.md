# Compact Flye

## Bases

This repository contains a compact version for the Third Generation Reads assembler called Flye. Basically, it allows to assemble long prone error reads keeping or in some cases improving (dependant on the k-mer size) the speed while reducing the memory usage.

## Usage

The improvements are due to the use of compact data structures through some of the stages of Flye's original algorithm. This allows us to achieve better results in terms of space and time. Changes are related only in how the internal data of the algorithm is stored and managed, but not in how Flye receives data. Therefore, the way to use it is exactly the same as Flye.

Directory docs has a wide explanation about how to build/install and use flye, thus compactflye aswell. Some easy examples could be:

* flye --pacbio-raw E.coli_PacBio_40x.fasta --out-dir out_pacbio --genome-size 5m --threads 4
* flye --nano-raw Loman_E.coli_MAP006-1_2D_50x.fasta --out-dir out_nano --genome-size 5m --threads 4

A great collection of datasets is available in PacBio's GitHub: https://github.com/PacificBiosciences/DevNet

#### Original GitHub:

https://github.com/fenderglass/Flye

