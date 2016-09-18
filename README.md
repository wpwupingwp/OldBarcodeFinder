# BarcodeFinder

## What is BarcodeFinder? 

BarcodeFinder is a BLAST based program written by python3, which could find
single-copy barcode sequence of given data in reasonable time.

The main features of BarcodeFinder contains:

1. **Output barcode sequence of given data**. Till now, it could find out barcode
   of given multiple genome sequences of two bacteria organism. Notice that
   the program was designed for higher taxonomic class and bigger organism
   number, which would be available in followed version after more test.

2. **Only generate single-copy barcode**. A nice barcode should be single-copy among
   genome. Herein, BarcodeFinder could automatically drop multi-copy candidate
   barcodes which lead to less result (but more accurate). 

3. **Appropriate resolution of barcode**. With given barcode, species come from different higher taxonomy grade could be neatly divided while remain same group in lower level. However, among lower grade, subgroups are possible if the barcode has enough mutation to cause distinguishable distinction which has no effect on higher level.

4. **Reasonable resource consumption**. Most of CPU time, memory and hard disk space were used by BLAST while the number is unexpectedly small. For instance, with **less than half minutes**, **~500mb memory** and **twice hard disk usage for given data size**, 4 barcodes were successfully generated from **~700mb input data**.

## Requirements

### Hardware

It depends on input data. If less than 1gb, the program could normally run on mainstream PC with less than 1 minute. Bigger data need more memory (not test).

### Software
1. [python3] (https://www.python.org/downloads/)

Be sure to install python3 rather than python 2.7.

2. [biopython] (http://biopython.org/wiki/Download)

3. [BLAST Suite] (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

4. [mafft] (http://mafft.cbrc.jp/alignment/software/)
　　
For Windows users, please make sure:

1. correct version and platform (x86/x64) for best performance. 
2. Biopython has lots of version,  its python version must be same with python3’s you installed before. For example, you installed python 3.5.1 with x86_64 architecture (depends on your CPU, most of them on PC now support x86_64 mode). Then you have to install biopython for python 3.5.1 with same architecture.
3. It is better to choose **all-in-one** package for mafft.

For Linux users, just run:

>sudo apt-get install python3-biopython ncbi-blast+ mafft

Or

>sudo yum install python3-biopython ncbi-blast+ mafft

### Operating system

Although technically BarcodeFinder and its dependent software could run on various operating system, till now, only Linux and Microsoft Windows were tested.

If you obtained a binary distribution you don't need to compile Vim.  If you
obtained a source distribution, all the stuff for compiling Vim is in the
`src` directory.  See `src/INSTALL` for instructions.


## Installation

1. Download BarcodeFinder.py
2. Finished

## Usage

Put BarcodeFinder.py into same directory of input fasta files. Then run the command:

For Microsoft Windows user:

>python barcode.py 

For Linux user:

>python3 barcode.py

For details, run **”python barcode.py –h”** or **”python barcode.py –h”**. You can change the default parameters of the program according to specific requirement.

There are other `README_*.txt` files, depending on the distribution you used.

## Copying ##

All rights reserved (**temporarily**). The license will change to GPL after the program becomes mature.

## Information ##

The latest version of BarcodeFinder can be found on [github] (https://github.com/wpwupingwp/201609).

If you have problems or any other questions, feel free to send email to the [author] (wuping@ibcas.ac.cn).

## Main author 

Send any other comments, patches, flowers and suggestions to:
	[Ping Wu] (wuping@ibcas.ac.cn).
