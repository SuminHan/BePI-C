# BePI-C

## Introduction

This program runs on top of [BePI](https://datalab.snu.ac.kr/bepi/) (SIGMOD, 2017), the work done by Jinhong Jung *et al*. I make it public for the researchers who are struggling for implementation of the basic algorithm. Although it requires further improvements such as faster matrix inversion, finer dynamic memory management, matrix visualization, *etc.*, I believe this work can be the cornerstone of the search engine based on network graph, edge prediction, or related fields.

## Setting
The code has been tested on the Ubuntu 18.04 environment, g++.

## Installation

### 1. Paralution
You need to install [Paralution](https://www.paralution.com/) module which enables easier implementation for the parallel computing on calculating matrix operations.

    $ cd /usr/local
    $ wget http://www.paralution.com/downloads/paralution-1.1.0.tar.gz
    $ tar zxvf paralution-1.1.0.tar.gz
    $ cd paralution-1.1.0/src
    $ make all
    $ make install

Then add the path under Linux in the *LD_LIBRARY_PATH* variable.

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/paralution-1.1.0/build/lib

Simple Test:

    $ cd /usr/local/paralution-1.1.0/build/bin
    $ wget ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/laplace/gr_30_30.mtx.gz
    $ gzip -d gr_30_30.mtx.gz 

### 2. Makefile

If you had installed paralution elsewhere, change the following paths as your environment.

Makefile:

    PARALUTION_LIB =/usr/local/paralution-1.1.0/build/lib/libparalution.so
    PARALUTION_INC =-I/usr/local/paralution-1.1.0/build/inc


### 3. Bugs

I may try to catch up the bugs as you issue on this repository.
