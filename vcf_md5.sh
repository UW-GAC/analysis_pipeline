#! /bin/bash

cd $1
md5deep -b *.vcf.gz > md5sum.txt
