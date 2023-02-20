#!/usr/bin/env bash

cd $CONDA_PREFIX/bin

git clone https://bitbucket.org/tguenther/read.git

cp read/READ.py read/READscript.R .

chmod +x READscript.R

rm -rf ./read
