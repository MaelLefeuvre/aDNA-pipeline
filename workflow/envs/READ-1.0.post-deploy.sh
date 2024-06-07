#!/usr/bin/env bash

#VERSION="ea376e2"
VERSION="v1.0"

cd $CONDA_PREFIX/bin

git clone https://bitbucket.org/tguenther/read.git

git -C read checkout "${VERSION}"

cp read/READ.py read/READscript.R .

chmod +x READscript.R

rm -rf ./read
