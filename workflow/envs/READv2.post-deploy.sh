#!/usr/bin/env bash


#eval "$(conda shell.bash hook)"
#conda activate $(basename "$(readlink -f ${BASH_SOURCE[0]%.post-deploy.sh})")

cd $CONDA_PREFIX/bin

git clone https://github.com/GuntherLab/READv2.git
cd READv2
git checkout 'READv2.00'

cd $CONDA_PREFIX/bin
cp READv2/READ2.py ./ 
chmod +x READ2.py

ln -s READ2.py READ2
rm -rf ./READv2
