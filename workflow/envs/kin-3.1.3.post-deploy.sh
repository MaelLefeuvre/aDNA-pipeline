#!/usr/bin/env bash

VERSION="3.1.3"
REPO_URL="https://github.com/DivyaratanPopli/Kinship_Inference/archive/refs/tags/v${VERSION}.tar.gz"
TMP_SOURCE="Kinship_Inference-${VERSION}"
CWD=`pwd`

wget -O - ${REPO_URL} | tar -xvzf - 
cd ${TMP_SOURCE} &&
pip3 install pypackage/kingaroo &&
pip3 install pypackage/kin

cd $CWD 
rm -r ${TMP_SOURCE}
