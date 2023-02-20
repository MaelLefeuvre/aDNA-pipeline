#!/usr/bin/env bash

TK_URL="https://github.com/MaelLefeuvre/tkgwv2/archive/refs/heads/develop.zip"


cd $CONDA_PREFIX/bin

wget "${TK_URL}"
unzip $(basename ${TK_URL}) && rm $(basename "${TK_URL}")

chmod +x tkgwv2-develop/TKGWV2.py
chmod +x tkgwv2-develop/TK-helpers.py
chmod +x tkgwv2-develop/scripts/*
chmod +x tkgwv2-develop/helpers/*.R

ln -s tkgwv2-develop/helpers/* .
ln -s tkgwv2-develop/TKGWV2.py
ln -s tkgwv2-develop/TK-helpers.py

# Permanently add to path
#conda develop `pwd`/tkgwv2-master/scripts
#conda develop `pwd`/tkgwv2-master/helpers
