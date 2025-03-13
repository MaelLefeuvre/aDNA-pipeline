#!/usr/bin/env bash

VERSION="v0.3.2"
REPO_URL="https://github.com/MaelLefeuvre/grups-rs"

REPO="$CONDA_PREFIX/$(basename ${REPO_URL%.git})"
BUILD_DIR="${REPO}/target"

cd $CONDA_PREFIX
git clone --recursive $REPO_URL
cd $REPO
git checkout "$VERSION"

export CARGO_INCREMENTAL=0
export CC="x86_64-conda-linux-gnu-gcc"
CMAKE_C_FLAGS="-I./ -L./ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" \
RUSTFLAGS="-Ctarget-cpu=native" cargo install --path . --root $CONDA_PREFIX

R --slave -e 'devtools::install("./grups.plots")'

