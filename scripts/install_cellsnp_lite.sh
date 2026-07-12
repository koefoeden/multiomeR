#!/usr/bin/env bash
set -euo pipefail

version=1.2.3
ref=8934b0d58b5d0cec5ab2d9576f0ebba689b619f1

if [[ $(cellsnp-lite --version 2>&1) == *"$version"* ]]; then
    echo "cellsnp-lite $version is already installed"
    exit 0
fi

build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT

curl -fsSL "https://github.com/single-cell-genetics/cellsnp-lite/archive/$ref.tar.gz" \
    | tar -xz -C "$build_dir"
cd "$build_dir/cellsnp-lite-$ref"
if [[ $(uname -s) == Darwin ]]; then
    sed -i.bak 's/#define _POSIX_C_SOURCE 200809L/#define _DARWIN_C_SOURCE/' src/thpool.c
    sed -i.bak 's/cellsnp_lite_LDADD = @HTSLIB_LIB@/cellsnp_lite_LDADD = @HTSLIB_LIB@ -ldeflate/' Makefile.am
fi
autoreconf -iv
./configure --prefix="$CONDA_PREFIX" --with-htslib="$CONDA_PREFIX"
make
make install
