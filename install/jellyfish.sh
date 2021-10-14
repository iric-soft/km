#!/usr/bin/env bash

set -euo pipefail

root=$(dirname $(dirname $(which python3)))

echo "### Installing Jellyfish in ${root}"

link='https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz'
tmpdir=`mktemp -d`
echo TEMPDIR=$tmpdir

pushd $tmpdir
curl -L $link --output jellyfish-2.2.6.tar.gz
tar zxvf jellyfish-2.2.6.tar.gz
cd jellyfish-2.2.6
# General case
./configure --prefix=$root --enable-python-binding
# If your python was installing with MacPort
# ./configure --prefix=$root --enable-python-binding PYTHON_EXTRA_LDFLAGS="-u _PyMac_Error" LDFLAGS="-L/opt/local/lib `python3-config --ldflags` `python-config --libs`"
make -j 4
make install
popd

echo "==> Jellyfish installed"
