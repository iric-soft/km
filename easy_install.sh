#!/usr/bin/env bash

set -e

# NEED:
# - pip : https://pip.pypa.io/en/stable/installing/
# - venv: https://virtualenv.pypa.io/en/stable/installation/

# Change INSTALL_DIR path, if you don't want to use $HOME
if [ -z "$INSTALL_DIR" ]
then
    INSTALL_DIR=$HOME
fi

INSTALL_DIR=$(realpath $INSTALL_DIR)

echo "### Installing relative to: $INSTALL_DIR ... ###"

## Create virtual environment
echo "### Creating virtual environment ... ###"
python -m venv $INSTALL_DIR/.virtualenvs/km
source $INSTALL_DIR/.virtualenvs/km/bin/activate
echo "==> Virtual environment created"

## Download and install km
pip install km-walk && echo "==> km installed"

## Execute km on a small example

# Need to reload the virtual environment each time you open a new terminal
# with: source $INSTALL_DIR/.virtualenvs/km/bin/activate
echo "### Testing km ... ###"

SCRATCH=$(mktemp -d)

if [ ! -d $SCRATCH ]
then
    echo "Directory creation was unsuccessful. Exiting..." >&2
    exit 1
fi

trap 'rm -fr "$SCRATCH"' EXIT

cd $SCRATCH

## Clone km in a temporary location
git clone https://github.com/iric-soft/km.git

km find_mutation km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa km/data/jf/02H025_NPM1.jf | km find_report -t km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa

km find_mutation km/data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa km/data/jf/03H116_ITD.jf | km find_report -t km/data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa

km find_mutation km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa km/data/jf/02H025_NPM1.jf -g
