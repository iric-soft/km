# NEED:
# - pip : https://pip.pypa.io/en/stable/installing/
# - virtualenv: https://virtualenv.pypa.io/en/stable/installation/

## Setup your software directory
mkdir -p $HOME/software
mkdir -p $HOME/archive
cd $HOME/archive

## Create virtual environment
virtualenv $HOME/.virtualenvs/km
source $HOME/.virtualenvs/km/bin/activate

## Download and install jellyfish
VERSION=2.2.10
curl -L https://github.com/gmarcais/Jellyfish/releases/download/v$VERSION/jellyfish-$VERSION.tar.gz --output jellyfish-$VERSION.tar.gz
tar zxvf jellyfish-$VERSION.tar.gz
cd jellyfish-$VERSION
# General case
./configure --prefix=$VIRTUAL_ENV --enable-python-binding
# For users who wish to install python with macport
# ./configure --prefix=$HOME/.virtualenvs/km --enable-python-binding PYTHON_EXTRA_LDFLAGS="-u _PyMac_Error"  LDFLAGS="-L/opt/local/lib `python-config --ldflags` `python-config --libs`"
make -j 4
make install

## Download and install km
cd $HOME/archive
# TODO: replace this part with a release when it's ready
git clone https://github.com/iric-soft/km.git
cd km
python setup.py install

## Execute km on a small example
# Need to reload the virtual environment each time you open a new terminal
# with: source $HOME/.virtualenvs/km/bin/activate
km find_mutation $HOME/archive/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $HOME/archive/km/data/jf/02H025_NPM1.jf | km find_report -t $HOME/archive/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa
# run with graphical option
# km find_mutation $HOME/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $HOME/software/km/data/jf/02H025_NPM1.jf -g
