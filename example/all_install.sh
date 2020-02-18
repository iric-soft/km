# NEED:
# - pip : https://pip.pypa.io/en/stable/installing/

## Setup your software directory
mkdir -p $HOME/software
cd $HOME/software

## Create virtual environment
python -m venv $HOME/.virtualenvs/km
source $HOME/.virtualenvs/km/bin/activate

## Download and install jellyfish
curl -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz --output jellyfish-2.2.6.tar.gz
tar zxvf jellyfish-2.2.6.tar.gz
cd jellyfish-2.2.6
# General case
./configure --prefix=$VIRTUAL_ENV --enable-python-binding
# For user which have some difficulty to install python binding on ios 
# ./configure --prefix=$HOME/.virtualenvs/km --enable-python-binding PYTHON_EXTRA_LDFLAGS="-u _PyMac_Error"  LDFLAGS="`python-config --ldflags` `python-config --libs`"
make -j 4
make install

## Download and install km
cd $HOME/software
# TODO: replace this part with a release when it's ready
git clone https://github.com/iric-soft/km.git
cd km
python setup.py install

## Execute km on a small example
# Need to reload the virtual environment each time you open a new terminal
# with: source $HOME/.virtualenvs/km/bin/activate
km find_mutation $HOME/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $HOME/software/km/data/jf/02H025_NPM1.jf | km find_report -t $HOME/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa
# run with graphical option
# km find_mutation $HOME/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $HOME/software/km/data/jf/02H025_NPM1.jf -g
