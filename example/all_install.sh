
## Setup your software directory
mkdir -p $HOME/software
cd $HOME/software

## Create virtual environment
virtualenv $HOME/.virtualenvs/km
source $HOME/.virtualenvs/km/bin/activate

## Download and install jellyfish
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz
tar zxvf jellyfish-2.2.6.tar.gz
cd jellyfish-2.2.6
./configure --prefix=$HOME/.virtualenvs/km --enable-python-binding
make -j 4
make install

## Download and install km
cd ..
# TODO: replace this part with a relase when it's ready
git clone https://github.com/iric-soft/km.git
cd km
python setup.py install

## Execute km on a small example
# Need to reload the virtual environment each time you open a new terminal
# with: source $HOME/.virtualenvs/km/bin/activate
km find_mutation ./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa ./data/jf/02H025_NPM1.jf | km find_report -t ./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa
