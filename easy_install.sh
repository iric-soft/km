# NEED:
# - pip : https://pip.pypa.io/en/stable/installing/
# - virtualenv: https://virtualenv.pypa.io/en/stable/installation/

# Change INSTALL_DIR path, if you don't want to use $HOME
INSTALL_DIR=$HOME

echo "### Setup install directory ... ###"
## Clean old version of previews try
rm -rf $INSTALL_DIR/software/km
rm -rf $INSTALL_DIR/software/jellyfish-2.2.6
## Setup your software directory
mkdir -p $INSTALL_DIR/software
cd $INSTALL_DIR/software
echo "==> Done!"

## Create virtual environment
echo "### Create virtual environment ... ###"
python -m venv $INSTALL_DIR/.virtualenvs/km
source $INSTALL_DIR/.virtualenvs/km/bin/activate
echo "==> Virtual environment created"

## Download and install jellyfish
echo "### Install Jellyfish ... ###"
curl -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz --output jellyfish-2.2.6.tar.gz
tar zxvf jellyfish-2.2.6.tar.gz
cd jellyfish-2.2.6
# General case
./configure --prefix=$VIRTUAL_ENV --enable-python-binding
# make clean
# For user which install python with macport
# ./configure --prefix=$INSTALL_DIR/.virtualenvs/km --enable-python-binding PYTHON_EXTRA_LDFLAGS="-u _PyMac_Error"  LDFLAGS="-L/opt/local/lib `python-config --ldflags` `python-config --libs`"
make -j 4 && make install && echo "==> Jellyfish installed"

## Download and install km
cd $INSTALL_DIR/software
# TODO: replace this part with a release when it's ready
echo "### Install km ... ###"
git clone https://github.com/iric-soft/km.git
cd km
# crucial for earlier versions of python (e.g. 3.5)
$VIRTUAL_ENV/bin/pip install pip setuptools wheel --upgrade
pip install . --use-feature=in-tree-build && echo "==> km installed"

## Execute km on a small example
# Need to reload the virtual environment each time you open a new terminal
# with: source $INSTALL_DIR/.virtualenvs/km/bin/activate
echo "### Run km test ... ###"
km find_mutation $INSTALL_DIR/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $INSTALL_DIR/software/km/data/jf/02H025_NPM1.jf | km find_report -t $INSTALL_DIR/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa
km find_mutation $INSTALL_DIR/software/km/data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa $INSTALL_DIR/software/km/data/jf/02H025_NPM1.jf -g
km find_mutation $INSTALL_DIR/software/km/data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa $INSTALL_DIR/software/km/data/jf/03H116_ITD.jf | km find_report -t $INSTALL_DIR/software/km/data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa
