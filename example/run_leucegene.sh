
# NEED:
# - sratoolkit: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software
# - jellyfish and km: https://github.com/iric-soft/km/blob/master/example/all_install.sh
# - 46GB available space (without --gzip in fastq-dump cmd)

## Move to km directory
# cd path_to_km

## Download fastq from GEO, see https://leucegene.ca/fr/pre-clinique/ressources/

# 03H041: FLT-ITD, NUP98-NSD1, MYC
sample="03H041"
mkdir -p data/leucegene/${sample}
cd  data/leucegene/${sample}
fastq-dump -I --split-files SRR949078
# less space but much longer
# fastq-dump -I --split-files --gzip SRR949078

## Execute Jellyfish and km
# -s = (0.50 * (8 * 1073741824 * [RAM]) / ([k_len] + [-c]))
jellyfish count -m 31 -o ./${sample}.jf -c 12 -s 799063683 -t 4 -C -L 2 '-Q+' <(cat ./*.fastq)
# execute jellyfish with fastq.gz files
# jellyfish count -m 31 -o ./${sample}.jf -c 12 -s 799063683 -t 8 -C -L 2 '-Q+' <(gunzip -c ./*.fastq.gz)

# Load the count table in RAM to improve the execution time of find_mutation
wc -l ./${sample}.jf

for file in ../../catalog/GRCh38/*.fa
do
  filename=$(basename "$file")
  filename="${filename%.*}"
  echo "Run on $filename"
  km find_mutation $file ./${sample}.jf | km find_report -t $file > ./${filename}.xls
done
