ReadDirectory="/home/minion/MinIONData/2017_Jun_6_1D2_ScrappieFirstLook/1dsq_analysis/More_1d2_reads_split"
ResultsOutputDirectory="/home/minion/MinIONData/2017_Jun_6_1D2_ScrappieFirstLook/1dsq_analysis/nitpick_results"

source /home/ben/anaconda2/bin/activate minionvironment

cd src

python nit_picker_main.py \
 --reads=$ReadDirectory \
 --outputdir=$ResultsOutputDirectory \
 --barcode


source /home/ben/anaconda2/bin/deactivate
