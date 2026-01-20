Qiime2 script used to process dataset

#!/bin/bash

#This command activates the Qiime2 environment. Depending on the computer and installation, the environment name (in this case, “qiime2-2019.4”) may differ. 
source activate qiime2-2019.4;
source tab-qiime;
#BEFORE DOING ANYTHING ELSE, cd into the folder you will be working in.

read -p "Input directory path where FASTQ files are contained (no subdirectories): " directory;
read -p "Input file path for map file (.txt): " mapfile;
read -p "Specify sampling depth for alpha and beta diversity statistics: " depth;

echo "FASTQ directory is $directory";
echo "Map file path is $mapfile";
echo "Sampling depth is $depth";
echo "Database is Greengenes and Silva";
#This command is the start of the official workflow. We are importing our raw FASTQ read1 and read2 files downloaded from basespace here. The –-type and –-input-format denoted is specific to Illumina files and to paired end sequencing (this would be different for single end read data). The –-input-path you feed should be a folder with all of your reads within and the reads should not be in sub-folders. You will get an error. I recommend that you search for “.fastq.gz” in the Basespace download directory of interest and copy the files from the search into a new folder, which in this case I named “FastQ”.
#the output is a .qza file which is specific to Qiime2 and cannot be directly used in other workflows like a .fna file might be.

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $directory  --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza;
echo "Sequences imported";
#The next command uses a dada2 plugin to filter out poor quality reads, determine the error rate from the specific Illumina run, and remove chimeras. You should note that if you are processing data from multiple sequencing runs, each run should be run through this step separately and then combined together after. This is because each run will have its own specific error rate and when combined, data can be lost or kept when it should not have been. The input is your file from above. The four parameters directly following the input are to inform dada2 on where to trim or truncate your sequences so they are all the same size. The “trim” parameter is where you want to cut on the left side and “trunc” is where you want to cut on the right. For this run, we would like to keep the full-length read so we pick “0” and “250” as cutting locations. 
#IMPORTANT: This command will take a long time to run depending on the size of your data set. I recommend that you run this overnight. Additionally, you will want to assign the number of processors that will be working on this job. The “-p--threads” command will let you choose how many processors. Typing “0” will, counter intuitively, assign all available processors to the job. I would recommend you do this to maximize the speed. The “–-verbose” flag will narrate what is happening with the process and should let you monitor things so that you will be able to judge whether the process is hung up on something or not.
qiime demux summarize --i-data demux-paired-end.qza --o-visualization demux-paired-end.qzv;

qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 250 --p-trunc-len-r 250 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 0 --verbose;
echo "dada2 filter and denoise complete";

qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file $mapfile;
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv;
echo "Table and Rep Seqs Generated";
#Align sequences and produce tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --p-n-threads 0 --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#This command will generate a large variety of alpha and beta diversity statistics, as well as a rarefied table of your samples. Make sure that you enter the correct sampling depth based on the previous step. All of your outputs will be put into the created directory. If you want to generate only specific metrics or a specific kind of metric, look on the Qiime2 tutorial pages for the desired command. 
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth $depth --m-metadata-file $mapfile --output-dir core-metrics-results;

#This next section will detail how you assign and visualize taxonomy of your data. This is done using your rep-seqs file generated from the dada2 step and a classifier. This classifier is basically a big library of sequences with a taxonomic ID assigned. This is a file you can download from the Qiime2 website and should be put into your working directory. In this case I am using a Greengenes reference with 99% OTUs and specifically the v4 region.
echo "Classification step";
qiime feature-classifier classify-sklearn --i-classifier $gg --i-reads rep-seqs.qza --p-n-jobs -1 --o-classification taxonomy-gg.qza;
qiime feature-classifier classify-sklearn --i-classifier $sv --i-reads rep-seqs.qza --p-n-jobs -1 --o-classification taxonomy-silva.qza;
qiime metadata tabulate --m-input-file taxonomy-gg.qza --o-visualization taxonomy-gg.qzv;
qiime metadata tabulate --m-input-file taxonomy-silva.qza --o-visualization taxonomy-silva.qzv;
#Here we take our generated taxonomy data and combine that with our original table and our metadata file to generate an interactive taxa bar plot. 
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy-gg.qza --m-metadata-file $mapfile --o-visualization taxa-bar-plots-gg.qzv;
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy-silva.qza --m-metadata-file $mapfile --o-visualization taxa-bar-plots-silva.qzv;
#The interactive bar plots can be downloaded as SVG files or a CSV file can be generated with the raw numbers. 

#Import Qiime2 artifacts into phyloseq for a more in-depth analysis


echo "Script complete";
#You can then go and view the plots generated from this. Make sure you specify .qzv files rather than .qza files. The desired ones will also be labeled as “emperor” plots.
