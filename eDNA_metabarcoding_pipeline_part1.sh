eDNA metabarcoding pipeline:

Rename the raw fastq.gz files by executing,
for file in *.gz; do newfile=$(echo "$file" | awk -F "_" '{print $1"_"$2"_"$3"_"$5"_"$6"_"$7}' | sed 's/ //g'); mv -- "$file" "$newfile"; done

Get QC metrics by executing script 'seqkit stats *.fastq.gz -a' inside the external drive folder containing raw data

Copy all raw renamed fastq.gz files from external drive to asax super computer by executing 'scp /Volumes/External/Metabarcoding_Raw_Data/*.gz uahjbb001@asax.asc.edu:/home/uahjbb001/Metabarcoding_Raw_Data'

Generate the QC reports for all raw fastq.gz files by executing the 'DADA2_QC_script.R' script in asax,
#!/bin/bash
module load R/4.1.0

R CMD BATCH DADA2_QC_script.R
printf "%s\n" "QC reports generated successfully!"

DADA2_QC_script.R
"library(dada2); packageVersion("dada2")
path <- "/home/uahjbb001/Metabarcoding_Raw_Data"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

pdffn = "QC1_R1.pdf"
pdf(file=pdffn, width=24, height=10)
plotQualityProfile(fnFs[1:49])
dev.off()

pdffn = "QC1_R2.pdf"
pdf(file=pdffn, width=24, height=10)
plotQualityProfile(fnRs[1:49])
dev.off()"
Similarly, execute a batch script to generate QC report for all samples.

Qiime2 (in asax):

Create the manifest file 'eDNA_pe-33-manifest.tsv' with three columns, 'Sample-id', 'forward-absolute-filepath', and 'reverse-absolute-filepath'. Sample-id contains sample names and other columns contains the file path, /home/uahjbb001/Metabarcoding_Raw_Data/*_R1/R2_001.fastq.gz. Copy the manifest file to asax super computer.

Execute the Qiime2 import data and cutadapat script as a batch script,
#!/bin/bash
source /home/uahjbb001/.bashrc
source /apps/profiles/modules_asax.sh.dyn
module load qiime2/2023.9 
conda activate qiime2-amplicon-2023.9  

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path eDNA_pe-33-manifest.tsv --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2
qiime cutadapt trim-paired --i-demultiplexed-sequences paired-end-demux.qza --p-front-f GGTCAACAAATCATAAAGATATTGG --p-front-r GGWACTAATCAATTTCCAAATCC --p-match-read-wildcards --p-match-adapter-wildcards --p-times 3 --p-overlap 10 --p-minimum-length 50 --p-discard-untrimmed --p-cores 9 --o-trimmed-sequences demux-paired-end-trimmed-fonly.qza --verbose
qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-trimmed-fonly.qza --p-adapter-f GGATTTGGAAATTGATTAGTWCC --p-adapter-r CCAATATCTTTATGATTTGTTGACC --p-match-read-wildcards --p-match-adapter-wildcards --p-times 3 --p-overlap 10 --p-minimum-length 50 --p-cores 9 --o-trimmed-sequences trimmed-seqs.qza --verbose
qiime demux summarize --i-data trimmed-seqs.qza --o-visualization trimmed-seqs.qzv

Execute the Qiime2 dada2 denoise script as a batch script,
qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trunc-len-f 120 --p-trunc-len-r 120 --p-trim-left-f 0 --p-trim-left-r 0 --p-n-threads 9 --p-n-reads-learn 1000000 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --verbose
qiime tools export --input-path table.qza --output-path Feature_table/
biom convert -i Feature_table/feature-table.biom -o Feature_table/feature-table.tsv --to-tsv

mkCOInr Classifier (in laptop): 

Commands for a quick installation of the mkCOInr conda environment and dependencies,
conda create --name mkcoinr python=3.9 -y
conda activate mkcoinr

python3 -m pip install cutadapt
conda install -c bioconda blast -y
conda install -c bioconda vsearch -y
pip install nsdpy

Download the mkCOInr working folder by executing git clone https://github.com/meglecz/mkCOInr

Then download the latest release of COInr database from Zenodo at https://zenodo.org/record/7898363, unzip, and rename the database folder. Execute,
wget https://zenodo.org/record/7898363/files/COInr_2023_05_03.tar.gz
tar -zxvf COInr_2023_05_03.tar.gz
rm COInr_2023_05_03.tar.gz
mv COInr_2023_05_03 COInr

Download a list of taxa occurring in the Hawaii islands,
1. Go to https://www.gbif.org/ and click Get data --> Occurrences
2. On the left panel under Occurrences, click Advanced
3. For Occurrence status, click Present
4. Scroll down to Location and click State province. Search 'state' and select to display a list of taxa on the right. 
5. Scroll down to Taxon and click Scientific name. Search "Animalia" and select to display a list of Animalia taxa from the chosen location on the right. 
6. Click on Map on the right panel to confirm the occurrence on the chosen location
7. Then, click on Download and download the taxa table using the Species list option
8. Once the table is downloaded as a csv file, rename the file as .tsv, and convert the table to a tab-separated-value table.
9. Then, copy and paste the entire column 'family' into a new worksheet and remove duplicates
10. Rename the header 'family' to 'taxon_name' and save the table as 'State_animalia_families.tsv' in the mkCOInr/data/example folder

Select sequences for a list of Animalia taxa with a minimum taxonomic rank (i.e., family). Execute,
perl scripts/select_taxa.pl -taxon_list data/example/State_animalia_families.tsv -tsv COInr/COInr.tsv -taxonomy COInr/taxonomy.tsv -min_taxlevel family -outdir tutorial/select_taxa -out State_animalia_COInr_selected.tsv

Select region using the bait_fas option. Execute,
perl scripts/select_region.pl -tsv tutorial/select_taxa/State_animalia_COInr_selected.tsv -outdir tutorial/select_region/bait_fas -e_pcr 0 -bait_fas data/one_seq_per_order_658-metazoa-trimmed.fasta -tcov 0.8 -identity 0.7

Dereplicate the sequences. Execute,
perl scripts/dereplicate.pl -tsv tutorial/select_region/bait_fas/trimmed.tsv -outdir tutorial/dereplicate -out dereplicated_trimmed_sequences.tsv

Format the database to qiime formats. Execute,
perl scripts/format_db.pl -tsv tutorial/dereplicate/dereplicated_trimmed_sequences.tsv -taxonomy COInr/taxonomy.tsv -outfmt qiime -outdir COInr/qiime_state -out COInr_qiime_state

Copy the qiime database files to asax super computer. Execute the Qiime2 mkCOInr classify script as a batch script (if classify using Qiime),
qiime tools import --type 'FeatureData[Sequence]' --input-path COInr_qiime_state_trainseq.fasta --output-path ref-fasta-state.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path COInr_qiime_state_taxon.txt --output-path ref-taxonomy-state.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-fasta-state.qza --i-reference-taxonomy ref-taxonomy-state.qza --o-classifier mkCOInr_state_classifier.qza --verbose
qiime feature-classifier classify-sklearn --i-classifier mkCOInr_state_classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy_mkCOInr_state.qza --verbose

Combine the taxonomy table with the ASV count table:
qiime feature-table transpose --i-table table.qza --o-transposed-feature-table transposed-table.qza
qiime metadata tabulate --m-input-file rep-seqs.qza --m-input-file taxonomy.qza --m-input-file transposed-table.qza --o-visualization merged-data.qzv
qiime tools export --input-path merged-data.qzv --output-path merged-data