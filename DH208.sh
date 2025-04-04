#!/bin/bash
#SBATCH --job-name=DH208 #job name
#SBATCH --output=/nobackup/vdvk82/DH208.out # output
#SBATCH --error=/nobackup/vdvk82/DH208.err # error
#SBATCH --time=48:00:00 # max runtime
#SBATCH --ntasks=1 #number of tasks
#SBATCH --cpus-per-task=8 # no of CPUs per tasks
#SBATCH --mem=100G

module load python/3.10.8
module load bioinformatics
module load filtlong/0.2.1
module load flye/2.9.2
module load prokka/1.14.6
module load samtools
module load minimap2
module load nanoplot

mkdir -p /nobackup/vdvk82/Nanoplot_DH208
mkdir -p /nobackup/vdvk82/Flye_DH208
mkdir -p /nobackup/vdvk82/Prokka_DH208
mkdir -p /nobackup/vdvk82/Filtlong_DH208

samtools view -b -f 4 /nobackup/vdvk82/dorado_demux_run1/48ff95e8b5ab9529617943432bf77443961cb71f_Sample5.bam > /nobackup/vdvk82/unmapped_DH208.bam
samtools fastq /nobackup/vdvk82/unmapped_DH208.bam > /nobackup/vdvk82/unmapped_DH208.fastq

NanoPlot --fastq /nobackup/vdvk82/unmapped_DH208.fastq --outdir /nobackup/vdvk82/Nanoplot_DH208

cd /nobackup/vdvk82/Filtlong_DH208
filtlong --min_length 1000 --min_mean_q 10 --target_bases 5000000000 /nobackup/vdvk82/unmapped_DH208.fastq > /nobackup/vdvk82/filtered_DH208.fastq

cd /nobackup/vdvk82
flye --meta --nano-raw /nobackup/vdvk82/filtered_DH208.fastq --out-dir /nobackup/vdvk82/Flye_DH208

prokka --outdir /nobackup/vdvk82/Prokka_DH208 --prefix 1 --force /nobackup/vdvk82/Flye_DH208/assembly.fasta

minimap2 -ax map-ont /nobackup/vdvk82/Flye_DH208/assembly.fasta /nobackup/vdvk82/filtered_DH208.fastq | samtools sort -o /nobackup/vdvk82/aligned_DH208.bam
samtools index /nobackup/vdvk82/aligned_DH208.bam

module load bedtools

bedtools bamtobed -i /nobackup/vdvk82/aligned_DH208.bam > /nobackup/vdvk82/mydata_DH208.bed

module load kraken2
module load krakentools
mkdir -p /nobackup/vdvk82/Kraken_DH208
kraken2 --db /nobackup/vdvk82/my_kraken_db --output /nobackup/vdvk82/Kraken_DH208/classification.txt --report /nobackup/vdvk82/Kraken_DH208/report.txt /nobackup/vdvk82/Flye_DH208/assembly.fasta

cd /nobackup/vdvk82
module load maxbin2
run_MaxBin.pl
mkdir -p /nobackup/vdvk82/maxbin_output_DH208

run_MaxBin.pl -contig /nobackup/vdvk82/Flye_DH208/assembly.fasta -abund /nobackup/vdvk82/depth_DH208.txt -out /nobackup/vdvk82/maxbin_output_DH208

mkdir -p /nobackup/vdvk82/prokka_bin_DH208

for bin in /nobackup/vdvk82/maxbin_output_DH208/*.fa /nobackup/vdvk82/maxbin_output_DH208/*.fna; do
    BIN_NAME=$(basename "$bin" | sed 's/\.[^.]*$//')
    prokka --outdir /nobackup/vdvk82/prokka_bin_DH208 --force --prefix "${BIN_NAME}_annotation" "$bin"
done

mkdir -p /nobackup/vdvk82/Kraken_bins_DH208
for bin in /nobackup/vdvk82/maxbin_output_DH208/*.fa /nobackup/vdvk82/maxbin_output_DH208/*.fna; do
    if [[ -f "$bin" ]]; then
        BIN_NAME=$(basename "$bin" | sed 's/\.[^.]*$//')

         kraken2 \
            --db /nobackup/vdvk82/my_kraken_db \
            --threads 8 \
            --report "/nobackup/vdvk82/Kraken_bins_DH208/${BIN_NAME}_report.txt" \
            --output "/nobackup/vdvk82/Kraken_bins_DH208/${BIN_NAME}_kraken_output.txt" \
            --use-names "$bin"

              echo "Finished processing $BIN_NAME!"
    fi
done

echo "All bins processed successfully!"
