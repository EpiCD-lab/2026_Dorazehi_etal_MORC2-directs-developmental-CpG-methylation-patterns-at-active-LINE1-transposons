# Dorazehi et al. (2026), MORC2 directs developmental CpG methylation patterns at active LINE-1 transposons

## Bulk RNA-seq data analysis

### FastQ file processing
Bulk RNA-seq FastQ files were processed using the Snakemake pipeline *Bulk-RNAseq_analysis.sn*. 
The pipeline was run using:
```bash
snakemake -s Bulk-RNAseq_analysis.sn -j 6 --latency-wait 120 --rerun-incomplete
```
### Read count
Reads were counted using FeatureCounts. 
Example of code:
```bash
featureCounts -s2 -F GTF -g gene_id -a <genome_annotation.gtf> -o <output.csv> <input.bam>
```
The subsequent differential expression analysis was conducted in R using the package DESeq2. 
An example of the DESeq2 analysis workflow can be found in *Bulk_RNAseq_DEA_DESeq2_example.Rmd*. 


## DNA methylation analysis

Basecalling and modified-base calling (for 5mC in CG context) was conducted using Dorado. 
Example of code:
```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.3.0 --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v --reference <hg38_ref.fasta> <input.pod5> > <output.bam>
```


Methylartist was used for the extraction of methylation calls from modbam files. The resulting files were visualised in R. 
Example of code:

```bash
methylartist segmeth -b <modbam.bam> --motif CG --ref <hg38_ref.fasta> -i <region.bed>
```


Modkit and BedGraphToBigWig were used to create BigWig files used for visualisation through heatmaps (DeepTools) and IGV. 
Example of code:

```bash
# Tabulation of CG methylation 
modkit pileup --cpg --ignore h --combine-strands --bedgraph --prefix <output> --ref <hg38_ref.fasta> <modbam.bam> <output.bedgraph>
# Sorting of the bedgraph by coordinates and retaining only relevant information
bedtools sort -i <output.bedgraph> | cut -f 1-4 > <sorted.bedgraph>
# creating the bigwig
bedGraphToBigWig <sorted.bedgraph> <chrom_sizes> <output.bw>
```


## CUT&Run data processing

### FastQ file processing

FastQ files were processed using the Snakemake pipeline *CUT_and_RUN_analysis.sn*. 
The pipeline was run using:
```bash
snakemake -s CUT_and_RUN_analysis.sn -j 6 --latency-wait 120 --rerun-incomplete
```

### Peak calling

Peak calling was conducted using HOMER. 
Example of code:

```bash
# Creation of tag directories for both control and target samples
makeTagDirectory <output_tag_dir> <input.bam>
```
```bash
# Peak calling over the control background
findPeaks <target_tagDir> -style histone -region -minDist 1000 -size 300 -o <output_peaks_list.txt> -i <control_tagDir>
```
Called peaks were manually inspected in IGV, and additional filtering was applied as necessary. 
For more details about the filtering, as well as specific parameters and samples used as control, please refer to the publication.

### Peak characterisation

Peaks were annotated using a priority order specific to each peak type. Each peak received a single annotation: if it matched multiple categories, it was assigned to the first matching category in the order below.
- MORC2 peaks:
```
1. All L1 elements
2. Peaks overlapping gene regions
4. Peaks not assigned to any of the categories above 
```

- H3K9me3-positive H3K4me3 peaks:
For these peaks, L1 annotation was restricted to full-length elements from the youngest subfamilies, which are more likely to be transcriptionally active and carry the H3K4me3 mark.
```
1. Selected young L1 elements (Full length L1HS, L1PA2, or L1PA3 - longer than 6 kb)
2. Transcript promoter windows
3. Other annotated repeats
4. Peaks not assigned to any of the categories above 
```

Peak characterisation and subsequent visualisation was conducted in R.

### Heatmaps 
Heatmaps of Cut&Run data, ATAC-seq data and DNA methylation data were created with DeepTools. 
Example of code:

```bash
computeMatrix reference-point -S <big.wig1> <big.wig2> ... --referencePoint center -R <target_genomic_region_coordinates.bed> --missingDataAsZero -b 10000 -a 10000 -o  <output_matrix>
# heatmaps were created without clustering
plotHeatmap -m <output_matrix> --sortRegions descend --plotFileFormat pdf -o <output.pdf>
# or with clustering
plotHeatmap -m <output_matrix> --sortRegions descend --plotFileFormat pdf --kmeans 2 --clusterUsingSamples <sample to base the clustering on> -o <output.pdf>
```

## Single nucleus RNA-seq data analysis
All single nucleus RNA-seq FastQ and subsequent files were processed using cellRanger as specified in *snRNAseq_Seurat_Processing.R*.

---

For more information on specific parameters and/or version of the bioinformatic tools, please refer to the publication.
All subsequent R analyses and visualisation were conduction in RStudio 2024.09.1+394 (R v4.4.2).
