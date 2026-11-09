# Dorazehi et al. (2026), MORC2 directs developmental CpG methylation patterns at active LINE-1 transposons


## DNA methylation analysis

Basecalling and modified-base calling (for 5mC in CG context) was conducted using Dorado. 
Example of code:
```bash
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.3.0 --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v --reference <hg38_ref.fasta> > <output.bam>
```


Methylartist was used for the extraction of methylation calls from modbam files. The resulting files were visualised in R. 
Example of code:

```bash
methylartist segmeth -b <modbam.bam> --motif CG --ref <hg38_ref.fasta> -i <region.bed>
```

Modkit and BedGraphToBigWig were used to create BigWig file used for visualisation through heatmaps and IGV. 
Example of code:

```bash
# Tabulation of CG methylation 
modkit pileup --cpg --ignore h --combine-strands --bedgraph --prefix <output> --ref <hg38_ref.fasta> <modbam.bam> <output.bedgraph>
# Sorting of the bedgraph by coordinates and retaining only relevant information
bedtools sort -i <output.bedgraph> | cut -f 1-4 > <sorted.bedgraph>
# creating the bigwig
bedGraphToBigWig <sorted.bedgraph> <chrom_sizes> <output.bw>
```


## CUTandRun data processing

### Peak calling

Peak calling was conducted using HOMER. 
Example of code:

```bash
# Creation of tag directories for both control and target samples
makeTagDirectory <output_tag_dir> <bam>
```
```bash
# Peak calling over the control background
findPeaks <target_tagDir> -style histone -region -minDist 1000 -size 300 -o <output_peaks_list.txt> -i <control_tagDir>
```
Called peaks were manually inspected in IGV, and additional filtering was applied as necessary 
For more details about the filtering, as well as specific parameters and samples used as control, please refer to the publication.


## Heatmaps 
Heatmaps of Cut&Run data, ATAC-seq data and DNA methylation data were created with DeepTools. 
Example of code:

```bash
computeMatrix reference-point -S <big.wig1> <big.wig2> ... --referencePoint center -R <target_genomic_region_coordinates.bed> --missingDataAsZero -b 10000 -a 10000 -o  <output_matrix>
# heatmaps were created without clustering
plotHeatmap -m <output_matrix> --sortRegions descend --plotFileFormat pdf -o <output.pdf>
# or with clustering
plotHeatmap -m <output_matrix> --sortRegions descend --plotFileFormat pdf --kmeans 2 --clusterUsingSamples <sample to base the clustering on> -o <output.pdf>
```

---

For more information on specific parameters and/or version of the bioinformatic tools, please refer to the publication.
