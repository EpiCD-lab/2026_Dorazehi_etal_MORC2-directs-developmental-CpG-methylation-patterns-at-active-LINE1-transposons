# Dorazehi et al. (2026), MORC2 directs developmental CpG methylation patterns at active LINE-1 transposons

---
### Methylation analysis

Methylartist was used for the extraction of methylation calls from modbam files. The resulting files were visualised in R. Example of code:

```bash
methylartist segmeth -b <modbam.bam> --motif CG --ref <ref.fasta> -i <region.bed>
```

Modkit and BedGraphToBigWig were used to create BigWig file used for visualisation through heatmaps and IGV. Example of code:

```bash
modkit pileup --cpg --ignore h --combine-strands --bedgraph --prefix <output> --ref <ref.fasta> <modbam.bam> <output.bedgraph>  
bedtools sort -i <output.bedgraph> | cut -f 1-4 > <sorted.bedgraph>
bedGraphToBigWig <sorted.bedgraph> <chrom_sizes> <output.bw>
```
