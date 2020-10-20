#### Gronau filters

We used the following filters to select regions representing neutral genetic diversity. Filters come from Gronau *et al.,* in [this paper](https://www.nature.com/articles/ng.937). 
Raw filters are available from [raw](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/Gronau_filters/raw).

Details of these regions:

>**filter_hotspot1000g.starch**		- 1000 Genomes recombination hotspots (avoid putting loci in these regions)

>**filter_Map20.starch**		- poor mapping quality

>**filter_rmsk20.starch**		- recent duplications (RepeatMasker score < 20)

>**filter_segDups.starch**		- recent segmental duplications

>**filter_selection_10000_100.starch**	- gene exons + 1000bp flanking + conderved elements + 100bp flanking

>**filter_simpleRepeat.starch**		- simple repeats 

>**filter_SysErrHCB.starch**		- positions with systematic seuqencing errors 

>**filter_SysErr.starch**		- positions with systematic sequencing errors


These files are in starch format which is a binary bed format. To work with these files, we used [bedops](https://bedops.readthedocs.io/en/latest/). 
First, we merged the files together:

```bash
bedops -m filter_hotspot1000g.starch filter_Map20.starch filter_rmsk20.starch filter_segDups.starch filter_selection_10000_100.starch filter_simpleRepeat.starch filter_SysErrHCB.starch filter_SysErr.starch | starch - > Gronau_filters_merge.starch
```

Gronau_filters_merge.starch contains positions to remove. Get list of positions to keep.
Make bed file containing full autosome coordinates: **autosomes.bed** (We manually generated a .idxstats file for this step):

```bash
awk '{print $1"\t"0"\t"$2}' MOS104_idxstats | head -n 22 > autosomes_sort.bed
```

Convert this bed file to starch format:

```bash
sort-bed autosomes.bed > autosomes_sort.bed
starch autosomes_sort.bed > autosomes.starch
```

Remove Gronau positions so we only have positions which we want to call left:

```
bedops -d autosomes.starch global_starch_filters/Gronau_filters_merge.starch | starch - > autosomes_with_Gronau_filters.starch
```

Split up Gronau filters by chromosome:

```bash
unstarch autosomes_with_Gronau_filters.starch | sed s/chr// - > autosomes_with_Gronau_filters.bed
for i in {1..22};do
awk '{if ($1=='"${i}"') print $0}' ../autosomes_with_Gronau_filters.bed > autosomes_with_Gronau_filters_chr${i}.bed &
done
```
These split files are available in [bed_files](https://github.com/EvolEcolGroup/data_paper_genetic_pipeline/tree/main/Gronau_filters/bed_files). Depending on the reference sequence, you would need to add "chr" in front of each chromosome number.
