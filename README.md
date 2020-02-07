![logo](https://github.com/gstarrett/oncovirus_tools/blob/master/oncovirus_tools.png)  
*In progress*  
Suite of tools for the analysis of viruses involved in cancer with a specific focus on understanding the attributes of integrated viruses.

## Workflow
1. Align reads to reference genome containing host and/or viral sequences using software such as bowtie2 or bwa that will generate a bam file
2. Input bam file and reference fasta into `get_int_loc.pl`
   * `get_int_loc.pl` will also output a fastq file of all viral sequence containing pairs and will assemble these using SPAdes
   * `normVirReads.py` can be used at this step to normalize the viral coverage or "copy number" relative to the remaining viral genome (in case of large deletions) and the number of human reads
3. Prepare SPAdes fastg file for plotting with R using `draw_fastg.pl`
4. Annotate SPAdes fastg using `annotate_fastg.pl`
5. Plot fastg file using R and ggraph with `plot_annotated_fastg.pl`
6. Based on these values and the output of `get_int_loc.pl` refine integration sites
7. Get host and virus sequences using the aforemented refined integration sites with `get_break_seqs.pl`

## Requirements
* SPAdes
* samtools
* bedtools
### Perl modules
* Bio::DB::FASTA
* Text::Levenshtein
* Bio::Perl
* Algorithm::BloomFilter
### R packages
* ggplot2
* ggraph
* tidygraph
* reshape2

## Script descriptions and usage
### bootstrap_locations.pl
Randomly get 1000 human chromosomal locations.
### bootstrap_microhomology.pl
Randomly get 1000  human and virus integration pairs and calculate the degree of microhomology
### annotate_fastg.pl
Annotate the fastg data in a tabular format based on given reference sequences using BLAST.
### draw_fastg.pl
Make a tabular text file that can be further annotated or plotted using R.
### get_break_seqs.pl
Get the viral and host sequences overlapping "break" or viral integration sites to assess microhomology and other features.
### get_int_loc.pl
Find the likely integration sites of a particular virus in a host genome while extracting all read pairs containing viral sequences identified by a bloom filter built on 32bp k-mers based on the viral reference sequence.
### integration_cnvs.pl
Coming soon.
### mutual_exclusivity.R
R script that calculated mutual exclusivity/co-occurence from an input matrix, such as for somatic alterations in genes (rows) by patient/sample (columns)
### normVirReads.py
Calculate the fraction of viral genome bases covered and normalize the coverage based on these values as well as number of human reads.
### plot_annotated_fastg.pl
Make the pdf file containing the ggraph network plot of fastg data
