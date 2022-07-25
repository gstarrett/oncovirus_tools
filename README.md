[![DOI](https://zenodo.org/badge/123922184.svg)](https://zenodo.org/badge/latestdoi/123922184)

![logo](https://github.com/gstarrett/oncovirus_tools/blob/master/oncovirus_tools.png)  

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
### Perl modules (Perl v5.18)
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
Randomly get 1000 human and virus integration pairs and calculate the degree of microhomology
### annotate_fastg.pl
`perl annotate_fastg.pl input.fastg.txt input.blastn.txt > input.fastg.ann.txt`

Annotate the fastg data in a tabular format based on given reference sequences using BLAST (blast output is generated from get_int_loc.pl).
### draw_fastg.pl
`perl draw_fastg.pl input.fastg > input.fastg.txt`

Make a tabular text file that can be further annotated or plotted using R.
### get_break_seqs.pl
`perl get_break_seqs.pl input.int.txt`

Get the viral and host sequences overlapping "break" or viral integration sites to assess microhomology and other features.
### get_int_loc.pl
`perl get_int_loc.pl input.bam "virus_contig_name" ref.fa`

Find the likely integration sites of a particular virus in a host genome while extracting all read pairs containing viral sequences identified by a bloom filter built on 32bp k-mers based on the viral reference sequence.
### integration_cnvs.pl
Coming soon.
### mutual_exclusivity.R
R script that calculated mutual exclusivity/co-occurence from an input matrix, such as for somatic alterations in genes (rows) by patient/sample (columns)
### normVirReads.py
`python normVirReads.py input.virus.bedgraph input.stats.txt`

Example commands to generate the needed input
```
java -jar picard.jar BamIndexStats I=input.bam > input.stats.txt
samtools view -b input.bam 'virus' | bedtools genomecov -bg -ibam stdin > input.virus.bedgraph`
```
Calculate the fraction of viral genome bases covered and normalize the coverage based on these values as well as number of human reads.
### plot_annotated_fastg.pl
`perl plot_annotated_fastg.pl input.fastg.ann.txt`

Make the pdf file containing the ggraph network plot of fastg data

## Specific example for MCPyV
<p>This assumes that the user has already generated a deduplicated and indexed bam alignment against a reference containing the human genome hg38 and MCPyV named hg38_MCPyV.fa. This command needs 4 threads and will need at least 60GB (up to 120GB) of RAM for the BLAST and SPAdes steps.</p>

```
perl get_int_loc.pl MCC000.bam MCPyV hg38_MCPyV.fa
```

<p>This command will generate the following files and directories:</p>

* <b>MCC000.bam.MCPyV.TIMESTAMP.fastq :</b> unpaired fastq file of all virus-aligned, unaligned, and 25bp virus k-mer containing reads
* <b>MCC000.bam.MCPyV.TIMESTAMP.bedpe :</b> paired bed file containing virus-host discordant pair positions and strands to determine location and orientation of integration event
* <b>MCC000.bam.MCPyV.TIMESTAMP.int.txt :</b> Integration summary file with the following columns:
 * <b>hostChrom :</b> Host chromosome with the integration event
 * <b>hostStart :</b> Start position of the integration event on the host chromosome
 * <b>hostEnd :</b> Start position of the integration event on the host chromosome
 * <b>virChrom :</b> Virus "chromosome" name
 * <b>virStart :</b> Start position of the integration event on the virus genome
 * <b>virEnd :</b> End position of the integration event on the virus genome
 * <b>numPairs :</b> Number of read pairs supporting this integration event
 * <b>hostSkew :</b> Skewdness of read coverage on the host chromosome over the integration event (only relevant for hybrid capture approaches like ViroPanel). Negative skew indicates that the junction is at the hostStart position, positive skew indicates its at the hostEnd position
 * <b>virSkew :</b> Skewdness of read coverage on the viral genome over the integration event (only relevant for hybrid capture approaches like ViroPanel). Negative skew indicates that the junction is at the virStart position, positive skew indicates its at the virEnd position
* <b>MCC000_spades :</b> Default SPAdes output directory generated from MCC000.bam.MCPyV.TIMESTAMP.fastq
* <b>MCC000_spades/assembly_graph.fastg.blastn.txt :</b> BLASTn output of assembly graph
* <b>MCC000_spades/mh.10.txt :</b> Counts of microhomology events up in a 10bp window around integration junctions
* <b>MCC000_spades/mh.10.seq.txt :</b> Summary of microhomology sequences in a 10bp window around integration junctions

<p>Annotate the MCPyV integration assemblies

```
perl draw_fastg.pl MCC000_spades/assembly_graph.fastg > MCC000_spades/assembly_graph.fastg.txt
perl annotate_fastg.pl MCC000_spades/assembly_graph.fastg.txt MCC000_spades/assembly_graph.blastn.txt > MCC000_spades/assembly_graph.fastg.ann.txt
perl plot_annotated_fastg.pl MCC000_spades/assembly_graph.fastg.ann.txt
```
<p>The final command will generate a pdf file with the annotated assembly graph of the integrated virus structure. This can be used to determine number of copies and copy structure.</p>
<p>Normalize viral read abundance by host reads and percent viral genome coverage. This is useful in estimating the number of copies per tumor cell.

```
java -jar picard.jar BamIndexStats I=MCC000.bam > MCC000.stats.txt
samtools view -b MCC000.bam 'MCPyV' | bedtools genomecov -bg -ibam stdin > MCC000.MCPyV.bedgraph
python normVirReads.py MCC000.MCPyV.bedgraph MCC000.stats.txt
```
<p>Summary file contains the following columns in this order:</p>

* <b>Virus name</b>
* <b>File name</b>
* <b>normCp :</b> Normalized viral genome copies by percent coverage per 1000 human reads
* <b>huReads :</b> Total number of human reads
* <b>bases :</b> Number viral genome bases covered by reads
* <b>avgCov :</b> Number of viral reads over number of viral bases covered
* <b>cov :</b> Fraction of the viral genome covered by reads
* <b>length :</b> Length of the complete viral genome
