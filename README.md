# CBB752 Final Project 2.6, R card, by Julian Q Zhou

## Objective

Given a FASTA file and a k-mer of interest, assess enrichment of the k-mer in the sequence(s) contained in the FASTA file.

## Source code

Available [here](https://github.com/jqz752/cbb752_2.6_R)

* `r_kmer_enrich.R`: main script

## Sample input
* A FASTA file	that looks like below (e.g. `sample_input_mapk3.fa`):

`>gi|568815582|ref|NC_000016.10|:30113184-30124230 Homo sapiens chromosome 16, GRCh38.p2 Primary Assembly`

`CCAGGCTGCGCCAGCATCTTCTTCCTCGTCCACACCCTGCCCTGCCACTTCGCTCTCCTTCTCTCTTGGT`

`CCCTGCCCCGTTTCTAGCATGCCCCCTTGGACCTACCCCTCTGTGCTGTCCACTTTGGCACCTGTTCTCA`

`CCCCTACCCGGCTCACCTCCTCGGTGGGCCCCCAGGCGGATGCGGAAGGTGGGAGCCCTGGGCGTGTGCA`

`GCAGATGAGGCCGGCGCAGGAAGAAGATGGAGAGCATGGCATAGCTGCCCAGGGCAGGGAGGGCATAGTA`

## Sample output
* A comma-separated .csv file that looks like below (e.g. `sample_output_cd47_brca2.csv`):

kmer.hits | kmer.p | count.q1 | count.mean | count.median | count.q3 | count.max | kmer.p.adj |significant
---|---|---|---|---|---|---|---|---
44|0.390625|26|56.2412109375|53|76|534|0.4921875|0
87|0.4921875|49|95.1396484375|87.5|128|1128|0.4921875|0

## Usage

`is.enriched.wrapper(input.fasta, input.alphabet, target.kmer, output.csv, significance = 0.05, multiple.testing = NULL)`

* `input.fasta`: FASTA file containing (possibly multiple) sequences; case-insensitive
* `input.alphabet`: a character vector of all possible letters; case-insensitive
* `target.kmer`: a string of k-mer of interest; case-insensitive
* `output.csv`: filename of output csv file
* `significance`: significance level; between 0 and 1
* `multiple.testing`: if `NULL`, no correction for multiple testing; otherwise, one of the possible arguments for `method` in `p.adjust()` (i.e. `bonferroni`, `holm`, `hochberg`, `hommel`, `BH`, `fdr`, or `BY`)

The output csv file contains the following columns for each input sequence:
* `kmer.hits`: number of times the k-mer of interest appears in the sequence
* `kmer.p`: empirical p-value for enrichment of k-mer of interest
* `kmer.p.adj`: if applicable, adjusted p-value for enrichment of the k-mer of interest
* `count.q1`/`mean`/`median`/`q3`/`max`: summary statistics of counts of all k-mers present
* `significnt`: `0`/`1`; `1` if `kmer.p`(`.adj`) < `significance`

For a given sequence, `r_kmer_enrich.R` computes the frequency of all the k-mers present in that sequence. Based on this distribution, the empirical (one-tailed) p-value for enrichment of the k-mer of interest is computed. If input contains multiple sequences, correction for multiple testing can be performed using a method of the user's choice (`multple.testing`). The empirical p-value (adjusted p-value in the case of correcting for multiple testing) is then compared with user-defined `significance` level, and a value of `1` indicates significance (p-value < `significance`).

#### Example 1 (Input FASTA file has only 1 sequence; no correction for multiple testing)

`is.enriched.wrapper(input.fasta = 'sample_input_mapk3.fa',
                    input.alphabet = c('A', 'T', 'G', 'C'),
                    target.kmer = 'GTC',
                    output.csv = 'sample_output_mapk3.csv',
                    significance = 0.01, multiple.testing = NULL)`

See `sample_output_mapk3.csv` for output.

#### Example 2 (Input FASTA file has >1 sequences; correction for multiple testing)
 
`is.enriched.wrapper(input.fasta = 'sample_input_cd47_brca1.fa',
                    input.alphabet = c('A', 'T', 'G', 'C'),
                    target.kmer = 'ATTGC',
                    output.csv = 'sample_output_cd47_brca1.csv',
                    significance = 0.01, multiple.testing = 'fdr')`

See `sample_output_cd47_brca1.csv` for output.
