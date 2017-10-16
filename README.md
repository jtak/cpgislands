# cpgislands
## A simple tool to find [CpG-islands](https://en.wikipedia.org/wiki/CpG_site#CpG_islands) in nucleotide sequences

#### Function

1. Read sequences from a file or directly from GenBank
2. Split the sequences to segments of a certain size
3. Decide if the segments are CpG-islands
4. Print all found islands

##### CpG-island definition
By default CpGislands uses the following definition for CpG-island:

* CG% at least 50%
* Observed-to-Expected CpG-site ratio at least 60%
* length at least 200 bp


#### Examples

Print help: `python3 cpgislands.py -h`

tiny.fasta contains a sequence for testing calculations with small segments:

`python3 cpgislands.py -f tiny.fasta -l 10`

Use `-v` to see the number of CpG-sites etc.

Use `python3 cpgislands.py -f SH1_genome.fasta` to test a longer sequence

Fetch a sequence from GenBank with ID NC\_012211:

`python3 cpgislands.py --genbank NC_012211 --email EMAIL-ADDRESS`

Entrez requires an email address to prevent abuse.


#### Requirements

* Biopython
* numpy



