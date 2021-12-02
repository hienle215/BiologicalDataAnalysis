### LECTURE 11: SURVIVAL ANALYSIS and OTHER TOPIC

### TRANSCRIPT DATABASE-Example
# let's load the org.Hs.eg.db and TxDb.Hsapiens.UCSC.hg19.knowGene package first

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#let's retrieve the transcript information for gene LYN
lyn_key = select(org.Hs.eg.db,
                 keys = "LYN",
                 columns = "ENTREZID",
                 keytype = "ALIAS")
lyn_key

lyn_key = lyn_key$ENTREZID
lyn_key

transcripts = select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                     keys = lyn_key,
                     columns = c("GENEID", "TXCHROM", "TXSTART", "TXEND", "TXID"),
                     keytype = "GENEID"
                     )
transcripts

# we can see that the gene have two transcript variants stored in the databased
# To find the exon structure, letäs retrieve some more columns of data
transcripts_2 = select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                       keys = lyn_key,
                       columns = c("GENEID", "TXCHROM", "TXSTART", "TXEND", "TXID", "EXONSTART", "EXONEND"),
                       keytype = "GENEID")
transcripts_2

# in face in this case, two transcripts both have 13 exons, with a small difference visible above in exon 2.

### HUMAN GENOME SEQUENCE
# The librar BSgenome.Hsapiens.UCSC.hg1 in Bioconductor contains the human genome sequence build 19. Let's use this for some examples.
# Functions of the Biostrings package provide the foundations for handling sequence data in R

library(BSgenome.Hsapiens.UCSC.hg19)

hsgenome = BSgenome.Hsapiens.UCSC.hg19
hsgenome

# taking an individual chromosome, we find that the data is stored as a DNA string objech
hsgenome$chr20
alphabetFrequency(hsgenome$chr20)

# these strings can also be e.g indexed or converted into regular character sequences with as.character(). See also translate, reverse, complement, reverseComplement. Other strings types also exist RNAStrings, AAString
as.character(hsgenome$chr20)
translate(hsgenome$chr20)
reverse(hsgenome$chr20)
complement(hsgenome$chr20)
reverseComplement(hsgenome$chr20)
RNAString(hsgenome$chr20)
AAString(hsgenome$chr20)

### MASKED GENOME
# a msked genome is available from BSgenome.Hspaiens.UCSC.hg19.masked. Looking at the data here, we see:
library(BSgenome.Hsapiens.UCSC.hg19.masked)
hsmasked = BSgenome.Hsapiens.UCSC.hg19.masked 
hsmasked$chr20
# difficult regions are now masked (# in the sequence)
# you can use e.g active(masks(chr20)) = c(TRUE, FALSE, TRUE, FALSE) to set the active masks
active(masks(chr20)) = c(TRUE, FALSE, TRUE, FALSE)

### FINDINF PATTERNS
# a common task is to find all matches to a given pattern from the genome
# a function for doing this is called Biosrings:matchPattern()
library(Biostrings)
plismatches = matchPattern(DNAString("ATTTTCGGGG"),
                           hsgenome$chr4)
plismatches
### notice that the output is something called a View object: The object presents windows into an underlying data sequence, both coordinates and the underlying data (which therefore must be kept in memory).
# in fact, in many cases, the view information can actually be replaced simply by the coordinates for processing. This part of the output can be obtained with ranges(plus_matches)
# by default, matchPattern() only reports perfect matches so that all the sequences shown through the view will be identical in our example.
# in case, we want to search both strand, we nee to also look for the reverse complement of our pattern. 
minusmatches = matchPattern(reverseComplement(DNAString("ATTTTCGGGG")), 
                            hsgenome$chr4)
minusmatches

# Using this basic matchPattern() is slow if we want to take a lot of searches. For a quicker approach, if the search patterns are of equal width you can use matchPDict after building a PDict object from the search patterns

#### OPTIMAL ALIGNMENTS IN R
# from Biostrings package, you can also find pairwiseAlignment, which can do both Smith_waterman (local) and Needlemand-Wunsch (global) optional algnments for any pair of DNA or protein sequence
pairwiseAlignment(DNAString("CGGCTAAACT"),
                  DNAString("CTAACC")) #Global
pairwiseAlignment(DNAString("CGGCTAAACT"),
                  DNAString("CTAACC"),
                  type="local") # local

### GRanges
# gene coordinates as GRangs
# as our first example of a GRanges object, let's extract the gene coordinates from our transcript database. Similarly, we ca extract exone or transcript
gn_coords = GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
print(gn_coords)

# accessing GRanges data
# different "columns" of GRanges can be accessed with e.g mcol(), strand(), and seqnames(). For example, to see how the genes are divided into different chromosome, we can do
table(seqnames(gn_coords))

#if we want to see all the genes on chromosome 4, we can do indexing
gn_coords[seqnames(gn_coords) == "chr4"]


# retriveing sequence using GRanges
# to find sequences matching with selected GRanges location, we can used getSeq()
seq_res = BSgenome::getSeq(hsgenome, gn_coords[1:5])
seq_res

# to access individual sequence strings, you can use indexing
seq_res = BSgenome::getSeq(hsgenome, gn_coords[1:5])
seq_res[[1]]

### GRangesList
#transcript information can also be extracted in a GRangesList format. This is a subsettable GRangesList giving all the exones belonging to each of the 23459 genes
exons_by_genes=exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
print(exons_by_genes)

#You can also do e.g. exonsBy(You can also do e.g. exonsBy(ts,"tx"), see help for more.,"tx")
?exonsBy


#### Retrieving sequences with GRangesList
# we can also retrieve sequences with the GRangesList
eseqs = BSgenome::getSeq(hsgenome, exons_by_genes)
eseqs
# individual elements for this list can be accessed to get DNAStringSets as before
eseqs = BSgenome::getSeq(hsgenome, exons_by_genes[1:5])
eseqs[["10"]]


### SNP DATA
# related to the genome data, there are also SNP packages, the list of which cab be seen by using available.SNP fro BSgenome
BSgenome::available.SNPs()

# let's use SNPlocs.Hsapients.dbSNP.20120608. Notice that this SNP library is compatible with hg19 (SNPs will be mapped correctly)
library(SNPlocs.Hsapiens.dbSNP.20120608)
dbsnp = SNPlocs.Hsapiens.dbSNP.20120608
snps_result = snplocs(dbsnp, "ch4")
dim(snps_result)
head(snps_result,4)

#SNP data can also obtained as GRanges object
snps_results_granges = snplocs(dbsnp, "ch4", as.GRanges = TRUE)
snps_results_granges

### BIOMART onliune database
# let's start by loading the biomaRt package
library(biomaRt)

# BimaRT databases
# we start by listing available web service with listMarts()
listMarts()

# accessing dataset 
# each sever (called a mart) contains a number of datasets, which can be listed with listDatasets(). Let's see what we can access from ensembl.
ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useMart("ensembl")
listDatasets(ensembl)

# the dataset can be selected with useDataset(), let's do
hs = useDataset("hsapiens_gene_ensembl", 
                mart = ensembl)
hs

### filter
# We can check the available filters with listFilters().
filters = listFilters(hs)
dim(filters)
head(filters, 6)

# attributes
# similarly, the availabe attributes can be checked with listAttributes()
attributes = listAttributes(hs)
dim(attributes)
head(attributes, 6)

### USING getBM
# for example, we can retrieve all items in a given stretch of chromosome 4 and lsit their transcription start sites as follows
# notice that single items can have several start sites and there are start sites which do not have symbol. We only have on protein-coding gene in this part of the output
getBM(mart = hs,
      filters = c("chromosome_name", "start", "end"),
      values = list("4", 1, 500000),
      attributes = c("hgnc_symbol", "transcript_start")
)

# genomeGraph
# genomeGraphs is a package that can be used for visualizing e.g gene models or some measurement data such a microarrays. The gene models for genomeGraph can be downloads from Ensembl using biomaRT
libray(GenomeGraph)
cd33.gene = makeGene("hgnc_symbol", 
                     id = "CD33", 
                     biomart = hs)
cd33.transcripts = makeTranscript("hgnc_symbol",
                                  id = "CD33",
                                  biomart = hs)
cd33.gene@ens

chrom.axis = makeGenomeAxis(add53 = TRUE)
cd33toplot=list(makeTitle("CD33"), chrom.axis, "gene" = cd33.gene, "transcripts" = cd33.transcripts)
gdPlot(cd33toplot)


#### SNP DATA
# as another example of the kinds of data that can be retrieved. Let's see how we can retrieve SNP information
snpdb = useMart("ENSEMBL_MART_SNP", 
                dataset = "hsapiens_snp")
snpdb

res = getBM(snpdb,
            filters = c("chr_name", "start", "end"),
            values = list("4", 46500, 47000),
            attributes = c("chrom_start", "chrom_end", "phenotype_description", "allele"))
res