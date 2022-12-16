#!/usr/bin/env Rscript
library("optparse")
library("tximport")


# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
option_list <- list( 
	make_option(c("-t","--tx2gene"), default=NULL, 
		help="Table with 2 columns: TXNAME and GENEID;",
		metavar="FILE.tsv"),
	make_option(c("-a","--kallisto_abundance"), default=NULL, 
		help="abundance.h5 file path (e.g. --kallisto_abundance path/to/abundance.h5);
		tximport require the name abundance.h5 :(",
		metavar="PATH/abundance.h5"),
	make_option(c("-s","--sample"), default=NULL, 
		help="Sample name (will be the .tsv column)",
		metavar="SAMPLE"),
	make_option(c("-o","--output"), default="STOUT", 
		help="filename for the output.
		output file: tsv file with 2 columns GENEID <sample>",
		metavar="PATH/OUTPUT")
	)

description <- "Convert abundance.h5 kallisto to a GENE count table.

	e.g. 
	%prog \\
		--kallisto_abundance a/abundance.h5 \\
		--tx2gene      gencode_tsx2gene.tsv \\
		--sample             sample1 \\
		--output             gene_tpm.tsv"

opt <- parse_args(OptionParser(option_list=option_list,description = description))

abundance_h5_file <- strsplit(opt$kallisto_abundance,"\\s+")[[1]]
tx2gene <- read.delim(opt$tx2gene, header = TRUE)
sample <- opt$sample

names(abundance_h5_file) <- sample

txi.kallisto.tsv <- tximport(abundance_h5_file, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

if(opt$output == "STOUT"){
	head(txi.kallisto.tsv$abundances)
} else {
	tpm <- txi.kallisto.tsv$abundance
	GENEID <- rownames(tpm)
	tpm_final <- cbind(GENEID = GENEID, tpm) 
	write.table(file = opt$output, x = tpm_final, row.names = FALSE, sep = "\t", quote = FALSE)
}
