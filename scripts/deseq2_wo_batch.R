#!/usr/bin/env Rscript
## COMMAND LINE PARSE ######
library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
option_list <- list( 
	make_option(c("-i","--input_dir"), default=NULL, 
		help="Path to gene count files (.tsv);
		input files must be tsv files with 2 columns: GENEID, <SAMPLE>",
		metavar="PATH"),
	make_option(c("-e","--input_extension"), default=NULL, 
		help="the difference between SAMPLE in <sample_table> and the input files

		e.g.
		files: A.tsv, B.tsv
		sampleTable$SAMPLE; A B
		in this case: --input_extension .tsv;",
		metavar=".EXTENSION"),
	make_option(c("-a","--group_a"), default=NULL, 
	            help="Group A in col_categ (e.g. KD); UP in A: log2FC > 0",
	            metavar="STR"),
	make_option(c("-b","--group_b"), default=NULL, 
	            help="Group B in col_categ (e.g. WT); UP in B: log2FC < 0",
	            metavar="STR"),
	make_option(c("-s","--sample_table"), default=NULL, 
		help="table with at least two columns: SAMPLE <col_categ>

		e.g.
		SAMPLE KO_STATUS
		A      KO
		B      KO
		C      WT
		D      KO
		E      WT",
		metavar="PATH/FILE.TSV"),
	make_option(c("-c","--col_categ"), default=NULL, 
		help="In the samples table the column header with the category to compair

		e.g. 
		KO_STATUS
		KO
		KO
		WT
		KO
		WT",
		metavar="COLNAME"),
	make_option(c("-o","--output"), default="STOUT", 
		help="filename for the output.
		output file: tsv file with 2 columns GENEID <sample>",
		metavar="PATH/OUTPUT")
	)

description <- ""

opt <- parse_args(OptionParser(option_list=option_list,description = description))

############################

suppressPackageStartupMessages({
  library('DESeq2')
  library('tidyverse')
  library('limma')
})

my_theme <- theme(axis.title = element_text(size=12,color='black'),
                  axis.text= element_text(size=10,color='grey60'),
                  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                  axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.ticks = element_line(color='grey80'),
                  strip.background =element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_line(color = "grey80",size=rel(1)),
                  panel.grid.minor = element_line(color = "#f3f3f3",size=rel(1)),
                  plot.title = element_text(size = 20, face = 'plain', hjust = 0,color='grey10'),
                  plot.subtitle = element_text(size = 12,color='grey40', hjust = 0),
                  legend.position = 'bottom',
                  plot.caption=element_text(hjust=0))

## FUNCTIONS ############

run_deseq2 <- function(colData,groups,group_a,group_b,countData){
	#e.g. groups: TNM; group1: I; group2: II;
  colData <- colData %>%
    arrange(sampleName) %>%
    select(sampleName,!!groups) %>%
    mutate(!!groups := if_else(.data[[groups]] == !!group_a, 
              factor("A",levels = c("A","B")),
              factor("B",levels = c("A","B")))) %>%
    column_to_rownames("sampleName")
	

	dds_data <- DESeqDataSetFromMatrix(colData   = colData,
	                                   countData = countData,
	                                   design    = as.formula(paste0("~ ",groups)))
	dds      <- DESeq(dds_data)
	keep     <- rowSums(counts(dds)) >= 10
	dds      <- dds[keep,]

	res      <- results(dds,contrast = c(groups,"A","B"))
	
	vsd <- vst(dds, blind=FALSE)
	
	group_translate <- data.frame(group = c("A","B"), contrast=c(group_a,group_b))
	
	
	write.table(assay(vsd) %>% as.data.frame %>% rownames_to_column("GENEID"),file = paste0(opt$output,"/",group_a,"_vs_",group_b,"_transformed_expression.tsv"),quote = FALSE,sep = "\t",row.names = FALSE)
	
	PCAdata <- plotPCA(vsd, intgroup=c(groups), returnData = TRUE)
	write.table(PCAdata %>% left_join(group_translate, by = "group") ,file = paste0(opt$output,"/",group_a,"_vs_",group_b,"_PCA.deseq2.tsv"),quote = FALSE,sep = "\t",row.names = FALSE)
	PCA_percent <- data.frame(PC = c("PC1","PC2"),percentVar = round(100 * attr(PCAdata, "percentVar")))
	write.table(PCA_percent,file = paste0(opt$output,"/",group_a,"_vs_",group_b,"_PCA_percent.deseq2.tsv"),quote = FALSE,sep = "\t",row.names = FALSE)
	res
}

do_write_file <- function(res,groups,group1,group2){
  x <- as.data.frame(res)
  x <- x %>% rownames_to_column(var="GENES") %>%
    mutate(UP_IN = if_else(log2FoldChange > 0, !!group1, !!group2))
  
     
  write.table(x,file = paste0(opt$output,"/",group1,"_vs_",group2,".deseq2.tsv"),quote = FALSE,sep = "\t",row.names = FALSE)
}

get_countData <- function(my_dir, samples_in_sampleTable) {
	
	all_files <- list.files(my_dir) %>% .[. %in% samples_in_sampleTable]
	countData = data.frame()
	for (my_file in all_files) {
		aux <- read.delim(file.path(my_dir,my_file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
		aux[,2] <- round(aux[,2])

		if(length(countData) == 0){
			countData <- aux
		} else {
			countData <- merge(countData,aux, by = "GENEID")
		}
	}

	rownames(countData) <- countData %>% select(GENEID) %>% pull()
	countData <- countData %>% select(-GENEID)

	as.matrix(countData)
}

#########################

# --

sampleTable <- read.delim(opt$sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sampleTable <- sampleTable %>%
   mutate(sampleName = paste0(SAMPLE)) %>%
   mutate(fileName = paste0(SAMPLE,opt$input_extension))

print(sampleTable)
countData  <- get_countData(opt$input_dir,sampleTable$fileName)

dir.create(opt$output, showWarnings = FALSE)
res <- run_deseq2(sampleTable,opt$col_categ,opt$group_a,opt$group_b,countData)
do_write_file(res,opt$col_categ,opt$group_a,opt$group_b)
