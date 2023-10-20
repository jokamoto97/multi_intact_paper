#Aggregate liver Multi-INTACT posterior results across lipid trait, then perform GSEA

library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
#library(biomaRt)
library(INTACT)
#library(org.Hs.eg.db)
#library(GO.db)


if(!dir.exists(paste0("results/ukbb_gtex_multi_intact/lipid_gsea"))){
        dir.create(paste0("results/ukbb_gtex_multi_intact/lipid_gsea"));
};


outfile <- "results/ukbb_gtex_multi_intact/lipid_gsea/lipid_gsea.txt"

outfile_input <- "results/ukbb_gtex_multi_intact/lipid_gsea/gsea_input.txt"

dir <- "results/ukbb_gtex_multi_intact/gppc_rst/"

#Metabolite annotation data

annot <- fread("data/metsim_merge_metabolite_annotation_kegg_hmdb.txt",h=T) %>% filter(SUPER_PATHWAY == "Lipid")

metab_ids <- read.table('data/metab_names.txt',h=F)$V1

#List of tested lipid metabolites

tested_lipids <- intersect(annot$CID,metab_ids)

rst_files <- paste0(dir,"Liver/",tested_lipids,"_Liver_gppc_rst.txt")

#GO annotations

annot_file <- "data/go_annotations.RData"

load(annot_file)

#Read in Multi-INTACT results

rst <- fread(cmd = paste0("awk FNR!=1 ",paste(rst_files,collapse=" ")),h=F)

#Format INTACT-GSE input

gse_input <- rst %>%
	group_by(V1) %>%
	summarise(aggreg_GPPC = 1-prod(1-V8)) %>% 
	arrange(desc(aggreg_GPPC)) %>%
	data.frame() 
	
gse_input$protein_id <- str_split_fixed(gse_input$V1,":",n = 3)[,2]
gse_input$gene <- gsub("\\_.*", "", gse_input$protein_id)

#Save GSEA input

gse_input %>% 
	dplyr::select(gene,aggreg_GPPC) %>%
	write.table(file = outfile_input,quote = F, sep = '\t',row.names = F, col.names = T)

#Gather GO BP data

#filter_genes <- gse_input %>% filter(aggreg_GPPC != 0) %>%
#                dplyr::select(gene) %>% unlist() %>% as.vector()

#ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#go.table <- getBM(attributes= c("ensembl_gene_id",
#				"hgnc_symbol",
#				"go_id","name_1006",
#				"namespace_1003"),
#		  mart= ensembl,filter = "hgnc_symbol",values = filter_genes)


#go.table <- go.table %>% filter(namespace_1003 == "biological_process")
#go.list <- go.table %>% dplyr::select(go_id) %>% unlist() %>% unique()

#Build list of annotations

#test_set<- list()

#for (j in 1:length(go.list)){

#	results <- AnnotationDbi::select(org.Hs.eg.db, keys=c(go.list[j]),
#					 columns = c('SYMBOL'), keytype = "GOALL")
#	gene_symbols <- unique(results$SYMBOL)

#	tmp <- list("go1" = gene_symbols)

#	test_set <- c(test_set,tmp)
#}

#names(test_set) <- go.list

#annot_file <- "data/go_annotations.RData"

#save(test_set,file = annot_file)

#Run INTACT-GSE

gene_data <- gse_input

rst <- NULL

.enrich_res <- utils::getFromNamespace(".enrich_res", "INTACT")

for (i in seq(1,length(test_set))){

    gene_data$d_vec <- 0

    gene_data$d_vec[gene_data$gene %in% test_set[[i]]] <- 1

    out <- .enrich_res(sig_lev = 0.05,
                      pprobs = gene_data$aggreg_GPPC,
                      d_vec = gene_data$d_vec,
                      SE_type = "NDS",
                      boot_rep = NULL)

    out$Gene_Set <- names(test_set)[i]

    out <- out[,c(8,seq(1,7))]

    rst <- rbind.data.frame(rst,out)
}


#save rst file

write.table(rst,file=outfile,row.names = F, col.names = T, sep = '\t',quote = F)

rst <- fread(outfile,h=T)

rst %>%
	filter(pval < 0.05) %>%
	arrange(desc(Estimate)) %>%
	dplyr::select(Gene_Set,Estimate,pval,CI_Leftlim,CI_Rightlim,CONVERGED) %>%
	mutate(Estimate = round(Estimate,digits = 3),
	       CI_Leftlim = round(CI_Leftlim,digits = 3),
	       CI_Rightlim = round(CI_Rightlim,digits = 3))
