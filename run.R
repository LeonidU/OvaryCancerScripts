#mirna_list<-c("hsa-mir-124-1","hsa-mir-125b-2","hsa-mir-129-1", "hsa-mir-127", "hsa-mir-132", "hsa-mir-137", "hsa-mir-193a" ,"hsa-mir-148a","hsa-mir-191","hsa-mir-203a","hsa-mir-339","hsa-mir-375")
lncrna_list<-c("MEG3","MAGI2-AS3", "ZEB1-AS1", "NBLA00301")
# NBLA00301 - HAND2-AS1
protein_list <- protein_gene_names_symbol[protein_gene_names_symbol %in% rownames(z1)]
for (n1 in protein_list) {
for (n2 in lncrna_list) {
    a1<-cor.test(as.double(z3[n1,]), as.double(z3[n2,]), method="spearman")
    min_r <- as.double(a1$estimate)
    min_cc <- a1
	if ((!is.na(a1$estimate))&&(min_r > 0.6)&&(min_cc$p.value < 0.01)) {
		print(paste0(n1, " ", n2, " ", as.character(min_r), " ", as.character(min_cc$p.value)))
	}
}
}
