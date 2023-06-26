library(Biostrings)
library(rentrez)
library(dplyr)
library(purrr)
library(text.alignment)

longestListIntersect <- function(a, b) {
    result <- 0
    max_len <- 0
    flag <- FALSE
    v <- (a == b)
    for (i in 1:length(a)) {
        if (v[i]) {
            if (!flag) {
                flag <- TRUE
                result <- 1
            } else {
                result <- result + 1
            }
        } else {
            flag <- FALSE
            if (max_len < result) {
                max_len <- result
            }
            result <- 0
        }
    }
    return(max_len)
}

# function to get complementary sequence
getComplement <- function(sequence) {
  dna_seq <- DNAString(sequence)
  return(reverseComplement(dna_seq))
}

# main function to find complementary regions between two gene lists
findComplementaryRegions <- function(gene1, gene2) {
  
  # calculate complements
  complements1 <- sapply(sequences1, get_complement)
  
  # find matching sequences
  matching_sequences <- intersect(complements1, complements2)
  
  return(matching_sequences)
}

getGeneSequence <- function(geneID){
    request <- paste0("Homo Sapiens[ORGN] AND ",geneID, "[GENE]")
    seq_id <- entrez_search(db="gene", term=request, retmax=1)
    linked_seq_ids <- entrez_link(dbfrom="gene", id=seq_id$ids, db="nuccore")
    linked_transripts <- linked_seq_ids$links$gene_nuccore_refseqrna
    result <- rentrez::entrez_fetch(db = "nuccore", id = linked_transripts, rettype = "fasta")
    return(result)
}

fastaToList <- function(fasta) {
    result <- list()
    fasta_splitted <- strsplit(fasta, "\n\n")[[1]]
    for (fasta_record in fasta_splitted) {
        result[strsplit(fasta_record, "\n")[[1]][1]] <- paste(strsplit(fasta_record, "\n")[[1]][-1], collapse="")
    }
    return(result) 
}

table <- read.csv("/home/leo/Braga/mRNA-lncRNA-pos.xls", sep="\t", stringsAsFactor=FALSE)[,1:2]
final_table <- data.frame(name_mRNA=character(0), name_lncRNA=character(0), intersect=numeric(0))
for (i in 1:nrow(table)) {
    a <- fastaToList(getGeneSequence(table[i,1]))
    b <- fastaToList(getGeneSequence(table[i,2]))
#    print(as.character(getComplement(b[[1]])))
    print(a)
    print(b)
    ss <- smith_waterman(a[[1]], as.character(getComplement(b[[1]])))
    level <- longestListIntersect(ss$a$alignment$tokens, ss$b$alignment$tokens)
    final_table <- rbind(final_table, data.frame(name_mRNA=table[i,1], name_lncRNA=table[i,2], intersect=level, fraction=level/nchar(a[[1]])))
    #print(fastaToList(getGeneSequence("SEMA3B-AS1")))
    #print(final_table)
    write.csv(final_table, file = "/home/leo/Braga/lncRNA-mRNA-binding-pos.txt", quote = FALSE, sep = "\t")
}
#write.csv(final_table, file = "/home/leo/Braga/lncRNA-mRNA-binding-pos.txt", quote = FALSE, sep = "\t")