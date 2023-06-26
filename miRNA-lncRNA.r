library(Biostrings)
library(rentrez)
library(dplyr)
library(purrr)
library(text.alignment)
library(scanMiR)

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

table_mRNA <- read.csv("/home/leo/Braga/mRNA-lncRNA-pos.xls", sep="\t", stringsAsFactor=FALSE)[,1]
mirnas <- data.frame(name = c("hsa-miR-135a", "hsa-miR-33b", "hsa-miR-377", "hsa-miR-500a", "hsa-miR-550a", "hsa-miR-8072", "hsa-miR-484", "hsa-miR-942", "hsa-miR-219a",  "hsa-miR-4708", "hsa-miR-3175", "hsa-miR-4513"),seed = c("tatggctttttattcctatgtga", "gtgcattgctgttgcattgc", "agaggttgcccttggtgaattc", "taatccttgctacctgggtgaga", "agtgcctgagggagtaagagccc", "ggcggcggggaggtaggcag", "tcaggctcagtcccctcccgat", "tcttctctgttttggccatgtg", "tgattgtccaaacgcaattct", "agagatgccgccttgctcctt", "cggggagagaacgcagtgacgt", "agactgacggctggaggccca"))
#data.frame(name = c("MIR124","MIR17","MIR191","MIR24-2","MIR148A","MIR137","MIR203","MIR375"), seed=c("cgtgttcacagcggaccttgat", "caaagtgcttacagtgcaggtag", "caacggaatcccaaaagcagctg", "tgcctactgagctgaaacacag", "aaagttctgagacactccgact", "acgggtattcttgggtggataat", "agtggttcttaacagttcaacagtt", "gcgacgagcccctcgcacaaacc"))
final_table <- data.frame(name_mRNA=character(0), name_miRNA=character(0), intersect=character(0), fraction=numeric(0))
for (a in table_mRNA) {
    mRNA <- fastaToList(getGeneSequence(a))[[1]]
#    b <- fastaToList(getGeneSequence(table[i,2]))
#    print(as.character(getComplement(b[[1]])))
    print(a)
    for (i in 1:nrow(mirnas)) {
        miRNA <- mirnas[i,]$name
        aa <- findSeedMatches(mRNA, mirnas[i,]$seed)
        level <- as.character(aa@elementMetadata$type)
        print(aa)
        final_table <- rbind(final_table, data.frame(name_mRNA=a, name_miRNA=miRNA, intersect=level, fraction=level/nchar(mRNA)))
        write.csv(final_table, file = "/home/leo/Braga/lncRNA-miRNA-binding-pos.txt", quote = FALSE, sep = "\t")
    }
}
#write.csv(final_table, file = "/home/leo/Braga/lncRNA-mRNA-binding-pos.txt", quote = FALSE, sep = "\t")