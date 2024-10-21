b01_l <- read.table("barcode01_l_results.txt")
b01_r <- read.table("barcode01_r_results.txt")


col_name <- c('query_id', 'subject_id', 'identity', 'length', 'mismatch', 'gapopen', 'pos_1_query', 'pos_2_query', 'pos_1_subject', 'pos_2_subkect', 'evalue', 'score')
colnames(b01_l) <- col_name
colnames(b01_r) <- col_name