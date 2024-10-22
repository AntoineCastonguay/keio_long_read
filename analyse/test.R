test <- read.table("../keio_long_read/res/out4/4_alignment/all_output/barcode02_aln_output.txt", header = T)

# Créer un data frame vide
df <- data.frame(query_id = character(), 
                 match = numeric(), 
                 miss_match = numeric(), 
                 gene = character(), 
                 pos_l_read = numeric(), 
                 pos_l_gene = numeric(),                  
                 diff_l = numeric(),
                 pos_r_read = numeric(), 
                 pos_r_gene = numeric(),
                 diff_r =numeric())

for (i in 1:nrow(test)) {
  if (test$pos_2_l_subject[i] != 'nd' & test$pos_1_r_subject[i] != 'nd') {
      if (as.numeric(test$pos_2_l_subject[i])>as.numeric(test$pos_1_r_subject[i])) {
        ens1 <- seq(as.numeric(test$pos_1_r_subject[i]),as.numeric(test$pos_2_l_subject[i]))
        var <- 'l'
      }else{
        ens1 <- seq(as.numeric(test$pos_2_l_subject[i]),as.numeric(test$pos_1_r_subject[i]))
        var <- 'r'
      }
      print(test$query_id[i])
    for (j in 1:(nrow(ecoli_positif))) {
      ens2 <- seq(ecoli_positif$first_pos[j], ecoli_positif$second_pos[j])
      if (length(ens1) > length(ens2)) {
        next
      }
      tab <- table(ens2 %in% ens1)
      if (length(tab) == 2) {
        print(tab)
        print(ecoli_positif$gene[j])
        print(var)
        
        if (var == 'l') {
          pos_r_l <- test$pos_2_l_subject[i]
          pos_g_l <- ecoli_positif$second_pos[j]
          pos_r_r <- test$pos_1_r_subject[i]
          pos_g_r <- ecoli_positif$first_pos[j]
          res_l <- as.numeric(test$pos_2_l_subject[i]) - ecoli_positif$second_pos[j]
          res_r <- as.numeric(test$pos_1_r_subject[i]) - ecoli_positif$first_pos[j]
        }else{
          pos_r_l <- test$pos_2_l_subject[i]
          pos_g_l <- ecoli_positif$first_pos[j]
          pos_r_r <- test$pos_1_r_subject[i]
          pos_g_r <- ecoli_positif$second_pos[j]
          res_l <- as.numeric(test$pos_2_l_subject[i]) - ecoli_positif$first_pos[j]
          res_r <- as.numeric(test$pos_1_r_subject[i]) - ecoli_positif$second_pos[j]
        }
        
        # Ajouter des lignes au data frame
        df <- rbind(df, data.frame(query_id = test$query_id[i], 
                                   match = tab[2], 
                                   miss_match = tab[1], 
                                   gene = ecoli_positif$gene[j], 
                                   pos_l_read = pos_r_l,
                                   pos_l_gene = pos_g_l,
                                   diff_l = res_l,
                                   pos_r_read = pos_r_r,
                                   pos_r_gene = pos_g_r,
                                   diff_r = res_r))
      }
    }
  }
}

write.csv(df,file = "../keio_long_read/analyse/b2.csv")
