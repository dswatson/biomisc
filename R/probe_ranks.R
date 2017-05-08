# Rank Gene Comp

probe_ranks <- function(x, y) {
  
  df_x <- data_frame(Probe = x, Rank_x = row_number())
  df_y <- data_frame(Probe = y, Rank_y = row_number())
  df <- df_x %>% inner_join(df_y, by = 'Probe')
  cor.test(df$Rank_x, df$Rank_y, method = 'spearman')
  
}