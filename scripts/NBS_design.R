design_matrix_nbs <- subj_info$BPD
design_matrix_nbs <- cbind(design_matrix_nbs, ifelse(design_matrix_nbs == 0, 1, 0))
write.table(design_matrix_nbs, file = "/Users/nth7/Documents/MATLAB/NBS1.2/schaefer422_aroma_ridge.net_partial/design_matrix_group.txt", row.names = FALSE, col.names = FALSE)
