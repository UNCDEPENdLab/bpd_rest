plot_group_mlm <- function(mlm_res, browse = FALSE){
  require(ggplot2); require(multcomp); require(wesanderson); require(ggpubr)
  plot_file <- paste0(basedir, "/Figures/bpd_mlm_plot_lmerCM", parcellation, "_", preproc_pipeline, "_", conn_method, "_", data.reduce, ".pdf")
  #get network names
  if(use.yeo == 1){ 
    
    community <- yeo7_community(agg.g)
    
    membership <- data.frame(membership = factor(community$membership), 
                             community.name=mapvalues(community$membership, from = c("1","2","3","4","5","6","7"), to = c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")),
                             stringsAsFactors = FALSE)
  }
  
  if (browse == TRUE) {browser()}
  
  pdf(file = plot_file, width=10, height=7)
  #pdf(file = paste0(plot_file, "_testing.pdf"), width=10, height=7)
  for (i in names(mlm_res)){
    cm <- lmerCellMeans(mlm_res[[i]])
    cm <- inner_join(data.table(cm), data.table(membership), by = "membership")
   
    pobj <- pairs(lsmeans(mlm_res[[i]], ~BPD|membership))
    mycon <- pobj@linfct
    comparisons <- glht(mlm_res[[i]], linfct=mycon)
    comparisons.table <- data.frame(table_glht(comparisons))
    comparisons.table[,"Pr...z.."] <- format(comparisons.table[,"Pr...z.."], scientific = F)
    # as.matrix(as.numeric(comparisons.table))
    community.name <- unique(cm$community.name)
    asterisks.box <- ifelse(comparisons.table$Pr...z.. <= .001, "***", 
                                          ifelse(comparisons.table$Pr...z.. <= .01, "**", 
                                                 ifelse(comparisons.table$Pr...z.. <= .05, "*",
                                                        ifelse(comparisons.table$Pr...z.. <= .1, "†", "")) ))
   max.y <- array(NA, length(community.name))
    for(net in 1:length(community.name)){
      comm <- community.name[net]
      subset <- filter(cm,community.name == comm)
      max.y[net] <- as.numeric(max(subset[[i]]))
      
    }
   
    
   annotation_df <- data.frame(community.name = community.name, asterisks.box = asterisks.box, max.y =max.y)
     
   pobj.age <- summary(pairs(lstrends(mlm_res[[i]], ~BPD|membership, var = "Age")))
   age.pvals <- pobj.age$p.value
   asterisks.age <- ifelse(age.pvals <= .001, "***", 
                           ifelse(age.pvals <= .01, "**", 
                                  ifelse(age.pvals <= .05, "*",
                                         ifelse(age.pvals <= .1, "†", "")) ))
   
   annot.y <- array(NA, length(community.name))
   annot.x <- array(NA, length(community.name))
   
   for(net in 1:length(community.name)){
     comm <- community.name[net]
     subset <- filter(cm,community.name == comm)
     # max.y[net] <- as.numeric(max(subset[[i]]))
     x.range <- range(subset$Age)
     annot.x[net] <- (x.range[1] + x.range[2])/2
     annot.y[net] <- as.numeric(max(subset[[i]]))
     
   }
   annotation_df_age <- data.frame(x = annot.x, y = annot.y, community.name = community.name, labs = asterisks.age)

    colnames(cm)[4] <- i
    
    # pdf(file = paste0(plot_file, "_testing.pdf"), width=10, height=7)
    a <- ggplot(cm, aes(x=Age)) +
      geom_smooth(aes_string(y = i, color = "BPD"), method = lm)  + 
      labs(title = i) + scale_color_manual(values=wes_palette(n=2, name="Darjeeling")) +  facet_wrap(~community.name) +
      geom_text(data = annotation_df_age, position = position_nudge(y = +0.1), aes(x = x, y =y, label = labs, group = NULL), size = 6) 
    
    print(a)
    
    main.box <- ggplot(data = cm, aes(x = community.name)) + 
      geom_boxplot(aes_string(y = i, fill = "BPD")) + scale_fill_manual(values=wes_palette(n=2, name="Darjeeling")) + 
      labs(title = i) + geom_text(data = annotation_df, position = position_nudge(y = +0.1), aes(x = community.name, y =max.y, label = asterisks.box), size = 6)
    print(main.box)
    
    }
  dev.off()

  message("wrote MLM group comparisons by network to: ", plot_file)
}
 
  