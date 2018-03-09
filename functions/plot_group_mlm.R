plot_group_mlm <- function(mlm_res, browse = FALSE, weighted = FALSE, pdf = TRUE){
  suppressMessages(require(ggplot2)); suppressMessages(require(multcomp)); suppressMessages(require(wesanderson)); suppressMessages(require(ggpubr))
  
  #thresh can be FALSE, fc, or prop
  #weighted is either TRUE or FALSE
  
  if(weighted) {plot_file <- paste0(basedir, "/Figures/bpd_mlm_plot_lmerCM", file_tag_nothresh,".pdf")} else{
    plot_file <- paste0(basedir, "/Figures/bpd_mlm_plot_lmerCM", file_tag,".pdf")
  }
  
  #get network names
  if(use.yeo == 1){ 
    
    community <- yeo7_community(agg.g)
    
    membership <- data.frame(membership = factor(community$membership), 
                             community.name=mapvalues(community$membership, from = c("1","2","3","4","5","6","7"), to = c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")),
                             stringsAsFactors = FALSE)
  }
  
  if (browse == TRUE) {browser()}
  
  #mlm_res_both <- mlm_res
  
  pdf(file = plot_file, width=10, height=7)
 
  # for(j in 1:2){ 
  #   mlm_res <- mlm_res[[j]]
  for (i in names(mlm_res)){
    
    
    # cm <- lmerCellMeans(mlm_res[[i]])
    # cm <- inner_join(data.table(cm), data.table(membership), by = "membership")
    # 
    
    mlm_res <- mlm_res[[1]]
    cm <- lmerCellMeans(mlm_res[[i]], divide="Age")
    ggplot(cm, aes(x=membership, y=central, ymin=central-se, ymax=central+se, color=BPD)) + geom_pointrange() + facet_wrap(~Age)
    lstrends(mlm_res[[i]], var="Age", ~BPD*membership)
    
    # vv <- toanalyze
    # vv$membership <- factor(vv$membership)
    # vv$BPD <- factor(vv$BPD)
    # 
    # mm_slow <- lmer(central ~ 1 + membership + (1|node) + (1 + membership|id), vv)
    # 
    # 
    # mm <- lmer(central ~ 1 + membership*BPD + (1|node) + (1 |id), vv)
    # cm <- lmerCellMeans(mm)
    # ggplot(cm, aes(x=membership, y=central, ymin=central-se, ymax=central+se, color=BPD)) + geom_pointrange(position=position_dodge(width=0.5))
    # 
    # lsmeans(mm, ~BPD | membership)
    # 
    # idmat <- ranef(mm)$id
    # #ggplot(toanalyze, aes(x=factor(membership), y=central, color=factor(BPD))) + geom_boxplot()
    # ggplot(cm, aes(x=membership, y=central, color=BPD)) + geom_boxplot()
    # 
    
    pobj <- pairs(lsmeans(mlm_res[[i]], ~BPD|membership))
    mycon <- pobj@linfct
    comparisons <- glht(mlm_res[[i]], linfct=mycon)
    comparisons.table <- data.frame(table_glht(comparisons))
    comparisons.table[,"Pr...z.."] <- format(comparisons.table[,"Pr...z.."], scientific = F) 
    
    # as.matrix(as.numeric(comparisons.table))
    community.name <- unique(cm$community.name)
    
    display_table <- cbind(community.name, comparisons.table) 
    colnames(display_table) <- c("Network", "Estimate", "Std.E", "Z.value", "p")
    
    asterisks.box <- ifelse(comparisons.table$Pr...z.. <= .001, "***", 
                            ifelse(comparisons.table$Pr...z.. <= .01, "**", 
                                   ifelse(comparisons.table$Pr...z.. <= .05, "*",
                                          ifelse(comparisons.table$Pr...z.. <= .1, ".", "")) ))
    max.y <- array(NA, length(community.name))
    for(net in 1:length(community.name)){
      comm <- community.name[net]
      subset <- filter(cm,community.name == comm)
      max.y[net] <- as.numeric(max(subset[[i]]))
      
      
      
      
      annotation_df <- data.frame(community.name = community.name, asterisks.box = asterisks.box, max.y =max.y)
      
      pobj.age <- summary(pairs(lstrends(mlm_res[[i]], ~BPD|membership, var = "Age")))
      display_table_age <- data.frame(Network = community.name, Estimate = pobj.age$estimate, Std.E = pobj.age$SE, t.ratio = pobj.age$t.ratio, p = pobj.age$p.value)
      
      age.pvals <- pobj.age$p.value
      asterisks.age <- ifelse(age.pvals <= .001, "***", 
                              ifelse(age.pvals <= .01, "**", 
                                     ifelse(age.pvals <= .05, "*",
                                            ifelse(age.pvals <= .1, ".", "")) ))
      
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
    }
      
      
      # pdf(file = paste0(plot_file, "_testing.pdf"), width=10, height=7)
      a <- ggplot(cm, aes(x=Age)) + 
        geom_smooth(aes_string(y = i, color = "BPD"), method = lm)  + 
        labs(title = i) + scale_fill_manual(values=wes_palette(n=2, name="Darjeeling")) +  facet_wrap(~community.name, scales = "free") +
        geom_text(data = annotation_df_age, position = position_nudge(y = +0.1), aes(x = x, y =y, label = labs, group = NULL), size = 6) + theme_bw()
      
      print(a)
      
      #print(display_table_age)
      
      
      main.box <- ggplot(data = cm, aes(x = community.name)) + 
        geom_boxplot(aes_string(y = i, fill = "BPD")) + scale_fill_manual(values=wes_palette(n=2, name="Darjeeling")) + 
        labs(title = i) #+ geom_text(data = annotation_df, position = position_nudge(y = +0.1), aes(x = community.name, y =max.y, label = asterisks.box), size = 6) + theme_bw()
      print(main.box)
      #print(display_table)
      
      
  }
#}
    dev.off()
    message("wrote MLM group comparisons by network to: ", plot_file)
} 
  


