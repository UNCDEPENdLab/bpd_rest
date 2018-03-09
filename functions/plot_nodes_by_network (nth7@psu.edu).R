plot_nodes_by_network <- function(toanalyze){
  require("taRifx")
  
  toanalyze <- data.frame(toanalyze, community.name=mapvalues(toanalyze$membership, from = c("1","2","3","4","5","6","7"), 
                                                              to = c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN"))) 
 
  
  toanalyze <- remove.factors(toanalyze)
  tomerge <- dplyr::select(atlas, name, anat_label) 
  colnames(tomerge) <- c("node", "anat_label")
  
  toanalyze <- left_join(toanalyze, tomerge, by = "node")
  
  toanalyze$BPD <- ifelse(toanalyze$BPD == 0, "Control", "BPD")
  
  pdf(file = paste0(basedir, "/figures/nodal_metrics_by_network_ixns_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", data.reduce,".pdf"), width = 11, height = 8)
  
  for(net in unique(toanalyze$community.name)) {
    this.metric <- toanalyze[which(toanalyze$community.name == net),]
    g <- ggplot(this.metric, aes(x = Age, y = central, color = factor(BPD))) + geom_point(alpha = .2) + stat_smooth(method="lm", se=FALSE) + 
      labs(x = "Age", y = "Centrality", title = paste0(as.character(this.metric$community.name), " Centrality Estimates by node")) + scale_color_brewer("Group", palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
      facet_wrap(~node )
    
    print(g) 
    
    g2 <- ggplot(this.metric, aes(x = Age, y = integration, color = factor(BPD))) + geom_point(alpha = .2) + stat_smooth(method="lm", se=FALSE) + 
      labs(x = "Age", y = "Integration", title = paste0(as.character(this.metric$community.name), " Integration Estimates by node")) + scale_color_brewer("Group", palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
      facet_wrap(~node )
    
    print(g2)
    
    g3 <- ggplot(this.metric, aes(x = Age, y = within.mod, color = factor(BPD))) + geom_point(alpha = .2) + stat_smooth(method="lm", se=FALSE) + 
      labs(x = "Age", y = "Within-Module Connectivity", title = paste0(as.character(this.metric$community.name), " Within-Module Conectivity Estimates by node")) + scale_color_brewer("Group", palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
      facet_wrap(~node )
    
    print(g3)
    
  }
  dev.off()
  
  pdf(file = paste0(basedir, "/figures/nodal_metrics_by_network_main_effects_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", data.reduce,".pdf"), width = 11, height = 8)
  
  for(net in unique(toanalyze$community.name)) {
    this.metric <- toanalyze[which(toanalyze$community.name == net),]
    g <- ggplot(this.metric, aes(x = factor(BPD), y = central, color = factor(BPD))) + geom_point(alpha = .02) + geom_boxplot() + scale_color_brewer("Group", palette="Set1") + theme_bw() +
      labs(x = "Group Membership", y = "Centrality", title = paste0(as.character(this.metric$community.name), " Centrality Estimates by node")) + facet_wrap(~node)
    print(g)  
    
      stat_smooth(method="lm", se=FALSE) + 
      labs(x = "Age", y = "Centrality", title = paste0(as.character(this.metric$community.name), " Centrality Estimates by node")) + scale_color_brewer("Group", palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
      facet_wrap(~node )
    
    print(g) 
    
    g2 <- ggplot(this.metric, aes(x = Age, y = integration, color = factor(BPD))) + geom_point(alpha = .2) + stat_smooth(method="lm", se=FALSE) + 
      labs(x = "Age", y = "Integration", title = paste0(as.character(this.metric$community.name), " Integration Estimates by node")) + scale_color_brewer("Group", palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
      facet_wrap(~node )
    
    print(g2)
    
    g3 <- ggplot(this.metric, aes(x = Age, y = within.mod, color = factor(BPD))) + geom_point(alpha = .2) + stat_smooth(method="lm", se=FALSE) + 
      labs(x = "Age", y = "Within-Module Connectivity", title = paste0(as.character(this.metric$community.name), " Within-Module Conectivity Estimates by node")) + scale_color_brewer("Group", palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5)) + theme_bw() +
      facet_wrap(~node )
    
    print(g3)
    
  }
  dev.off()
  
}


