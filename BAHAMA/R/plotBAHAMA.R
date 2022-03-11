plotBAHAMA <- function(BAHAMA, PT_labels=0, HLT_labels=0, HLGT_labels=0, SOC_labels=0){
  if(PT_labels == 0){
    data_plot <- data.frame(RR_PT_x = BAHAMA@tidy_rr_pt$estimate,
                            RR_PT_y = abs(BAHAMA@tidy_rr_pt$estimate/BAHAMA@tidy_rr_pt$std.error),
                            RR_PT_size = BAHAMA@tidy_rr_pt$abs_diff,
                            PT_label = BAHAMA@tidy_rr_pt$label,
                            SOC = BAHAMA@tidy_rr_pt$SOC)
  }
  else{
    data_plot <- data.frame(RR_PT_x = BAHAMA@tidy_rr_pt$estimate,
                            RR_PT_y = abs(BAHAMA@tidy_rr_pt$estimate/BAHAMA@tidy_rr_pt$std.error),
                            RR_PT_size = BAHAMA@tidy_rr_pt$abs_diff,
                            PT_label = PT_labels,
                            SOC = BAHAMA@tidy_rr_pt$SOC)
  }

  p_pt <- ggplot(data_plot, aes(x = RR_PT_x, y = RR_PT_y, color=SOC, size=RR_PT_size))+
    geom_point() + theme_bw() +
    labs(x = "log(RR)", y = "log(RR)/SD", size="Difference in number of events", tag = "PT-level")+
    #geom_hline(yintercept=-log(0.05, base=10), linetype="dashed", alpha=0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
    geom_text_repel(data = subset(data_plot, y>1.96 | x<(-0.8) | x>1),
                    aes(x,y,label=PT_label), vjust=2, hjust=0.5, size=3, color = "black") +
    guides(alpha = FALSE)

  #HLT -level
  if(HLT_labels == 0){
    data_plot <- data.frame(RR_HLT_x = BAHAMA@tidy_rr_hlt$estimate,
                            RR_HLT_y = abs(BAHAMA@tidy_rr_hlt$estimate/BAHAMA@tidy_rr_hlt$std.error),
                            RR_HLT_size = BAHAMA@tidy_rr_hlt$abs_diff,
                            HLT_label = BAHAMA@tidy_rr_hlt$label,
                            SOC = BAHAMA@tidy_rr_hlt$SOC)
  }
  else{
    data_plot <- data.frame(RR_HLT_x = BAHAMA@tidy_rr_pt$estimate,
                            RR_HLT_y = abs(BAHAMA@tidy_rr_pt$estimate/BAHAMA@tidy_rr_pt$std.error),
                            RR_HLT_size = BAHAMA@tidy_rr_pt$abs_diff,
                            HLT_label = HLT_label,
                            SOC = BAHAMA@tidy_rr_hlt$SOC)
  }

  p_hlt <- ggplot(data_plot, aes(x = RR_HLT_x, y = RR_HLT_y, color=SOC, size=RR_HLT_size))+
    geom_point() + theme_bw() +
    labs(x = "log(RR)", y = "log(RR)/SD", size="Difference in number of events", tag = "HLT-level")+
    #geom_hline(yintercept=-log(0.05, base=10), linetype="dashed", alpha=0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
    geom_text_repel(data = subset(data_plot, y>1.96 | x<(-0.8) | x>1),
                    aes(x,y,label=HLT_label), vjust=2, hjust=0.5, size=3, color = "black") +
    guides(alpha = FALSE)

  #HLGT -level
  if(HLGT_labels == 0){
    data_plot <- data.frame(RR_HLGT_x = BAHAMA@tidy_rr_hlgt$estimate,
                            RR_HLGT_y = abs(BAHAMA@tidy_rr_hlgt$estimate/BAHAMA@tidy_rr_hlgt$std.error),
                            RR_HLGT_size = BAHAMA@tidy_rr_hlgt$abs_diff,
                            HLGT_labels = BAHAMA@tidy_rr_hlgt$label,
                            SOC = BAHAMA@tidy_rr_hlgt$SOC)
  }
  else{
    data_plot <- data.frame(RR_HLGT_x = BAHAMA@tidy_rr_pt$estimate,
                            RR_HLGT_y = abs(BAHAMA@tidy_rr_pt$estimate/BAHAMA@tidy_rr_pt$std.error),
                            RR_HLGT_size = BAHAMA@tidy_rr_pt$abs_diff,
                            HLGT_labels = HLGT_labels,
                            SOC = BAHAMA@tidy_rr_hlgt$SOC)
  }

  p_hlgt <- ggplot(data_plot, aes(x = RR_HLGT_x, y = RR_HLGT_y, color=SOC, size=RR_HLGT_size))+
    geom_point() + theme_bw() +
    labs(x = "log(RR)", y = "log(RR)/SD", size="Difference in number of events", tag = "HLT-level")+
    #geom_hline(yintercept=-log(0.05, base=10), linetype="dashed", alpha=0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
    geom_text_repel(data = subset(data_plot, y>1.96 | x<(-0.8) | x>1),
                    aes(x,y,label=HLGT_labels), vjust=2, hjust=0.5, size=3, color = "black") +
    guides(alpha = FALSE)

  #SOC -level
  if(SOC_labels == 0){
    data_plot <- data.frame(RR_SOC_x = BAHAMA@tidy_rr_soc$estimate,
                            RR_SOC_y = abs(BAHAMA@tidy_rr_soc$estimate/BAHAMA@tidy_rr_soc$std.error),
                            RR_SOC_size = BAHAMA@tidy_rr_soc$abs_diff,
                            SOC_labels = BAHAMA@tidy_rr_soc$label)
  }
  else{
    data_plot <- data.frame(RR_SOC_x = BAHAMA@tidy_rr_soc$estimate,
                            RR_SOC_y = abs(BAHAMA@tidy_rr_soc$estimate/BAHAMA@tidy_rr_soc$std.error),
                            RR_SOC_size = BAHAMA@tidy_rr_soc$abs_diff,
                            soc_labels = SOC_labels)
  }

  p_soc <- ggplot(data_plot, aes(x = RR_SOC_x, y = RR_SOC_y, color=SOC_labels, size=RR_SOC_size))+
    geom_point() + theme_bw() +
    labs(x = "log(RR)", y = "log(RR)/SD", size="Difference in number of events", tag = "HLT-level")+
    #geom_hline(yintercept=-log(0.05, base=10), linetype="dashed", alpha=0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
    geom_text_repel(data = subset(data_plot, y>1.96 | x<(-0.8) | x>1),
                    aes(x,y,label=SOC_labels), vjust=2, hjust=0.5, size=3, color = "black") +
    guides(alpha = FALSE)

  return(list(p_pt = p_pt,
              p_hlt = p_hlt,
              p_hlgt = p_hlgt,
              p_soc = p_soc))
}
