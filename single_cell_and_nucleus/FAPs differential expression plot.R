# load in datasets
WT_vs_mdx <- openxlsx::read.xlsx("differential expression FAPs/mdxvsWT_cell.xlsx", sheet = 1)
WT_vs_TG <- openxlsx::read.xlsx("differential expression FAPs/TGvsWT_cell.xlsx", sheet = 1)
mdx_vs_TG <- openxlsx::read.xlsx("differential expression FAPs/TGvsmdx_cell.xlsx", sheet = 1)


require(ggplot2)
require(tidyverse)
require(ggrepel)
require(patchwork)

data <- list(WT_vs_mdx, WT_vs_TG, mdx_vs_TG)

for (df in seq(data)) {
  data[[df]]$threshold <- data[[df]]$pvals_adj < 0.05
  data[[df]]$color <- ifelse(data[[df]]$adj_fix < 0.05, "black", "grey80")
  data[[df]]$names <- str_to_title(data[[df]]$names)
}


scatterPlot <- function(df, 
                        scatter = TRUE,
                        label = TRUE) {
  
  ## start
  p <- ggplot(data = df, 
              aes(x = diff, y = metric)) +
    
    ## line for (0,0)
    geom_hline(yintercept = 0, 
               linetype = 2, 
               color = "black",
               linewidth = 0.35) +
    
    geom_vline(xintercept = 0, 
               linetype = 2, 
               color = "black",
               linewidth = 0.35) +
    
    ## fill colors by p value cutoff
    scale_fill_identity() +
    
    ## make ns more transparent
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), 
                       guide = "none") +
    
    labs(y = "Score", x = "Expressing cells (%)") +
    
    #  geom_density_2d(data = subset(WT_mdx_total, 
    #                                WT_mdx_total$pval_bulk < 0.05 & 
    #                                  WT_mdx_total$pval_decay < 0.05)) +
    
    theme(
      axis.title = element_text(colour = "black", size = 8),
      axis.text = element_text(colour = "black", size = 8),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      #axis.line = element_line(linewidth = 0.35, colour = "black"),
      panel.border = element_rect(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(-0.1, "cm"),
      aspect.ratio = 1 # make them all square
    )
  
  ## scatter plot
  if (scatter == TRUE) {
    p <- p + ggrastr::rasterise(
      geom_point(shape = 21, 
                 size = 0.5, 
                 stroke = 0,
                 aes(
                   # color by threshold
                   fill = color,
                   # transparency by threshold
                   alpha = threshold)),
      dpi = 1200) }
  
  if (label == TRUE) {
    
    p <- p + geom_text_repel(
      
      ## top up and down genes to label
      ## filter and bind to a single dataframe
      data = rbind(
        df %>% 
          slice_max(order_by = metric, n = 5),
        
        df %>% 
          slice_min(order_by = metric, n = 5)
        ),
      
      ## names column contains genes
      aes(label = names),
      
      ## make 6 pt font by convsersion
      size = 6 / ggplot2::.pt
      
    )
    
    }
  
  
  return(p)
  
} 

plot_list <- lapply(data, scatterPlot)

plot_list[[1]] | plot_list[[2]] | plot_list[[3]]


ggsave("Figures/FAPs differnetial expression.svg",
       height = 5, width = 5*3, units = "cm",
       dpi = 600)
