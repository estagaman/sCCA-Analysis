library(tidyverse)
library(dplyr)

#visualize the p-values for each method
#use the p-values column of each file 
unsup <- read.csv("/Users/elise/Downloads/CCA_results_collapsed.csv", row.names = 1)
sup_sub <- read.csv("/Users/elise/Downloads/CCA_results_supervised_sub.csv", row.names = 1)
sup_surv <- read.csv("/Users/elise/Downloads/CCA_results_supervised_surv.csv", row.names = 1)
unsup <- unsup[2:nrow(unsup), ]
sup_sub <- sup_sub[2:nrow(sup_sub), ]

library(ggplot2)

#create a dataframe of the p-values for each method and the associated chromosome
p_val_df <- data.frame(
  "sparse sCCA-subtype" = sup_sub$Subtype,
  "sparse CCA-subtype" = unsup$Subtype,
  "sparse sCCA-survival" = sup_surv$Survival,
  "sparse CCA-survival" = unsup$Survival,
  "Chromosome" = c(rep(unsup$Chromosome, 4))
)

#reshape the p-values dataframe for easier plotting
p_val_long <- p_val_df %>%
  pivot_longer(
    cols = -Chromosome,
    names_to = "Method",
    values_to = "p_value"
  )

#convert chromosome column to a factor
p_val_long$Chromosome <- factor(p_val_long$Chromosome,
                                levels = as.character(sort(as.numeric(unique(p_val_long$Chromosome)))))

#use ggplot to create a line plot, with each method denoted by color and line type
ggplot(p_val_long, aes(x = Chromosome, y = p_value, color = Method, group = Method, linetype = Method)) +
  geom_line(size = 0.7) +
  geom_point(size = 1) +
  scale_y_continuous(trans = "reverse") +
  scale_color_manual(values = c("#D95F02","#66C2A5", "#D95F02", "#66C2A5")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed")) +
  theme_minimal() +
  labs(title = "P-values by Chromosome",
       y = "P-value",
       x = "Chromosome")

#save the plot
ggsave("sCCA_results_plot.png", width = 8, height = 5, dpi = 300)
