

library(dplyr)


## cis positive direction
cis_Positive=df_filtered[ (df_filtered$`CIS.TRANS`=="CIS" &df_filtered$`Direction`=="+"), ] # TRANS ,CIS

df_defined<-cis_Positive
out_prefix<-deparse(substitute(df_defined))
pqtl_type="CIS" # CIS, TRANS
Direction="Positive" # "Positive", "Negative"
# Load your data and preprocess it

# Define the shapes and their order
shapes <- c("nominally significant" = 21, "FDR" = 22, "strictly significant" = 24)

# Add a new column "shape" with the shapes based on "shape_vale" values
df_defined <- df_defined %>% 
  mutate(shape = factor(shape_vale, levels = c("nominally significant", "FDR", "strictly significant")))

p <- ggplot(df_defined, aes(x = Phenotype, y = Gene)) +
  geom_point(size = 3, aes(fill = P.value, shape = shape), stroke = 0.2,color="black") +
    labs(title = glue("{pqtl_type} PQTLs Beta {Direction}"))+
  scale_shape_manual(values = shapes) +
  scale_fill_gradient(low = "#f1f1f6", high = "#0707f8") +  # Use a gradient color scale
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 8, face = "bold")) +
  theme(legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size = 8))

# Move the legend to the bottom in two rows
p <- p + theme(legend.position = "bottom")
p <- p + theme(legend.direction = "horizontal",
               legend.margin = margin(1, 1, 1, 1),
               panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
               panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "#9d9999"),
               panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "#9d9999"),
               legend.box = "vertical",  # Split the legend into two columns
               legend.box.margin = margin(t = 0, b =0))  # Margin to adjust spacing

# Adjust the legend title size
p <- p + theme(legend.title = element_text(size = 7),plot.title = element_text(hjust = 0.5))

# Save the plot to a PDF file
ggsave(glue("All_phenotype_{pqtl_type}_{Direction}Direction_shapePvalue_ColorPvalue.pdf"), plot = p, device = "pdf", width = 5, height = 12)
