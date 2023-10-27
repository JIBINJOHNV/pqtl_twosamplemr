


df_filtered$shape_value <- ifelse(is.na(df_filtered$P.value), "no_value",
        ifelse(df_filtered$P.value > 3 & df_filtered$P.value < 6, "FDR",
        ifelse(df_filtered$P.value >= 6, "Strictly significant", "not significant")))

library(ggplot2)

# Your existing code
p <- ggplot(df_filtered, aes(x = Phenotype, y = Grne, color = P.value)) +
  geom_point(aes(size = P.value, shape = shape_value)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(0, 3)) +
  facet_wrap(~CIS.TRANS) +
  theme_minimal()

# Modify the font size of y-axis tick labels
p <- p + theme(axis.text.y = element_text(size = 6))  # Adjust the size to 6

# Add a slightly darker background for facet headers and adjust the panel border
p <- p +
  theme(strip.background = element_rect(fill = "lightgray", color = "gray"),  # Adjust colors as needed
        panel.border = element_rect(color = "lightgray", fill = NA, size = 1))  # Adjust border color and size

# Save the plot to a PDF file
ggsave("test_color_and_shape_gradient.pdf", plot = p, device = "pdf")





library(ggplot2)

# Your existing code
p <- ggplot(df_filtered, aes(x = Phenotype, y = Grne, color = P.value)) +
  geom_point(aes(shape = shape_value), size = 2.5) +  # Set the size to 4 (or your desired fixed size)
  scale_color_gradient(low = "blue", high = "red") +
  facet_wrap(~CIS.TRANS) +
  theme_minimal()

# Modify the font size of y-axis tick labels
p <- p + theme(axis.text.y = element_text(size = 6))  # Adjust the size to 6

# Add a slightly darker background for facet headers and adjust the panel border
p <- p +
  theme(strip.background = element_rect(fill = "lightgray", color = "gray"),  # Adjust colors as needed
        panel.border = element_rect(color = "lightgray", fill = NA, size = 1))  # Adjust border color and size

# Save the plot to a PDF file
ggsave("Test_color_gradient_shape_fixed.pdf", plot = p, device = "pdf")
