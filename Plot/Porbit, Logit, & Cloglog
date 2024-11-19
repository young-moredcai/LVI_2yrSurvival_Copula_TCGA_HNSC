# Load required library
library(ggplot2)

# Define a sequence of x values
x <- seq(-5, 5, length.out = 100)

# Compute the CDFs for probit, cloglog, and logit
data <- data.frame(
  x = rep(x, 3),
  CDF = c(pnorm(x), plogis(x), 1 - exp(-exp(x))),
  Type = factor(rep(c("Probit (Normal CDF)", "Logit (Logistic CDF)", "Cloglog"), each = length(x)))
)

# Create the plot
ggplot(data, aes(x = x, y = CDF, color = Type)) +
  geom_line(size = 1) +
  labs(
    title = "Comparison of CDFs: Probit, Logit, and Cloglog",
    x = "x",
    y = "CDF",
    color = "Link Function"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )
