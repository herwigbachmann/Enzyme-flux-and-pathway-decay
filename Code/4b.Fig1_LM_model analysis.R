library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)


models <- rate_mean %>%
  dplyr::filter(lacY > 0.006) %>%
  group_by(preculture, addition) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rate ~ lacY, data = .)))


models_with_residuals <- models %>%
  mutate(results = map(model, augment)) %>%
  unnest(results)

#Fig.Suppl_to_1b----
pdf("Figures/Fig.Suppl_to_1B_Residual_Analysis.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(preculture ~ addition) +
  labs(
    title = "Residuals vs. Fitted Plot",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()
dev.off()


#Fig.Suppl_to_1b----
pdf("Figures/Fig.Suppl_to_1B_Q-Q.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(sample = .resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_grid(preculture ~ addition) +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw()
dev.off()
