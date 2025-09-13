##analyze combined data from Fig.3


models <- rate_mean %>%
  dplyr::filter(lacY > 0.008) %>%
  group_by(Bgal, manganese) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rate ~ lacY, data = .)))

# 1. Take the models data frame
# 2. Augment each model and create a new column 'results'
# 3. Unnest the 'results' to get a single, tidy data frame
models_with_residuals <- models %>%
  mutate(results = map(model, augment)) %>%
  unnest(results)

#Fig.Suppl_to_3----
pdf("Figures/Fig.Suppl_to_3_Residual_Analysis.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(Bgal ~ manganese) +
  labs(
    title = "Residuals vs. Fitted Plot",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()
dev.off()


#Fig.Suppl_to_3----
pdf("Figures/Fig.Suppl_to_3_Q-Q.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(sample = .resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_grid(Bgal ~ manganese) +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw()
dev.off()
