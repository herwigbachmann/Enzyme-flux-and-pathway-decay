library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)


##analyze Mn data model

models <- rate_mean_Mn %>%
  dplyr::filter(lacY > 0.008) %>%
  group_by(strain, manganese) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rate ~ lacY, data = .)))


models_with_residuals <- models %>%
  mutate(results = map(model, augment)) %>%
  unnest(results)

#Fig.Suppl_to_2a----
pdf("Figures/Fig.Suppl_to_2a_Residual_Analysis.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(strain ~ manganese) +
  labs(
    title = "Residuals vs. Fitted Plot",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()
dev.off()


#Fig.Suppl_to_2a----
pdf("Figures/Fig.Suppl_to_2a_Q-Q.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(sample = .resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_grid(strain ~ manganese) +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw()
dev.off()








#######analyze ATPase data model


models <- rate_mean_ATPase %>%
  dplyr::filter(lacY > 0.008) %>%
  group_by(strain, Bgal) %>%
  nest() %>%
  mutate(model = map(data, ~lm(rate ~ lacY, data = .)))


models_with_residuals <- models %>%
  mutate(results = map(model, augment)) %>%
  unnest(results)

#Fig.Suppl_to_2c----
pdf("Figures/Fig.Suppl_to_2c_Residual_Analysis.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(. ~ Bgal) +
  labs(
    title = "Residuals vs. Fitted Plot",
    x = "Fitted Values",
    y = "Residuals"
  ) +
  theme_bw()
dev.off()


#Fig.Suppl_to_2a----
pdf("Figures/Fig.Suppl_to_2c_Q-Q.pdf",width=5,height=5)
ggplot(models_with_residuals, aes(sample = .resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_grid(. ~ Bgal) +
  labs(
    title = "Normal Q-Q Plot",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_bw()
dev.off()

