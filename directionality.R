top20_loss <- imp_df %>%
  arrange(desc(loss)) %>%
  slice(1:20) %>%
  mutate(
    Importance    = loss,
    ImportanceLab = round(loss, 2)
  )
direction_df <- df_rf_ma %>%
  filter(pct_bmi_outcome_WM %in% c("loss", "gain")) %>%
  group_by(pct_bmi_outcome_WM) %>%
  summarise(across(all_of(top_loss_vars), mean, na.rm = TRUE)) %>%
  pivot_longer(-pct_bmi_outcome_WM, names_to = "Predictor", values_to = "Mean") %>%
  pivot_wider(names_from = pct_bmi_outcome_WM, values_from = Mean) %>%
  mutate(
    Direction     = loss - gain,  # positive = higher in "loss"
    DirectionSign = if_else(Direction > 0, "+", "-")
  )

direction_df %>% select(Predictor, Direction, DirectionSign)

top20_loss <- left_join(top20_loss,direction_df[, c("Predictor", "Direction", "DirectionSign")],by = "Predictor")

write.csv(top20_loss, "top20_loss_directions.csv", row.names = FALSE)
write.csv(df_rf_ma, "rf_model_input_data.csv", row.names = FALSE)