evaluate_model <- function(optimal_ranger, data, reference_class, plt_title) {
  library(ranger)
  library(readr)
  source("scripts/ggplotConfusionMatrix.R")
  
  # calculate prediction accuracy
  pred <- predict(optimal_ranger, data)
  pred_tab <- table(observed = reference_class, predicted = pred$predictions)
  print(pred_tab)
  
  # calculate performance
  performance <- sum(diag(pred_tab)) / sum(pred_tab)
  print(paste0("ACCURACY = ", performance))

  # plot pretty confusion matrix
  cm <- caret::confusionMatrix(data = pred$predictions, reference = reference_class)
  ggplotConfusionMatrix(cm, plot_title = plt_title)
  # ggsave(filename = confusion_pdf, plot = plt, scale = 1, width = 6, height = 4, dpi = 300)
}
