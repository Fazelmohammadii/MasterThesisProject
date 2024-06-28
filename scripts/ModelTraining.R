library(caret)
library(e1071)
library(ranger)
library(randomForest)
library(pROC)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


##################################################################################################################
################################################ DATA PREPARATION ################################################
##################################################################################################################

##DNMT3A##

#labeling based on the status _ 2 classes 
DN_2 <- dout
DN_2$LABEL <- ifelse(DN_2$status=="CH_only","CH","AML")
#all one hot encoded
DN_2$AA_position <- as.character(DN_2$AA_position) #AA position as categorical variable
DN_2_OH <-DN_2[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
DN_2_OH<- encFunc(DN_2_OH) #one-hot encoding process
DN_2_OH$LABEL <- as.factor(DN_2$LABEL)
#AA position as integer:
DN_2$AA_position <- as.integer(DN_2$AA_position) #AA position as categorical variable
DN_2_int <-DN_2[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
DN_2_int<- encFunc(DN_2_int) #one-hot encoding process
DN_2_int$LABEL <- as.factor(DN_2$LABEL)

#labeling based on the status _ 3 classes 
DN_3 <- dout
#all one hot encoded
DN_3$AA_position <- as.character(DN_3$AA_position) #AA position as categorical variable
DN_3_OH <-DN_3[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
DN_3_OH<- encFunc(DN_3_OH) #one-hot encoding process
DN_3_OH$LABEL <- as.factor(DN_3$status)
#AA position as integer:
DN_3$AA_position <- as.integer(DN_3$AA_position) #AA position as categorical variable
DN_3_int <-DN_3[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
DN_3_int- encFunc(DN_3_int) #one-hot encoding process
DN_3_intLABEL <- as.factor(DN_3$status)


##TET2##

#labeling based on the status _ 2 classes 
TT_2 <- tout
TT_2$LABEL <- ifelse(TT_2$status=="CH_only","CH","AML")
#all one hot encoded
TT_2$AA_position <- as.character(TT_2$AA_position) #AA position as categorical variable
TT_2_OH <-TT_2[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
TT_2_OH<- encFunc(TT_2_OH) #one-hot encoding process
TT_2_OH$LABEL <- as.factor(TT_2$LABEL)
#AA position as integer:
TT_2$AA_position <- as.integer(TT_2$AA_position) #AA position as categorical variable
TT_2_int <-TT_2[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
TT_2_int<- encFunc(TT_2_int) #one-hot encoding process
TT_2_int$LABEL <- as.factor(TT_2$LABEL)

#labeling based on the status _ 3 classes 
TT_3 <- tout
#all one hot encoded
TT_3$AA_position <- as.character(TT_3$AA_position) #AA position as categorical variable
TT_3_OH <-TT_3[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
TT_3_OH<- encFunc(TT_3_OH) #one-hot encoding process
TT_3_OH$LABEL <- as.factor(TT_3$status)
#AA position as integer:
TT_3$AA_position <- as.integer(TT_3$AA_position) #AA position as categorical variable
TT_3_int <-TT_3[,c(-7,-10)] # filter out labels column_avoiding one hot encoding
TT_3_int<- encFunc(TT_3_int) #one-hot encoding process
TT_3_int$LABEL <- as.factor(TT_3$status)


##################################################################################################################
############################################### FEATURE EVALUATION ###############################################
##################################################################################################################


#ROC Plot function
ggrocs <- function(rocs, breaks = seq(0,1,0.1), legendTitel = "Legend") {
  if (length(rocs) == 0) {
    stop("No ROC objects available in param rocs.")
  } else {
    require(plyr)
    # Store all sensitivities and specifivities in a data frame
    # which an be used in ggplot
    RocVals <- plyr::ldply(names(rocs), function(rocName) {
      if(class(rocs[[rocName]]) != "roc") {
        stop("Please provide roc object from pROC package")
      }
      data.frame(
        fpr = rev(rocs[[rocName]]$specificities),
        tpr = rev(rocs[[rocName]]$sensitivities),
        names = rep(rocName, length(rocs[[rocName]]$sensitivities)),
        stringAsFactors = T
      )
    })
    
    aucAvg <- (sapply(rocs, "[[", "auc"))
    
    rocPlot <- ggplot(RocVals, aes(x = fpr, y = tpr, colour = names)) +
      geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5, colour = "gray") + 
      geom_step() +
      scale_x_reverse(name = "False Positive Rate (1 - Specificity)",limits = c(1,0), breaks = breaks) + 
      scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,1), breaks = breaks) +
      theme_bw() + 
      coord_equal() + 
      #annotate("text", x = 0.1, y = 0.1, vjust = 0, label = paste("AUC =",sprintf("%.3f",aucAvg))) +
      guides(colour = guide_legend(legendTitel)) +
      theme(axis.ticks = element_line(color = "grey80"))
    
    rocPlot
  }
}

#Evaluating feature contribution with the model training based on each feature separately
#we select one hot encoded data set with two labels due to more rigid quality of the data 
#using GLM model due to its simplicity comparing to the other algorithms
# Function to extract features based on pattern
extract_features <- function(data, patterns) {
  features <- lapply(patterns, function(pattern) {
    select(data, matches(pattern))
  })
  names(features) <- patterns
  return(features)
}

# Function to train model and calculate ROC curve
train_and_roc <- function(x, y, outer_ctrl, inner_ctrl, tunegrid) {
  model <- train(x = x, y = y, method = "glmnet", trControl = outer_ctrl,
                 tuneGrid = tunegrid, metric = "ROC", tuneControl = inner_ctrl)
  roc_curve <- roc(y, predict(model, newdata = x, type = "prob")[, 2], plot = FALSE)
  return(list(model = model, roc_curve = roc_curve))
}

# Patterns to match columns
patterns <- c("domain", "trimers", "Variant_Classification", "Variant_Type", 
              "AA_position", "subs", "AAref", "AAalt", "polarity_ref", 
              "charge_ref", "sidechain_ref", "pI", "MW", "hydrophobicity", 
              "polarity_alt", "charge_alt", "sidechain_alt")

# Extract features for DNohpcwgtw and TTohpcwgtw datasets
dn_features <- extract_features(DN_2_OH, patterns)
tt_features <- extract_features(TT_2_OH, patterns)

# Add all features excluding the last column
dn_features$all <- DN_2_OH[,-ncol(DN_2_OH)]
tt_features$all <- TT_2_OH[,-ncol(TT_2_OH)]

# Define binary outcome
dn_y <- DN_2_OH$LABEL
tt_y <- TT_2_OH$LABEL

# Define training controls
outer_ctrl <- trainControl(method = "cv", number = 5, allowParallel = TRUE, 
                           verboseIter = FALSE, classProbs = TRUE, summaryFunction = twoClassSummary)
inner_ctrl <- trainControl(method = "cv", number = 5, allowParallel = TRUE, 
                           verboseIter = FALSE, classProbs = TRUE, summaryFunction = twoClassSummary)

# Define hyperparameters to tune
tunegrid <- expand.grid(alpha = 0:10/10, lambda = 10^seq(-4, -1, by = 0.1))

# Train models and calculate ROC curves for DNMT3A dataset
dn_results <- lapply(dn_features, function(x) {
  train_and_roc(x, dn_y, outer_ctrl, inner_ctrl, tunegrid)
})

# Train models and calculate ROC curves for TET2 dataset
tt_results <- lapply(tt_features, function(x) {
  train_and_roc(x, tt_y, outer_ctrl, inner_ctrl, tunegrid)
})

# Extract ROC curves for plotting
dn_rocs <- lapply(dn_results, function(result) result$roc_curve)
tt_rocs <- lapply(tt_results, function(result) result$roc_curve)

# Plot ROC curves
ggrocs(dn_rocs)
ggrocs(tt_rocs)

# Print AUCs for reference
dn_aucs <- sapply(dn_rocs, function(roc_curve) roc_curve$auc)
tt_aucs <- sapply(tt_rocs, function(roc_curve) roc_curve$auc)

print(dn_aucs)
print(tt_aucs)


##################################################################################################################
################################################# MODEL TRAINING #################################################
##################################################################################################################


################                  
#####DNMT3A#####
################                  
# Set seed for reproducibility
set.seed(123)

# Define a function to calculate class weights
calc_weights <- function(y) {
  freq <- table(y)
  weights <- unlist(lapply(y, function(x) {
    weight <- freq[as.character(x)]
    weight <- 1 / weight
    return(weight)
  }))
  return(weights)
}

# Split data into train and test sets
split_data <- function(data, label_col, split_ratio = 0.9) {
  index <- createDataPartition(data[[label_col]], p = split_ratio, list = FALSE)
  train_data <- data[index,]
  test_data <- data[-index,]
  test_data[[label_col]] <- as.factor(test_data[[label_col]])
  return(list(train_data = train_data, test_data = test_data))
}

# Train and evaluate models
train_and_evaluate <- function(train_data, test_data, label_col, tunegrid_RF, tunegrid_GLM, outerCtrl, innerCtrl) {
  x_train <- train_data[,-ncol(train_data)]
  y_train <- train_data[,label_col]
  x_test <- test_data[,-ncol(test_data)]
  y_test <- test_data[,label_col]

  # Calculate class weights
  weights <- calc_weights(y_train)

  # Train Random Forest model
  RF_model <- train(
    x = x_train, y = as.vector(y_train), method = "rf",
    trControl = outerCtrl, tuneGrid = tunegrid_RF,
    trainControl = innerCtrl, weights = weights
  )
  RF_results <- RF_model$results[RF_model$results$Kappa == max(RF_model$results$Kappa),]
  RF_results <- RF_results[1, c("Accuracy", "Kappa", "AccuracySD", "KappaSD")]

  # Train GLM model
  GLM_model <- train(
    x = x_train, y = y_train, method = "glmnet",
    trControl = outerCtrl, tuneGrid = tunegrid_GLM,
    tuneControl = innerCtrl, weights = weights
  )
  GLM_results <- GLM_model$results[GLM_model$results$Kappa == max(GLM_model$results$Kappa),]
  GLM_results <- GLM_results[1, c("Accuracy", "Kappa", "AccuracySD", "KappaSD")]

  # Test models on the unseen test data
  pred_RF <- predict(RF_model, newdata = test_data)
  cm_RF <- confusionMatrix(pred_RF, test_data[[label_col]])
  pred_GLM <- predict(GLM_model, newdata = test_data)
  cm_GLM <- confusionMatrix(pred_GLM, test_data[[label_col]])

  return(list(
    RF_model = RF_model, GLM_model = GLM_model,
    RF_results = RF_results, GLM_results = GLM_results,
    cm_RF = cm_RF, cm_GLM = cm_GLM
  ))
}

# Plot confusion matrices
plot_confusion_matrix <- function(cm, title) {
  df <- data.frame(cm$table)
  ggplot(df, aes(x = Prediction, y = Reference)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    scale_fill_gradient(low = "white", high = "red") +
    geom_text(aes(label = Freq), colour = "black") +
    theme_minimal() +
    ggtitle(title)
}

# Define training controls
innerCtrl <- trainControl(method = "cv", number = 5)
outerCtrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
tunegrid_RF <- expand.grid(mtry = c(sqrt(ncol(train_data))))
tunegrid_GLM <- expand.grid(alpha = 0:10/10, lambda = 10^seq(-4, -1, by = 0.1))

# Load datasets and split into train and test sets
data_DNstchpc <- DN_3_OH
data_DNaastchpc <- DN_3_int
data_tw <- DN_2_OH
data_twaa <- DN_2_int

splits_DNstchpc <- split_data(data_DNstchpc, "LABEL")
splits_DNaastchpc <- split_data(data_DNaastchpc, "LABEL")
splits_tw <- split_data(data_tw, "LABEL")
splits_twaa <- split_data(data_twaa, "LABEL")

# Train and evaluate models for 3-class datasets
results_DNstchpc <- train_and_evaluate(
  splits_DNstchpc$train_data, splits_DNstchpc$test_data, "LABEL",
  tunegrid_RF, tunegrid_GLM, outerCtrl, innerCtrl
)
results_DNaastchpc <- train_and_evaluate(
  splits_DNaastchpc$train_data, splits_DNaastchpc$test_data, "LABEL",
  tunegrid_RF, tunegrid_GLM, outerCtrl, innerCtrl
)

# Combine results
results_3class <- rbind(
  results_DNstchpc$RF_results, results_DNstchpc$GLM_results,
  results_DNaastchpc$RF_results, results_DNaastchpc$GLM_results
)
results_3class$models <- c("RF_OH", "GLM_OH", "RF_AAint", "GLM_AAint")

# Plot results
df_melted_DNstchpcw <- melt(results_3class, id.vars = "models")
ggplot(df_melted_DNstchpcw, aes(x = variable, y = models, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 6) +
  scale_fill_gradient(low = "white", high = "green") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  ggtitle("Model Metrics DNMT3A 3-class") +
  labs(x = "", y = "", fill = "Values") +
  theme(legend.position = "none")

# Plot confusion matrices
p1stchpcw <- plot_confusion_matrix(results_DNstchpc$cm_RF, "RF_all_OH")
p2stchpcw <- plot_confusion_matrix(results_DNstchpc$cm_GLM, "GLM_all_OH")
p1aastchpcw <- plot_confusion_matrix(results_DNaastchpc$cm_RF, "RF_AA_Num")
p2aastchpcw <- plot_confusion_matrix(results_DNaastchpc$cm_GLM, "GLM_AA_Num")
grid.arrange(p1stchpcw, p2stchpcw, p1aastchpcw, p2aastchpcw, ncol = 2, top = "Confusion Matrices of DNMT3A models (3-class)")

# Train and evaluate models for 2-class datasets
results_tw <- train_and_evaluate(
  splits_tw$train_data, splits_tw$test_data, "LABEL",
  tunegrid_RF, tunegrid_GLM, outerCtrl, innerCtrl
)
results_twaa <- train_and_evaluate(
  splits_twaa$train_data, splits_twaa$test_data, "LABEL",
  tunegrid_RF, tunegrid_GLM, outerCtrl, innerCtrl
)

# Combine results
results_2class <- rbind(
  results_tw$RF_results, results_tw$GLM_results,
  results_twaa$RF_results, results_twaa$GLM_results
)
results_2class$models <- c("RF_OH", "GLM_OH", "RF_AAint", "GLM_AAint")

# Plot results
df_melted_DNtw <- melt(results_2class, id.vars = "models")
ggplot(df_melted_DNtw, aes(x = variable, y = models, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 6) +
  scale_fill_gradient(low = "white", high = "green") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  ggtitle("Model Metrics DNMT3A 2-class") +
  labs(x = "", y = "", fill = "Values")

# Plot confusion matrices
p1tw <- plot_confusion_matrix(results_tw$cm_RF, "RF_OH")
p2tw <- plot_confusion_matrix(results_tw$cm_GLM, "GLM_OH")
p1twaa <- plot_confusion_matrix(results_twaa$cm_RF, "RF_aaInt")
p2twaa <- plot_confusion_matrix(results_twaa$cm_GLM, "GLM_aaInt")
grid.arrange(p1tw, p2tw, p1twaa, p2twaa, ncol = 2, top = "Confusion Matrices of DNMT3A models (2-class)")


################                  
######TET2######
################  
# Function to split data into training and testing sets
split_data <- function(data, label_col, split_ratio = 0.9) {
  index <- createDataPartition(data[[label_col]], p = split_ratio, list = FALSE, times = 1)
  train_data <- data[index, ]
  test_data <- data[-index, ]
  test_data[[label_col]] <- as.factor(test_data[[label_col]])
  
  list(
    x_train = train_data[,-ncol(train_data)],
    y_train = train_data[,ncol(train_data)],
    x_test = test_data[,-ncol(test_data)],
    y_test = test_data[,ncol(test_data)]
  )
}

# Function to calculate weights
calc_weights <- function(y) {
  freq <- table(y)
  weights <- unlist(lapply(y, function(x) {
    weight <- freq[as.character(x)]
    1 / weight
  }))
  return(weights)
}

# Function to train models and extract results
train_model <- function(x_train, y_train, model_method, tune_grid, outer_ctrl, inner_ctrl, weights) {
  set.seed(123)
  model <- train(
    x = x_train,
    y = as.vector(y_train),
    method = model_method,
    trControl = outer_ctrl,
    tuneGrid = tune_grid,
    tuneControl = inner_ctrl,
    weights = weights
  )
  best_model <- model$results[model$results$Kappa == max(model$results$Kappa), ]
  best_model[1, c("Accuracy", "Kappa", "AccuracySD", "KappaSD")]
}

# Function to create plots
create_plot <- function(results, title) {
  df_melted <- melt(results, id.vars = "models")
  ggplot(df_melted, aes(x = variable, y = models, fill = value)) + 
    geom_tile() +
    geom_text(aes(label = round(value, digits = 2)), color = "black", size = 6) +
    geom_rect(aes(xmin = as.numeric(variable), xmax = as.numeric(variable),
                  ymin = as.numeric(models), ymax = as.numeric(models)),
              fill = NA, color = "black") +
    scale_fill_gradient(low = "white", high = "green") +
    scale_x_discrete(position = "top") + 
    theme_minimal() +
    ggtitle(title) +
    labs(x = "", y = "", fill = "Values") +
    theme(legend.position = "none")
}

# Function to get confusion matrix metrics
get_cm_metrics <- function(model, x_test, y_test) {
  predictions <- predict(model, newdata = x_test)
  cm <- confusionMatrix(predictions, y_test)
  cm$overall[c('Accuracy', 'Kappa')]
}

# Fully one-hot encoded dataset: (TET2)
data_TTstchpc <- TT_3_OH
split_TTstchpc <- split_data(data_TTstchpc, "LABEL")
TT_weights <- calc_weights(split_TTstchpc$y_train)

RF_TTstchpcw <- train_model(split_TTstchpc$x_train, split_TTstchpc$y_train, "rf", tunegrid_RF, outerCtrl, innerCtrl, TT_weights)
GLM_TTstchpcw <- train_model(split_TTstchpc$x_train, split_TTstchpc$y_train, "glmnet", tunegrid_GLMst, outerCtrl, innerCtrl, TT_weights)

# Numeric value encoded dataset: (TET2)
data_TTaastchpc <- TT_3_int
split_TTaastchpc <- split_data(data_TTaastchpc, "LABEL")

RF_TTaastchpcw <- train_model(split_TTaastchpc$x_train, split_TTaastchpc$y_train, "rf", tunegrid_RF, outerCtrl, innerCtrl, TT_weights)
GLM_TTaastchpcw <- train_model(split_TTaastchpc$x_train, split_TTaastchpc$y_train, "glmnet", tunegrid_GLMst, outerCtrl, innerCtrl, TT_weights)

# Combine results
results_TTstchpcw <- rbind(RF_TTstchpcw, GLM_TTstchpcw, RF_TTaastchpcw, GLM_TTaastchpcw)
rownames(results_TTstchpcw) <- NULL
results_TTstchpcw$models <- c("RF_OH", "GLM_OH", "RF_AAint", "GLM_AAint")

# Plot results
accuracies_TTstchpcw <- create_plot(results_TTstchpcw, "Model_Metrics_TET2_3classes")

# Test confusion matrix and metrics
cm_RF_TTstchpcw <- get_cm_metrics(RF_model_TTstchpcw, split_TTstchpc$x_test, split_TTstchpc$y_test)
cm_GLM_TTstchpcw <- get_cm_metrics(GLM_model_TTstchpcw, split_TTstchpc$x_test, split_TTstchpc$y_test)
cm_RF_TTaastchpcw <- get_cm_metrics(RF_model_TTaastchpcw, split_TTaastchpc$x_test, split_TTaastchpc$y_test)
cm_GLM_TTaastchpcw <- get_cm_metrics(GLM_model_TTaastchpcw, split_TTaastchpc$x_test, split_TTaastchpc$y_test)

cm_test_TTstchpcw <- t(data.frame(
  GLM_TToh = cm_GLM_TTstchpcw,
  GLM_TTaaint = cm_GLM_TTaastchpcw,
  RF_TToh = cm_RF_TTstchpcw,
  RF_TTaaint = cm_RF_TTaastchpcw
))

# Two Class, Weighted
data_tw <- TT_2_OH
split_tw <- split_data(data_tw, "LABEL")
TT_weights_tw <- calc_weights(split_tw$y_train)

RF_tw <- train_model(split_tw$x_train, split_tw$y_train, "rf", tunegrid_RF, outerCtrl, innerCtrl, TT_weights_tw)
GLM_tw <- train_model(split_tw$x_train, split_tw$y_train, "glmnet", tunegrid_GLMst, outerCtrl, innerCtrl, TT_weights_tw)

data_twaa <- TT_2_int
split_twaa <- split_data(data_twaa, "LABEL")

RF_twaa <- train_model(split_twaa$x_train, split_twaa$y_train, "rf", tunegrid_RF, outerCtrl, innerCtrl, TT_weights_tw)
GLM_twaa <- train_model(split_twaa$x_train, split_twaa$y_train, "glmnet", tunegrid_GLMst, outerCtrl, innerCtrl, TT_weights_tw)

results_TTtw <- rbind(RF_tw, GLM_tw, RF_twaa, GLM_twaa)
rownames(results_TTtw) <- NULL
results_TTtw$models <- c("RF_OH", "GLM_OH", "RF_AAint", "GLM_AAint")

# Plot results
accuracies_TTtw <- create_plot(results_TTtw, "Model_Metrics_TET2_2classes")

# Test confusion matrix and metrics
cm_RF_tw <- get_cm_metrics(RF_model_tw, split_tw$x_test, split_tw$y_test)
cm_GLM_tw <- get_cm_metrics(GLM_model_tw, split_tw$x_test, split_tw$y_test)
cm_RF_twaa <- get_cm_metrics(RF_model_twaa, split_twaa$x_test, split_twaa$y_test)
cm_GLM_twaa <- get_cm_metrics(GLM_model_twaa, split_twaa$x_test, split_twaa$y_test)

cm_TT_tw <- t(data.frame(
  RF_OH = c(cm_RF_tw, cm_RF_tw$byClass),
  RF_AA = c(cm_RF_twaa, cm_RF_twaa$byClass),
  GLM_OH = c(cm_GLM_tw, cm_GLM_tw$byClass),
  GLM_AA = c(cm_GLM_twaa, cm_GLM_twaa$byClass)
))

# Create confusion matrix plots
create_cm_plot <- function(cm_data, title) {
  ggplot(data.frame(cm_data$table), aes(x = Prediction, y = Reference)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    scale_fill_gradient(low = "white", high = "red") +
    geom_text(aes(label = Freq), colour = "black") +
    theme_minimal() +
    ggtitle(title)
}

P3stchpcw <- create_cm_plot(cm_RF_TTstchpcw, "RF_all_OH")
P4stchpcw <- create_cm_plot(cm_GLM_TTstchpcw, "GLM_all_OH")
P3aastchpcw <- create_cm_plot(cm_RF_TTaastchpcw, "RF_AA_Num")
P4aastchpcw <- create_cm_plot(cm_GLM_TTaastchpcw, "GLM_AA_Num")

cms_TTstchpcwplot <- grid.arrange(P3stchpcw, P4stchpcw, P3aastchpcw, P4aastchpcw, ncol = 2, top = "Confusion Matrices of TET2 models (3class)")

# Return results
list(
  accuracies_TTstchpcw = accuracies_TTstchpcw,
  cm_test_TTstchpcw = cm_test_TTstchpcw,
  accuracies_TTtw = accuracies_TTtw,
  cm_TT_tw = cm_TT_tw,
  cms_TTstchpcwplot = cms_TTstchpcwplot
)
