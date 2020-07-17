#libraries required 
#y ~ f(x)
#x=vector de variables (predictores)
#y=vector de numeros (Tran)
library(tidyverse)
library(caret)
library(leaps)
#linear model con regsubsets
 models <- regsubsets(Tran~., data = soil1a, nvmax = 10)
#summary of the model 
summary(models)
#plot models
plot(models, scale='r2')
#visualize result
res.sum <- summary(models)
data.frame(
Adj.R2 = which.max(res.sum$adjr2),
CP = which.min(res.sum$cp),
BIC = which.min(res.sum$bic)
)
# id: model id
# object: regsubsets object
# data: data used to fit regsubsets
# outcome: outcome variable
get_model_formula <- function(id, object, outcome){
  # get models data
  models <- summary(object)$which[id,-1]
  # Get outcome variable
  #form <- as.formula(object$call[[2]])
  #outcome <- all.vars(form)[1]
  # Get model predictors
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  # Build model formula
  as.formula(paste0(outcome, "~", predictors))
}

#visualize model formulas 
get_model_formula(3, models, "Tran")

get_cv_error <- function(model.formula, data){
  set.seed(1)
  train.control <- trainControl(method = "repeatedcv", number = 5, repeats=5)
  cv <- train(model.formula, data = data, method = "lm",
              trControl = train.control)
  cv$results$RMSE
}

model.ids <- 1:10
cv.errors <-  map(model.ids, get_model_formula, models, "Tran") %>%
  map(get_cv_error, data = soil1a) %>%
  unlist()
cv.errors
print(paste0('model ', which.min(cv.errors), ' yields lowest error'))
# Select the model that minimize the CV error
coef(models, which.min(cv.errors))

get_model_formula(5, models, "Tran")

