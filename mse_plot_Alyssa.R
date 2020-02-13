library(tidyverse)

### load data
load("~/Desktop/SPRING 2020/Computing/group_project1_computing/equalweak.RData")

#### updated metrics function
metrics = function(models, p = 50, n.strong, n.weak.ind, n.weak.cor){
  truth.df = lapply(models, `[[`, 1)
  forward = lapply(models, `[[`, 2)
  lasso = lapply(models, `[[`, 3)
  beta.true = lapply(models, `[[`, 4)
  beta.forward = lapply(models, `[[`, 5)
  beta.lasso = lapply(models, `[[`, 6)
  weak.preds = lapply(models, `[[`, 7)
  
  n.true = sum(n.strong, n.weak.ind, n.weak.cor)
  
  condition = 4*sqrt(log(p)/1000)
  
  
  
  ## Size of models
  forward.size = do.call(rbind,
                         lapply(1:N, function(i){
                           sum(forward[i][[1]])
                         })
  )
  lasso.size = do.call(rbind,
                       lapply(1:N, function(i){
                         sum(lasso[i][[1]])
                       })
  )
  
  ## Number of strong predictors in model
  forward.strong.n = do.call(rbind,
                             lapply(1:N, function(i){
                               strong.index = which(beta.true[i][[1]] > condition)
                               model.index = which(forward[i][[1]] == TRUE)
                               
                               sum(model.index %in% strong.index)
                             })
  )
  lasso.strong.n = do.call(rbind,
                           lapply(1:N, function(i){
                             strong.index = which(beta.true[i][[1]] > condition)
                             model.index = which(lasso[i][[1]] == TRUE)
                             
                             sum(model.index %in% strong.index)
                           })
  ) 
  
  ## Average percentage of true predictors chosen by each model
  forward.true = do.call(rbind,
                         lapply(1:N, function(i){
                           ind = which(truth.df[[i]] == 1) # index of which parameters are true
                           forward.nonnull = forward[[i]][ind] # did the forward method select them?
                           
                           return(forward.nonnull)
                         })
  )
  forward.pcttrue = mean(rowSums(forward.true)) / n.true
  
  lasso.true = do.call(rbind,
                       lapply(1:N, function(i){
                         ind = which(truth.df[[i]] == 1) # index of which parameters are true
                         lasso.nonnull = lasso[[i]][ind] # did the lasso method select them?
                         
                         return(lasso.nonnull)
                       })
  )
  lasso.pcttrue = mean(rowSums(lasso.true)) / n.true
  
  ## Average percentage of true predictors chosen by each model
  forward.false = do.call(rbind,
                          lapply(1:N, function(i){
                            ind = which(truth.df[[i]] == 0) # index of which parameters are null
                            forward.nonnull = forward[[i]][ind] # did the forward method select them?
                            
                            return(forward.nonnull)
                          })
  )
  forward.pctnull = mean(rowSums(forward.false)) / n.true
  
  lasso.false = do.call(rbind,
                        lapply(1:N, function(i){
                          ind = which(truth.df[[i]] == 0) # index of which parameters are null
                          lasso.nonnull = lasso[[i]][ind] # did the lasso method select them?
                          
                          return(lasso.nonnull)
                        })
  )
  lasso.pctnull = mean(rowSums(lasso.false)) / n.true
  
  
  
  
  ### percentage of weak predictors
  forward.weak = do.call(rbind,
                         lapply(1:N, function(i){
                           ind = weak.preds[[i]] # index of which parameters are null
                           forward.nonnull = forward[[i]][ind] # did the forward method select them?
                           
                           return(forward.nonnull)
                         })
  )
  forward.pctweak = rowSums(forward.weak) / (n.weak.cor + n.weak.ind)
  
  lasso.weak = do.call(rbind,
                       lapply(1:N, function(i){
                         ind = weak.preds[[i]]  # index of which parameters are null
                         lasso.nonnull = lasso[[i]][ind] # did the lasso method select them?
                         
                         return(lasso.nonnull)
                       })
  )
  lasso.pctweak = rowSums(lasso.weak) / (n.weak.cor + n.weak.ind)
  
  ## MSE for strong vs weak vs null predictors
  forward. = do.call(rbind,
                             lapply(1:N, function(i){
                               strong.index = which(beta.true[i][[1]] > condition)
                               model.index = which(forward[i][[1]] == TRUE)
                               
                               sum(model.index %in% strong.index)
                             })
  )
  lasso.strong.n = do.call(rbind,
                           lapply(1:N, function(i){
                             strong.index = which(beta.true[i][[1]] > condition)
                             model.index = which(lasso[i][[1]] == TRUE)
                             
                             sum(model.index %in% strong.index)
                           })
  ) 
  
  MSE.df = do.call(rbind, lapply(1:N, function(i){
    params.forward = beta.forward[[i]]
    params.lasso = beta.lasso[[i]]
    b.true = beta.true[[i]]
    strong.index = which(beta.true[i][[1]] > condition)
    weak.index = which(beta.true[i][[1]] <= condition & beta.true[i][[1]] > 0)
    null.index = which(beta.true[i][[1]] == 0)
    
    # get beta estimates from each model
    names(params.forward) = as.numeric(gsub("X", "", names(params.forward)))
    forward.betas = data.frame(X = as.numeric(names(params.forward)),
                               b.forward = params.forward)
    names(params.lasso) = as.numeric(gsub("V", "", names(params.lasso)))
    lasso.betas = data.frame(X = as.numeric(names(params.lasso)),
                             b.lasso = params.lasso)
    
    # calculated MSE across betas for each model
    beta.estimates = data.frame(X = 1:p,
                                b.true = b.true) %>% # true betas
      left_join(forward.betas, by = "X") %>% # forward selection
      left_join(lasso.betas, by = "X") %>% # lasso 
      mutate_all(funs(replace_na(., 0))) %>% # let missing params be 0
      mutate(MSE.forward = (b.forward - b.true)^2, # calculate MSE of each beta
             MSE.lasso = (b.lasso - b.true)^2) %>%
      summarise(avg.MSE.forward = mean(MSE.forward, na.rm = TRUE), # average MSE across betas
                avg.MSE.lasso = mean(MSE.lasso, na.rm = TRUE),
                # calculate MSE for strong and weak predictors
                avg.forward.strong = mean((b.forward[strong.index] - b.true[strong.index])^2, na.rm = TRUE),
                avg.forward.weak = mean((b.forward[weak.index] - b.true[weak.index])^2, na.rm = TRUE),
                avg.lasso.strong = mean((b.lasso[strong.index] - b.true[strong.index])^2, na.rm = TRUE),
                avg.lasso.weak = mean((b.lasso[weak.index] - b.true[weak.index])^2, na.rm = TRUE),
                avg.forward.null = mean((b.forward[null.index] - b.true[null.index])^2, na.rm = TRUE),
                avg.lasso.null = mean((b.lasso[null.index] - b.true[null.index])^2, na.rm = TRUE)) %>%
      unlist
    
    return(beta.estimates)
  }) 
  ) %>%
    cbind(., forward.pctweak) %>%
    cbind(., lasso.pctweak)
  
  return(list(
    model.size = c(forward = forward.size, lasso = lasso.size),
    n.weak = c(forward = rowSums(forward.weak), lasso = rowSums(lasso.weak)),
    n.strong = c(forward = forward.strong.n, lasso = lasso.strong.n),
    pcttrue = c(forward = forward.pcttrue, lasso = lasso.pcttrue),
    pctnull = c(forward = forward.pctnull, lasso = lasso.pctnull),
    MSE.df = MSE.df
  ))
}



### recalculate metrics

metrics1 = metrics(models = models.equ.weak1,
                    p = 50,
                    n.strong = 8,
                    n.weak.ind = 1,
                    n.weak.cor = 1)
metrics4 = metrics(models = models.equ.weak4,
                   p = 50,
                   n.strong = 8,
                   n.weak.ind = 4,
                   n.weak.cor = 4)
metrics7 = metrics(models = models.equ.weak7,
                     p = 50,
                     n.strong = 8,
                     n.weak.ind = 7,
                     n.weak.cor = 7)
metrics10 = metrics(models = models.equ.weak10,
                     p = 50,
                     n.strong = 8,
                     n.weak.ind = 10,
                     n.weak.cor = 10)



# df for model size
model.size = rbind(
  matrix(data = metrics1$model.size, nrow = 2000, ncol = 1),
  matrix(data = metrics4$model.size, nrow = 2000, ncol = 1),
  matrix(data = metrics7$model.size, nrow = 2000, ncol = 1),
  matrix(data = metrics10$model.size, nrow = 2000, ncol = 1)
) %>%
  as.data.frame %>%
  rename(model.size = "V1") %>%
  mutate(model = rep(c(rep("forward", 1000), rep("lasso", 1000)), 4))

# df for number of weak predictors included
n.weak.model = rbind(
  matrix(data = metrics1$n.weak, nrow = 2000, ncol = 1),
  matrix(data = metrics4$n.weak, nrow = 2000, ncol = 1),
  matrix(data = metrics7$n.weak, nrow = 2000, ncol = 1),
  matrix(data = metrics10$n.weak, nrow = 2000, ncol = 1)
) %>%
  as.data.frame %>%
  rename(n.weak = "V1") 

# df for number of strong predictors included
n.strong.model = rbind(
  matrix(data = metrics1$n.strong, nrow = 2000, ncol = 1),
  matrix(data = metrics4$n.strong, nrow = 2000, ncol = 1),
  matrix(data = metrics7$n.strong, nrow = 2000, ncol = 1),
  matrix(data = metrics10$n.strong, nrow = 2000, ncol = 1)
) %>%
  as.data.frame %>%
  rename(n.strong = "V1") 

# combine above df's 
missing_weak_df = model.size %>%
  bind_cols(n.weak.model) %>%
  bind_cols(n.strong.model) %>%
  mutate(n.weak.total = c(rep(2, 2000),
                          rep(8, 2000),
                          rep(14, 2000),
                          rep(20, 2000))) %>%
  mutate(n.missing = n.weak.total - n.weak)

# retrieve MSE from scenarios
mse1 = metrics1$MSE.df[,c(1,2)] %>%
  as.data.frame() %>%
  rename(forward = "avg.MSE.forward", lasso = "avg.MSE.lasso") %>%
  gather(key = "model", value = "MSE")
mse4 = metrics4$MSE.df[,c(1,2)] %>%
  as.data.frame() %>%
  rename(forward = "avg.MSE.forward", lasso = "avg.MSE.lasso") %>%
  gather(key = "model", value = "MSE")
mse7 = metrics7$MSE.df[,c(1,2)] %>%
  as.data.frame() %>%
  rename(forward = "avg.MSE.forward", lasso = "avg.MSE.lasso") %>%
  gather(key = "model", value = "MSE")
mse10 = metrics10$MSE.df[,c(1,2)] %>%
  as.data.frame() %>%
  rename(forward = "avg.MSE.forward", lasso = "avg.MSE.lasso") %>%
  gather(key = "model", value = "MSE")

mse.df = bind_rows(mse1, mse4, mse7, mse10)

# df of MSE to plot with
plot.df = bind_cols(missing_weak_df, mse.df) %>%
  select(-model1)

## PLOT proportion missing versus MSE
plot.df %>%
  mutate(prop.missing = n.missing/n.weak.total) %>% # proportion of weak predictors absent from model
  group_by(prop.missing, n.weak.total, model) %>%
  summarise(mse = mean(MSE),
            n = n(),
            size = mean(model.size),
            strong = mean(n.strong)) %>%
  ggplot(aes(x = prop.missing, y = mse, size = strong, color = model)) + 
  geom_point() +
  #scale_size(breaks = c(5:8), range = c(1, 7)) +
  facet_grid(. ~ n.weak.total) +
  labs(
    x = "Proportion of all weak predictors that are absent",
    y = "MSE",
    size = "# strong predictors in model"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


# MSE of strong predictors by each model
mse.strong1 = metrics1$MSE.df[,c(3, 5)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.strong", lasso = "avg.lasso.strong") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "strong",
         n.weak = 2)
mse.strong4 = metrics4$MSE.df[,c(3, 5)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.strong", lasso = "avg.lasso.strong") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "strong",
         n.weak = 8)
mse.strong7 = metrics7$MSE.df[,c(3, 5)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.strong", lasso = "avg.lasso.strong") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "strong",
         n.weak = 14)  
mse.strong10 = metrics10$MSE.df[,c(3, 5)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.strong", lasso = "avg.lasso.strong") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "strong",
         n.weak = 20)

# MSE of weak predictors by each model
mse.weak1 = metrics1$MSE.df[,c(4, 6)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.weak", lasso = "avg.lasso.weak") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "weak",
         n.weak = 2)
mse.weak4 = metrics4$MSE.df[,c(4, 6)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.weak", lasso = "avg.lasso.weak") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "weak",
         n.weak = 8)
mse.weak7 = metrics7$MSE.df[,c(4, 6)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.weak", lasso = "avg.lasso.weak") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "weak",
         n.weak = 14)  
mse.weak10 = metrics10$MSE.df[,c(4, 6)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.weak", lasso = "avg.lasso.weak") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "weak",
         n.weak = 20)

# MSE of null predictors by each model
mse.null1 = metrics1$MSE.df[,c(7, 8)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.null", lasso = "avg.lasso.null") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "null",
         n.weak = 2)
mse.null4 = metrics4$MSE.df[,c(7, 8)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.null", lasso = "avg.lasso.null") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "null",
         n.weak = 8)
mse.null7 = metrics7$MSE.df[,c(7, 8)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.null", lasso = "avg.lasso.null") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "null",
         n.weak = 14)  
mse.null10 = metrics10$MSE.df[,c(7, 8)] %>%
  as.data.frame() %>%
  rename(forward = "avg.forward.null", lasso = "avg.lasso.null") %>%
  gather(key = "model", value = "MSE") %>%
  mutate(predictor = "null",
         n.weak = 20)


mse_predtype_df = bind_rows(mse.strong1, mse.strong4, mse.strong7, mse.strong10,
                            mse.weak1, mse.weak4, mse.weak7, mse.weak10,
                            mse.null1, mse.null4, mse.null7, mse.null10)

## PLOT of MSE for strong and weak predictor by model across scenarios
mse_predtype_df %>%
  group_by(predictor, model, n.weak) %>%
  summarise(mse = mean(MSE)) %>%
  ggplot(aes(x = n.weak, y = mse, color = predictor)) +
  geom_point() +
  geom_line() +
  facet_grid(.~model) +
  labs(
    x = "Total number of weak predictors",
    y = "MSE of beta coefficients estimates"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
  
