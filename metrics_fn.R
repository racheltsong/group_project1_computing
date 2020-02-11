metrics = function(models, p = 50, n.strong, n.weak.ind, n.weak.cor){
  truth.df = lapply(models, `[[`, 1)
  forward = lapply(models, `[[`, 2)
  lasso = lapply(models, `[[`, 3)
  beta.true = lapply(models, `[[`, 4)
  beta.forward = lapply(models, `[[`, 5)
  beta.lasso = lapply(models, `[[`, 6)
  weak.preds = lapply(models, `[[`, 7)
  
  n.true = sum(n.strong, n.weak.ind, n.weak.cor)
  
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
  
  # 
  # ### percentage of weak predictors
  # forward.strong = do.call(rbind,
  #                         lapply(1:N, function(i){
  #                           ind = strong.preds[[i]] # index of which parameters are null
  #                           forward.nonnull = forward[[i]][ind] # did the forward method select them?
  #                           
  #                           return(forward.nonnull)
  #                         })
  # )
  # forward.pctstrong = rowSums(forward.strong) / (n.strong.cor + n.strong.ind)
  # 
  # lasso.strong = do.call(rbind,
  #                       lapply(1:N, function(i){
  #                         ind = strong.preds[[i]]  # index of which parameters are null
  #                         lasso.nonnull = lasso[[i]][ind] # did the lasso method select them?
  #                         
  #                         return(lasso.nonnull)
  #                       })
  # )
  # lasso.pctstrong = rowSums(lasso.strong) / (n.strong.cor + n.strong.ind)
  # 
  # 
  
  MSE.df = do.call(rbind, lapply(1:N, function(i){
    params.forward = beta.forward[[i]]
    params.lasso = beta.lasso[[i]]
    b.true = beta.true[[i]]
    
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
                avg.MSE.lasso = mean(MSE.lasso, na.rm = TRUE)) %>%
      unlist
    
    return(beta.estimates)
  }) 
  ) %>%
    cbind(., forward.pctweak) %>%
    cbind(., lasso.pctweak)
  
  return(list(
    model.size = c(forward = forward.size, lasso = lasso.size),
    n.weak = c(forward = rowSums(forward.weak), lasso = rowSums(lasso.weak)),
    pcttrue = c(forward = forward.pcttrue, lasso = lasso.pcttrue),
    pctnull = c(forward = forward.pctnull, lasso = lasso.pctnull),
    MSE.df = MSE.df
  ))
}