
library(lme4)
library(lmerTest)

library(gtsummary)

library(tidyverse)
library(broom)
library(caret)
library(pROC)

set.seed(1219)

get_logistic_pred = function(mod, data, res = "y", pos = 1, neg = 0, cut = 0.5) {
  probs = predict(mod, newdata = data, type = "response")
  ifelse(probs > cut, pos, neg)
}

# set number of regions as 2
nregions = 5
nfamilies = 100
randintvar = 3
randintvar2 = 3
intercept = 0
groupeffect = 1
residvar = 1
nperregion = 500
bx1 = 0.25
bfem = 0.01



replications = 1000
npatients = rep(nperregion, nregions)
region = seq(1,nregions)
sampsize <- sum(npatients)
train_index <- as.logical(rbinom(sampsize, 1, .75))

# using purrr:::map to generate datasets
res1 <- data.frame(is = c(1:replications)) %>%
  mutate(
    dat = map(is, function(.x){
      set.seed(.x)
      # make random intercept for each region
      randint_reg = rnorm(nregions, sd = sqrt(randintvar))
      # keep same number of patients per region for now
      
      data1 <- as.data.frame(cbind(region, randint_reg, npatients))
      data2 <- slice(data1, rep(seq(nregions), npatients))
      # make groups of correlated data within region called "families"
      # keeping just 2 families for now for simplicity
      # make random intercept for each family
      randint2 = data.frame(
        familyid = c(1:(nregions*nfamilies)),
        randintfam = rnorm(nregions*nfamilies, sd = sqrt(randintvar2)))
      # now collect the data together
      simple <- data2 %>%
        mutate(familyid = rep(seq(nregions*nfamilies), c(rep(sampsize/(nregions*nfamilies),nregions*nfamilies)))) %>%
        left_join(randint2, by = "familyid") %>%
        mutate(familyid = as.factor(familyid))
      finaldat <- simple %>%
        # make normal predictor
        mutate(x1 = (rnorm(sampsize)),
               # make binary demographic variable (male/female)
               female = rep(c(0,1), times = sampsize/2)) %>%
        mutate(
          # now make linear predictor model
          z = randint_reg + randintfam + bx1*x1 + bfem*female + + (rnorm(sampsize, sd = sqrt(residvar))), # add random effects, fixed effects, and random noise
          # trying out 2 ways to make binary outcome
          pr1 = exp(z)/(1+exp(z)),
          pr2 = 1/(1+exp(-z)),
          runis = runif(sampsize, 0, 1),
          y1 = ifelse(runis < pr1, 1, 0), # pulled from Ken's blog
          y2 = rbinom(sampsize,1,pr2))
      finaldat
    }))

start_time_glm <- Sys.time()
# now add glmer model and pull out comptime, accuracy, auc etc
res2 <- res1 %>%
  mutate(
    mod_glmer = map(dat, possibly(function(.x){
      # make a glmer model 
      start_time <- Sys.time()
      mod1 <- glmer(y1~ x1 + female + (1|region) + (1|familyid), .x[train_index,], family = "binomial")
      end_time <- Sys.time()
      comptime = end_time - start_time
      
      test_pred_50 = get_logistic_pred(mod1, data = .x[!train_index,], res = "default", cut = 0.5)
      test_tab_50 = table(predicted = test_pred_50, actual = .x$y1[!train_index])
      
      test_con_mat_50 = confusionMatrix(test_tab_50)
      
      test_prob = predict(mod1, newdata = .x, type = "response")
      test_roc = roc(.x$y1 ~ test_prob, plot = FALSE)
      #test_roc$auc
      
      considerations = c(test_con_mat_50$overall["Accuracy"],
                         test_con_mat_50$byClass["Sensitivity"],
                         test_con_mat_50$byClass["Specificity"],
                         comptime = comptime,
                         auc = test_roc$auc
      )
      considerations
    }
    )
    ) # possibly
  )
end_time_glm = Sys.time()

time_glm = end_time_glm - start_time_glm

save(res2, file = "scenario1glm.Rdata")