install.packages("glmnet")
set.seed(2021)
library(glmnet)
data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- data$PF
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

for(j in 1:nrow(indep_var)){
  for(k in 1:ncol(indep_var))
  {
    mean_col <- mean(indep_var[,k])
    indep_var[j,k] <- (indep_var[j,k]-mean(indep_var[,k]))/sd(indep_var[,k])
  }
}

indep_var_mat = as.matrix(indep_var)
lambdas <- 10^seq(-4, 2, by = 0.01)

# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(indep_var_mat, dep_var, alpha = 1, lambda = lambdas, standardize = FALSE)

# Best 
lambda_best <- lasso_reg$lambda.min 
lambda_best

lasso_model <- glmnet(indep_var_mat, dep_var, alpha = 1, lambda = lambda_best, standardize = FALSE)
coef(lasso_model)

##MAKE PREDICTIONS
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')

bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data$TS_GP)/bracket_indep_var$TS_GP[j]
  bracket_indep_var$TS_GP[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$TS_GP[j]
  bracket_indep_var$STA_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_W[j]
  bracket_indep_var$STA_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_L[j]
  bracket_indep_var$STA_Conf_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Conf_W[j]
  bracket_indep_var$STA_Conf_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Conf_L[j]
  bracket_indep_var$STA_Home_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Home_W[j]
  bracket_indep_var$STA_Home_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Home_L[j]
  bracket_indep_var$STA_Away_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Away_W[j]
  bracket_indep_var$STA_Away_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Away_L[j]
  bracket_indep_var$STA_Neutral_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Neutral_W[j]
  bracket_indep_var$STA_Neutral_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Neutral_L[j]
  bracket_indep_var$STA_Top25_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Top25_W[j]
  bracket_indep_var$STA_Top25_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Top25_L[j]
  bracket_indep_var$STA_Conf_PF[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Conf_PF[j]
  bracket_indep_var$STA_Conf_PA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_Conf_PA[j]
  bracket_indep_var$STA_PF[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_PF[j]
  bracket_indep_var$STA_PA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$STA_PA[j]
  bracket_indep_var$TS_FGM[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$TS_FGM[j]
  bracket_indep_var$TS_FGA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$TS_FGA[j]
  bracket_indep_var$TS_3PM[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$TS_3PM[j]
  bracket_indep_var$TS_3PA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$TS_3PA[j]
  bracket_indep_var$TS_FTM[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$TS_FTM[j]
  bracket_indep_var$TS_FTA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$TS_FTA[j]
  bracket_indep_var$Opp_TS_GP[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_TS_GP[j]
  bracket_indep_var$Opp_STA_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_W[j]
  bracket_indep_var$Opp_STA_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_L[j]
  bracket_indep_var$Opp_STA_Conf_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Conf_W[j]
  bracket_indep_var$Opp_STA_Conf_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Conf_L[j]
  bracket_indep_var$Opp_STA_Home_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Home_W[j]
  bracket_indep_var$Opp_STA_Away_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Away_L[j]
  bracket_indep_var$Opp_STA_Neutral_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Neutral_W[j]
  bracket_indep_var$Opp_STA_Neutral_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Neutral_L[j]
  bracket_indep_var$Opp_STA_Top25_W[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Top25_W[j]
  bracket_indep_var$Opp_STA_Top25_L[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Top25_L[j]
  bracket_indep_var$Opp_STA_Conf_PF[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Conf_PF[j]
  bracket_indep_var$Opp_STA_Conf_PA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_Conf_PA[j]
  bracket_indep_var$Opp_STA_PF[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_PF[j]
  bracket_indep_var$Opp_STA_PA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_STA_PA[j]
  bracket_indep_var$Opp_TS_FGM[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_TS_FGM[j]
  bracket_indep_var$Opp_TS_FGA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_TS_FGA[j]
  bracket_indep_var$Opp_TS_3PM[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_TS_3PM[j]
  bracket_indep_var$Opp_TS_3PA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_TS_3PA[j]
  bracket_indep_var$Opp_TS_FTM[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_TS_FTM[j]
  bracket_indep_var$Opp_TS_FTA[j] <- bracket_indep_var$Adjustment[j]*bracket_indep_var$Opp_TS_FTA[j]
}

bracket_indep_var <- subset(bracket_indep_var, select = -c(Adjustment))

bracket_indep_var <- bracket_indep_var[names(indep_var)]

for(j in 1:nrow(bracket_indep_var)){
  for(k in 1:ncol(bracket_indep_var))
  {
    bracket_indep_var[j,k] <- (bracket_indep_var[j,k]-mean(indep_var1[,k]))/sd(indep_var1[,k])
  }
}

#bracket_pred<-predict(MLR,MLR_A_indep_var)
#MLR_A_pred <- as.data.frame(MLR_A_pred)
#MLR_A_pred$Teams <- MLR_A$Team
#MLR_A_pred

bracket_indep_var_mat = as.matrix(bracket_indep_var)
pred <- predict(lasso_model, s = lambda_best, newx =bracket_indep_var_mat)
pred <- as.data.frame(pred)
pred$Teams <- bracket$Team
pred

##Calculate MSE AOE Rank OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
pred <- pred[-c(21,22), ] 
rownames(pred) <- NULL

mse_lasso = 0
for(i in 1:nrow(act_scores)){
  mse_lasso = mse_lasso + (act_scores[i,2]-pred[i,1])^2
}
mse_lasso

