set.seed(2021)
data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- subset(data, select = PF)
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

for(j in 1:nrow(indep_var)){
  for(k in 1:ncol(indep_var))
  {
    mean_col <- mean(indep_var[,k])
    indep_var[j,k] <- (indep_var[j,k]-mean(indep_var[,k]))/sd(indep_var[,k])
  }
}

indep_var <- as.data.frame(indep_var)
MLR<-lm(data$PF ~ STA_PF + Opp_KP_ADE_Rk + Opp_KP_Non.Conference_Strength_of_Schedule + KP_Ajusted_Efficiency_Margin + TS_3PA + KP_Non.Conference_Strength_of_Schedule + AP_Rnk + KP_AOOE_Rk + Opp_KP_AOOE_Rk + STA_Top25_L + Opp_AP_Rnk + Opp_STA_Home_L + Opp_STA_Top25_W + Opp_KP_AOW_Rk + Opp_KP_Luck_Rk + KP_AODE_Rk + KP_AT_Rk + Opp_KP_AT_Rk + KP_AOE_Rk + Opp_STA_Neutral_W,data=indep_var)
summary(MLR)

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

bracket_indep_var <- bracket_indep_var[c("STA_PF","Opp_KP_ADE_Rk","Opp_KP_Non.Conference_Strength_of_Schedule","KP_Ajusted_Efficiency_Margin","TS_3PA","KP_Non.Conference_Strength_of_Schedule","AP_Rnk","KP_AOOE_Rk","Opp_KP_AOOE_Rk","STA_Top25_L","Opp_AP_Rnk","Opp_STA_Home_L","Opp_STA_Top25_W","Opp_KP_AOW_Rk","Opp_KP_Luck_Rk","KP_AODE_Rk","KP_AT_Rk","Opp_KP_AT_Rk","KP_AOE_Rk","Opp_STA_Neutral_W")]

indep_var1 <- indep_var1[c("STA_PF","Opp_KP_ADE_Rk","Opp_KP_Non.Conference_Strength_of_Schedule","KP_Ajusted_Efficiency_Margin","TS_3PA","KP_Non.Conference_Strength_of_Schedule","AP_Rnk","KP_AOOE_Rk","Opp_KP_AOOE_Rk","STA_Top25_L","Opp_AP_Rnk","Opp_STA_Home_L","Opp_STA_Top25_W","Opp_KP_AOW_Rk","Opp_KP_Luck_Rk","KP_AODE_Rk","KP_AT_Rk","Opp_KP_AT_Rk","KP_AOE_Rk","Opp_STA_Neutral_W")]

for(j in 1:nrow(bracket_indep_var)){
  for(k in 1:ncol(bracket_indep_var))
  {
    bracket_indep_var[j,k] <- (bracket_indep_var[j,k]-mean(indep_var1[,k]))/sd(indep_var1[,k])
  }
}

pred<-predict(MLR,bracket_indep_var)
pred <- as.data.frame(pred)
pred$Teams <- bracket$Team
pred

##Calculate MSE AOE Rank OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
pred <- pred[-c(21,22), ] 
rownames(pred) <- NULL

mse_linear = 0
for(i in 1:nrow(act_scores)){
  mse_linear = mse_linear + (act_scores[i,2]-pred[i,1])^2
}
mse_linear

