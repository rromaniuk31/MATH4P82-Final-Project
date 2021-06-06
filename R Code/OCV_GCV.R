set.seed(2021)
library(psych)
library(fda)

##################
##Getting STA_PF##
##################

data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- data$PF
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))
indep_var <- scale(indep_var)

summary(indep_var[,24])
sta_pf = indep_var[,24]

##GCV
domain = c(-2.5,3)
k = 50
norder = 4
nknots = k + 1
knots = seq(domain[1], domain[2], length.out = nknots)
nbasis = nknots + norder - 2
basis = create.bspline.basis(knots,nbasis,norder)
plot(basis)

basismat  = eval.basis(sta_pf,basis)
V = eval.penalty(basis,int2Lfd(2))

lambda = exp(seq(-0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
gcv = rep(NA,length(lambda))
n = length(dep_var)

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  gcv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-v/n))^2)
}

plot(log(lambda), gcv, type="l", col=2,lty=1, lwd = 2,ylab="GCV")
gcv.lam = lambda[which.min(gcv)]
gcv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var

plot(sta_pf,dep_var)
data = data.frame(sta_pf,dep_var,yhatgcv)
data = data[order(data$sta_pf),]
lines(data$sta_pf,data$yhatgcv,lwd=2,col=4)

##Creating Confidence Bands For GCV
L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatgcv)^2)/(n - 2*v + vtilde)
xdens = seq(-2.5,3,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)%*%dep_var
plot(sta_pf,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make GCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_ppg <- bracket_indep_var$STA_PF
ppg_gcv_pred <- vector()
for(i in 1:length(test_ppg)){
  x <- 2.5 + test_ppg[i]
  z <- x/(5.5/10000)
  ppg_gcv_pred[i] <- fhat[z]
}

ppg_gcv_pred <- as.data.frame(ppg_gcv_pred)
ppg_gcv_pred$Teams <- bracket$Team
ppg_gcv_pred

##Calculate MSE PPG GCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
ppg_gcv_pred <- ppg_gcv_pred[-c(21,22), ] 
rownames(ppg_gcv_pred) <- NULL

mse_ppg_gcv = 0
for(i in 1:nrow(act_scores)){
  mse_ppg_gcv = mse_ppg_gcv + (act_scores[i,2]-ppg_gcv_pred[i,1])^2
}
mse_ppg_gcv

##OCV
lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
ocv = rep(NA,length(lambda))

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  Lii = diag(L)
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  ocv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-Lii))^2)
}

plot(log(lambda), ocv, type="l", col=2,lty=1, lwd = 2,ylab="OCV")
ocv.lam = lambda[which.min(ocv)]
ocv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var

plot(sta_pf,dep_var)
data = data.frame(sta_pf,dep_var,yhatocv)
data = data[order(data$sta_pf),]
lines(data$sta_pf,data$yhatocv,lwd=2,col=4)

# confidence band
# ocv
L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatocv)^2)/(n - 2*v + vtilde)
xdens = seq(-2.5,3,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)%*%dep_var
plot(sta_pf,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make OCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_ppg <- bracket_indep_var$STA_PF
ppg_ocv_pred <- vector()
for(i in 1:length(test_ppg)){
  x <- 2.5 + test_ppg[i]
  z <- x/(5.5/10000)
  ppg_ocv_pred[i] <- fhat[z]
}

ppg_ocv_pred <- as.data.frame(ppg_ocv_pred)
ppg_ocv_pred$Teams <- bracket$Team
ppg_ocv_pred

##Calculate MSE PPG OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
ppg_ocv_pred <- ppg_ocv_pred[-c(21,22), ] 
rownames(ppg_ocv_pred) <- NULL

mse_ppg_ocv = 0
for(i in 1:nrow(act_scores)){
  mse_ppg_ocv = mse_ppg_ocv + (act_scores[i,2]-ppg_ocv_pred[i,1])^2
}
mse_ppg_ocv



## K-Fold Cross Validation (k=10)
#K = 10 
#cv.err.k = array(NA,c(K, length(lambda)))
#cv.lam = rep(NA, length(lambda))

# randomly shuffle 268 observations
#random.sample = sample(1:268, 268, replace = FALSE) 
#lab = c(rep(1:10,rep(26)),1,2,3,4,5,6,7,8)

#for(iterlambda in 1:length(lambda))
#{
#  for(itercv in 1:K)
#  {
#    test.ind = random.sample[which(lab==itercv)]
#    xtest = sta_pf[test.ind]
#    ytest = dep_var[test.ind]
    
    # use training data to fit model
#    basismat.train = basismat[-test.ind,]
#    betahat = ginv(t(basismat.train)%*%basismat.train+lambda[iterlambda]*V)%*%t(basismat.train)%*%dep_var[-test.ind]
    
    # use test data to calculate residuals
#    yhat.test = basismat[test.ind,]%*%betahat
    
#    cv.err.k[itercv, iterlambda] = sum((dep_var[test.ind]-yhat.test)^2)
#  }
#  cv.lam[iterlambda] = mean(cv.err.k[, iterlambda])
#}


#plot(log(lambda), cv.lam, type="l", col=2,lty=1, lwd=2,ylab="CV")

# select lambda
#cv.lam.opt = lambda[which.min(cv.lam)]
#cv.lam.opt

# use the selected lambda to estimate f
#L = basismat%*%ginv(t(basismat)%*%basismat+cv.lam.opt*V)%*%t(basismat)
#yhatcv = L%*%dep_var

#n  = 10000
#x  = seq(domain[1],domain[2],length.out = n)
#y0 = apply(as.matrix(x),1,funf)

#domain = c(-2.5,3)
#funf = function(x){
#  nknots   = 51
#  knots    = seq(domain[1],domain[2],length.out=nknots)
#  norder   = 4
#  nbasis   = nknots + norder - 2 
#  basis    = create.bspline.basis(knots,nbasis,norder)
#  basismat = eval.basis(x, basis)
#  set.seed(2021)
#  b  = rnorm(nbasis)
#  return(basismat%*%b)
#}
  
#plot(sta_pf,dep_var) # plot data
#lines(x,y0,col=1,lwd=2) # add true function f(x)
#lines(sta_pf,yhatcv,col=4, lty=1, lwd=2) # add estimated function with lambda chosen by 10-fold CV

# confidence band
# k-fold
#v = tr(L)
#L = basismat%*%ginv(t(basismat)%*%basismat+cv.lam*V)%*%t(basismat)
#yhatocv = L%*%dep_var
#v = tr(L)
#vtilde = tr(L%*%t(L))
#sigma2 = sum((dep_var-yhatocv)^2)/(n - 2*v + vtilde)
#xdens = seq(-2.5,3,length.out = 10000)
#basismatdens = eval.basis(xdens,basis)
#mtc = ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
#off = rep(NA,length(xdens))
#for(i in 1:length(xdens))
#{
#  lx = basismatdens[i,]%*%mtc
#  sehat = sqrt(sigma2*lx%*%t(lx))
#  off[i] = 1.96*sehat
#}

# plot data together with fhat and its confidence band
#fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)%*%dep_var
#plot(sta_pf,dep_var,col=1)
#lines(xdens,fhat,col=2,lwd=2)
#lines(xdens,fhat+off, lty=2, col=4,lwd=2)
#lines(xdens,fhat-off, lty=2, col=4,lwd=2)

#sigma2


#############################
##Getting Opp_STA_Neutral_W##
#############################

data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- data$PF
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))
indep_var <- scale(indep_var)

summary(indep_var[,70])
opp_neu_w = indep_var[,70]

##GCV
domain = c(-2,3.5)
k = 50
norder = 4
nknots = k + 1
knots = seq(domain[1], domain[2], length.out = nknots)
nbasis = nknots + norder - 2
basis = create.bspline.basis(knots,nbasis,norder)
plot(basis)

basismat  = eval.basis(opp_neu_w,basis)
V = eval.penalty(basis,int2Lfd(2))

lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
gcv = rep(NA,length(lambda))
n = length(dep_var)

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  gcv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-v/n))^2)
}

plot(log(lambda), gcv, type="l", col=2,lty=1, lwd = 2,ylab="GCV")
gcv.lam = lambda[which.min(gcv)]
gcv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var

plot(opp_neu_w,dep_var)
data = data.frame(opp_neu_w,dep_var,yhatgcv)
data = data[order(data$opp_neu_w),]
lines(data$opp_neu_w,data$yhatgcv,lwd=2,col=4)

##Creating Confidence Bands For GCV
L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatgcv)^2)/(n - 2*v + vtilde)
xdens = seq(-2,3.5,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_neu_w,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make GCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_neu_w <- bracket_indep_var$Opp_STA_Neutral_W
opp_neu_w_gcv_pred <- vector()
for(i in 1:length(test_opp_neu_w)){
  x <- 2 + test_opp_neu_w[i]
  z <- x/(5.5/10000)
  opp_neu_w_gcv_pred[i] <- fhat[z]
}

opp_neu_w_gcv_pred <- as.data.frame(opp_neu_w_gcv_pred)
opp_neu_w_gcv_pred$Teams <- bracket$Team
opp_neu_w_gcv_pred

##Calculate MSE Opp Neutral W GCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_neu_w_gcv_pred <- opp_neu_w_gcv_pred[-c(21,22), ] 
rownames(opp_neu_w_gcv_pred) <- NULL

mse_opp_neu_w_gcv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_neu_w_gcv = mse_opp_neu_w_gcv + (act_scores[i,2]-opp_neu_w_gcv_pred[i,1])^2
}
mse_opp_neu_w_gcv

##OCV
lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
ocv = rep(NA,length(lambda))

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  Lii = diag(L)
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  ocv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-Lii))^2)
}

plot(log(lambda), ocv, type="l", col=2,lty=1, lwd = 2,ylab="OCV")
ocv.lam = lambda[which.min(ocv)]
ocv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var

plot(opp_neu_w,dep_var)
data = data.frame(opp_neu_w,dep_var,yhatocv)
data = data[order(data$opp_neu_w),]
lines(data$opp_neu_w,data$yhatocv,lwd=2,col=4)

# confidence band
# ocv
L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatocv)^2)/(n - 2*v + vtilde)
xdens = seq(-2,3.5,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_neu_w,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make OCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_neu_w <- bracket_indep_var$Opp_STA_Neutral_W
opp_neu_w_ocv_pred <- vector()
for(i in 1:length(test_opp_neu_w)){
  x <- 2 + test_opp_neu_w[i]
  z <- x/(5.5/10000)
  opp_neu_w_ocv_pred[i] <- fhat[z]
}

opp_neu_w_ocv_pred <- as.data.frame(opp_neu_w_ocv_pred)
opp_neu_w_ocv_pred$Teams <- bracket$Team
opp_neu_w_ocv_pred

##Calculate MSE Opp Neutral W OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_neu_w_ocv_pred <- opp_neu_w_ocv_pred[-c(21,22), ] 
rownames(opp_neu_w_ocv_pred) <- NULL

mse_opp_neu_w_ocv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_neu_w_ocv = mse_opp_neu_w_ocv + (act_scores[i,2]-opp_neu_w_ocv_pred[i,1])^2
}
mse_opp_neu_w_ocv

#########################
##Getting Opp_KP_ADE_Rk##
#########################

data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- data$PF
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))
indep_var <- scale(indep_var)

summary(indep_var[,95])
opp_ade_rk = indep_var[,95]

##GCV
domain = c(-1,4)
k = 50
norder = 4
nknots = k + 1
knots = seq(domain[1], domain[2], length.out = nknots)
nbasis = nknots + norder - 2
basis = create.bspline.basis(knots,nbasis,norder)
plot(basis)

basismat  = eval.basis(opp_ade_rk,basis)
V = eval.penalty(basis,int2Lfd(2))

lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
gcv = rep(NA,length(lambda))
n = length(dep_var)

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  gcv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-v/n))^2)
}

plot(log(lambda), gcv, type="l", col=2,lty=1, lwd = 2,ylab="GCV")
gcv.lam = lambda[which.min(gcv)]
gcv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var

plot(opp_ade_rk,dep_var)
data = data.frame(opp_ade_rk,dep_var,yhatgcv)
data = data[order(data$opp_ade_rk),]
lines(data$opp_ade_rk,data$yhatgcv,lwd=2,col=4)

##Creating Confidence Bands For GCV
L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatgcv)^2)/(n - 2*v + vtilde)
xdens = seq(-1,4,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_ade_rk,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make GCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_ade_rk <- bracket_indep_var$Opp_KP_ADE_Rk
opp_ade_rk_gcv_pred <- vector()
for(i in 1:length(test_opp_ade_rk)){
  x <- 1 + test_opp_ade_rk[i]
  z <- x/(5/10000)
  opp_ade_rk_gcv_pred[i] <- fhat[z]
}

opp_ade_rk_gcv_pred <- as.data.frame(opp_ade_rk_gcv_pred)
opp_ade_rk_gcv_pred$Teams <- bracket$Team
opp_ade_rk_gcv_pred

##Calculate MSE Opp ADE Rank GCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_ade_rk_gcv_pred <- opp_ade_rk_gcv_pred[-c(21,22), ] 
rownames(opp_ade_rk_gcv_pred) <- NULL

mse_opp_ade_rk_gcv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_ade_rk_gcv = mse_opp_neu_w_gcv + (act_scores[i,2]-opp_ade_rk_gcv_pred[i,1])^2
}
mse_opp_ade_rk_gcv

##OCV
lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
ocv = rep(NA,length(lambda))

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  Lii = diag(L)
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  ocv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-Lii))^2)
}

plot(log(lambda), ocv, type="l", col=2,lty=1, lwd = 2,ylab="OCV")
ocv.lam = lambda[which.min(ocv)]
ocv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var

plot(opp_ade_rk,dep_var)
data = data.frame(opp_ade_rk,dep_var,yhatocv)
data = data[order(data$opp_ade_rk),]
lines(data$opp_ade_rk,data$yhatocv,lwd=2,col=4)

# confidence band
# ocv
L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatocv)^2)/(n - 2*v + vtilde)
xdens = seq(-1,4,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_ade_rk,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make OCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_ade_rk <- bracket_indep_var$Opp_KP_ADE_Rk
opp_ade_rk_ocv_pred <- vector()
for(i in 1:length(test_opp_ade_rk)){
  x <- 1 + test_opp_ade_rk[i]
  z <- x/(5/10000)
  opp_ade_rk_ocv_pred[i] <- fhat[z]
}

opp_ade_rk_ocv_pred <- as.data.frame(opp_ade_rk_ocv_pred)
opp_ade_rk_ocv_pred$Teams <- bracket$Team
opp_ade_rk_ocv_pred

##Calculate MSE Opp ADE Rank OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_ade_rk_ocv_pred <- opp_ade_rk_ocv_pred[-c(21,22), ] 
rownames(opp_ade_rk_ocv_pred) <- NULL

mse_opp_ade_rk_ocv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_ade_rk_ocv = mse_opp_ade_rk_ocv + (act_scores[i,2]-opp_ade_rk_ocv_pred[i,1])^2
}
mse_opp_ade_rk_ocv

#####################
##Getting KP_AOE_Rk##
#####################

data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- data$PF
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))
indep_var <- scale(indep_var)

summary(indep_var[,39])
aoe_rk = indep_var[,39]

##GCV
domain = c(-1,4.5)
k = 50
norder = 4
nknots = k + 1
knots = seq(domain[1], domain[2], length.out = nknots)
nbasis = nknots + norder - 2
basis = create.bspline.basis(knots,nbasis,norder)
plot(basis)

basismat  = eval.basis(aoe_rk,basis)
V = eval.penalty(basis,int2Lfd(2))

lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
gcv = rep(NA,length(lambda))
n = length(dep_var)

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  gcv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-v/n))^2)
}

plot(log(lambda), gcv, type="l", col=2,lty=1, lwd = 2,ylab="GCV")
gcv.lam = lambda[which.min(gcv)]
gcv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var

plot(aoe_rk,dep_var)
data = data.frame(aoe_rk,dep_var,yhatgcv)
data = data[order(data$aoe_rk),]
lines(data$aoe_rk,data$yhatgcv,lwd=2,col=4)

##Creating Confidence Bands For GCV
L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatgcv)^2)/(n - 2*v + vtilde)
xdens = seq(-1,4.5,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)%*%dep_var
plot(aoe_rk,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make GCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_aoe_rk <- bracket_indep_var$KP_AOE_Rk
aoe_rk_gcv_pred <- vector()
for(i in 1:length(test_aoe_rk)){
  x <- 1 + test_aoe_rk[i]
  z <- x/(5.5/10000)
  aoe_rk_gcv_pred[i] <- fhat[z]
}

aoe_rk_gcv_pred <- as.data.frame(aoe_rk_gcv_pred)
aoe_rk_gcv_pred$Teams <- bracket$Team
aoe_rk_gcv_pred

##Calculate MSE AOE Rank GCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
aoe_rk_gcv_pred <- aoe_rk_gcv_pred[-c(21,22), ] 
rownames(aoe_rk_gcv_pred) <- NULL

mse_aoe_rk_gcv = 0
for(i in 1:nrow(act_scores)){
  mse_aoe_rk_gcv = mse_aoe_rk_gcv + (act_scores[i,2]-aoe_rk_gcv_pred[i,1])^2
}
mse_aoe_rk_gcv

##OCV
lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
ocv = rep(NA,length(lambda))

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  Lii = diag(L)
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  ocv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-Lii))^2)
}

plot(log(lambda), ocv, type="l", col=2,lty=1, lwd = 2,ylab="OCV")
ocv.lam = lambda[which.min(ocv)]
ocv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var

plot(aoe_rk,dep_var)
data = data.frame(aoe_rk,dep_var,yhatocv)
data = data[order(data$aoe_rk),]
lines(data$aoe_rk,data$yhatocv,lwd=2,col=4)

# confidence band
# ocv
L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatocv)^2)/(n - 2*v + vtilde)
xdens = seq(-1,4.5,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)%*%dep_var
plot(aoe_rk,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make OCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_aoe_rk <- bracket_indep_var$KP_AOE_Rk
aoe_rk_ocv_pred <- vector()
for(i in 1:length(test_aoe_rk)){
  x <- 1 + test_aoe_rk[i]
  z <- x/(5.5/10000)
  aoe_rk_ocv_pred[i] <- fhat[z]
}

aoe_rk_ocv_pred <- as.data.frame(aoe_rk_ocv_pred)
aoe_rk_ocv_pred$Teams <- bracket$Team
aoe_rk_ocv_pred

##Calculate MSE AOE Rank OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
aoe_rk_ocv_pred <- aoe_rk_ocv_pred[-c(21,22), ] 
rownames(aoe_rk_ocv_pred) <- NULL

mse_aoe_rk_ocv = 0
for(i in 1:nrow(act_scores)){
  mse_aoe_rk_ocv = mse_aoe_rk_ocv + (act_scores[i,2]-aoe_rk_ocv_pred[i,1])^2
}
mse_aoe_rk_ocv





###################################################
##Getting Opp Non-Conference Strength of Schedule##
###################################################

data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- data$PF
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))
indep_var <- scale(indep_var)

summary(indep_var[,102])
opp_ncsos = indep_var[,102]

##GCV
domain = c(-2.5,5)
k = 50
norder = 4
nknots = k + 1
knots = seq(domain[1], domain[2], length.out = nknots)
nbasis = nknots + norder - 2
basis = create.bspline.basis(knots,nbasis,norder)
plot(basis)

basismat  = eval.basis(opp_ncsos,basis)
V = eval.penalty(basis,int2Lfd(2))

lambda = exp(seq(-5,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
gcv = rep(NA,length(lambda))
n = length(dep_var)

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  gcv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-v/n))^2)
}

plot(log(lambda), gcv, type="l", col=2,lty=1, lwd = 2,ylab="GCV")
gcv.lam = lambda[which.min(gcv)]
gcv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var

plot(opp_ncsos,dep_var)
data = data.frame(opp_ncsos,dep_var,yhatgcv)
data = data[order(data$opp_ncsos),]
lines(data$opp_ncsos,data$yhatgcv,lwd=2,col=4)

##Creating Confidence Bands For GCV
L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatgcv)^2)/(n - 2*v + vtilde)
xdens = seq(-2.5,4,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_ncsos,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make GCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_ncsos <- bracket_indep_var$Opp_KP_Non.Conference_Strength_of_Schedule
opp_ncsos_gcv_pred <- vector()
for(i in 1:length(test_opp_ncsos)){
  x <- 2.5 + test_opp_ncsos[i]
  z <- x/(7.5/10000)
  opp_ncsos_gcv_pred[i] <- fhat[z]
}

opp_ncsos_gcv_pred <- as.data.frame(opp_ncsos_gcv_pred)
opp_ncsos_gcv_pred$Teams <- bracket$Team
opp_ncsos_gcv_pred

##Calculate MSE AOE Rank GCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_ncsos_gcv_pred <- opp_ncsos_gcv_pred[-c(21,22), ] 
rownames(opp_ncsos_gcv_pred) <- NULL

mse_opp_ncsos_gcv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_ncsos_gcv = mse_opp_ncsos_gcv + (act_scores[i,2]-opp_ncsos_gcv_pred[i,1])^2
}
mse_opp_ncsos_gcv

##OCV
lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
ocv = rep(NA,length(lambda))

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  Lii = diag(L)
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  ocv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-Lii))^2)
}

plot(log(lambda), ocv, type="l", col=2,lty=1, lwd = 2,ylab="OCV")
ocv.lam = lambda[which.min(ocv)]
ocv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var

plot(opp_ncsos,dep_var)
data = data.frame(opp_ncsos,dep_var,yhatocv)
data = data[order(data$opp_ncsos),]
lines(data$opp_ncsos,data$yhatocv,lwd=2,col=4)

# confidence band
# ocv
L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatocv)^2)/(n - 2*v + vtilde)
xdens = seq(-2.5,5,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_ncsos,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make OCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_ncsos <- bracket_indep_var$KP_Non.Conference_Strength_of_Schedule
opp_ncsos_ocv_pred <- vector()
for(i in 1:length(test_opp_ncsos)){
  x <- 2.5 + test_opp_ncsos[i]
  z <- x/(7.5/10000)
  opp_ncsos_ocv_pred[i] <- fhat[z]
}

opp_ncsos_ocv_pred <- as.data.frame(opp_ncsos_ocv_pred)
opp_ncsos_ocv_pred$Teams <- bracket$Team
opp_ncsos_ocv_pred

##Calculate MSE AOE Rank OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_ncsos_ocv_pred <- opp_ncsos_ocv_pred[-c(21,22), ] 
rownames(opp_ncsos_ocv_pred) <- NULL

mse_opp_ncsos_ocv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_ncsos_ocv = mse_opp_ncsos_ocv + (act_scores[i,2]-opp_ncsos_ocv_pred[i,1])^2
}
mse_opp_ncsos_ocv








###################################
##Getting Opp Adjusted Tempo Rank##
###################################

data <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data[is.na(data)] <- 100

dep_var <- data$PF
indep_var <- subset(data, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))
indep_var <- scale(indep_var)

summary(indep_var[,97])
opp_at_rk = indep_var[,97]

##GCV
domain = c(-2,2)
k = 50
norder = 4
nknots = k + 1
knots = seq(domain[1], domain[2], length.out = nknots)
nbasis = nknots + norder - 2
basis = create.bspline.basis(knots,nbasis,norder)
plot(basis)

basismat  = eval.basis(opp_at_rk,basis)
V = eval.penalty(basis,int2Lfd(2))

lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
gcv = rep(NA,length(lambda))
n = length(dep_var)

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  gcv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-v/n))^2)
}

plot(log(lambda), gcv, type="l", col=2,lty=1, lwd = 2,ylab="GCV")
gcv.lam = lambda[which.min(gcv)]
gcv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var

plot(opp_at_rk,dep_var)
data = data.frame(opp_at_rk,dep_var,yhatgcv)
data = data[order(data$opp_at_rk),]
lines(data$opp_at_rk,data$yhatgcv,lwd=2,col=4)

##Creating Confidence Bands For GCV
L = basismat%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
yhatgcv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatgcv)^2)/(n - 2*v + vtilde)
xdens = seq(-2,2,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+gcv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_at_rk,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make GCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_at_rk <- bracket_indep_var$Opp_KP_AT_Rk
opp_at_rk_gcv_pred <- vector()
for(i in 1:length(test_opp_at_rk)){
  x <- 2 + test_opp_at_rk[i]
  z <- x/(4/10000)
  opp_at_rk_gcv_pred[i] <- fhat[z]
}

opp_at_rk_gcv_pred <- as.data.frame(opp_at_rk_gcv_pred)
opp_at_rk_gcv_pred$Teams <- bracket$Team
opp_at_rk_gcv_pred

##Calculate MSE AOE Rank GCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_at_rk_gcv_pred <- opp_at_rk_gcv_pred[-c(21,22), ] 
rownames(opp_at_rk_gcv_pred) <- NULL

mse_opp_at_rk_gcv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_at_rk_gcv = mse_opp_at_rk_gcv + (act_scores[i,2]-opp_at_rk_gcv_pred[i,1])^2
}
mse_opp_at_rk_gcv

##OCV
lambda = exp(seq(0,10,length.out = 2000))
df = RSS = rep(NA,length(lambda))
ocv = rep(NA,length(lambda))

for(iterlambda in 1:length(lambda))
{
  L = basismat%*%ginv(t(basismat)%*%basismat+lambda[iterlambda]*V)%*%t(basismat)
  v = tr(L)
  
  yhat = L%*%dep_var
  Lii = diag(L)
  df[iterlambda]  = v
  RSS[iterlambda] = sum((dep_var-yhat)^2)
  ocv[iterlambda] = 1/n*sum(((dep_var-yhat)/(1-Lii))^2)
}

plot(log(lambda), ocv, type="l", col=2,lty=1, lwd = 2,ylab="OCV")
ocv.lam = lambda[which.min(ocv)]
ocv.lam

L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var

plot(opp_at_rk,dep_var)
data = data.frame(opp_at_rk,dep_var,yhatocv)
data = data[order(data$opp_at_rk),]
lines(data$opp_at_rk,data$yhatocv,lwd=2,col=4)

# confidence band
# ocv
L = basismat%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
yhatocv = L%*%dep_var
v = tr(L)
vtilde = tr(L%*%t(L))
sigma2 = sum((dep_var-yhatocv)^2)/(n - 2*v + vtilde)
xdens = seq(-2,2,length.out = 10000)
basismatdens = eval.basis(xdens,basis)
mtc = ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)
off = rep(NA,length(xdens))
for(i in 1:length(xdens))
{
  lx = basismatdens[i,]%*%mtc
  sehat = sqrt(sigma2*lx%*%t(lx))
  off[i] = 1.96*sehat
}

# plot data together with fhat and its confidence band
fhat  = basismatdens%*%ginv(t(basismat)%*%basismat+ocv.lam*V)%*%t(basismat)%*%dep_var
plot(opp_at_rk,dep_var,col=1)
lines(xdens,fhat,col=2,lwd=2)
lines(xdens,fhat+off, lty=2, col=4,lwd=2)
lines(xdens,fhat-off, lty=2, col=4,lwd=2)

sigma2

#########################
## Make OCV Predictions##
#########################
data1 <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Dataset.csv')
data1[is.na(data1)] <- 100

dep_var1 <- data1$PF
indep_var1 <- subset(data1, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,Opponent,PF,PA,WL,Won_Game,OppTeamID,Opp_TeamID,Opp_TS_Team,Opp_KP_Conference,GameID))

bracket <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Bracket.csv')
bracket[is.na(bracket)] <- 100

bracket_indep_var <- subset(bracket, select = -c(TournamentYear,TeamID,TS_Team,KP_Conference,Round,Round_Number,Team,OppTeamID,Opp_KP_Conference,OppTeam,Opp_TS_Team))

for(j in 1:nrow(bracket_indep_var)){
  bracket_indep_var$Adjustment[j] <- mean(data1$TS_GP)/bracket_indep_var$TS_GP[j]
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

test_opp_at_rk <- bracket_indep_var$Opp_KP_AT_Rk
opp_at_rk_ocv_pred <- vector()
for(i in 1:length(test_opp_at_rk)){
  x <- 2 + test_opp_at_rk[i]
  z <- x/(4/10000)
  opp_at_rk_ocv_pred[i] <- fhat[z]
}

opp_at_rk_ocv_pred <- as.data.frame(opp_at_rk_ocv_pred)
opp_at_rk_ocv_pred$Teams <- bracket$Team
opp_at_rk_ocv_pred

##Calculate MSE AOE Rank OCV
act_scores <- read.csv(file = 'D:\\Fourth Year\\Winter\\MATH 4P82 (Nonparametric Statistics)\\Final Project\\Actual Scores.csv')
opp_at_rk_ocv_pred <- opp_at_rk_ocv_pred[-c(21,22), ] 
rownames(opp_at_rk_ocv_pred) <- NULL

mse_opp_at_rk_ocv = 0
for(i in 1:nrow(act_scores)){
  mse_opp_at_rk_ocv = mse_opp_at_rk_ocv + (act_scores[i,2]-opp_at_rk_ocv_pred[i,1])^2
}
mse_opp_at_rk_ocv
