# processing of simulation results 
# model: Poisson - Poisson


library(dplyr)

tab<-readRDS("P-P_symul.rds")

table(tab$n)
table(tab$p)
table(tab$lambda1)
table(tab$lambda2)

summary(tab)


d <- filter(tab,n==250,p==0.3,lambda1==2,lambda2==2)
boxplot(d[,-c(1:4)])




quality2<-function(data,x_name,param_name)
# quality indicators of the 'param' parameter estimator
{
  tab<-data[[x_name]]
  param<-data[[param_name]][1]
  MEAN=mean(tab)
  BIAS=MEAN-param
  BIASp=100*BIAS/param
  ABSMAX=max(tab-param)
  ABSMAXp=100*ABSMAX/param
  ABSBIAS=mean(abs(tab-param))
  ABSBIASp=100*ABSBIAS/param
  MSE=mean((tab-param)^2)
  RMSE=sqrt(MSE)
  RMSEp=100*RMSE/param
      
  return(data.frame(MEAN=MEAN,BIAS=BIAS,BIASp=BIASp,ABSMAX=ABSMAX,ABSMAXp=ABSMAXp,
                    ABSBIAS=ABSBIAS,ABSBIASp=ABSBIASp,
                    MSE=MSE,RMSE=RMSE,RMSEp=RMSEp))
  
}



model<-"Poisson"

grp<-group_by(tab,n,p,lambda1,lambda2)
summarize(grp,n1=n())

qq1<-do(grp,quality2(.,"mm","p"))  %>% mutate(EST="MM")
qq2<-do(grp,quality2(.,"ml","p"))  %>% mutate(EST="ML")
qq3<-do(grp,quality2(.,"mm0","p"))  %>% mutate(EST="MM0")
qq4<-do(grp,quality2(.,"ml0","p"))  %>% mutate(EST="ML0")

qq<-rbind(qq1,qq2,qq3,qq4)
qq$EST<-factor(qq$EST,levels=c("MM","ML","MM0","ML0"))
  
table(qq$EST)

# quality measures for tables in publication
write.csv2(qq,"P-P_qq.csv")




