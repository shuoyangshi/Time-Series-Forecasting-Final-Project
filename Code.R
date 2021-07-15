library(readxl)
library(astsa)
library(tseries)
library(forecast)
library(stats)
library(fGarch)
library(TSA)
library(dynlm)
library(xts)
library(modelr)
library(dplyr)

setwd("/Users/Hedy/Desktop/Final Project")

data_carbon <- read.csv("Carbon dioxide emissions in Hawaii.csv")
data_econ <- read_excel("economic-indicators-time-series.xlsx")
#data_legal <- read_excel("legal-reserve-ratio-for-germany.xlsx")
data_mfr <- read_excel("manufacturers-index-of-new-order.xlsx")
data_nzb <- read.csv("NZBirths.csv")
data_pub <- read_excel("pub-prinv.xlsx")
#data_public <- read_excel("public-construction-for-united-s.xlsx")
#data_total <- read_excel("total-hours-worked-per-month-per.xlsx")
data_wheather <- read_excel("wheatherPr.xlsx")

### Data 1 ###

data_carbon <- read.csv("Carbon dioxide emissions in Hawaii.csv")
carbon <- data_carbon$Carbondioxide
carbon <- ts(carbon, start=c(1959,1), frequency=12)
ts.plot(carbon)

### (1) use an (S)ARIMA approach
# this data looks like have constant variance
# first, detrending and differencing
fit_c <- lm(carbon ~ time(carbon), na.action=NULL)
d_c <- diff(carbon)
dd_c <- diff(diff(carbon, 12))
par(mfrow=c(3,1))
plot(resid(fit_c), main="detrended")
plot(d_c, main="first difference")
plot(dd_c, main="12+1 difference")
# We go for differencing
# plot ACFs
acf2(carbon, 48, main="carbon")
acf2(d_c, 48, main="first difference")
acf2(dd_c, 48, main="12+1 difference")
# The resulting series seems now stationary: at nonseasonal level, its acf cutting off after lag 1 and 
# pacf cutting off after lag 1 suggesting an ARMA(1,1).
# At seasonal level: acf dying down and pacf dying down.
# Here I will try three models as following
c_model1 <- sarima(carbon, 1, 1, 1, 1, 1, 1, 12)
c_model2 <- sarima(carbon, 1, 1, 1, 1, 1, 0, 12)
c_model3 <- sarima(carbon, 1, 1, 1, 0, 1, 1, 12)
c_model1$AIC; c_model2$AIC; c_model3$AIC
# model3: ARIMA(1,1,1)(0,1,1)[12] is the best model.
c_model4 <- sarima(carbon, 0, 1, 1, 0, 1, 1, 12)
c_model5 <- sarima(carbon, 1, 1, 0, 0, 1, 1, 12)
c_model1$AIC; c_model2$AIC; c_model3$AIC; c_model4$AIC; c_model5$AIC
c_model1$BIC; c_model2$BIC; c_model3$BIC; c_model4$BIC; c_model5$BIC
# The AIC prefer the ARIMA(1,1,1)(0,1,1)[12], whereas the BIC prefers the simpler ARIMA(0,1,1)(0,1,1)[12] model.
# It is often the case that the BIC will select a model of smaller order than the AIC or AICc.
# For now, I would say ARIMA(1,1,1)(0,1,1)[12] is better after comparing the diagnostic check.
sarima.for(carbon, 4, 1, 1, 1, 0, 1, 1, 12)

### (2) unit root test, de-trending the data if suitable, and the possibility of finding a suitable ARMA model

adf.test(resid(fit_c), k=0) # DF test
adf.test(resid(fit_c)) # ADF test
pp.test(resid(fit_c)) # PP test
# The Dickey-Fuller test confirms that the data has a unit root (P-value=.4)

#carbongr <- diff(log(carbon)) # growth rate
#plot(carbongr)
#acf2(carbongr)

### (3) spectral analysis

# method 1
n <- length(carbon)
periodogram(d_c, xlab="Frequency of diff(carbon)")
omega <- order(periodogram(d_c, plot=FALSE)$spec, decreasing=T)[1:2]/n
(1/omega)/12 # years
# method 2
mvspec(carbon, log='n')
# There are two other peaks: one at 1 simply referring at the seasonal nature of the data (as they are monthly). 
# Another small one at 2 which indicates a possible 1/2=0.5 year cycle in the Carbon dioxide emissions in Hawaii.
# All mehtods agree the same result.

wave=cbind(cos((1:n)*2*pi*omega[1]),sin((1:n)*2*pi*omega[1]),
            cos((1:n)*2*pi*omega[2]),sin((1:n)*2*pi*omega[2]))
fit1=lm(carbon~wave[,1]+wave[,2]+wave[,3]+wave[,4])
summary(fit1)
plot(ts(fit1$residuals))

fitted_lm=fit1$fitted.values
fit0=sarima(carbon,1,1,1,xreg=ts(wave))
fitted_xreg=carbon-fit0$fit$residuals
fitted_signal=SigExtract(carbon)

plot(carbon)
lines(fitted_xreg,col="red") #best
lines(fitted_signal,col="blue") #worst
lines(fitted_lm,col="green")  #moderate

### Data 2 ###

data_mfr <- read_excel("manufacturers-index-of-new-order.xlsx")
mfr <- data_mfr$`Manufacturers' Index of New Orders of Durable Goods for United States`
mfr <- ts(mfr, start=c(1929,1), frequency=12)
ts.plot(mfr)

### (1) use an (S)ARIMA approach
# first, detrending and differencing
fit_m <- lm(mfr ~ time(mfr), na.action=NULL)
d_m <- diff(mfr)
par(mfrow=c(2,1))
plot(resid(fit_m), main="detrended")
plot(d_m, main="first difference")
# plot ACFs
acf2(mfr, main="mfr")
acf2(resid(fit_m), main="detrended")
acf2(d_m, main="first difference")
# We go for differencing
# its acf cutting off after lag 2 and pacf cutting off after lag 2 suggesting an ARMA(2,2).
# We will try following models
m_model1 <- sarima(mfr, p=2, d=1, q=2)
m_model2 <- sarima(mfr, p=1, d=1, q=2)
m_model3 <- sarima(mfr, p=2, d=1, q=1)
m_model4 <- sarima(mfr, p=1, d=1, q=1)
m_model1$AIC; m_model2$AIC; m_model3$AIC; m_model4$AIC
m_model1$BIC; m_model2$BIC; m_model3$BIC; m_model4$BIC
# The AIC prefer the ARIMA(2,1,1) fit, whereas the BIC prefers the simpler ARIMA(1,1,1) model. It is often
# the case that the BIC will select a model of smaller order than the AIC or AICc.
# For now, I would say ARIMA(2,1,1) is better after comparing the diagnostic check.
sarima.for(mfr, 4, p=2, d=1, q=1)

### (2) unit root test, de-trending the data if suitable, and the possibility of finding a suitable ARMA model

adf.test(mfr, k=0) # DF test
adf.test(mfr) # ADF test
pp.test(mfr) # PP test
mvspec(mfr, log='n')
# The Dickey-Fuller test confirms that the data has a unit root (P-value=.48)

#mfrgr <- diff(log(mfr)) # growth rate
#plot(mfrgr)
#acf2(mfrgr)

### (3) spectral analysis

# method 1
n <- length(mfr)
omega=order(periodogram(mfr)$spec,decreasing = T)[1:2]/n
(1/omega)/12 # years
# method 2
mvspec(mfr, log='n')
# Method 2 is hard to find the result. From method 1, we could get there are 11-year and 5.5-year cycles in 
# the Manufacturers' Index of New Orders of Durable Goods for United States.

wave=cbind(cos((1:n)*2*pi*omega[1]),sin((1:n)*2*pi*omega[1]),
           cos((1:n)*2*pi*omega[2]),sin((1:n)*2*pi*omega[2]))
fit1=lm(mfr~wave[,1]+wave[,2]+wave[,3]+wave[,4])
summary(fit1)
plot(ts(fit1$residuals))

fitted_lm=fit1$fitted.values
fit0=sarima(mfr,2,1,1,xreg=ts(wave))
fitted_xreg=mfr-fit0$fit$residuals
fitted_signal=SigExtract(mfr)

plot(mfr)
lines(fitted_xreg,col="red") #best
lines(fitted_signal,col="blue") #worst
lines(fitted_lm,col="green")  #moderate


### Data 3 ###

data_econ <- read_excel("economic-indicators-time-series.xlsx")
private <- data_econ$`Total Private Construction (Millions of Dollars)`
public <- data_econ$`Total Public Construction (Millions of Dollars)`
private <- ts(private, start=c(2002,1), frequency=12)
public <- ts(public, start=c(2002,1), frequency=12)
par(mfrow=c(2,1))
ts.plot(private)
ts.plot(public)

### (1) linear regression with dependent errors

plot(x=data_econ$`Total Public Construction (Millions of Dollars)`,
     y=data_econ$`Total Private Construction (Millions of Dollars)`)
# It seems have a convex, we can have a try
trend <- time(public)
pub <- public-mean(public)
pub2 <- pub^2
fit1 <- lm(private ~ trend + pub, na.action=NULL)
fit2 <- lm(private ~ trend + pub + pub2, na.action=NULL)
summary(aov(lm(private ~ cbind(trend, pub, pub2))))

num <- length(private) # sample size
AIC(fit1)/num - log(2*pi) # AIC
BIC(fit1)/num - log(2*pi) # BIC
AIC(fit2)/num - log(2*pi) # AIC
BIC(fit2)/num - log(2*pi) # BIC
# fit 2's AIC and BIC are smaller, so fit 2 is a better model.
summary(fit2)
# The linear regression model would be private = 4538000 + -2229*trend + 1.283*pub + -0.0001614*pub2
# Residual standard error is 10100, which means it is a really bad model.

### (2) lagged regression and TFM

# method 1:
lag2.plot(public, private, 41, corr=FALSE)

# It shows fairly strong linear relationships between Private, PRIt, and the Public series at many lags,
# such as PUBt-7, PUBt-8, PUBt-15, PUBt-16, PUBt-17, PUBt-19, etc. Indicating the coefficients are negative,
# implying that increases in the Public lead to decreases in the Private. However, some plots have very
# dispersive points. In order to find a good fit, I would try PUBt-18, PUBt-19, PUBt-30 and PUBt-31, which
# have more agminated points.

fit3 <- dynlm(private ~ L(public,18))
summary(fit3)
fit4 <- dynlm(private ~ L(public,19))
summary(fit4)
fit5 <- dynlm(private ~ L(public,30))
summary(fit5)
fit6 <- dynlm(private ~ L(public,31))
summary(fit6)
fit7 <- dynlm(private ~ L(public,18) + L(public,19) + L(public,30) + L(public,31))
summary(fit7)
# fit7 is best for now.

# method 2
par(mfrow=c(3,1))
acf(private, 48, main="Total Private Construction (Millions of Dollars)")
acf(public, 48, main="Total Public Construction (Millions of Dollars)")
ccf(public, private, 48, main="Private vs Public", ylab="CCF")
# Both of the ACFs exhibit periodicities corresponding to the correlation between values separated by 12 units.
# Observations separated by six months are negatively correlated, showing that positive excursions tend to be 
# associated with negative excursions six months removed.
# however, shows some departure from the cyclic component of each series and there is an obvious peak at h = −6,
# h = -18, h = -30, h = -42. This result implies that public measured at time t−6, t-18, t-30, t-42 months
# are associated with the private series at time t.
fit8 <- dynlm(private ~ L(public,6) + L(public,18) + L(public,30) + L(public,42))
summary(fit8)
# fit8 is best for now.
fit9 <- dynlm(private ~ L(public,6) + L(public,18) + L(public,30) + L(public,42) +
                + L(public,19) + L(public,30) + L(public,31))
summary(fit9)
# fit9 is best for now.

# method 3
summary(LagReg(public, private, L=15, M=32, threshold=1))
summary(LagReg(private, public, L=15, M=32, inverse=TRUE, threshold=0.1))
mse(fit9, as.data.frame(cbind(private, public))) # 271361832
# The second LagReg is best.

### (3) coherence analysis (cross periodogram)

econ_data<-as.data.frame(cbind(private, public))
econ=mvspec(econ_data,spans=c(3,3),taper=.5)
plot(econ,plot.type="coh",ci=-1)

econ$df
f = qf(.999, 2, econ$df-2)
C = f/(18+f) 

### Data 4 ###

data_wheather <- read_excel("wheatherPr.xlsx")
temp <- data_wheather$Temp
dp <- data_wheather$DewPt
cc <- data_wheather$CldCvr
ws <- data_wheather$WndSpd
pr <- data_wheather$Precip

### (1) linear regression with dependent errors

pairs(cbind(TEMP=temp, DP=dp, CC=cc, WS=ws, P=sqrt(pr)))

trend <- time(temp)
P <- sqrt(pr)
fit1 <- lm(temp ~ trend + dp + cc + ws + pr, na.action=NULL)
fit2 <- lm(temp ~ trend + dp + cc + ws + pr + P, na.action=NULL)
summary(aov(lm(temp ~ cbind(trend, dp, cc, ws, pr, P))))

num <- length(temp) # sample size
AIC(fit1)/num - log(2*pi) # AIC
BIC(fit1)/num - log(2*pi) # BIC
AIC(fit2)/num - log(2*pi) # AIC
BIC(fit2)/num - log(2*pi) # BIC

summary(fit2)
# fit 2 is a good fit.

### (2) lagged regression and TFM

lag2.plot(dp, temp, 41, corr=FALSE)
lag2.plot(cc, temp, 41, corr=FALSE)
lag2.plot(ws, temp, 41, corr=FALSE)
lag2.plot(P, temp, 41, corr=FALSE)

# It shows fairly strong linear relationships between Temp(Tt) and the DewPt, CldCvr, WndSpd, sqrt(Precip)
# series at some lags, such as DPt-38, CCt-1 , WSt-9, P-8.

fit3 <- dynlm(temp ~ L(dp,38) + L(cc,1) + L(ws,9) + L(P,8))
summary(fit3)
fit4 <- dynlm(temp ~ L(dp,38) + L(cc,1) + L(P,8))
summary(fit4)
# fit4 is better, but fit2 is best for now.

### (3) coherence analysis (cross periodogram)

wheather_trans = rename(as.data.frame(cbind(temp, dp, cc, ws, pr)),
                        TEMP=temp, DP=dp, CC=cc, WS=ws, P=pr)
wheather_trans[,"P"]=sqrt(wheather_trans[,"P"])
weather=mvspec(wheather_trans,spans=c(6,6),taper=.5)
plot(weather,plot.type="coh",ci=-1)
weather$df
f = qf(.999, 2, weather$df-2)
C = f/(18+f)


### Data 5 ###

data_nzb <- read.csv("NZBirths.csv")
# I would study two variables
mm <- ts(data_nzb$MaoriMale, start=c(2000,1), frequency=4)
tm <- ts(data_nzb$TotalMale, start=c(2000,1), frequency=4)
mf <- ts(data_nzb$MaoriFemale, start=c(2000,1), frequency=4)
tf <- ts(data_nzb$TotalFemale, start=c(2000,1), frequency=4)
ts.plot(mm, tm, gpars=list(col=c("black","red")), main="MaoriMale vs TotalMale")
ts.plot(mf, tf, gpars=list(col=c("black","red")), main="MaoriFemale vs TotalFemale")

### (1) linear regression with dependent errors

pairs(cbind(MM=data_nzb$MaoriMale, TM=data_nzb$TotalMale, MF=data_nzb$MaoriFemale, TF=data_nzb$TotalFemale))

fit1 <- lm(tm ~ mm + mf + tf, na.action=NULL)
summary(fit1)
summary(aov(fit1))

### (2) lagged regression and TFM

summary(LagReg(mm, tm, L=15, M=32, threshold=0.2))
summary(LagReg(mf, tm, L=15, M=32, threshold=0.2))
summary(LagReg(tf, tm, L=15, M=32, threshold=0.2))

fit2 <- dynlm(tm ~ mm + L(mm,1) + mf + L(mf,1) + L(mf,2) + tf)
summary(fit2)

summary(LagReg(mm, tm, L=15, M=32, threshold=0.1))
summary(LagReg(mf, tm, L=15, M=32, threshold=0.1))
summary(LagReg(tf, tm, L=15, M=32, threshold=0.1))

fit3 <- dynlm(tm ~ mm + L(mm,1) + L(mm,5) + L(mm,7) + L(mm,9) + L(mm,11) + L(mm,13) + L(mm,15) +
                mf + L(mf,1) + L(mf,2) + L(mf,3) + L(mf,5) + L(mf,7) + L(mf,9) + L(mf,11) + 
                L(mf,13) + L(mf,15) + tf)
summary(fit3)

fit4 <- dynlm(tm ~ mm + tf)
summary(fit4)

### (3) coherence analysis (cross periodogram)

nzb=mvspec(data_nzb[-1],spans=c(5,5),taper=.5)
plot(nzb,plot.type="coh",ci=-1)
nzb$df
f = qf(.999, 2, nzb$df-2)
C = f/(18+f)



### Data 6 ###

data_pub <- read_excel("pub-prinv.xlsx")
govinv <- ts(data_pub$govinv, start=c(1948,3), frequency=4)
prinv <- ts(data_pub$prinv, start=c(1948,3), frequency=4)
ts.plot(govinv, prinv, gpars=list(col=c("black","red")), main="govinv vs prinv")

### (1) linear regression with dependent errors

plot(x=data_pub$prinv, y=data_pub$govinv)
fit1 <- lm(govinv ~ prinv)
summary(fit1)

### (2) lagged regression and TFM

# method 1
summary(LagReg(prinv, govinv, L=15, M=32, threshold=0.02))
summary(LagReg(govinv, prinv, L=15, M=32, inverse=TRUE, threshold=0.02))
# The first one is more simple and mse is smaller.
fit2 <- dynlm(govinv ~ prinv + L(prinv,1) +  L(prinv,2) +  L(prinv,3) +  L(prinv,5) +  L(prinv,8))
summary(fit2)
fit3 <- dynlm(govinv ~ prinv + L(prinv,8))
summary(fit3)

# method 2, to be continued
lag2.plot(prinv, govinv, 19, corr=FALSE)
# due to method 1, we focus on prinv(t-8), there are two solpes we need to find out (when x>400 or <400)
dummy = ifelse(prinv<400, 0, 1)
priL8 <- stats::lag(prinv,-8)
dL8 = stats::lag(dummy,-8)
inv = ts.intersect(govinv, prinv, priL8, dL8, dframe=TRUE)
fit4 <- lm(govinv ~ priL8*dL8, data=inv, na.action=NULL)
summary(fit4)
fit5 <- lm(govinv ~ prinv + priL8*dL8, data=inv, na.action=NULL)
summary(fit5)
AIC(fit5)

# TFM, input:prinv, output:govinv

pri_d = resid(lm(prinv~time(prinv), na.action=NULL)) # detrended prinv 
acf2(pri_d)
fit = arima(pri_d, order=c(1,0,0))
ar1 = as.numeric(coef(fit)[1]) # = 0.8912
pri_pw = resid(fit)
govi_fil = stats::filter(govinv, filter=c(1, -ar1), sides=1)
ccf(pri_pw, govi_fil, ylab="CCF", na.action=na.omit, panel.first=grid(), 60)
# Noting the apparent shift of d = 8 quarters
gov_pr = ts.intersect(govinv, govL1=stats::lag(govinv,-1), priL8=stats::lag(pri_d,-8))
u = lm(gov_pr[,1]~gov_pr[,2:3], na.action=NULL)
acf2(resid(u)) # suggests ar1
arx = sarima(gov_pr[,1], 1, 0, 0, xreg=gov_pr[,2:3]) # final model
arx$AIC
pred = govinv + resid(arx$fit) # 1-step-ahead predictions
ts.plot(pred, govinv, col=c('gray90',1), lwd=c(7,1))

### (3) coherence analysis (cross periodogram)

pub=mvspec(data_pub[-1],spans=c(3,3),taper=.5)
plot(pub,plot.type="coh",ci=-1)

pub$df
f = qf(.999, 2, pub$df-2)
C = f/(18+f)

