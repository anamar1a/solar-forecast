library(forecast)

flux<-read.csv("/data/sw/Flux_2010_2017_allY.csv")
sixmonth_flux<-flux[which(as.POSIXct(flux$Time) >= "2014-01-01" & as.POSIXct(flux$Time) < "2014-07-01"),]

# fluxA<-flux[1:(length(flux[,1]) - 30),]
# fluxB<-tail(flux, 30)

fluxA<-sixmonth_flux
fluxB<-flux[which(as.POSIXct(flux$Time) >= "2014-07-01" & as.POSIXct(flux$Time) < "2014-07-02"),]

## transform the data into a time series; the frequency has a lot of influence on the results

## frequency values can be: 262800 for 2 min periodicity; 12 for monthly periodicity; 365 for yearly periodicity; 8760 for hourly periodicity
freq=30
myts<-ts(fluxA$Flux[!is.na(fluxA$Flux)], frequency = freq)
fit <- stl(myts, s.window=freq)

# plot(fit)
# seasonplot(myts)

acf(myts, plot = F)
acf(myts, type = "covariance", plot = F)$acf

doubledriftmyts<-diff(diff(myts))
# plot(doubledriftmyts)
# plot(HoltWinters(myts))

forecastmyts<-forecast(HoltWinters(myts), h = length(fluxB[,1]))
# plot(predict(HoltWinters(myts), 30))
# plot(HoltWinters(myts))

predictedHW<-predict(HoltWinters(myts), length(fluxB[,1]))

mseHW<-sum(fluxB$Flux - predictedHW)^2/length(fluxB[,1])

plot(fluxB$Flux, col = "red", type = "l", ylim = c(min(predictedHW, fluxB$Flux), max(predictedHW, fluxB$Flux)), main = "Forecasting Accuracy from HW Model")
points(as.numeric(predictedHW), col = "blue", type = "l")
legend("topleft", col = c("red", "blue"), c("real data", "forecasted data"), pch = 19)

## random walk model
predictedRW<-rwf(myts, drift = T, h=length(fluxB[,1]))
plot(rwf(myts, drift = T, h=length(fluxB[,1])), col = "blue")

mseRW<-sum(fluxB$Flux - predictedRW$mean)^2/length(fluxB[,1])

plot(fluxB$Flux, col = "red", type = "l", ylim = c(min(predictedRW$mean, fluxB$Flux), max(predictedRW$mean, fluxB$Flux)), main = "Forecasting Accuracy from RW Model")
points(as.numeric(predictedRW$mean), col = "blue", type = "l")
legend("topleft", col = c("red", "blue"), c("real data", "forecasted data"), pch = 19)

## simple exponential smoothing
predictedSE<-ses(myts, alpha = .2, h=length(fluxB[,1]))

mseSE<-sum(fluxB$Flux - predictedSE$mean)^2/length(fluxB[,1])

plot(fluxB$Flux, col = "red", type = "l", ylim = c(min(predictedSE$mean, fluxB$Flux), max(predictedSE$mean, fluxB$Flux)), main = "Forecasting Accuracy from SE Model")
points(as.numeric(predictedSE$mean), col = "blue", type = "l")
legend("topleft", col = c("red", "blue"), c("real data", "forecasted data"), pch = 19)

## moving average 
movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

predictedMA<-tail(movingAverage(myts, n=1), length(fluxB[,1]))

mseMA<-sum(fluxB$Flux - predictedMA)^2/length(fluxB[,1])

plot(fluxB$Flux, col = "red", type = "l", ylim = c(min(predictedMA, fluxB$Flux), max(predictedMA, fluxB$Flux)), main = "Forecasting Accuracy from MA Model")
points(as.numeric(predictedMA), col = "blue", type = "l")
legend("topleft", col = c("red", "blue"), c("real data", "forecasted data"), pch = 19)

predictedValues<-as.data.frame(cbind(fluxB$Flux, as.numeric(predictedMA), predictedSE$mean, predictedRW$mean, as.numeric(predictedHW)))

write.csv(predictedValues, "PredictedFromBaselines.csv")


## autoregressive processes
ar(myts)
plot(ar(myts))


## power law fit
library(sm)
library(e1071)
n=5000
bins<-binning(flux$Flux[!is.na(flux$Flux)], nbins=n)
plot(log(bins$x), log(bins$x.freq), pch=18, col = "red")
abline(lm(log(bins$x.freq)~log(bins$x)))

power.law.fit(flux$Flux[!is.na(flux$Flux)])
power.law.fit(bins$x)

## Lorenz curves
library(ineq)
plot(Lc(flux$Flux[!is.na(flux$Flux)]))

gini_flux<-ineq(flux$Flux[!is.na(flux$Flux)],type="Gini")


### use the flare classification for the test data:
subsetflux<-sixmonth_flux #(this can be changed with the fluxB when needed to validate with the first day of July)

flareA<-subsetflux[which(subsetflux$Flux <= 10^(-7)),]
flareB<-subsetflux[which(subsetflux$Flux > 10^(-7) & subsetflux$Flux <= 10^(-6)),]
flareC<-subsetflux[which(subsetflux$Flux > 10^(-6) & subsetflux$Flux <= 10^(-5)),]
flareM<-subsetflux[which(subsetflux$Flux > 10^(-5) & subsetflux$Flux <= 10^(-4)),]
flareX<-subsetflux[which(subsetflux$Flux > 10^(-4)),]

flareA$Flare<-c("A")
flareB$Flare<-c("B")
flareC$Flare<-c("C")
flareM$Flare<-c("M")
flareX$Flare<-c("X")

newsubsetflux<-as.data.frame(rbind(flareA, flareB, flareC, flareM, flareX))

### time clasification (based on maximum flux values)
hourlyflux<-aggregate(newsubsetflux[2], by = list(format(as.POSIXct(newsubsetflux$Time), "%m-%d %H")), FUN = max, na.rm=T)
dailyflux<-aggregate(newsubsetflux[2], by = list(as.Date(newsubsetflux$Time)), FUN = max, na.rm=T)
monthlyflux<-aggregate(newsubsetflux[2], by = list(format(as.Date(newsubsetflux$Time), "%Y-%m")), FUN = max, na.rm=T)

### differenced time series
newflux<-newsubsetflux[order(as.POSIXct(newsubsetflux$Time)),]
newflux$diff_flux<-c(0, diff(newflux$Flux))

poz_flux<-newflux[which(newflux$diff_flux > 0),]
neg_flux<-newflux[which(newflux$diff_flux < 0),]


### predictions from NOAA data

setwd("~/NOAA_2014_sixmonths")
noaa_files<-list.files()
noaa_estimates<-lapply(noaa_files, readLines)

idM<-charmatch("Class M",noaa_estimates[1][[1]])
idX<-charmatch("Class X",noaa_estimates[1][[1]])
noaa_M<-as.numeric(strsplit(strsplit(noaa_estimates[1][[1]][idM], " ")[[1]][6], "/")[[1]][1])
noaa_X<-as.numeric(strsplit(strsplit(noaa_estimates[1][[1]][idX], " ")[[1]][6], "/")[[1]][1])
for (i in 1:length(noaa_files)){
  idM<-charmatch("Class M",noaa_estimates[i][[1]])
  idX<-charmatch("Class X",noaa_estimates[i][[1]])
  noaa_M[i]<-as.numeric(strsplit(strsplit(noaa_estimates[i][[1]][idM], " ")[[1]][6], "/")[[1]][1])
  noaa_X[i]<-as.numeric(strsplit(strsplit(noaa_estimates[i][[1]][idX], " ")[[1]][6], "/")[[1]][1])
}

flareA<-dailyflux[which(dailyflux$Flux <= 10^(-7)),]
flareB<-dailyflux[which(dailyflux$Flux > 10^(-7) & dailyflux$Flux <= 10^(-6)),]
flareC<-dailyflux[which(dailyflux$Flux > 10^(-6) & dailyflux$Flux <= 10^(-5)),]
flareM<-dailyflux[which(dailyflux$Flux > 10^(-5) & dailyflux$Flux <= 10^(-4)),]
flareX<-dailyflux[which(dailyflux$Flux > 10^(-4)),]

flareA$Flare<-c("A")
flareB$Flare<-c("B")
flareC$Flare<-c("C")
flareM$Flare<-c("M")
flareX$Flare<-c("X")

newdailyflux<-as.data.frame(rbind(flareA, flareB, flareC, flareM, flareX))
newdailyflux<-newdailyflux[sort(newdailyflux, by = newdailyflux$Time),]

cbind(newdailyflux[-1,], noaa_M, noaa_X)

write.csv(cbind(newdailyflux[-1,], noaa_M, noaa_X), "NOAA 2014 predictions.csv")
