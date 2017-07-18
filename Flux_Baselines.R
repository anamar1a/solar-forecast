library(forecast)

flux<-read.csv("~/Flux_2010_2015.csv")

fluxA<-flux[1:(length(flux[,1]) - 30),]
fluxB<-tail(flux, 30)

## transform the data into a time series; the frequency has a lot of influence on the results

## frequency values can be: 262800 for 2 min periodicity; 12 for monthly periodicity; 365 for yearly periodicity; 8760 for hourly periodicity
freq=262800
myts<-ts(fluxA$Flux[!is.na(fluxA$Flux)], start=c(2010, 1, 1), end = c(2015, 12, 30), frequency = freq)

fit <- stl(myts, s.window="period")

# plot(fit)
# seasonplot(myts)

acf(myts, plot = F)
acf(myts, type = "covariance", plot = F)$acf

doubledriftmyts<-diff(diff(myts))
# plot(doubledriftmyts)
# plot(HoltWinters(myts))

forecastmyts<-forecast(HoltWinters(myts), h = 30)
# plot(predict(HoltWinters(myts), 30))
# plot(HoltWinters(myts))

predictedHW<-predict(HoltWinters(myts), 30)

mseHW<-sum(fluxB$Flux - predictedHW)^2/30

plot(fluxB$Flux, col = "red", type = "l", ylim = c(min(predictedHW, fluxB$Flux), max(predictedHW, fluxB$Flux)), main = "Forecasting Accuracy from HW Model")
points(as.numeric(predictedHW), col = "blue", type = "l")
legend("topleft", col = c("red", "blue"), c("real data", "forecasted data"), pch = 19)

## random walk model
predictedRW<-rwf(myts, drift = T, h=30)
plot(rwf(myts, drift = T, h=30), col = "blue")

mseRW<-sum(fluxB$Flux - predictedRW$mean)^2/30

plot(fluxB$Flux, col = "red", type = "l", ylim = c(min(predictedRW$mean, fluxB$Flux), max(predictedRW$mean, fluxB$Flux)), main = "Forecasting Accuracy from RW Model")
points(as.numeric(predictedRW$mean), col = "blue", type = "l")
legend("topleft", col = c("red", "blue"), c("real data", "forecasted data"), pch = 19)

## simple exponential smoothing
predictedSE<-ses(myts, alpha = .2, h=30)

mseSE<-sum(fluxB$Flux - predictedSE$mean)^2/30

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

predictedMA<-tail(movingAverage(myts, n=1), 30)

mseMA<-sum(fluxB$Flux - predictedMA)^2/30

plot(fluxB$Flux, col = "red", type = "l", ylim = c(min(predictedMA, fluxB$Flux), max(predictedMA, fluxB$Flux)), main = "Forecasting Accuracy from MA Model")
points(as.numeric(predictedMA), col = "blue", type = "l")
legend("topleft", col = c("red", "blue"), c("real data", "forecasted data"), pch = 19)

predictedValues<-as.data.frame(c(fluxB$Flux, predictedMA, predictedSE$mean, predictedRW$mean, predictedHW))

write.csv(predictedValues, "PredictedFromBaselines.csv")
