
remove(list = ls())

install.packages("readxl")
library(readxl)

# Load the data from excel file, weekly spread 

data <- read_excel("Spread-weekly.xlsx")

# Remove first row, as it contains NA
data <- data[-1,]

# The data is inverse, we need to change sign
data$`IT10DE10=TWEB (BID)` <- as.numeric(data$`IT10DE10=TWEB (BID)`) * -1

# Inspect data
head(data)
str(data)

# parse the data

library(zoo)

data$Date <- as.Date(data$Date)
spread <- zoo(data$`IT10DE10=TWEB (BID)`, order.by = data$Date)

# Plot the data

plot(spread, main="BTP-BUND Spread (Weekly", ylab="Spread (bps)", xlab= "Date")



# Forecast #1, using Naive Forecast (Mean)

mean <- mean(data$`IT10DE10=TWEB (BID)`, na.rm=TRUE)

# Plot spread values with the mean value

plot(spread, main="BTP-BUND Spread (Weekly", ylab="Spread (bps)", xlab= "Date")
abline(h = mean, col="red", lwd=2, lty=2)
text(x = start(spread), y= mean,
     labels=paste("Mean =", round(mean, 2)),
     col="red", adj= c(0, -0.5))





