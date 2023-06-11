###
### This is a companion to the paper which introduces the complexlm package. 
### It demonstrates the functions in said package by using them to explore and 
### analyze two sets of complex data. It produces all the figures featured in the
### Illustrations section of the paper.
### 
### William Ryan
###
### 6 June, 2023
###

### First we must load up some libraries that will help us with generating, manipulating, and visualizing the data.
library(stats)
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggforce)
### Now load complexlm. It is important to do so after loading stats and MASS, as complexlm masks a few functions from those packages.
library(complexlm)

### Now begin the first example problem. Consider n = m^2 points arranged in a square grid on the complex plane, centered on the origin. These are our predictor variables, x.
m <- 7
boundary <- m %/% 2
x <- expand.grid(c(-boundary:boundary), c(-boundary:boundary)) # Create a 7 by 7 grid of regular points.
x <- complex(real = x[[1]], imaginary = x[[2]]) # Make complex numbers
nn <- length(x) # Number of elements in x

### Next establish the true transformation vector, call it beta. This is a length two complex vector, the 2nd element of which describes the scaling and rotation that x will undergo, the first describes the translation.
beta <- c(complex(real = -3.5, imaginary = -4), complex(real = 7, imaginary = -4)) # Arbitrarily chosen.

### Now we can multiply beta[[2]] by x, then add beta[[1]], in order to generate a clean, transformed response vector, y.clean.
y.clean <- beta[[2]] * x + rep(beta[[1]], nn)
### We'd never encounter a clean, noiseless signal in the wild though, so to make this a good example we need to add some noise. 
set.seed(4242)
err <- mvrnorm(n = nn, mu = c(0,0), Sigma = sigma / 2) # Generate realizations of a symmetrical two dimensional normal random variable.
err <- complex(real = err[,1], imaginary = err[,2]) # Convert them to complex numbers.
y <- y.clean + err # Add the noise to the clean transformed vector to get a more realistic vector of responses, called y.

### Next we will collect our example data into a dataframe
exonedf <- data.frame(x = x, y.clean = y.clean, err = err, y = y)
### And now we are ready to try out the complex least squares fitting of 'complexlm'.
fitone.clean.ols <- lm(y ~ x, exonedf)

### Now that we have a fit, we can try out the fit diagnostics and summary methods included in 'complexlm'.
### First we will check the classes that this fit object falls into.
class(fitone.clean.ols) # 'zlm' and 'lm'
### Next summarize and print te summary.
summary.fitone.clean.ols <- summary(fitone.clean.ols)
print(summary.fitone.clean.ols)
### Then generate an analysis of variance (ANOVA)
anova.lm(fitone.clean.ols)
### We can also calculate the standardized residuals, hat matrix or influence scores, and the Cook's distances.
rstandard(fitone.clean.ols)
zhatvalues(fitone.clean.ols)
cooks.distance(fitone.clean.ols)
### These three things are used to generate the traditional diagnostic plots for linear models. We can draw them like so,
plot(fitone.clean.ols, which = c(1,3,4,5,6))

### Now we'll make a nicer plot of the 'measured' and fitted response values. This will serve the role of the classic x-y plot with a trend line.
### First add the fitted response to the dataframe.
exonedf$y.fit <- fitted(fitone.clean.ols)
exonedf$r <- residuals(fitone.clean.ols)
### Convert to long form data frame
melt.exonedf <- melt(exonedf, id = "x")
labels <- c('y' = 'measured response', 'y.clean' = 'relationship without noise', 'y.fit' = 'fitted response', 'r' = 'residuals')
### Then make the plot.
ggplot(melt.exonedf[grepl('y', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y'=19, 'y.clean'=0, 'y.fit'=18), labels = labels) +
  scale_color_manual(values = c('y'="cyan3", 'y.clean'='forestgreen', 'y.fit'='red'), labels = labels) +
  geom_line(data = melt.exonedf[melt.exonedf$variable == "y" | melt.exonedf$variable == "y.fit",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  scale_linetype_manual("Residuals", values = 'solid', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated and OLS Fit Values Without Outliers") 
### In this plot the measured (simulated) responses are represented by cyan circles, the noiseless responses are shown as empty green squares, 
### and the fitted responses are red diamonds. Fitted values are connected to thier corresponding measured points by solid coral lines, which represent the residuals.
### All red diamonds lie mostly within green squares, which indicates that the complex ordinary least squares fit successfully extracted the relationship between the
### the predictor and response variables in this case.

### We can also plot the residuals in the complex plane.
ggplot(melt.exonedf[grepl('r', melt.exonedf$variable),], aes(x = Re(value), y = Im(value))) +
  geom_point(aes(xend = 0, yend = 0)) +
labs(y = "Imaginary", x = "Real", title = "Residuals of OLS Fit Without Outliers") +
  coord_fixed()

### That worked quite well, but what if we add some outliiers?
### Say that there is a probability theta < 0.5 that each element of the data was not transformed by beta, but by betap.
### We will arbitrarily set,
set.seed(4242)
theta <- 0.3
isoutl <- as.logical(rbinom(nn, 1, theta)) # Logical vector that decides which elements of y should be replaced with an outlier.
sum(isoutl) / nn # Check the actual proportion of outliers.
### Now we will arbitrarily choose a value for betap.
betap <- c(complex(real = 27, imaginary = 29), complex(real = -5, imaginary = 2))
err.outl <- mvrnorm(n = nn, mu = c(0,0), Sigma = sigma) # Generate realizations of a symmetrical two dimensional normal random variable, to represent the uncertainty / error in the outliers.
err.outl <- complex(real = err.outl[,1], imaginary = err.outl[,2])
### Now add another column to exonedf populated with x * beta or x * betap, depending on the value of isoutl.
exonedf <- exonedf %>% mutate(y.outl = if_else(isoutl, betap[1] + x * betap[2] + err.outl, y))

### Let's see how well complex least squares does with these outliers.
fitone.outl.ols <- lm(y.outl ~ x, exonedf)
summary(fitone.outl.ols)
plot(fitone.outl.ols)
### Repeat the same proceedure as before to plot the response data.
exonedf$y.fit.outl.ols <- fitted(fitone.clean.ols)
exonedf$r.outl.ols <- residuals(fitone.clean.ols)
labels <- c(labels, 'y.outl' = 'measured reesponse with outliers', 'y.fit.outl.ols' = 'OLS fitted response', 'r.outl.ols' = 'OLS residuals')
melt.exonedf <- melt(exonedf, id = 'x')

ggplot(melt.exonedf[grepl('y|outl|y.clean', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y.outl'=19, 'y.clean'=0, 'y.fit.outl.ols'=18), labels = labels) +
  scale_color_manual(values = c('y.outl'="cyan3", 'y.clean'='forestgreen', 'y.fit.outl'='red'), labels = labels) +
  geom_line(data = melt.exonedf[melt.exonedf$variable == "y.outl" | melt.exonedf$variable == "y.fit.outl.ols",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  scale_linetype_manual("Residuals", values = 'solid', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated and OLS Fit Values With Outliers") 
