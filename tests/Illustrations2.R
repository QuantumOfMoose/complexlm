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
err <- mvrnorm(n = nn, mu = c(0,0), Sigma = sigma / 6) # Generate realizations of a symmetrical two dimensional normal random variable.
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
plot(fitone.clean.ols)
