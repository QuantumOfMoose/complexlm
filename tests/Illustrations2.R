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
m <- 5
boundary <- m %/% 2
x <- expand.grid(c(-boundary:boundary), c(-boundary:boundary)) # Create a 5 by 5 grid of regular points.
x <- complex(real = x[[1]], imaginary = x[[2]]) # Make complex numbers
nn <- length(x) # Number of elements in x

### Plot our predictors as points in the complex plane. They are a simple square grid.
ggplot(data = data.frame(x), aes(x = Re(x), y = Im(x))) +
  geom_point(size = 2) +
  labs(x = "Real", y = "Imaginary", title = "Complex Predictor Variable")

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
#print(fitone.clean.ols)

### Now that we have a fit, we can try out the fit diagnostics and summary methods included in 'complexlm'.
### First we will check the classes that this fit object falls into.
class(fitone.clean.ols) # 'zlm' and 'lm'
### Many S3 methods written for objects of class 'lm' also work for 'zlm' objects.
coefficients(fitone.clean.ols)
formula(fitone.clean.ols)
kappa(fitone.clean.ols)
### Others, do not.
#summary.lm(fitone.clean.ols)
### Next summarize and print the summary.
summary.fitone.clean.ols <- summary(fitone.clean.ols)
print(summary.fitone.clean.ols)
### Then generate an analysis of variance (ANOVA)
anova.zlm(fitone.clean.ols)
### We can also calculate the standardized residuals, hat matrix or influence scores, and the Cook's distances.
head(rstandard(fitone.clean.ols))
exonedf$hats.clean.ols <- zhatvalues(fitone.clean.ols)
#head(exonedf$hats.clean.ols)
head(cooks.distance(fitone.clean.ols))
### These three things are used to generate the traditional diagnostic plots for linear models. We can draw them like so,
plot(fitone.clean.ols, which = c(1,3,4,5,6))
### The following lines generate the same plots, but save them to individual .pdf files.
# setwd() # Change the working directory, if desired.
# pdf("resvfit.pdf") # Draw Residuals vs. Fitted plot to this file.
# plot(fitone.clean.ols, which = 1) # Plot Residuals vs. Fitted
# dev.off() # Close the file / graphics device.
# pdf("scalloc.pdf") # Draw the Scale-Location plot in this file.
# plot(fitone.clean.ols, which = 3) # Draw the Scale-Location plot. Note that which = 2 is not available for complex input.
# dev.off()
# pdf("cookdist.pdf") # Draw the Cook's Distance plot here.
# plot(fitone.clean.ols, which = 4) # Draw the Cook's Distance plot.
# dev.off()
# pdf("reslev.pdf") # Put the Residuals vs. Leverage plot in this file.
# plot(fitone.clean.ols, which = 5)
# dev.off()
# pdf("cooklev.pdf") # Draw the Cook's Distance vs. Scaled Leverage plot in this file.
# plot(fitone.clean.ols, which = 6) # Plot the Cook's Distance vs. Scaled Leverage.
# dev.off()

### Now we'll make a nicer plot of the 'measured' and fitted response values. This will serve the role of the classic x-y plot with a trend line.
### First add the fitted response to the dataframe.
exonedf$y.fit <- fitted(fitone.clean.ols)
exonedf$r <- residuals(fitone.clean.ols)
### Convert to long form data frame
melt.exonedf <- melt(exonedf, id = c("x", "hats.clean.ols"))
labels <- c('y' = 'measured response', 'y.clean' = 'relationship without noise', 'y.fit' = 'fitted response', 'r' = 'residuals') # Nice labels for plots.
### Define some consistent colors to make the plots more aesthetically pleasing.
clreen <- 'darkgreen'
brout <- 'brown'
respring <- 'springgreen2'
olsalmon <- 'salmon'
hubsky <- 'deepskyblue2'
hampurp <- 'mediumorchid3'
bislate <- 'slateblue2'

### Then make the plot.
ggplot(melt.exonedf[grepl('y', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y'=19, 'y.clean'=0, 'y.fit'=18), labels = labels) +
  scale_color_manual(values = c('y'= respring, 'y.clean'= clreen, 'y.fit'= olsalmon), labels = labels) +
  geom_line(data = melt.exonedf[melt.exonedf$variable == "y" | melt.exonedf$variable == "y.fit",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  scale_linetype_manual("Residuals", values = 'solid', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated and OLS Fit Values Without Outliers") 
### In this plot the measured (simulated) responses are represented by cyan circles, the noiseless responses are shown as empty green squares, 
### and the fitted responses are red diamonds. Fitted values are connected to their corresponding measured points by solid coral lines, which represent the residuals.
### All red diamonds lie mostly within green squares, which indicates that the complex ordinary least squares fit successfully extracted the relationship between the
### the predictor and response variables in this case.

### We can also plot the residuals in the complex plane.
ssrlab <- function(name){ # A function that calculates the sum of the squared elements of a vector-like object named 'name', then formats it into a string.
  paste("Sum of Squared Residuals:\n", round(sum(Mod(name)^2), digits = 4)) # A string to write on the residual plot that gives the sum of squared residuals.
} 
ggplot(melt.exonedf[grepl('r', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), alpha = hats.clean.ols)) +
  geom_segment(aes(xend = 0, yend = 0)) +
  labs(y = "Imaginary", x = "Real", title = "Residuals of OLS Fit Without Outliers", alpha = "Influence Score") +
  coord_fixed() +
  #geom_label(x = 20, y = -10, label = ssrlab(value)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

### That worked quite well, but what if we add some outliers?
### Say that there is a probability theta < 0.5 that each element of the data was not transformed by beta, but by betap.
### We will arbitrarily set,
set.seed(4242)
theta <- 0.2
isoutl <- as.logical(rbinom(nn, 1, theta)) # Logical vector that decides which elements of y should be replaced with an outlier.
sum(isoutl) / nn # Check the actual proportion of outliers.
### Now we will arbitrarily choose a value for betap.
betap <- c(complex(real = 42, imaginary = 31), complex(real = .5, imaginary = .82))
err.outl <- mvrnorm(n = nn, mu = c(0,0), Sigma = sigma) # Generate realizations of a symmetrical two dimensional normal random variable, to represent the uncertainty / error in the outliers.
err.outl <- complex(real = err.outl[,1], imaginary = err.outl[,2])
### Now add another column to exonedf populated with x * beta or x * betap, depending on the value of isoutl.
exonedf <- exonedf %>% mutate(y.outl = if_else(isoutl, betap[1] + x * betap[2] + err.outl, y))
### And add a column that represents the phenomenon responsible for the outliers. i.e. the predictor variables transformed by betap, without noise.
exonedf$clout <- exonedf$x * betap[2] + betap[1]

### Let's see how well complex least squares does with these outliers.
fitone.outl.ols <- lm(y.outl ~ x, exonedf)
summary(fitone.outl.ols)
plot(fitone.outl.ols, which = c(1,3,4,5,6))
### Repeat the same procedure as before to plot the response data.
exonedf$y.fit.outl.ols <- fitted(fitone.outl.ols)
exonedf$r.outl.ols <- residuals(fitone.outl.ols)
exonedf$hats.outl.ols <- zhatvalues(fitone.outl.ols)
labels <- c(labels, 'y.outl' = 'measured response with outliers', 'y.fit.outl.ols' = 'OLS fitted response', 'r.outl.ols' = 'OLS residuals', 'clout' = 'outlier relationship')
melt.exonedf <- melt(exonedf, id = c('x', 'hats.outl.ols'))

ggplot(melt.exonedf[grepl('y|outl|y.clean|clout', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y.outl'=19, 'y.clean'=0, 'y.fit.outl.ols'=18, 'clout' = 5), labels = labels) +
  scale_color_manual(values = c('y.outl' = respring, 'y.clean' = clreen, 'y.fit.outl.ols' = olsalmon, 'clout' = brout), labels = labels) +
  geom_line(data = melt.exonedf[melt.exonedf$variable == "y.outl" | melt.exonedf$variable == "y.fit.outl.ols",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  #scale_linetype_manual("Residuals", values = 'dashed', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated and OLS Fit Values With Outliers") +
  coord_fixed()
### And plot the residuals.
ggplot(melt.exonedf[grepl('r.outl', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), size = hats.outl.ols)) +
  geom_point(color = olsalmon) +
  geom_label(aes(x = 20, y = -10, label = ssrlab(value)), size = 3.5) +
  labs(y = "Imaginary", x = "Real", title = "Residuals of OLS Fit With Outliers", size = "Influence Score") +
  coord_fixed() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

### That didn't work very well. While the sum of the squared residuals was minimized, the fit coefficients are very far from the actual relationship between predictor and response.
### Now let's try the complex robust M-estimator. We'll try the Huber loss objective function, which is the default.
### But first we need to establish some parameters for the IWLS algorithm.
accept = 1e-20 # Shared acceptance criterion
iterations = 70 # Shared max iterations
fitone.outl.hub <- rlm(y.outl ~ x, exonedf, maxit = iterations, acc = accept)
summary(fitone.outl.hub)
plot(fitone.outl.hub)
### Repeat the same procedure as before to plot the response data.
exonedf$y.fit.outl.hub <- fitted(fitone.outl.hub)
exonedf$r.outl.hub <- residuals(fitone.outl.hub)
exonedf$hats.outl.hub <- zhatvalues(fitone.outl.hub)
labels <- c(labels, 'y.fit.outl.hub' = 'Huber rlm fitted response', 'r.outl.hub' = 'Huber rlm residuals')
melt.exonedf <- melt(exonedf, id = c('x', 'hats.outl.hub'))

ggplot(melt.exonedf[grepl('y|outl|y.clean|clout', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y.outl'=19, 'y.clean'=0, 'y.fit.outl.ols'=18, 'y.fit.outl.hub' = 13, 'clout' = 5), labels = labels) +
  scale_color_manual(values = c('y.outl'= respring, 'y.clean'= clreen, 'y.fit.outl.ols'= olsalmon, 'y.fit.outl.hub' = hubsky, 'clout' = brout), labels = labels) +
  #geom_line(data = melt.exonedf[melt.exonedf$variable == "y.outl" | melt.exonedf$variable == "y.fit.outl.ols",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  #scale_linetype_manual("Residuals", values = 'dashed', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated, OLS Fit, and Huber Robust Fit Values With Outliers") +
  coord_fixed()
### And plot the residuals.
ggplot(melt.exonedf[grepl('r.outl', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), size = )) +
  #geom_point(aes(color = variable)) +
  geom_segment(aes(x = 0, y = 0, xend = Re(value), yend = Im(value), color = variable)) +
  labs(y = "Imaginary", x = "Real", title = "Residuals of OLS Fit With Outliers") +
  scale_color_manual("Residuals", values = c('r.outl.ols' = olsalmon, 'r.outl.hub' = hubsky)) +
  coord_fixed() +
  geom_label(x = 20, y = -10, label = ssrlab(value)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
### It's not perfect, but significantly better than what the ordinary least squares fit found. While the overall sum of squared residuals for the Huber M-estimation 
### is larger than that of the ols fit, the residuals due to the non-outliers are much smaller.

### Now let's try some re-descending M-estimators. First with the Hampel objective function.
fitone.outl.ham <- rlm(y.outl ~ x, exonedf, maxit = iterations, acc = accept/100, psi = psi.hampel, a = 1.345)
exonedf$y.fit.outl.ham <- fitted(fitone.outl.ham)
exonedf$r.outl.ham <- residuals(fitone.outl.ham)
exonedf$hats.outl.ham <- zhatvalues(fitone.outl.ham)
### Then with Tukey's bisquare
fitone.outl.bis <- rlm(y.outl ~ x, exonedf, maxit = iterations + 90, acc = 10000 * accept, psi = psi.bisquare)
exonedf$y.fit.outl.bis <- fitted(fitone.outl.bis)
exonedf$r.outl.bis <- residuals(fitone.outl.bis)
exonedf$hats.outl.bis <- zhatvalues(fitone.outl.bis)
### Update the labels.
labels <- c(labels, 'y.fit.outl.ham' = 'Hampel rlm fitted response', 'r.outl.ham' = 'Hampel rlm residuals', 'y.fit.outl.bis' = 'Tukey Bisquare rlm fitted response', 'r.outl.bis' = 'Tukey Bisquare rlm residuals')
melt.exonedf <- melt(exonedf, id = 'x')

ggplot(melt.exonedf[grepl('y|outl|y.clean|clout', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y.outl'=19, 'y.clean'=0, 'y.fit.outl.ols'=18, 'y.fit.outl.hub' = 13, 'y.fit.outl.ham' = 2, 'y.fit.outl.bis' = 8, 'clout' = 5), labels = labels) +
  scale_color_manual(values = c('y.outl'= respring, 'y.clean'= clreen, 'y.fit.outl.ols'= olsalmon, 'y.fit.outl.hub' = hubsky, 'y.fit.outl.ham' = hampurp, 'y.fit.outl.bis' = bislate, 'clout' = brout), labels = labels) +
  #geom_line(data = melt.exonedf[melt.exonedf$variable == "y.outl" | melt.exonedf$variable == "y.fit.outl.ols",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  #scale_linetype_manual("Residuals", values = 'dashed', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated, OLS Fit, and Huber Robust Fit Values With Outliers") +
  coord_fixed()
### And plot the residuals.
rlablframe <- melt.exonedf[grepl('r.outl', melt.exonedf$variable),] %>% group_by(variable) %>% summarize(label = ssrlab(value))
ggplot(melt.exonedf[grepl('r.outl', melt.exonedf$variable),], aes(x = Re(value), y = Im(value))) +
  #geom_point(aes(color = variable)) +
  geom_segment(aes(x = 0, y = 0, xend = Re(value), yend = Im(value), color = variable)) +
  facet_wrap(vars(variable), nrow = 3, ncol = 2, labeller = as_labeller(labels)) +
  labs(y = "Imaginary", x = "Real", title = "Residuals of OLS Fit With Outliers") +
  scale_color_manual("Residuals", values = c('r.outl.ols' = olsalmon, 'r.outl.hub' = hubsky, 'r.outl.ham' = hampurp, 'r.outl.bis' = bislate), labels = labels) +
  coord_fixed() +
  #geom_label(data = x = 20, y = -10, aes(label = ssrlab(value))) + # Unfortunately, this doesn't work with the facets. :(
  geom_label(data = rlablframe, x = 27, y = -10, aes(label = label)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

### Notice that the rlm() calls that used redescending influence functions achieved estimates of beta closer to the true value. 
### This is because the algorithm downweighted them into oblivion, so that they have no influence upon the final estimate. Thus they 
### are exceedingly powerful tools of discrimination, but care must be taken when using them as they can end up ignoring useful data as well.

### Finally, we will draw the estimated coefficients from each of these fits, along with the true relationship and the outlier relationship, on the complex plane.
fitnames <- paste('fitone.outl', c('ols', 'hub', 'ham', 'bis'), sep = '.')
coefgetter <- function(name) { # A function that collects the coefficients and their variances and pseudo-variances of a fit with name name, and puts them in a named vector.
  fit <- as.name(name)
  fitsum <- summary(eval(fit))
  xstuf <- fitsum$coefficients[1,][0:3] # We don't need the t value.
  interstuf <- fitsum$coefficients[2,][0:3]
  names(xstuf) <- make.names(paste("x", names(xstuf))) # Add "x" to front of element names, to differentiate them from the values corresponding to the intercept when we combine the vectors.
  names(interstuf) <- make.names(paste('interc', names(interstuf))) # Same here, add prefix to names to disambiguate them.
  return(c(xstuf, interstuf)) # combine and return the vectors of values and errors.
}
coefframe <- do.call(rbind, lapply(fitnames, coefgetter))
coefframe <- as.data.frame(coefframe)
coefframe <- rbind(coefframe, c(beta[2], NA_complex_, NA_complex_, beta[1], NA_complex_, NA_complex_)) # Add the true coefficient values.
coefframe <- rbind(coefframe, c(betap[2], NA_complex_, NA_complex_, betap[1], NA_complex_, NA_complex_)) # Add the true outlier coefficient values.
coefframe$variable <- c(fitnames, "true", "true.outlll")

### Plot the scale rotation (slope) estimates
ggplot(coefframe, aes(color = variable)) +
  geom_point(aes(x = Re(x.Estimate), y = Im(x.Estimate))) +
  geom_ellipse(aes(x0 = Re(x.Estimate), y0 = Im(x.Estimate), a = Mod(x.Pseudo.Std..Error), b = Re(x.Std..Error)^2 / Mod(x.Pseudo.Std..Error), angle = Arg(x.Pseudo.Std..Error))) +
  scale_color_manual(values = c('true'= clreen, 'fitone.outl.ols'= olsalmon, 'fitone.outl.hub' = hubsky, 'fitone.outl.ham' = hampurp, 'fitone.outl.bis' = bislate, 'true.outlll' = brout), labels = c('true'= "True Value", 'fitone.outl.ols'= "Least-Squares", 'fitone.outl.hub' = "M-Estimation w/ Huber Objective", 'fitone.outl.ham' = 'M-Estimation w/ Hampel Objective', 'fitone.outl.bis' = 'M-Estimation w/ Bisquare Objective', 'true.outlll' = 'Outlier Value')) +
  labs(title = "Scale Rotation Coefficient Estimates", color = 'Regression', x = 'Real', y = 'Imaginary')
### This plot shows that the M-estimator regressions using Tukey's bisquare objective function and Hampel's objective function produce estimates of the scale rotation coefficient that are very close to the true value,
### despite the presence of outliers pulling them toward the marroon dot. The bisquare regression results in a markedly lower variance though. The Huber M-estimator, on the other hand, did not fair much better than the least-squares fit.

### Now we plot the translation (intercept) coefficient estimates.
ggplot(coefframe, aes(color = variable)) +
  geom_point(aes(x = Re(interc.Estimate), y = Im(interc.Estimate))) +
  geom_ellipse(aes(x0 = Re(interc.Estimate), y0 = Im(interc.Estimate), a = Mod(interc.Pseudo.Std..Error), b = Re(interc.Std..Error)^2 / Mod(interc.Pseudo.Std..Error), angle = Arg(interc.Pseudo.Std..Error))) +
  scale_color_manual(values = c('true'= clreen, 'fitone.outl.ols'= olsalmon, 'fitone.outl.hub' = hubsky, 'fitone.outl.ham' = hampurp, 'fitone.outl.bis' = bislate, 'true.outlll' = brout), labels = c('true'= "True Value", 'fitone.outl.ols'= "Least-Squares", 'fitone.outl.hub' = "M-Estimation w/ Huber Objective", 'fitone.outl.ham' = 'M-Estimation w/ Hampel Objective', 'fitone.outl.bis' = 'M-Estimation w/ Bisquare Objective', 'true.outlll' = 'Outlier Value')) +
  labs(title = "Scale Rotation Coefficient Estimates", color = 'Regression', x = 'Real', y = 'Imaginary')
### In estimating the translation, the robust M-estimators all outperformed the least-squares estiimate by a considerable degree. The true value is well within the first standard deviation ellipse of all three of them.

###
### This concludes the first example of ordinary and robust fitting using 'complexlm'.
### Next we will take on some actual, measured data.
###

###
### Example 2
### AC Hall Effect Data
###

### In this example we will analyze the results of a room temperature Hall effect measurement of a copper sample,
### with the goal of extracting its carrier mobility and concentration.
### Since the complex voltage output of an AC Hall measurement is proportional to the current input, the core of this task is complex linear fitting.

### First we must load the Hall effect data contained within CuHallData.rda, which is included with 'complexlm'.
data(CuHallData)
summary(CuHallData)
### Notice that the current columns are all real numbers, but 'complexlm' requires that both the response and predictors be complex. Luckily, the reals are a subset
### of the complex numbers, so we can convert our currents to complex representation by adding 0*i to them.
CuHallData$zCurrent <- complex(real = CuHallData$Curent.In.meas.A., imaginary = 0) # Use the input current.
#CuHallData$zCurrent <- (CuHallData$Curent.In.meas.A. - CuHallData$Curent.Out.meas.A.)/2 # We take the average of the measured input and output in an attempt to neutralize the effect of current leaks

### First we'll have a look at all the output voltage data in the complex plane. Since the input currents are a arranged in a simple line along the real axis
### We would expect the output voltages to be arranged in rough lines as well, though not just on the real axis.
### Since the data in CuHallData is from several different current-voltage measurements, taken at different magnetic field frequencies and one of two possible contact arrangements,
### we will group the data by the corresponding columns. We'll use 'dplyr' and tibble for convenience.
Halltibble <- CuHallData %>% group_by(Contact.Arangement, Magnet.Field.Frequency.Hz.) # First group by the contact arrangement and frequency.
### And now we can draw our plot.
ggplot(Halltibble, aes(x = Re(OutputV), y = Im(OutputV), shape = Contact.Arangement, color = as.factor(Magnet.Field.Frequency.Hz.))) +
  geom_point() 
### Well that looks like a mess.
### Let's see how it looks when we plot real and imaginary output voltage vs input current.
ggplot(Halltibble, aes(x = Re(zCurrent), y = Re(OutputV))) +
  geom_point(aes(color = 'real')) +
  geom_point(aes(y = Im(OutputV), color = 'imaginary')) +
  scale_color_manual(values = c('real' = "forestgreen", 'imaginary' = "coral")) +
  facet_grid(rows = vars(Contact.Arangement), cols = vars(Magnet.Field.Frequency.Hz.)) +
  labs(x = "Current", y = "Output Voltage", color = "Voltage Component")
### Here we have organized the data into a grid based on the magnetic field frequency and contact orientation with which they were taken.
### The label at the top of each column gives the frequency in Hertz, and the label at the right of the rows give the contact orientation.
### These plots are certainly neater, but they also show that the slopes we need to extract are exceedingly small. 
### And, especially for the measurement taken with f_B = 0.57 Hz and contact orientation "D", filled with outliers.
### This is exactly what we should expect for a metal.
  
### Now we can apply the linear fitting routines of 'complexlm'. We should note that this dataframe contains information from multiple different current-voltage measurements. 
### We can use summarize from 'dplyr' to efficiently find the linear relationship between current and voltage for each of these. 
fitter <- function(data, formula, ...) { # A function that will perform the fits then collect the coefficients and uncertainties into a tibble.
  thisfit <- rlm(formula, data, ...)
  thissumary <- summary(thisfit)
  as_tibble(thissumary$coefficients[1:2,1:3], rownames = "coefficient")
}
Hallcoeftibble <- Halltibble %>% summarize(fitter(OutputV ~ zCurrent, data = .))

fitterlm <- function(data, formula, ...) { # A function that will perform ols fits then collect the coefficients and uncertainties into a tibble.
  thisfit <- lm(formula, data, ...)
  thissumary <- summary(thisfit)
  as_tibble(thissumary$coefficients[1:2,1:3], rownames = "coefficient")
}
lmHallcoeftibble <- Halltibble %>% summarize(fitterlm(OutputV ~ zCurrent, data = .))

### Add a column for the average magnetic field amplitudes for I-V sweep, and columns for the resitivity and thickness.
otherdf <- Halltibble %>% summarize(Magnet.Field.T = mean(Magnetic.Field.T), Resistivity.Ohm.cm = mean(Resistivity.Ohm.cm), thickness.cm. = mean(thickness.cm.))
Hallcoeftibble <- inner_join(x = Hallcoeftibble, y = otherdf, by = c("Contact.Arangement", "Magnet.Field.Frequency.Hz."))
lmHallcoeftibble <- inner_join(x = lmHallcoeftibble, y = otherdf, by = c("Contact.Arangement", "Magnet.Field.Frequency.Hz."))

### The mobility in cm^2/(V s) can be calculated (absent phase shift correction) from the slope (zCurrent) coefficients [which are the Hall resistances] of the fits like so,
### mobility mu = Mod([slope coefficient]) * [thickness] / [Magnet.Field.T] * [10^-4]
### First filter the coefficient tibbles to remove the (intercept) rows, then add a column of mobility values.

mobilitytibble <- Hallcoeftibble %>% filter(coefficient == "zCurrent") %>% mutate(mobility = (Mod(Value) * thickness.cm.) / (Magnet.Field.T * 10^-4 * Resistivity.Ohm.cm))
lmmobilitytibble <- lmHallcoeftibble %>% filter(coefficient == "zCurrent") %>% mutate(mobility = (Mod(Estimate) * thickness.cm.) / (Magnet.Field.T * 10^-4 * Resistivity.Ohm.cm))
summary(mobilitytibble$mobility)
### Which is pretty close (on the same order of magnitude) to that found in the literature. Yay!
### With this information, along with a value for the elementary charge, we can also calculate the carrier concentration. 
elemcharge <- 1.602176e-19 # in units of Coulombs
mobilitytibble$density <- 1 / (elemcharge * mobilitytibble$Resistivity.Ohm.cm * mobilitytibble$mobility) # This has units of 1/cm^3
### Now let's compare the 
