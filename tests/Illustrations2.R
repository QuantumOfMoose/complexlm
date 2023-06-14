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
### and the fitted responses are red diamonds. Fitted values are connected to thier corresponding measured points by solid coral lines, which represent the residuals.
### All red diamonds lie mostly within green squares, which indicates that the complex ordinary least squares fit successfully extracted the relationship between the
### the predictor and response variables in this case.

### We can also plot the residuals in the complex plane.
ssrlab <- function(name){ # A function that calculates the sum of the squared elements of a vector-like object named 'name', then formats it into a string.
  paste("Sum of Squared Residuals:\n", round(sum(Mod(name)^2), digits = 4)) # A string too write on the residual plot that gives the sum of squared residuals.
} 
ggplot(melt.exonedf[grepl('r', melt.exonedf$variable),], aes(x = Re(value), y = Im(value))) +
  geom_point(aes(xend = 0, yend = 0)) +
  labs(y = "Imaginary", x = "Real", title = "Residuals of OLS Fit Without Outliers") +
  coord_fixed() +
  geom_label(x = 20, y = -10, label = ssrlab(value)) +
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
plot(fitone.outl.ols)
### Repeat the same proceedure as before to plot the response data.
exonedf$y.fit.outl.ols <- fitted(fitone.outl.ols)
exonedf$r.outl.ols <- residuals(fitone.outl.ols)
labels <- c(labels, 'y.outl' = 'measured reesponse with outliers', 'y.fit.outl.ols' = 'OLS fitted response', 'r.outl.ols' = 'OLS residuals', 'clout' = 'outlier relationship')
melt.exonedf <- melt(exonedf, id = 'x')

ggplot(melt.exonedf[grepl('y|outl|y.clean|clout', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y.outl'=19, 'y.clean'=0, 'y.fit.outl.ols'=18, 'clout' = 5), labels = labels) +
  scale_color_manual(values = c('y.outl' = respring, 'y.clean' = clreen, 'y.fit.outl.ols' = olsalmon, 'clout' = brout), labels = labels) +
  geom_line(data = melt.exonedf[melt.exonedf$variable == "y.outl" | melt.exonedf$variable == "y.fit.outl.ols",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  #scale_linetype_manual("Residuals", values = 'dashed', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated and OLS Fit Values With Outliers") +
  coord_fixed()
### And plot the residuals.
ggplot(melt.exonedf[grepl('r.outl', melt.exonedf$variable),], aes(x = Re(value), y = Im(value))) +
  geom_point(color = olsalmon) +
  labs(y = "Imaginary", x = "Real", title = "Residuals of OLS Fit With Outliers") +
  coord_fixed() +
  geom_label(x = 20, y = -10, label = ssrlab(value)) +
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
labels <- c(labels, 'y.fit.outl.hub' = 'Huber rlm fitted response', 'r.outl.hub' = 'Huber rlm residuals')
melt.exonedf <- melt(exonedf, id = 'x')

ggplot(melt.exonedf[grepl('y|outl|y.clean|clout', melt.exonedf$variable),], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('y.outl'=19, 'y.clean'=0, 'y.fit.outl.ols'=18, 'y.fit.outl.hub' = 13, 'clout' = 5), labels = labels) +
  scale_color_manual(values = c('y.outl'= respring, 'y.clean'= clreen, 'y.fit.outl.ols'= olsalmon, 'y.fit.outl.hub' = hubsky, 'clout' = brout), labels = labels) +
  #geom_line(data = melt.exonedf[melt.exonedf$variable == "y.outl" | melt.exonedf$variable == "y.fit.outl.ols",], aes(group = x, lty = "residual"), size = 0.4, color = "lightcoral") +
  #scale_linetype_manual("Residuals", values = 'dashed', guide=guide_legend(override.aes = list(size = 0.5, color = "lightcoral"))) +
  labs(y = "Imaginary", x = "Real", shape = "Response Values", color = "Response Values", title = "Generated, OLS Fit, and Huber Robust Fit Values With Outliers") +
  coord_fixed()
### And plot the residuals.
ggplot(melt.exonedf[grepl('r.outl', melt.exonedf$variable),], aes(x = Re(value), y = Im(value))) +
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
### Then with Tukey's bisquare
fitone.outl.bis <- rlm(y.outl ~ x, exonedf, maxit = iterations + 90, acc = 10000 * accept, psi = psi.bisquare)
exonedf$y.fit.outl.bis <- fitted(fitone.outl.bis)
exonedf$r.outl.bis <- residuals(fitone.outl.bis)
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