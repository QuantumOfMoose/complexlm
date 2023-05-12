###
### This is a companion to the paper which introduces the complexlm package. 
### It demonstrates the functions in said package by using them to explore and 
### analyze two sets of complex data. It produces all the figures featured in the
### Illustrations section of the paper.
### 
### William Ryan
###
### 4 May, 2023
###

### First we must load up some libraries that will help us with generating, manipulating, and visualizing the data.
library(stats)
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggforce)
### Now load complexlm. It is important to do so after loading stats and MASS, as complexlm masks several functions from those packages.
library(complexlm)

### With that we are ready for our first example problem. It will be a fully two dimensional problem, based on randomly generated data. 
### To ensure repeatability we set the random seed.
set.seed(4242)
### Consider a set of points in two-dimensional space. They could represent stars in the sky, a herd of elephants viewed from above by a drone,
### cities on a map, or something else entirely. We will generate such data from a symmetrical two dimensional normal distribution with mvrnorm from MASS.
nn <- 40
sigma <- matrix(c(3,0,0,3), 2, 2)
locs <- mvrnorm(n = nn, mu = c(9.3, 6.2), Sigma = sigma) # Arbitrary mean.

### We need this to be complex, so convert the numeric matrix locs into the complex vector clocs.
clocs <- complex(real = locs[,1], imaginary = locs[,2])
### And put them in a dataframe for ease of use.
ex.one <- data.frame(clocs = clocs)

### Now let's plot this to see what we have to work with.
ggplot(ex.one, aes(x = Re(clocs), y = Im(clocs))) +
  geom_point() +
  labs(x = "Real Position", y = "Imaginary Position", title = "Initial Positions")

### Next we linearly transform these initial points. If these points represent the locations of objects in a field of view,
### such a mathematical operation could represent the objects moving as a group, or the imaging device rotating, zooming, and/or panning.
### We once again call on mvrnorm to generate random transformation coefficients.
transformer <- mvrnorm(n = 2, mu = c(0,0), Sigma = 2.5*sigma)
beta <- complex(real = transformer[,1], imaginary = transformer[,2]) # A length two complex vector, the first element will translate, the second will rotate and scale.

### Of course there isn't much point in regression if there is no error or variance in the transformed data; so we randomly generate another complex vector, err,
### which represents the uncertainty / variance / error in our response data. Thus the transformed / measured / response data is tclocs = beta*clocs + err.
err <- mvrnorm(n = nn, mu = c(0,0), Sigma = sigma / 6)
ex.one$err <- complex(real = err[,1], imaginary = err[,2])

### Calculate the transformed / measured / response data.
ex.one$tclocs <- as.vector(cbind(rep(1, nn), ex.one$clocs) %*% beta) + ex.one$err

### And then visualize them on the same plot as clocs.
ggplot(ex.one, aes(x = Re(clocs), y = Im(clocs))) +
  geom_point() +
  geom_point(aes(x = Re(tclocs), y = Im(tclocs)), color = "red") +
  labs(x = "Real Position", y = "Imaginary Position", title = "Transformed Positions")

### Since we simulated these transformed points ourselves, it's very easy to find the transformation coefficients. Just type "print(beta)", and there you go.
### Obviously, that doesn't work for most sets of data we encounter in the wild. So we'll use the lm function of complexlm instead.
lm.xone <- lm(tclocs ~ clocs, ex.one)

### And we can check how well it works by comparing the estimated beta to the one used to generate the data.
print(betahat <- lm.xone$coefficients)
print(beta)
### Those look pretty similar. Let's calculate the percent error as a rudimentary way to quantify the accuracy.
print(perror <- 100 * (betahat - beta) / beta )
### Not bad, especially the rotation and scaling coefficient.

### Take the fitted values and residuals and add them to ex.one
ex.one$fitted <- as.vector(lm.xone$fitted.values)
ex.one$resids <- as.vector(lm.xone$residuals)

### Plotting complex linear fits like this one is more difficult and complex than plotting linear models of real variables. To facilitate
### this task we will build a seperate dataframe that contains the same information as ex.one, but with the real and imaginary components 
### of each variable seperated into different columns.
rex.one <- mutate(ex.one, across(.cols = everything(), .fns = list(Re = Re, Im = Im), .names = "{.fn}[{.col}]"), .keep = "unused")

### Now use this to make a big plot matrix of the predictor, response, and predicted variables in rex.one.
plot(select(rex.one, !contains(c("resids", "err"))))
### That looks pretty good. 

### One last test before moving on to robust regression. As we saw in section [theory], the residuals should be equal almost surely. Let's
### take their difference.
summarize(erres.one <- ex.one$err - ex.one$resids)

### Were we performing real valued regression, this is the point at which we would plot a line on a scatter plot with the predictor variables on the horizontal
### axis and the response variables on the vertical. But, since complex numbers are two dimensional, we can't draw such a graph which shows the entire relationship
### between the predictor and response variables. Instead we will draw three graphs .....
### To make plotting easier, we'll melt ex.one to make a long dataframe.
#meltex.one <- melt(ex.one, id.vars = "clocs", measuure.vars = c("tclocs", "err", "fitted", "resids")) # This doesn't make plotting much easier. Skip it.
ggplot(ex.one, aes(x = Re(clocs), y = Re(tclocs))) +
  geom_point(aes(color = 'Real'), show.legend = T) +
  geom_point(aes(x = Im(clocs), y = Im(tclocs), color = 'Imaginary'), show.legend = T) +
  geom_smooth(aes(x = Re(clocs), y = Re(fitted)), color = "#00BFC4", method = "lm", se = F) +
  geom_smooth(aes(x = Im(clocs), y = Im(fitted)), color = "#F8766D", method = "lm", se = F) +
  labs(x = "Predictor", y = "Response") +
  theme(legend.position = "right")
### Ok, that plot looks pretty good.

### Now we will add some outliers. Say that there is a probability theta < 0.5 that each element of the data was not transformed by beta, but by beta.prime.
### We will arbitrarily set,
theta <- 0.1
### And randomly generate beta.prime in much the same way as we did for beta, except we will increase the variance.
set.seed(623284)
transformer <- mvrnorm(n = 2, mu = c(3,-2), Sigma = 5*sigma)
beta.prime <- complex(real = transformer[,1], imaginary = transformer[,2])
### We will make some new errors for the outliers as well.
err.outl <- mvrnorm(n = nn, mu = c(0,0), Sigma = sigma / 2)
ex.one$err.outl <- complex(real = err[,1], imaginary = err[,2])
### Add a column to ex.one of logicals generated by a single trial binomial distribution with probability theta.
ex.one$outlier <- as.logical(rbinom(n, 1, theta))
### Check the actual proportion of outliers.
print(outlier.proportion.exone <- sum(ex.one$outlier) / length(ex.one$outlier))
### Now add another column to ex.one populated with clocs * beta or clocs * beta.prime, depending on the value of ex.one$outlier. "oat" is intentional, it lets us order things in legends of future plots.
ex.one <- ex.one %>% mutate(oatl.tclocs = if_else(outlier, beta.prime[1] + clocs * beta.prime[2] + err.outl, beta[1] + clocs * beta[2] + err))

### Plot oatl.tclocs.
ggplot(ex.one, aes(x = Re(oatl.tclocs), y = Im(oatl.tclocs))) +
  geom_point() +
  geom_point(aes(x = Re(tclocs), y = Im(tclocs)), color ='red', shape = 1) +
  labs(x = "Real", y = 'Imaginary')

### Next we will perform four complex regressions on oatl.tclocs; an ordinary least-squares fit, and robust regressions using each of the three psi functions included in complexlm.
accept = 1e-6 # Shared acceptance criterion
iterations =70 # Shared max iterations
### First, the ols.
ex.one.ols <- lm(oatl.tclocs ~ clocs, ex.one)
ex.one$ols.fit <- as.vector(ex.one.ols$fitted.values)
ex.one$ols.resid <- as.vector(ex.one.ols$residuals)

### Next, the rlm with the default Huber loss psi function.
ex.one.hub <- rlm(oatl.tclocs ~ clocs, ex.one, maxit = iterations, acc = accept)
ex.one$hub.fit <- as.vector(ex.one.hub$fitted.values)
ex.one$hub.resid <- as.vector(ex.one.hub$residuals)

### Now we use the Hampel psi function for M-estimation.
ex.one.ham <- rlm(oatl.tclocs ~ clocs, ex.one, psi = psi.hampel, maxit = iterations, acc = accept)
ex.one$ham.fit <- as.vector(ex.one.ham$fitted.values)
ex.one$ham.resid <- as.vector(ex.one.ham$residuals)

### Finally we turn to Tukey's bisquare psi function.
ex.one.bis <- rlm(oatl.tclocs ~ clocs, ex.one, psi = psi.bisquare, maxit = iterations, acc = accept)
ex.one$bis.fit <- as.vector(ex.one.bis$fitted.values)
ex.one$bis.resid <- as.vector(ex.one.bis$residuals)

fitnames <- c('ols', 'hub', 'ham', 'bis') # Collect the name component of each kind of fit in a vector. For use with future call() functions.
### Plotting all of these will be tricky, but essential to seeing the benefits of the robust fits, so we will make a long form dataframe with only the colums that we plan to plot.
ex.one.plotdf <- melt(select(ex.one, !c(clocs, err, err.outl, outlier, tclocs))) %>% arrange(variable)
### This character vector will be used to produce a column of human readable labels.
assignLabels <- function(variable){
  if(variable == 'fitted') return("Least-Squares Fit to Uncontaminated Data")
  if(variable == 'resids') return("Residuals of Uncontaminated OLS")
  if(variable == 'oatl.tclocs') return('Contaminated Response')
  if(variable == 'ols.fit') return('Least-Squares Fit')
  if(variable == 'ols.resid') return('Residuals of Least-Squares Fit')
  if(variable == 'hub.fit') return('Robust Fit with Huber Psi')
  if(variable == 'hub.resid') return('Residuals of Robust Fit with Huber Psi')
  if(variable == 'ham.fit') return('Robust Fit with Hampel Psi')
  if(variable == 'ham.resid') return('Residuals of Robust Fit with Hampel Psi')
  if(variable == 'bis.fit') return('Robust Fit with Tukey Bisquare Psi')
  if(variable == 'bis.resid') return('Residuals of Robust Fit with Tukey Bisquare Psi')
}
ex.one.plotdf$label <- as.factor(unlist(lapply(ex.one.plotdf$variable, assignLabels)))

#### PPPPPPPPLLLLLLOOOOOOTTTTTTTTSSSSSSS (don't put that in the paper)
ggplot(filter(ex.one.plotdf, (!grepl('resid', variable, ignore.case = TRUE) & variable != 'fitted')), aes(x = Re(value), y = Im(value), color = label, shape = label)) +
  geom_point() + # Formerly called exoneplot
  labs(x = 'Real', y = 'Imaginarry', shape = 'Data Origin', color = 'Data Origin') +
  geom_ellipse(aes(x0 = 70, y0 = 19, a = 10, b = 19, angle = -0.84), linewidth = .03, color = 'red') +
  annotate('text', x = 70, y = 19, label = 'Outliers')
  #scale_shape_manual(labels = fitlabels[!grepl('resid', fitlabels)])
### As can be seen in this plot, the fitted data obtained by robust M-estimation hews much closer to the non-outlier response data.

### Now we'll look at the residuals of each fitting method.
residplot <- ggplot(filter(ex.one.plotdf, (grepl('resid', variable, ignore.case = TRUE) & variable != 'resids')), aes(x = Re(value), y = Im(value), color = as.factor(variable))) + 
  geom_segment(aes(x = 0, y = 0, xend = Re(value), yend = Im(value)), arrow = arrow(length = unit(0.03, "npc")), size = .4) +
  facet_wrap(vars(variable), nrow = 3, ncol = 2) +
  labs(y = "imaginary", x = "real", title = "Complex Residuals") +
  theme(legend.position = "none") +
  coord_fixed()
residplot
### The set of residuals from all these regression techniques are dominated by those of the 5 outliers. However, ols residuals of the other, uncontaminated
### data points are much larger than their contemporaries from any of the M-estimator fits. 

### Finally, we will have a look at the coefficient estimates produced by these regressions, and their uncertainties. 
### We shall see which best extracted the relationship between the predictor and response variables.
summarydf.exone <- data.frame()

plot(ex.one)
### Woops, not what I wanted...
ggplot(meltex.one, aes(x = Re(value), y = Im(value), color = variable)) +
  geom_point()
