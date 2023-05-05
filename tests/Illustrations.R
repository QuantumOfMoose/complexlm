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
ggplot(ex.one, aes(x =                                                                                                                                                                                                                   ))

plot(ex.one)
### Woops, not what I wanted...
ggplot(meltex.one, aes(x = Re(value), y = Im(value), color = variable)) +
  geom_point()
