####
#### William Ryan
#### 12 September, 2021
#### Test script for complex robust linear fit (rlm).
#### Modified from rlm in MASS.
#### Test script based on Wuber's answer on stackexchange https://stats.stackexchange.com/questions/66088/analysis-with-complex-data-anything-different#
####


library(complexlm)
library(MASS, include.only = "mvrnorm")
library(ggplot2)
# Synthesize data.
# (1) the independent variable `w`.
#
w.max <- 6 # Max extent of the independent values
w <- expand.grid(seq(-w.max,w.max, 3), seq(-w.max,w.max, 3))
w <- complex(real=w[[1]], imaginary=w[[2]])
w <- w[Mod(w) <= w.max] ### This drops any element of w that has a modulus greater than w.max
n <- length(w)
#
# (2) the dependent variable `z`.
#
beta <- c(-20+5i, complex(argument=2*pi/3, modulus=3/2)) ### The fist is the intercept, the 2nd the slope.
print(beta)
sigma <- 1.5; rho <- 0.7 # Parameters of the error distribution
set.seed(17)
e <- mvrnorm(n, c(0,0), matrix(c(1,rho,rho,1)*sigma^2, 2)) ### A bunch of random numbers to add errors to the example data.
e <- complex(real=e[,1], imaginary=e[,2]) ### Make those random numbers complex.
z <- as.vector((desX <- cbind(rep(1,n), w)) %*% beta + e) ### The design matrix desX is defined in this line.
z.pure <- as.vector((desX <- cbind(rep(1,n), w)) %*% beta) # z data without the noise.
# 
# Add m outliers to z. Draw the outliers from a 
#
z.clean <- z # Preserve outlier free z.
m <- 3
outlier.pos <- sample(1:n, m) # Positions of outliers in z.
outliers <- mvrnorm(m, c(5.8, -7.32),  matrix(c(1,rho,rho,1)*sigma^2, 2)) # The real and imaginary values of the outliers.
outliers <- complex(real=outliers[,1], imaginary=outliers[,2]) # Make the outliers comlex numbers.
z <- replace(z, outlier.pos, outliers)
#
# Collect everything into a dataframe.
#
fitframe <- data.frame(w, z.pure, z.clean, z)
#
# Whuber's ordinary complex linear fit.
#
print(beta.hat <- solve(Conj(t(desX)) %*% desX, Conj(t(desX)) %*% z), digits=4)
z.whuber <- beta.hat[2] * w + beta.hat[1]
fitframe$z.whuber <- z.whuber
fitframe$res.whuber <- as.vector(z - desX %*% beta.hat)
#
# Robust complex linear fit.
#
print(mytestfit <- rlm(x = w, y = z, interc = TRUE)) # Uses default psi=psi.huber, the Huber objective function.
rbeta.hat <- mytestfit$coefficients
fitframe$z.robust <- mytestfit$coefficients[2] * w + mytestfit$coefficients[1]
fitframe$res.robust <- mytestfit$residuals
# Robust complex linear fit, with Hampel objective function.
#
print(mytestfitHam <- rlm(x = w, y = z, psi = psi.hampel, interc = TRUE)) # Uses psi=psi.hampel, the Hampel objective function.
rHambeta.hat <- mytestfitHam$coefficients
fitframe$z.robustHam <- mytestfitHam$coefficients[2] * w + mytestfitHam$coefficients[1]
fitframe$res.robustHam <- mytestfitHam$residuals
# Robust complex linear fit, with Tukey's bisquare objective function.
#
print(mytestfitBi <- rlm(x = w, y = z, psi = psi.bisquare, interc = TRUE)) # Uses psi=psi.bisquare, Tukey's bisquare objective function.
rBibeta.hat <- mytestfitBi$coefficients
fitframe$z.robustBi <- mytestfitBi$coefficients[2] * w + mytestfitBi$coefficients[1]
fitframe$res.robustBi <- mytestfitBi$residuals

library(reshape2)
meltedfitframe <- melt(fitframe, id = "w")
meltedfitframe$variable <- factor(meltedfitframe$variable, ordered = TRUE)

betterplot2 <- ggplot(meltedfitframe[meltedfitframe$variable != "res.whuber" & meltedfitframe$variable != "res.robust" & meltedfitframe$variable != "res.robustMM" & meltedfitframe$variable != "z.clean",], aes(x = Re(value), y = Im(value), color = as.factor(variable), shape = as.factor(variable))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c('z.pure'=19, 'z.clean'=0, 'z'=20, 'z.whuber'=18, 'z.robust'=18), labels = c('z.pure'='pure relationsip', 'z.clean'='random noise added', 'z'='noise and outliers', 'z.whuber'='from ordinary fit', 'z.robust'='from robust fit')) +
  scale_color_manual(values = c('z.pure'="cyan3", 'z.clean'='forestgreen', 'z'='black', 'z.whuber'='red', 'z.robust'='blue'), labels = c('z.pure'='pure relationsip', 'z.clean'='random noise added', 'z'='noise and outliers', 'z.whuber'='from ordinary fit', 'z.robust'='from robust fit')) +
  geom_line(data = meltedfitframe[meltedfitframe$variable == "z" | meltedfitframe$variable == "z.whuber",], aes(group = w, lty = "ordinary lm"), size = 0.4, color = "lightcoral") +
  geom_line(data = meltedfitframe[meltedfitframe$variable == "z" | meltedfitframe$variable == "z.robust",], aes(group = w, lty = "robust lm"), size = 0.4, color = "royalblue") +
  scale_linetype_manual("Residuals", values = c('dotted', 'dotted'), guide=guide_legend(override.aes = list(size = c(0.5, 0.5), color = c("lightcoral", 'royalblue')))) +
  labs(y = "Imaginary", x = "Real", shape = "z values", color = "z values", title = "Generated and Fit Values of z") 
betterplot2 # Need to decide if I like cyan2 or cyan3 better.

#### Now plot the residuals.
residframe <- data.frame(w = w, pure.whub = fitframe$z.pure - fitframe$z.whuber, pure.rob = fitframe$z.pure - fitframe$z.robust,
                         clean.whub = fitframe$z.clean - fitframe$z.whuber, clean.rob = fitframe$z.clean - fitframe$z.robust,
                         zz.whuber = fitframe$res.whuber, zz.robust = fitframe$res.robust)
meltresidframe <- melt(residframe, id = "w")
assign_label <- function(varr) {
  if (varr == 'pure.whub') return('z.pure - z.olm (Ordinary residuals from pure relationsip.)')
  if (varr == 'pure.rob') return('z.pure - z.rlm (Robust residuals from pure relationship.)')
  if (varr == 'clean.whub') return('z.clean - z.olm (Ordinary residuals from noisy data.)')
  if (varr == 'clean.rob') return('z.clean - z.rlm (Robust residuals from noisy data.)')
  if (varr == 'zz.whuber') return('z - z.olm (Ordinary residuals from noisy data with outliers.)')
  if (varr == 'zz.robust') return('z - z.rlm (Robust residuals from noisy data with outliers.)')
}
meltresidframe$label <- as.factor(unlist(lapply(meltresidframe$variable, assign_label))) ### This gives us more readable labels for the facets of the next plot.

residplot <- ggplot(meltresidframe, aes(x = Re(value), y = Im(value), color = as.factor(variable))) + 
  geom_segment(aes(x = 0, y = 0, xend = Re(value), yend = Im(value)), arrow = arrow(length = unit(0.03, "npc")), size = .4) +
  facet_wrap(vars(label), nrow = 3, ncol = 2) +
  labs(y = "imaginary", x = "real", title = "Complex Residuals") +
  theme(legend.position = "none") +
  coord_fixed()
residplot

absresidplot <- ggplot(meltresidframe, aes(y = Mod(value)^2, x = label, color = Arg(value))) +
  geom_jitter(width = 0.1) + 
  scale_colour_gradientn(colours=rainbow(6, start = 0, end = 1, s = 1, v = 0.8)) +
  labs(y = "Modulus Squared", x = "Residual", color = 'Argument', title = "Squared Modulus of Residuals")
absresidplot

## A bunch of plots showing the transformation from w to z
ggplot(fitframe, aes(x = Re(w), y = Im(w))) +
  geom_point(size = 3) +
  labs(y = "Imaginary", x = "Real", title = "Independant Variable (Current)") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(size = 18)) +
  coord_fixed()
ggplot(fitframe, aes(x = Re(z.pure), y = Im(z.pure))) +
  geom_point(size = 3, shape = 19, color = 'cyan2') +
  labs(y = "Imaginary", x = "Real", title = "Dependant Variable (Voltage)") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(size = 18)) +
  coord_fixed() +
  ylim(-15,15) + xlim(-28, 2)
ggplot(fitframe, aes(x = Re(z.clean), y = Im(z.clean))) +
  geom_point(size = 3, shape = 0, color = "forestgreen") +
  labs(y = "Imaginary", x = "Real", title = "Noisy Dependant Variable (Voltage)") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(size = 18)) +
  coord_fixed() +
  ylim(-15,15) + xlim(-28, 2)
ggplot(fitframe, aes(x = Re(z), y = Im(z))) +
  geom_point(size = 3, shape = 20, color = "black") +
  labs(y = "Imaginary", x = "Real", title = "Noisy Dependant Variable with Outliers (Voltage)") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(size = 18)) +
  coord_fixed() +
  ylim(-19,19) + xlim(-28, 8)

####
#### A test to compare the performance of zmodel.matrix to model.matrix
####
library(profvis)
set.seed(4242)
slop <- 4.23
interc <- 1.4
Xt <- -20:20
tframe <- data.frame(Xt=Xt, Yt= Xt * slop + interc + rnorm(length(Xt)))
testerms <- terms(Yt ~ Xt)
zmodel.matrix(testerms, tframe)
model.matrix(testerms, tframe)
