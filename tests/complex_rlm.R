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
s.whuber <- sqrt(Re(mean(Conj(fitframe$res.whuber)*fitframe$res.whuber))) ### This might be standard deviation or something similar.
#
# Robust complex linear fit.
#
print(mytestfit <- rlm(x = w, y = z, interc = TRUE)) # Uses default psi=psi.huber, the Huber objective function.
rbeta.hat <- mytestfit$coefficients
fitframe$z.robust <- mytestfit$coefficients[2] * w + mytestfit$coefficients[1]
fitframe$res.robust <- mytestfit$residuals
s.robust <- sqrt(Re(mean(Conj(fitframe$res.robust)*fitframe$res.robust))) ### This might be standard deviation or something similar.
# Robust complex linear fit, with Hampel objective function.
#
print(mytestfitHam <- rlm(x = w, y = z, psi = psi.bisquare, interc = TRUE)) # Uses default psi=psi.huber, the Huber objective function.
rHambeta.hat <- mytestfit$coefficients
fitframe$z.robustHam <- mytestfit$coefficients[2] * w + mytestfit$coefficients[1]
fitframe$res.robustHam <- mytestfit$residuals
s.robustHam <- sqrt(Re(mean(Conj(fitframe$res.robust)*fitframe$res.robust))) ### This might be standard deviation or something similar.
#
# Robust complex linear fit, MM method.
#
print(mytestfitMM <- rlm(x = w, y = z, interc = TRUE, method = "MM"))
rMMbeta.hat <- mytestfitMM$coefficients
fitframe$z.robustMM <- mytestfitMM$coefficients[2] * w + mytestfitMM$coefficients[1]
fitframe$res.robustMM <- mytestfitMM$residuals
s.robustMM <- sqrt(Re(mean(Conj(fitframe$res.robustMM)*fitframe$res.robustMM))) ### This might be standard deviation or something similar.
#
# Show some diagnostics.
#
# par(mfrow=c(1,2))
# library(grid)
# col <- hsv((Arg(res)/pi + 1)/2, .8, .9)
# size <- Mod(res) / s
# resplot.whuber <- ggplot(fitframe, aes(x = Re(res.whuber), y = Im(res.whuber), color = Arg(res.whuber), size = Mod(res.whuber))) +
#   geom_point() +
#   theme(legend.position="none")
# resplot.robust <- ggplot(fitframe, aes(x = Re(res.robust), y = Im(res.robust), color = Arg(res.robust), size = Mod(res.robust))) +
#   geom_point() +
#   theme(legend.position="none")
# fitplot.whuber <- ggplot(fitframe, aes(x = Re(z.whuber), y = Im(z.whuber), color = Arg(z.whuber), size = Mod(z.whuber))) +
#   geom_point() +
#   theme(legend.position="none")
# fitplot.robust <- ggplot(fitframe, aes(x = Re(z.robust), y = Im(z.robust), color = Arg(z.robust), size = Mod(z.robust))) +
#   geom_point() +
#   theme(legend.position="none")
# 
# grid.draw(cbind(rbind(ggplotGrob(resplot.whuber), ggplotGrob(fitplot.whuber)),rbind(ggplotGrob(resplot.robust), ggplotGrob(fitplot.robust))))
#plot(res, pch=16, cex=size, col=col, main="Residuals")
#plot(Re(fit), Im(fit), pch=16, cex = size, col=col,
#     main="Residuals vs. Fitted")
# 
# plot(Re(c(z, fit)), Im(c(z, fit)), type="n",
#      main="Residuals as Fit --> Data", xlab="Real", ylab="Imaginary")
# points(Re(fit), Im(fit), col="Blue")
# points(Re(z), Im(z), pch=16, col="Red")
# arrows(Re(fit), Im(fit), Re(z), Im(z), col="Gray", length=0.1)
# 
# col.w <-  hsv((Arg(w)/pi + 1)/2, .8, .9)
# plot(Re(c(w, z)), Im(c(w, z)), type="n",
#      main="Fit as a Transformation", xlab="Real", ylab="Imaginary")
# points(Re(w), Im(w), pch=16, col=col.w)
# points(Re(w), Im(w))
# points(Re(z), Im(z), pch=16, col=col.w)
# arrows(Re(w), Im(w), Re(z), Im(z), col="#00000030", length=0.1)

#
# Visualize the data and different fits.
#
# ReRePlot <- ggplot(fitframe, aes(x = Re(w), y = Re(z))) +
#   geom_point() +
#   stat_function(fun = function(x) , color = "blue") +
#   geom_line(aes(x = Re(w), y = Re(z.whuber)), color = "red")
# ReRePlot
# ImImPlot <- ggplot(fitframe, aes(Im(w), Im(z))) +
#   geom_point() +
#   geom_line(aes(x = Im(w), y = Im(z.robust)), color = "blue") +
#   geom_line(aes(x = Im(w), y = Im(z.whuber)), color = "red")
# ImImPlot
# ReRePlot <- ggplot(fittestdf, aes(Re(w), Re(z))) +
#   geom_point() +
#   geom_line(aes(x = Re(w), y = Re(z.robust)), color = "blue") +
#   geom_line(aes(x = Re(w), y = Re(z.whuber)), color = "red")
# 
# ImImPlot <- ggplot(fittestdf, aes(Re(w), Im(z))) +
#   geom_point() +
#   geom_line(aes(x = Im(w), y = Im(z.robust)), color = "blue") +
#   geom_line(aes(x = Im(w), y = Im(z.whuber)), color = "red")
# 
# grid.draw(rbind(ggplotGrob(RealPlot),ggplotGrob(ImagPlot)))
# 
# par(mfrow=c(1,1))
# pairs(cbind(w.Re=Re(w), w.Im=Im(w), z.Re=Re(z), z.Im=Im(z),
#             ordfit.Re=Re(fitframe$z.whuber), ordfit.Im=Im(fitframe$z.whuber), robfit.Re = Re(fitframe$z.robust), robfit.Im = Im(fitframe$z.robust)), cex=1/2)
# 
# pairs(fitframe)

# betterplot <- ggplot(fitframe, aes(x = Re(w), y = Im(w))) + 
#   #geom_point(aes(color = Arg(z), size = Mod(z))) + scale_fill_viridis_c() +
#   geom_point(aes(x = Re(z.pure), y = Im(z.pure), color = Arg(z.pure)), size = 3) +
#   geom_point(aes(x = Re(z.pure), y = Im(z.pure)), color = "cyan2", size = 3, label = "pure z") +
#   geom_point(aes(x = Re(z.clean), y = Im(z.clean)), color = "forestgreen", size = 3, shape = 0, label = "clean z") +
#   geom_point(aes(x = Re(z), y = Im(z)), label = "z") +
#   geom_point(aes(x = Re(z.whuber), y = Im(z.whuber)), color = "red", shape = 17, size = 2, label = "olm z") +
#   geom_point(aes(x = Re(z.robust), y = Im(z.robust)), color = "blue", shape = 18, size = 3, label = "rlm z") +
#   geom_segment(aes(x = Re(z), y = Im(z), xend = Re(z.whuber), yend = Im(z.whuber)), linetype = "dotted", size = 0.4, color = "lightcoral", group = "olm resid") + 
#   geom_segment(aes(x = Re(z), y = Im(z), xend = Re(z.robust), yend = Im(z.robust), lty = "rlm resid"), linetype = "dotted", size = 0.4, color = "royalblue1", group = "rlm") + 
#   theme(legend.position = "right")
# betterplot


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
