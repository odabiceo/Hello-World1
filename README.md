# Hello-World1
I need help with an extra credit assignment from my Micro Biology Class

I'm new to the world of coding and I have to figure out the problem below to pass my class, I'm a senior in college and if I don't get this assingment done I won't pass my classes. :-( I'm so lost and confused please help?

Here is the assignment, I can't make sense of it:

s## Example of fitting prevalence data with observation error only.
## Steve Bellan, Meaningful Modeling of Epidemiologic Data 2012
## AIMS, Muizenberg

## First we need to create a system of ODEs.  We could also have used
## discrete time formulations.
library(deSolve)

## This is the model from Lab 1
sir <- function(t,y,parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    beta <- gamma*R0
    dSdt <- -beta*S*I/N
    dIdt <- beta*S*I/N - gamma*I
    # Note: Population size is constant, so don't need to specify dRdt
    list(c(dSdt,dIdt))
  })
}


time <- 0                               # set initial time
N0 <- 7781984                           # total population size
initial.SI <- c(S = 0.1*N0, I = 20.5) # Initially 10% of the population is susceptible

values <- c(R0 = 20,                # set R0
            gamma = 1/8,            # 1 / infectious period = 1/8 days
            N = N0)                 # population size (constant)

## Run ODE for a year with daily time steps output
time.out <- seq(0,365,by = .1)

## This is how we run ODE's in R (see Lab 1 for futher details)
ts.sir <- data.frame(lsoda(
  y = initial.SI,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = sir,                   # Function to evaluate
  parms = values                # Vector of parameters
  ))

## What is the prevalence of disease over time?
###################################################################### 
## Task 1: Replace ?? with the appropriate values to get our prevalence curve.
ts.sir$P <- ts.sir$?? / ??
  
## Look at our new data frame
head(ts.sir)

## Plot prevalence
plot(ts.sir$time, ts.sir$P, type="l", col="red")

######################################################################
## Generating the data from our simulation
######################################################################
## Say we sampled the population weekly during the epidemic, at the
## following time points:
time.samp <- seq(40,180, by = 7)

## What is the prevalence at those times?
prev.samp <- ts.sir$P[ts.sir$time %in% time.samp]
## Task 2: Figure out what %in% does and explain how it is used in prev.samp above.
c(-3,2,5,100,8,15) %in% 1:10

## And at those times we sampled 1000 people from the population
samp.sizes <- rep(1000, length(prev.samp))
## Sample from 1000 people at each time point based on the prevalence
## at that time point.
num.inf.samp <- rbinom(length(time.samp), size = samp.sizes, prob = prev.samp)
num.inf.samp

## Make a data frame of time sampled & number positive
mydata <- data.frame(time = time.samp,
                     num.inf.samp = num.inf.samp,
                     samp.sizes = samp.sizes)
mydata

## Calculate binomial confidence intervals for our samples so we can
## plot them. We could use the methods used in Lab 5, but it turns out
## that R has a function binom.test() that does this for us
## automatically.
?binom.test

## What is the confidence interval for HIV prevalence if we sampled
## 100 people from a population and 28 were HIV+?
binom.test(28, 100)
binom.test(28, 100)$conf.int[1:2]

## For each of our samples get CIs
for(ii in 1:length(num.inf.samp))
  {
    mydata$ci.l[ii] <- binom.test(mydata$num.inf.samp[ii], 1000)$conf.int[1]
    mydata$ci.u[ii] <- binom.test(mydata$num.inf.samp[ii], 1000)$conf.int[2]    
  }

## Plot data vs the true underlying epidemic.
plot(ts.sir$time, ts.sir$P, type="l", col="red", bty = "n",
     ylim = c(0, max(mydata$ci.u)), xlab = "day", ylab = "prevalence")
## Add data
points(mydata$time, mydata$num.inf.samp/1000, col = "black", pch = 19)
## Add CIs
arrows(mydata$time, mydata$ci.l, mydata$time, mydata$ci.u, length = .01,
       angle = 90, code = 3)
## Add legends
legend("topright", c("true prevalence", "sample prevalence estimate"),
       pch = 19, col=c("red","black"))

## As a start, we need to write a -log likelihood function that gives
## the probability that a given (beta,repp) would generate the
## observed data. Remember that if we are assuming that there is some
## true underlying epidemic curve that is deterministic and the data
## we observe are only noisy because of sampling/observation error
## (not because the underying curve is also noisy--i.e. process
## error--which is particularly likely for small epidemics).

## So we have binomial sampling errors. For instance the negative log
## likelihood of the true data is:

nll.true <- - sum(dbinom(mydata$num.inf.samp, samp.sizes, prev.samp, log = TRUE))

## But we want to try to estimate the best fit beta's and sigma's
## given our observed data. So let's write a function that takes beta
## & sigma as input and gives the negative log likelihood as output.

nll.fn <- function(logpars,                # log(R0) & log(gamma) in a vector, must be named
                   data = num.inf.samp,             # number positive in sample
                   n = samp.sizes,                        # sample sizes
                   pop.size = N0,
                   times = time.samp) # the times they were observed at
  {
    pars <- c(exp(logpars), N = pop.size)
        if(pars["R0"]<1) stop("R0<1 and epidemic doesn't take off, change your parameters")
    ts.sir.temp <- data.frame(lsoda(y = initial.SI, times = time.out,
                                    func = sir, parms = pars, atol = 1e-15))
    ts.sir.temp$P <- ts.sir.temp$I / pop.size
    prev.samp.temp <- ts.sir.temp$P[ts.sir.temp$time %in% times]
    nll <- - sum(dbinom(data, n, prev.samp.temp, log = TRUE))
    nll
  }


## Check that our function does everything we did above.
nll.fn(log(values[1:2]))
nll.true

## First remind yourself how optim() works. The more you read through
## the help file the easier this will be!!! In particular make sure
## you understand that the first argument of optim must be the initial
## values of the parameters to be fitted (i.e. beta & repp) and that
## any other parameters to be fixed are given as additional arguments
## (in the helpfile under "...")
?optim

## Select initial values for fitted parameters from which optimization
## routine will start. If you select bad initial values the algorithm
## can get stuck on a bad set of parameters. You can always try the
## true values as a starting point for this problem, although that's
## rarely possible in real problems.

init.pars <- c(R0 = log(10), gamma = log(.05))
## We will start with SANN optimization since it is stochastic and
## therefore less likely to get stuck in a local minima. But then
## finish with Nelder-Mead optimization which is much faster.

###  NOTE: for trace >0 you see more progress report, bigger numbers show more update
trace <- 3

## SANN: This is stochastic, be CAREFUL sometimes it gets stuck at
## local minima for unreasonble parameters. If you see this happen,
## run it again!
optim.vals <- optim(par = init.pars,
                    nll.fn,
                    control = list(trace = trace, maxit = 150),
                    method = "SANN")
    
optim.vals <- optim(par = optim.vals$par,
                    nll.fn,
                    control = list(trace = trace, maxit = 500, reltol = 10^-7),
                    method = "Nelder-Mead")
## While we use Nelder-Mead, normally you would use SANN first and
## then follow with Nelder-Mead since SANN is stochastic and will make
## sure to help you be sure that you aren't at local minima.

## Task 3: Look at the output of optim. Understand what it means. Did
## the algorithm converge? Look at ?optim to understand it.
optim.vals

## Compare the optimized values to the true values
exp(optim.vals$par)                     # fitted values
values[c("R0","gamma")]                 # True values

## Let's plot our true data vs the fitted plot. To do this we first
## must run a simulation with our optimized parameter values.

## This is how we run ODE's in R (see Lab 1 for futher details)
fitted.values <- c(exp(optim.vals$par), N = N0)
ts.sir.fitted <- data.frame(lsoda(y = initial.SI, times = time.out,
                                  func = sir, parms = fitted.values))
ts.sir.fitted$P <- ts.sir.fitted$I / N0
##  Now add a dashed red line to our plot showing the fitted prevalence
##  curve.
lines(ts.sir.fitted$time, ts.sir.fitted$P, col = "red", lty = 2)


######################################################################
## Contour plots

## With all other parameters fixed to their initial values, lets look
## at a contour likelihood plot over sigma and R0.  To do this we
## write wrapper functions of xx (R0) and yy (sigma) to feed to
## outer() and then contour(). This is confusing so make sure you
## understand every function.

## Thi function simply takes values xx and yy and feeds them into
## nll.fn above as a single variable called logpars.
nll2 <- function(xx,yy, browse=F)
    nll.fn(logpars = c(R0 = log(xx), gamma = log(yy)))
  }

## Now instead of giving a single argument on the log scale we give 2
## on the untransformed scale.
nll.fn(c(R0 = log(25), gamma = log(1/5)))
nll2(25, 1/5)

## If we try to give this function multiple values of R0 or gamma,
## however, it gets confused and gives an error because it only wants
## one value for each.
nll2(c(25:26), c(1/5,1/6))

## So we "Vectorize" this function so it can take in vectors of the
## parameters and return the output. nll3 hen calls on nll2() to take
## xx,yy as pairwise values and nll2() for all pairs
nll3 <- Vectorize(nll2, list("xx","yy"))

## Task 4: Explain how the following three lines are related.
nll2(25, 1/5)
nll2(26, 1/6)
nll3(c(25:26), c(1/5,1/6))

## Now we use the R function outer() to evaluate nll3() over a grid of
## {R0, gamma} combinations. This can take a long time because we have
## to do res^2 evaluations of nll.fn(), and recall that each time we
## do this we are running lsoda() inside nll.fn()
res <- 60                               # Grid resolution 60^2 = 3600

## Now create a sequence of R0 values for the grid
logR0.fit <- optim.vals$par["R0"]
## let's have the sequence be
R0.seq <- exp(seq(logR0.fit-.5, logR0.fit + .2, l = res))
R0.seq


## Now create a sequence of gamma values for the grid
loggamma.fit <- optim.vals$par["gamma"]
gamma.seq <- exp(seq(loggamma.fit-.4, loggamma.fit+.4, l = res))
gamma.seq
1/gamma.seq                             # more intuitively the range
                                        # of infectious periods we are
                                        # examining.

## The function outer() now evaluates nll3 on this grid.
?outer
mat <- outer(R0.seq, gamma.seq, nll3) # this can take a long time

## Make a contour plot that shows the confidence intervals in red.
## Likelihood Ratio Test confidence intervals uses chi squared
## distribution cutoff with degrees of freedom 2 (2 parameters)
ml.val <- optim.vals$value
conf.cutoff <- ml.val + qchisq(.95,2)/2

## Show likelihood contours
contour(R0.seq, gamma.seq,mat, xlab = "R0", ylab = "gamma",
        main = "NLL contours", bty = "n")
## Add red contour for 95% CI
contour(R0.seq, gamma.seq,mat, xlab = "R0", ylab = "gamma", levels = c(conf.cutoff),
        col = "red", lwd = 2, labels = "", labcex = .2, add = T)
## Plot MLE
points(exp(logR0.fit),exp(loggamma.fit), pch = 19)
text(exp(logR0.fit),exp(loggamma.fit), "MLE", pos = 4)
legend("topright", "95% Confidence Interval", col = "red", lty = 1, bg = "white")


## Challenge Question: Try to fit this same model while estimating S0
## (i.e. the susceptible proportion at the start of the outbreak.
