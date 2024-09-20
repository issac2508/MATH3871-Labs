#
# MATH3871 Lab 1
#

######### Question 1: Melbourne water

# a)

# In the tutorial, we computed by integration that the posterior mean of theta
# is given by a/b, where a = alpha + sum(x) and b = beta + n.

n <- 65
sum <- 24890
alpha <- 1
beta <- 0.01
theta_mean <- (alpha + sum) / (beta + n)
print(sprintf("Posterior mean of theta: %f", theta_mean))

# The 95% credible interval has the same idea as a confidence interval - we find
# the 2.5% and 97.5% percentiles, and there will be a 95% probability that
# theta will lie in this interval. We can interpret it as a direct probability
# (unlike the confidence intervals) because theta itself is a random variable.

shape <- alpha + sum
rate <- beta + n

left <- qgamma(0.025, shape=shape, rate=rate)
right <- qgamma(0.975, shape=shape, rate=rate)
print(sprintf("95%% central credible interval: (%f, %f)", left, right))

# b)

# Function to estimate the 95% credible interval by sampling 'n' points from
# the posterior and getting their 2.5% and 97.5% quantiles.
credible_int_est <- function(n) {
    points <- rgamma(n, shape=shape, rate=rate) # simulate points from posterior
    credible_int_estimate <- quantile(points, probs=c(0.025, 0.975))
    return(credible_int_estimate)
}

# Given a list of upper and lower bound estimates for the credible interval,
# plots a histogram.
plot_histogram <- function(upper, lower) {
    domain <- range(c(upper, lower))

    hist(
        lower, c=rgb(1, 0, 0, 0.5),
        main="Monte Carlo Estimates of 95% Credible Interval",
        xlim=domain,
        xlab="Estimate", ylab="Frequency"
    )
    hist(
        upper, col = rgb(0, 0, 1, 0.5), add = TRUE # add second histogram
    )
    legend(
        "top", legend = c("2.5% Quantile", "97.5% Quantile"),
        fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5))
    )
    points(x=left, y = 0, col = "black", pch = 19, cex = 1.5) # add true bounds
    points(x=right, y = 0, col = "black", pch = 19, cex = 1.5)
}

# Function to estimate the 95% credible interval using 'n_samples' from the
# posterior distribution and repeat the process 'n_runs' times.
simulate_credible_int <- function(n_samples, n_runs) {
    lower <- numeric(0)
    upper <- numeric(0)

    for (i in 1:n_runs) {
        estimate <- credible_int_est(n_samples)
        lower <- c(lower, estimate[[1]])
        upper <- c(upper, estimate[[2]])
    }
    return(list(upper=upper, lower=lower))
}

# Plot the result after estimating the 95% credible interval using 5000 samples
# and repeating this 250 times.
values <- simulate_credible_int(5000, 250)
upper <- values$upper
lower <- values$lower
plot_histogram(upper, lower)

# Find out what sample size is required to get the spread for each estimate
# below 0.15
spread <- 0.15
num <- 55000
while (TRUE) {
    values <- simulate_credible_int(num, 250)
    range_upper <- max(values$upper) - min(values$upper)
    range_lower <- max(values$lower) - min(values$lower)

    if (range_upper < spread && range_lower < spread) {
        print(sprintf("The spread for each estimate is less than %.2f when n is %d", spread, num))
        break
    }
    num <- num + 5000
}

# c)

# Each new observation follows a Poisson distribution with the parameter theta.
# So, we draw samples of theta from the posterior to draw samples of new
# observations from the Poisson distribution using the theta value.

num_obs <- 1000
new_observations <- numeric(0)

theta_est <- rgamma(num_obs, shape=shape, rate=rate)
new_observations <- rpois(num_obs, theta_est)

hist(new_observations, prob=TRUE)
xaxis <- 300:450
pred_dist <- dnbinom(xaxis, alpha+sum, 1 - 1/(beta+n+1))
# ^ since beta in the Gamma distribution is given as the rate parameter, the
#   corresponding probability parameter in the Negative Binomial needs to be
#   1 - (...) to make it a probability
#
# if we just give it as a scale parameter in Gamma, all is well
lines(xaxis, pred_dist, col = "red", lwd=3)



######## Question 2: Tuberculosis

# a)

samples <- read.table("tuberculosis.txt")
samples

# Marginal posterior histogram of each parameter
par(mfrow = c(2,2))
hist(samples$V1, main="Transmission rate (alpha)", xlab="Values")
hist(samples$V2, main="Transmission cost due to resistance (c)", xlab="Values")
hist(samples$V3, main="Mutation rate of drug resistance (p)", xlab="Values")
hist(samples$V4, main="Marker mutation rate (mu)", xlab="Values")
par(mfrow = c(1,1))

# Bivariate distributions
par(mfrow = c(2,3))
plot(x=samples$V1, y=samples$V2, xlab="alpha", ylab="c")
plot(x=samples$V1, y=samples$V3, xlab="alpha", ylab="p")
plot(x=samples$V1, y=samples$V4, xlab="alpha", ylab="mu")
plot(x=samples$V2, y=samples$V3, xlab="c", ylab="p")
plot(x=samples$V2, y=samples$V4, xlab="c", ylab="mu")
plot(x=samples$V3, y=samples$V4, xlab="p", ylab="mu")
par(mfrow = c(1,1))

# b)

# Function to compute the relative fitness of the drug-resistant strains,
# dependent on c and p and assuming all other variables are constant.
rel_fitness <- function(c, p) {
    delta <- 0.52
    tau <- 0.52
    eps_S <- 0.52
    eps_R <- 0.202
    return (1 - c) * (1/tau + 1/(delta + eps_R)) / (1/tau + 1/(delta + eps_S + p))
}

rel_fitness_values <- numeric(0)
for (i in 1:nrow(samples)) {
    c <- samples[i,]$V2
    p <- samples[i,]$V3
    rel_fitness_values <- c(rel_fitness_values, rel_fitness(c, p))
}

hist(rel_fitness_values, main="Posterior distribution of relative fitness")
print(sprintf("Mean relative fitness: %.4f", mean(rel_fitness_values)))
print(sprintf("Probability that relative fitness < 1: %.4f",
              sum(rel_fitness_values < 1) / length(rel_fitness_values)))

# c)

central <- quantile(rel_fitness_values, c(0.025, 0.975))

# The central 95% credible interval (0.893442, 1.837833) is very wide
# considering the small magnitude of the individual values. According to this
# credible interval, the relative fitness could be less than or greater than 1.

# d)

increment <- 0.05/250
shortest_interval_width <- 100
shortest_interval <- NULL

for (i in 0:250) {
    lower_quantile <- i * increment
    upper_quantile <- lower_quantile + 0.95
    credible_int <- quantile(rel_fitness_values, c(lower_quantile, upper_quantile))

    if (credible_int[2] - credible_int[1] < shortest_interval_width) {
        shortest_interval <- credible_int
        shortest_interval_width <- credible_int[2] - credible_int[1]
    }
}

shortest_interval

# The shortest credible interval is slightly shorter (0.937 compared to 0.944)
# than the central credible interval. The shortest interval is to the left of
# the central interval.



######## Question 3: Buffon's Needle

# a)
#
# bufffon <-function(n, d=1, ell=1, make.plot=TRUE) {
#   x <- runif(n, 0, d)
#   theta <- runif(n, 0, pi)
#   L <- sin(theta)
#   inters <- (ell*L) > x   # a vector of bools; whether each needle intersected or not
#   # ...
# }
#
