# --------------------------------------------------- # 
# bootstrap confidence intervals
# determine the 95% confidence interval for the population mean of seeded clouds
# --------------------------------------------------- # 
clouds <- read.table(file="clouds.txt", header=TRUE)
T_val <- mean(clouds$seeded)

B <- 1000
T_vals <- numeric(B)
for (i in 1:B) {
    # generate data sample X*_i, ..., X*_N by sampling N values from the original dataset X_i, ..., X_N with replacement
    sampled_seeded <- sample(clouds$seeded, replace=TRUE)
    # compute the test statistic T*_i = T(X*_i, ..., X*_N) for the sample
    T_vals[i] <- mean(sampled_seeded) # compute Ti*
}

# compute bootstrap confidence interval with confidence 1 - 2α
confidence <- 0.95
alpha <- (1 - confidence) / 2
# sum(T_vals < 2 * T_val - quantile(T_vals, alpha))
c(2 * T_val - quantile(T_vals, 1 - alpha), 2 * T_val - quantile(T_vals, alpha))

# --------------------------------------------------- # 
# bootstrap tests
# H_0: the data are a random sample from a population with some distribution
# e.g.: H_0: X_i ~ exp(1), i.i.d. i = 1, ..., N
# --------------------------------------------------- # 
data <- rexp(1000, 1)
# hist(data, prob=TRUE)
hist(data, prob=TRUE, ylim=c(0, 0.7))
x <- seq(0, max(data), length=1000)
lines(x, dexp(x), type="l", col="blue", lwd=2)

T_val <- max(data)
N <- length(data)
B <- 1000
T_vals = numeric(B)

for (i in 1:B) {
    # generate data sample according to H_0
    sampled_data <- rexp(N, 1)
    # compute the test statistic T*_i = T(X*_i, ..., X*_N) for the sample
    T_vals[i] <- max(sampled_data)
}

# hist(T_vals, prob=TRUE)
dens_max_exp = function(x, n) n * exp(-x) * (1 - exp(-x))^(n-1)

hist(T_vals, prob=TRUE, ylim=c(0, 0.4), main="histogram of tstar & true density curve of T")
lines(rep(T_val, 2), seq(0, 2 * dens_max_exp(T_val, N), length=2), type="l", col="red", lwd=3)
axis(1, T_val, expression(paste("t") ) )
u <- seq(0, max(T_vals), length=1000)
lines(u, dens_max_exp(u, N), type="l", col="blue")

# compare the T-value of the original data to the T*-values and determine a p-value
p_val_left <- sum(T_vals < T_val) / B
p_val_right <- sum(T_vals > T_val) / B
p_val <- 2 * min(p_val_left, p_val_right)

# --------------------------------------------------- # 
# one sample, normal, t-test
# assumes the data are a random sample from a normal population
# H_0: the mean of this population is μ_0
# --------------------------------------------------- # 
mu <- 0.2
x <- rnorm(50, mu, 1) # creating artificial data
par(mfrow=c(1, 3))
hist(x)
boxplot(x)
qqnorm(x)

t.test(x)

# --------------------------------------------------- # 
# one sample, not normal, sign test
# assumes the data are a random sample from a population with a certain median m
# H_0: m = m_0 that the median of this population is m_0
# T: #(X_i > m_0) which has the Bin(N, 0.5) distribution under H_0
# --------------------------------------------------- # 
exam_results <- c(
    3.7, 5.2, 6.9, 7.2, 6.4, 9.3, 4.3, 8.4, 6.5, 8.1, 7.3, 6.1, 5.8
)
binom.test(sum(exam_results > 6), length(exam_results), p=0.5)

# --------------------------------------------------- # 
# one sample, not normal, wilcoxon test
# assumes the data are a random sample from a symmetric population with a certain median m
# H_0: the median of this population is m_0
# --------------------------------------------------- # 
# exam_results - 6
# rank(abs(exam_results - 6))
# rank(abs(exam_results - 6))[exam_results - 6 > 0]
# sum(rank(abs(exam_results - 6))[exam_results - 6 > 0])
wilcox.test(exam_results, mu=6)

# --------------------------------------------------- # 
# two paired samples, paired t-test
# H_0: the mean of the population is 0
# --------------------------------------------------- # 
ashina <- read.table("ashina.txt", header=TRUE)
plot(vas.active ~ vas.plac, pch=grp, col=grp, data=ashina); abline(0, 1)
boxplot(ashina$vas.active, ashina$vas.plac)
boxplot(ashina$vas.active - ashina$vas.plac)

t.test(ashina$vas.active, ashina$vas.plac, paired=TRUE)
t.test(ashina$vas.active - ashina$vas.plac)

# --------------------------------------------------- # 
# two paired samples, pearson's correlation test
# assumes normality of the both Xi and Yi
# H_0: the correlation between the two popluation is 0
# --------------------------------------------------- # 
peruvians <- read.table("peruvians.txt", header=TRUE)
plot(systolic ~ weight, data=peruvians)

cor.test(peruvians$systolic, peruvians$weight)

# --------------------------------------------------- # 
# two paired samples, pearson's correlation test
# H_0: the rank correlation between the two popluation is 0
# --------------------------------------------------- # 
qqnorm(peruvians$weight, main="Q-Q Plot weight")
qqnorm(peruvians$systolic, main="Q-Q Plot systolic")

cor.test(peruvians$systolic, peruvians$weight, method="spearman")

# --------------------------------------------------- # 
# two paired samples, permutation tests
# H_0: no difference between the distribution of X_i and that of Y_i within samples
# --------------------------------------------------- # 
boxplot(ashina$vas.active, ashina$vas.plac, names=c("active","placebo"))
plot(ashina$vas.active, ashina$vas.plac); abline(0, 1)

T_val <- mean(ashina$vas.active - ashina$vas.plac)
B <- 1000
T_vals <- numeric(B)
for (i in 1:B) {
    permutated_ashina <- t(apply(cbind(ashina$vas.active, ashina$vas.plac), 1, sample))
    T_vals[i] <- mean(permutated_ashina[, 1] - permutated_ashina[, 2])
}

hist(T_vals)
p_val_left <- sum(T_vals < T_val) / B
p_val_right <- sum(T_vals > T_val) / B
p_val <- 2 * min(p_val_left, p_val_right)

# --------------------------------------------------- # 
# two independent samples, two samples t-test
# H_0: means of the populations are the same
# --------------------------------------------------- # 
light1 <- scan("light1.txt")
light2 <- scan("light2.txt")
par(mfrow=c(1, 3))
hist(light1); hist(light2); boxplot(light1, light2)

# --------------------------------------------------- # 
# two independent samples, Mann-Whitney test
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# two independent samples, Kolmogorov-Smirnov test
# --------------------------------------------------- # 
