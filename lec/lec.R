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
# H0: the data are a random sample from a population with some distribution
# e.g.: H0: X_i ~ exp(1), i.i.d. i = 1, ..., N
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
    # generate data sample according to H0
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
# H0: the mean of this population is μ_0
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
# H0: m = m_0 that the median of this population is m_0
# T: #(X_i > m_0) which has the Bin(N, 0.5) distribution under H0
# --------------------------------------------------- # 
exam_results <- c(
    3.7, 5.2, 6.9, 7.2, 6.4, 9.3, 4.3, 8.4, 6.5, 8.1, 7.3, 6.1, 5.8
)
binom.test(sum(exam_results > 6), length(exam_results), p=0.5)
# --------------------------------------------------- # 
# one sample, not normal, wilcoxon test
# assumes the data are a random sample from a symmetric population with a certain median m
# H0: the median of this population is m_0
# --------------------------------------------------- # 
# exam_results - 6
# rank(abs(exam_results - 6))
# rank(abs(exam_results - 6))[exam_results - 6 > 0]
# sum(rank(abs(exam_results - 6))[exam_results - 6 > 0])
wilcox.test(exam_results, mu=6)
# --------------------------------------------------- # 
# two paired samples, paired t-test
# H0: the mean of the population is 0
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
# H0: the correlation between the two popluation is 0
# --------------------------------------------------- # 
peruvians <- read.table("peruvians.txt", header=TRUE)
plot(systolic ~ weight, data=peruvians)

cor.test(peruvians$systolic, peruvians$weight)
# --------------------------------------------------- # 
# two paired samples, pearson's correlation test
# H0: the rank correlation between the two popluation is 0
# --------------------------------------------------- # 
qqnorm(peruvians$weight, main="Q-Q Plot weight")
qqnorm(peruvians$systolic, main="Q-Q Plot systolic")

cor.test(peruvians$systolic, peruvians$weight, method="spearman")
# --------------------------------------------------- # 
# two paired samples, permutation tests
# H0: no difference between the distribution of X_i and that of Y_i within samples
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
# H0: means of the populations are the same
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

# --------------------------------------------------- # 
# one-way ANOVA
# H0: μ1 = ... = μl (or, α1 = ... = αl = 0, μi = μ + αi)
# --------------------------------------------------- # 
n_level_of_factor <- 4 # l
n_observations_per_group <- 5 # ni (N if balanced design)
rbind(
    sample(rep(1:n_level_of_factor, n_observations_per_group)) # assigned level
    rep(1:n_level_of_factor, n_observations_per_group), # experimental units
)

melon <- read.table("melon.txt", header=TRUE)
par(mfrow=c(1, 2))
boxplot(melon); stripchart(melon, vertical=TRUE)

melon <- data.frame(
    yield = as.vector(as.matrix(melon)),
    variety = factor(rep(1:4, each=6))
)
is.factor(melon$variety)
is.numeric(melon$variety)

melon_lm <- lm(yield~variety, data=melon)
anova(melon_lm)

summary(melon_lm)
fitted(melon_lm)

confint(melon_lm)

contrasts(melon$variety) <- contr.sum
melon_lm <- lm(yield~variety, data=melon)
summary(melon_lm)

par(mfrow=c(2, 2))
for (i in 1:4) 
    qqnorm(melon[, i])

par(mfrow=c(1, 1))
qqnorm(residuals(melonaov))
# --------------------------------------------------- # 
# Kruskal-Wallis test
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# permutation test
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# two-way ANOVA
# HAB:  γij = 0 for every (i, j) (no interactions between factor A and B)
# HA: αi = 0 for every i (no main effect of factor A)
# HA: βi = 0 for every j (no main effect of factor B)
# --------------------------------------------------- # 
n_level_of_factor_a <- 4 # I
n_level_of_factor_b <- 2 # J
n_observations_per_group <- 3 # N
rbind(
    rep(1:n_level_of_factor_a, each=n_observations_per_group * n_level_of_factor_b), # assigned level of factor a
    rep(1:n_level_of_factor_b, times=n_observations_per_group * n_level_of_factor_a), # assigned level of factor b
    sample(1:(n_observations_per_group * n_level_of_factor_a * n_level_of_factor_b)) # experimental units
)

pvc <- read.table(file="pvc.txt", header=TRUE)
par(mfrow=c(1, 2))
boxplot(psize~operator, data=pvc); boxplot(psize~resin, data=pvc)

with(pvc, {
    par(mfrow=c(1, 2))
    interaction.plot(operator, resin, psize)
    interaction.plot(resin, operator, psize)
})

pvc$operator=as.factor(pvc$operator)
pvc$resin=as.factor(pvc$resin)

pvc_lm <- lm(psize ~ operator * resin, data=pvc)
anova(pvc_lm)

summary(pvc_lm)

contrasts(pvc$operator) = contr.sum
contrasts(pvc$resin) = contr.sum
pvc_lm = lm(psize ~ operator * resin, data=pvc)
summary(pvc_lm)

par(mfrow=c(1, 2))
qqnorm(residuals(pvc_lm))
plot(fitted(pvc_lm), residuals(pvc_lm))

pvc <- read.table(file="pvc.txt", header=TRUE)
pvc$operator=as.factor(pvc$operator)
pvc$resin=as.factor(pvc$resin)
pvc_lm = lm(psize ~ operator + resin, data=pvc)
anova(pvc_lm)

composite <- read.table("composite.txt", head=TRUE)
anova(lm(strength ~ laser * tape, data=composite))
anova(lm(strength ~ laser + tape, data=composite))

# --------------------------------------------------- # 
# randomized block design
# H0: α1 = ... = αl = 0
# --------------------------------------------------- # 
n_level_of_factor_a <- 4 # I
n_blocks <- 5 # J
n_observations_per_group <- 3 # N
for (i in 1:n_blocks) 
    sample(1:(n_observations_per_group * n_level_of_factor_a))

penicillin <- read.table("penicillin.txt", head=TRUE)
xtabs(yield ~ treat + blend, data=penicillin)

par(mfrow=c(1, 2))
boxplot(yield ~ treat, data=penicillin)
boxplot(yield ~ blend, data=penicillin)

par(mfrow=c(1, 2))
with(pvc, {
    par(mfrow=c(1, 2))
    interaction.plot(treat, blend, yield)
    interaction.plot(blend, treat, yield)
})

penicillin_lm <- lm(yield ~ treat + blend, data=penicillin)
anova(penicillin_lm)
summary(penicillin_lm)
# --------------------------------------------------- # 
# repeated measures
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# Friedman test
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# incomplete block designs
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# random effects
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# crossover design
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# split-plot design
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# contingency tables
# --------------------------------------------------- # 
grades <- matrix(
    c(8, 15, 13, 14, 19, 15, 15, 4, 7, 3, 1, 4), 
    byrow=TRUE, 
    ncol=3, nrow=4,
    dimnames=list(
        c("A", "B", "C", "D-F"), 
        c("Psychology", "Biology", "Other")
    )
)
row_sum <- apply(grades, 1, sum)
col_sum <- apply(grades, 2, sum)
total_sum <- sum(grades)
expected <- row_sum %*% t(col_sum) / total_sum
round(expected, 0)
T <- sum((grades - expected)^2 / expected)
p_val <- 1 - pchisq(T, 6)

z <- chisq.test(grades, simulate.p.value=TRUE)
z$expected
z$observed
z$residuals

handed <- matrix(
    c(2780,3281,311,300), 
    nrow=2, ncol=2, 
    byrow=TRUE, 
    dimnames=list(
        c("right-handed","other"), 
        c("men","women")
    )
)
fisher.test(handed)
# --------------------------------------------------- # 
# simple linear regression
# H0: β1 = 0
# --------------------------------------------------- # 
sat <- read.table("sat.txt", header=TRUE)
sat_lm <- lm(total~expend, data=sat)
summary(sat_lm)

plot(total~expend, data=sat)
abline(sat_lm)

par(mfrow=c(1, 2))
qqnorm(residuals(sat_lm))
plot(fitted(sat_lm), residuals(sat_lm))
# --------------------------------------------------- # 
# multiple linear regression
# H0: βi = 0
# --------------------------------------------------- # 
sat <- read.table("sat.txt", header=TRUE)
plot(sat[c("expend", "takers", "total")])
par(mfrow=c(1, 3))
for (col in c("expend", "takers", "total"))
    hist(sat[col][, 1], main=col)

sat_lm <- lm(total ~ expend + takers, data=sat)
summary(sat_lm)
confint(sat_lm)

par(mfrow=c(1, 2))
qqnorm(residuals(sat_lm))
plot(fitted(sat_lm), residuals(sat_lm))

sat$takers2 <- sat$takers ^ 2
sat_lm <- lm(total ~ expend + takers + takers2, data=sat)
summary(sat_lm)

par(mfrow=c(1, 2))
qqnorm(residuals(sat_lm))
plot(fitted(sat_lm), residuals(sat_lm))
# --------------------------------------------------- # 
# strategies to choose the variables
# --------------------------------------------------- # 

# --------------------------------------------------- # 
# diagnostics in linear regression
# --------------------------------------------------- # 
bodyfat <- read.table("bodyfat.txt", header=TRUE)
pairs(bodyfat)

bodyfat_lm <- lm(Fat ~ Thigh, data=bodyfat)
plot(residuals(bodyfat_lm), bodyfat$Thigh)

x <- residuals(lm(Thigh ~ Midarm + Triceps, data=bodyfat))
y <- residuals(lm(Fat ~ Midarm + Triceps, data=bodyfat))
plot(
    x, y, 
    main="Added variable plot for Thigh", 
    xlab="residual of Thigh", 
    ylab="residual of Fat"
)

plot(residuals(bodyfat_lm), bodyfat$Triceps)
plot(residuals(bodyfat_lm), bodyfat$Midarm)

plot(residuals(bodyfat_lm), bodyfat$Fat)
plot(res(bodyfat_lm), fitted(bodyfat_lm))

qqnorm(residuals(bodyfat_lm))
# --------------------------------------------------- # 
# outliers and influence points
# --------------------------------------------------- # 
forbes <- read.table("forbes.txt", head=TRUE)
forbes_lm <- lm(X100xLog10.Pressure. ~ BoilingPoint, data=forbes)
round(residuals(forbes_lm), 2)
order(abs(residuals(forbes_lm)))

u11 <- rep(0, 16)
u11[11] <- 1
summary(lm(X100xLog10.Pressure. ~ BoilingPoint + u11, data=forbes))

huber <- read.table("huber.txt", head=TRUE)
huber_lm <- lm(y ~ x, data=huber)
round(cooks.distance(huber_lm), 2)
plot(1:6, cooks.distance(huber_lm), type="b")
# --------------------------------------------------- # 
# collinearity
# --------------------------------------------------- # 
round(cor(bodyfat), 2)
pairs(bodyfat)

bodyfat_lm <- lm(Fat ~ Thigh + Triceps + Midarm, data=bodyfat)
vif(bodyfat_lm)
bodyfat_lm2 <- lm(Fat ~ Triceps + Midarm, data=bodyfat)
vif(bodyfat_lm2)
bodyfat_lm2 <- lm(Fat ~ Thigh, data=bodyfat)
vif(bodyfat_lm2)
