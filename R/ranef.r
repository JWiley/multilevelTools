if (FALSE) {
library(knitr)
library(data.table)
library(brms)
library(cmdstanr)
library(mice)
library(miceadds)
library(ggplot2)
library(bayesplot)

dmixed <- withr::with_seed(
  seed = 12345, code = {
    nGroups <- 100
    nObs <- 20
    theta.location <- matrix(rnorm(nGroups * 2), nrow = nGroups, ncol = 2)
    theta.location[, 1] <- theta.location[, 1] - mean(theta.location[, 1])
    theta.location[, 2] <- theta.location[, 2] - mean(theta.location[, 2])
    theta.location[, 1] <- theta.location[, 1] / sd(theta.location[, 1])
    theta.location[, 2] <- theta.location[, 2] / sd(theta.location[, 2])
    theta.location <- theta.location %*% chol(matrix(c(1.5, -.25, -.25, .5^2), 2))
    theta.location[, 1] <- theta.location[, 1] - 2.5
    theta.location[, 2] <- theta.location[, 2] + 1
    dmixed <- data.table(
      x = rep(rep(0:1, each = nObs / 2), times = nGroups))
    dmixed[, ID := rep(seq_len(nGroups), each = nObs)]

    for (i in seq_len(nGroups)) {
      dmixed[ID == i, y := rnorm(
        n = nObs,
        mean = theta.location[i, 1] + theta.location[i, 2] * x,
        sd = exp(1 + theta.location[i, 1] + theta.location[i, 2] * x))
        ]
    }
    copy(dmixed)
  })

ls.me <- brm(bf(
  y ~ 1 + x + (1 + x | p | ID),
  sigma ~ 1 + x + (1 + x | p | ID)),
  family = "gaussian",
  data = dmixed, seed = 1234,
  silent = 2, refresh = 0, iter = 4000, warmup = 1000, thin = 3,
  chains = 4L, cores = 4L, backend = "cmdstanr")

.re.data <- function(x, i) {
  xw <- as.data.table(t(x[, , i]))
  xw[, ID := dimnames(x)[[2]]]
  xlong <- melt(xw, id.vars = "ID", 
    value.name = dimnames(x)[[3]][i], 
    variable.name = ".imp")
  xlong[, .imp := as.integer(.imp)]
  return(xlong)
}
re.data <- function(x) {
    for (i in seq_along(dimnames(x)[[3]])) {
      if (i == 1) {
        out <- .re.data(x, i)
      } else {
        tmp <- .re.data(x, i)
        out <- merge(out, tmp, by = c("ID", ".imp"))
      }
    }
    return(out)
}

x <- ranef(ls.me, summary = FALSE)$ID

xlong <- re.data(x)


ggplot(xlong, aes(x = Intercept, y = x)) +
  geom_hex(show.legend = FALSE) +
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "white", linewidth = 2) + 
  scale_fill_continuous(type = "viridis") + 
  theme_minimal() + 
  xlab("Random Intercept") + ylab("Random 'x' Slope")

i <- 3
names(x[1, , i])
xw <- as.data.table(t(x[, , i]))
xw[, ID := dimnames(x)[[2]]]
xlong <- melt(xw, id.vars = "ID", value.name = dimnames(x)[[3]][i], variable.name = ".imp")

dmixed2 <- dmixed[, .(MX = sample(x, 1)), by = ID]
dmixed2[, ID := as.character(ID)]

xlong2 <- merge(xlong, dmixed2, by = "ID")
xlong2[, .imp := as.integer(.imp) - 1]

xlong3 <- as.mids(xlong2, .imp = ".imp", .id = "ID")

micombine.cor(xlong3, variables = 1:2)
xlong2[, .(r = cor(Intercept, MX)), by = .imp][, atanh(mean(tanh(r)))]

dim(t(x[, , 1]))



i <- 4
yw <- as.data.table(t(x[, , i]))
yw[, ID := dimnames(x)[[2]]]
ylong <- melt(yw, id.vars = "ID", value.name = dimnames(x)[[3]][i], variable.name = ".imp")

ylong[, .imp := as.integer(.imp) - 1]

ylong2 <- merge(ylong, xlong2, by = c("ID", ".imp"))

ggplot(ylong2, aes(x = sigma_x, y = sigma_Intercept)) +
  geom_hex() +
  stat_smooth(method = "lm", se = FALSE, colour = "white", linewidth = 2) + 
  theme_minimal()

ylong2[, .(r = cor(sigma_x, sigma_Intercept))]

ylong2[, .(r = cor(sigma_x, sigma_Intercept)), by = ID][, tanh(mean(atanh(r)))]

ylong2[, .(r = cor(x, Intercept)), by = ID][, mean(r)]


summary(ls.me)

mcmc_hex(ls.me, pars = c("sd_ID__Intercept", "sd_ID__x")) + 
  xlab("Random Intercept") + 
  ylab("Random Slope")
}