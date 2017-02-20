data("longley")
g <- lm(Employed ~ GNP + Population, data = longley)
summary(g, cor = T)

cor(g$res[-1], g$res[-16])

x <- model.matrix(g)
Sigma <- diag(16)
Sigma <- 0.31041 ^ abs(row(Sigma) - col(Sigma))
Sigi <- solve(Sigma)
xtxi <- solve(t(x) %*% Sigi %*% x)
beta <- xtxi %*% t(x) %*% Sigi %*% longley$Employed
beta

res <- longley$Employed - x %*% beta
sig <- sqrt(sum(res^2)/g$df)
sqrt(diag(xtxi)) * sig

sm <- chol(Sigma)
smi <- solve(t(sm))
sx <- smi %*% x
sy <- smi %*% longley$Employed
lm(sy ~ sx-1)$coef

# initial estimate of residual correlation was 0.31
    # now seeing what it is
cor(res[-1], res[-16])

# This is then iterated until convergence

rho_hat <- cor(g$res[-1], g$res[-16])
iter <- 1
repeat{
  rho_n <- rho_hat
  
  x <- model.matrix(g)
  Sigma <- diag(16)
  Sigma <- rho_n ^ abs(row(Sigma) - col(Sigma))
  Sigi <- solve(Sigma)
  xtxi <- solve(t(x) %*% Sigi %*% x)
  beta <- xtxi %*% t(x) %*% Sigi %*% longley$Employed
  
  res <- longley$Employed - x %*% beta
  # sig <- sqrt(sum(res^2)/g$df)
  # sqrt(diag(xtxi)) * sig
  # 
  # sm <- chol(Sigma)
  # smi <- solve(t(sm))
  # sx <- smi %*% x
  # sy <- smi %*% longley$Employed
  # lm(sy ~ sx-1)$coef
  
  # initial estimate of residual correlation was 0.31
  # now seeing what it is
  rho_hat <- cor(res[-1], res[-16])
  
  if(abs(rho_hat - rho_n) < 1e-8) break
  iter <- iter + 1
}

c(iter, rho_hat, rho_n)


library(nlme)
g <- gls(Employed ~ GNP + Population,
           correlation=corAR1(form = ~ Year), data = longley)
summary(g)
intervals(g)


