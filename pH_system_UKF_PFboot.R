library(mvtnorm) # multivariate normal distribution
#library(psych) # trace
library(pracma) # eye
#library("future.apply")
#plan("multisession", workers = 4)

x1i <- 1.2e-3
x2i <- 2e-3
x3i <- 2.5e-3
Kx <- 1e-7
Kw <- 1e-14
V <- 2.5
qA <- function(t) { 1+0.06*sin(0.01*t) }
qB <- 0.265
theta <- function(t) { V/qA(t) }
x0 <- c(9.484e-4, 4.194e-4, 5.242e-4)
var_state <- 2e-11
var_obs <- 1e-4
a <- 1
d <- -Kw^2/Kx
T <- 2000

at <- function(x, t)
{
  return(c(x[1] + (x1i - x[1])/theta(t) - x[1]*qB/V,
           x[2] - x[2]/theta(t) + (x2i - x[2])*qB/V,
           x[3] - x[3]/theta(t) + (x3i - x[3])*qB/V))
}
At <- function(x, t)
{
  b <- Kw/Kx + x[3] + x[2] - x[1]
  c <- (x[2] - x[1] - Kx)*Kw/Kx
  p <- -b/(3*a)
  q <- p^3 + (b*c - 3*a*d)/(6*a^2)
  r <- c/(3*a)
  return(Re(-log10((q + sqrt(as.complex(q^2 + (r-p^2)^3)))^(1/3) + (q - sqrt(as.complex(q^2 + (r-p^2)^3)))^(1/3) + p)))
}
bt <- diag(rep(sqrt(var_state), length(x0)))
Bt <- sqrt(var_obs)

pH_system <- function(x0, at, At, bt, Bt, T)
{
  sys <- matrix(NA, nrow = T+1, ncol = 4)
  colnames(sys) <- c('x1', 'x2', 'x3', 'pH')
  sys[1,1:3] <- x0
  U <- cbind(rnorm(T+1), rnorm(T+1), rnorm(T+1))
  W <- rnorm(T+1)
  for(i in 1:(T+1)){
    sys[i,4] <- At(sys[i,1:3], t) + Bt*W[i]
    if(i != T+1){
      sys[i+1,1:3] <- at(sys[i,1:3], i+1) + bt%*%U[i+1, ]
    }
  }
  return(sys)
}

set.seed(1)
sys <- pH_system(x0, at, At, bt, Bt, T)

# State Vector Plot
plot(sys[ ,'x1']*10^3, type = 'l', ylim = c(0, max(sys[ ,1:3])*10^3), col = 'blue', 
     xlab = 'Time, min', ylab = 'x*10^3, mol/L', main = 'State Vector')
lines(sys[ ,'x2']*10^3, col = 'red3')
lines(sys[ ,'x3']*10^3, col = 'magenta')
abline(h = seq(0, 1, by = 0.2), col = 'darkgray', lty = 'dotted')
abline(v = seq(0, 2000, by = 500), col = 'darkgray', lty = 'dotted')
legend('bottomleft', c('x1', 'x2', 'x3'), col = c('blue', 'red3', 'magenta'), lwd = 2)

# pH Plot
plot(sys[,'pH'], type = 'l', ylim = c(3.5, 6.5), xlab = 'Time, min', ylab = 'pH, -',
     main = 'Observations')
abline(h = seq(3.5, 6.5, by = 0.5), col = 'darkgray', lty = 'dotted')
abline(v = seq(0, 2000, by = 500), col = 'darkgray', lty = 'dotted')

# qA Plot
plot(0:2000, qA(0:2000), type = 'l', ylim = c(0.9, 1.1), lwd = 2, 
     xlab = 'Time, min', ylab = 'qA, L/min', main = 'qA, Acid Inlet Flow Rate')
abline(h = seq(0.9, 1.1, by = 0.05), col = 'darkgray', lty = 'dotted')
abline(v = seq(0, 2000, by = 500), col = 'darkgray', lty = 'dotted')

discrete <- function(Value, Prob)
{
  r <- runif(1)
  for(i in 1:length(Prob)){
    if(i == 1){
      if(r < Prob[1])
        return(Value[1, ])
    }else{
      if(r > sum(Prob[1:(i-1)]) & r < sum(Prob[1:i]))
        return(Value[i, ])
    }
  }
}

#Trivial Estimator
trivialEstimator <- function(at, x0, T)
{
  res <- matrix(NA, nrow = T+1, ncol = length(x0))
  res[1, ] <- x0
  for(i in 2:(T+1)){
    res[i, ] <- at(res[i-1, ], i)
  }
  return(res)
}

# Unscented Kalman Filter
sigmaPoints <- function(k, x, w0)
{
  N <- length(x)
  res <- matrix(NA, nrow = 2*N + 1, ncol = N + 1)
  res[2*N+1, ] <- c(x, w0)
  res[1:(2*N),N+1] <- rep((1-w0)/(2*N), 2*N)
  ch <- t(chol(k*N/(1-w0)))
  for(i in 1:N){
    res[i,1:N] <- x + ch[ ,i]
    res[N+i,1:N] <- x - ch[ ,i]
  }
  return(res)
}

UKF <- function(at, bt, At, Bt, m0, Y, w0)
{
  T <- length(Y)
  N <- length(m0)
  X_us <- matrix(NA, nrow = T, ncol = N)
  X_us[1, ] <- m0
  k <- eye(N)*var_state
  for(i in 1:(T-1)){
    sp <- sigmaPoints(k, X_us[i,], w0)
    Q <- dim(sp)[1]
    X_moon <- apply(sapply(1:Q, function(j) sp[j,N+1]*at(sp[j,1:N],i+1)), 1, FUN = sum)
    k_moon <- zeros(N)
    for(j in 1:Q)
      k_moon <- k_moon + sp[j,N+1]*(at(sp[j,1:N],i+1) - X_moon)%*%t(at(sp[j,1:N],i+1) - X_moon) 
    k_moon <- k_moon + bt%*%t(bt)
    sp <- sigmaPoints(k_moon, X_moon, w0)
    Y_moon <- sum(sapply(1:Q, function(j) sp[j,N+1]*At(sp[j,1:N],i+1)))
    kappa_moon <- sum(sapply(1:Q, function(j) sp[j,N+1]*(At(sp[j,1:N],i+1) - Y_moon)^2))
    kappa_moon <- kappa_moon + Bt^2
    mu_moon <- apply(sapply(1:Q, function(j) sp[j,N+1]*(sp[j,1:N] - X_moon)*(At(sp[j,1:N],i+1) - Y_moon)), 1, FUN = sum)
    X_us[i+1, ] <- X_moon + mu_moon*kappa_moon^(-1)*(Y[i+1] - Y_moon)
    k <- k_moon - kappa_moon^(-1)*mu_moon%*%t(mu_moon)
  }
  return(X_us)
}

#Bootstrap Particle Filter
bootstrapFilter <- function(at, bt, At, Bt, m0, Y, N)
{
  T <- length(Y)
  X_hat <- matrix(NA, nrow = T, ncol = 3)
  X <- matrix(rep(m0, N), nrow = N, byrow = TRUE)
  w <- rep(1/N, N)
  X_hat[1, ] <- m0
  for(i in 2:T){
    X <- t(sapply(1:N, function(j) rmvnorm(1, mean = at(X[j, ],j), sigma = bt%*%bt)))
    w_prev <- w
    w <- sapply(1:N, function(j) (dnorm(Y[i], mean = At(X[j, ],j), sd = Bt)/Bt)*w_prev[j])
    if(sum(w^2) == 0){
      w <- rep(1/N, N)
    }
    w <- w/sum(w)
    if(1/sum(w^2) < N/10){
      X <- t(sapply(1:N, function(j) discrete(X,w)))
      w <- rep(1/N, N)
    }
    X_hat[i, ] <- apply(sapply(1:N, function(j) w[j]*X[j, ]), 1, FUN = sum)
  }
  return(X_hat)
}

m0 <- x0
N <- 100
set.seed(1)
ukf <- UKF(at, bt, At, Bt, m0, Y = sys[,'pH'], w0 = 1/4)
bf <- bootstrapFilter(at, bt, At, Bt, m0, Y = sys[ ,'pH'], N)
triv <- trivialEstimator(at, x0, T)

# Results
par(mfrow = c(3,1))
name = 'Точные значения компонент и их оценки'
for(i in 1:3){
  plot(sys[ ,i]*1e3, type = 'l', lwd = 2, xlab = 'Time, min', ylab = paste0('x',i,'*10^3, mol/L'), main = name)
  abline(h = seq(min(sys[ ,i]*1e3), max(sys[ ,i]*1e3), length.out = 6), col = 'darkgray', lty = 'dotted')
  abline(v = seq(0, 2000, by = 250), col = 'darkgray', lty = 'dotted')
  name = ''
  lines(triv[ ,i]*1e3, lwd = 2, col = 'red3')
  lines(bf[ ,i]*1e3, lwd = 2, col = 'orange')
  lines(ukf[ ,i]*1e3, lwd = 2, col = 'darkviolet')
  #legend('bottomleft', c('State', 'UKF', 'PF(boot)', 'Trivial'), col = c('black', 'darkviolet', 'orange', 'red3'),
  #       lwd = 2, horiz = TRUE, bty = 'n')
}
par(mfrow = c(1,1))

calculateError <- function(at, bt, At, Bt, x0, T, w0, N)
{
  sys <- pH_system(x0, at, At, bt, Bt, T)
  ukf <- UKF(at, bt, At, Bt, m0 = x0, Y = sys[,'pH'], w0 = 1/4)
  bf <- bootstrapFilter(at, bt, At, Bt, m0 = x0, Y = sys[ ,'pH'], N)
  return(cbind(sys[ ,1:3] - ukf, sys[ ,1:3] - bf, sys[ ,1:3]))
}

w0 <- 1/4
N <- 100
n_samp <- 100
x1_u <- matrix(NA, nrow = T+1, ncol = n_samp)
x2_u <- matrix(NA, nrow = T+1, ncol = n_samp)
x3_u <- matrix(NA, nrow = T+1, ncol = n_samp)
x1_b <- matrix(NA, nrow = T+1, ncol = n_samp)
x2_b <- matrix(NA, nrow = T+1, ncol = n_samp)
x3_b <- matrix(NA, nrow = T+1, ncol = n_samp)
x1 <- matrix(NA, nrow = T+1, ncol = n_samp)
x2 <- matrix(NA, nrow = T+1, ncol = n_samp)
x3 <- matrix(NA, nrow = T+1, ncol = n_samp)
makeTables <- function(i)
{
  err <- calculateError(at, bt, At, Bt, x0, T, w0, N)
  x1_u[ ,i] <- err[ ,1]
  x2_u[ ,i] <- err[ ,2]
  x3_u[ ,i] <- err[ ,3]
  x1_b[ ,i] <- err[ ,4]
  x2_b[ ,i] <- err[ ,5]
  x3_b[ ,i] <- err[ ,6]
  x1[ ,i] <- err[ ,7]
  x2[ ,i] <- err[ ,8]
  x3[ ,i] <- err[ ,9]
  print(i)
}
start_time <- Sys.time()
sapply(1:n_samp, FUN = makeTables)
timediff <- difftime(Sys.time(),start_time)
cat('Расчёт занял: ', timediff, units(timediff))

start_time <- Sys.time()
for(i in 1:n_samp){
  err <- calculateError(at, bt, At, Bt, x0, T, w0, N)
  x1_u[ ,i] <- err[ ,1]
  x2_u[ ,i] <- err[ ,2]
  x3_u[ ,i] <- err[ ,3]
  x1_b[ ,i] <- err[ ,4]
  x2_b[ ,i] <- err[ ,5]
  x3_b[ ,i] <- err[ ,6]
  x1[ ,i] <- err[ ,7]
  x2[ ,i] <- err[ ,8]
  x3[ ,i] <- err[ ,9]
  print(i)
}
timediff <- difftime(Sys.time(),start_time)
cat('Расчёт занял: ', timediff, units(timediff))

comp_rmse <- cbind(apply(x1, 1, FUN = sd), apply(x2, 1, FUN = sd), apply(x3, 1, FUN = sd))
ukf_rmse <- cbind(apply(x1_u, 1, FUN = sd), apply(x2_u, 1, FUN = sd), apply(x3_u, 1, FUN = sd))
bf_rmse <- cbind(apply(x1_b, 1, FUN = sd), apply(x2_b, 1, FUN = sd), apply(x3_b, 1, FUN = sd))
colnames(comp_rmse) <- c('x1', 'x2', 'x3')
colnames(ukf_rmse) <- c('x1', 'x2', 'x3')
colnames(bf_rmse) <- c('x1', 'x2', 'x3')

comp_rmse1 <- read.csv('comp_rmse.csv')
bf_rmse1 <- read.csv('bf_rmse.csv')
ukf_rmse1 <- read.csv('ukf_rmse.csv')

par(mfrow = c(3,1))
name = 'СКО компонент и ошибок оценки'
for(i in 1:3){
  plot(bf_rmse1[ ,i], type = 'l', col = 'orange', xlab = 'Time, min', ylab = paste0('RMSE x', i, ', mol/L'), main = name)
  abline(h = seq(0, max(bf_rmse1[ ,i]), length.out = 6), col = 'darkgray', lty = 'dotted')
  abline(v = seq(0, 2000, by = 250), col = 'darkgray', lty = 'dotted')
  name = ''
  lines(ukf_rmse1[ ,i], col = 'darkviolet')
  lines(comp_rmse1[ ,i], col = 'darkgreen')
}
par(mfrow = c(1,1))

legend('bottomleft', c('x1', 'x2', 'x3','System State', 'Trivial Estimator'), 
       col = c('blue', 'red3', 'magenta', 'black', 'black'), lwd = 2, 
       x.intersp = 0.3, text.width = 250, lty = c(1,1,1,1,2), ncol = 2)
bad_index <- c(370:581, 970:1240, 1592:1859)
ukf_cov <- cov(sys[20:2001 ,1:3] - ukf[20:2001, ])
bf_cov <- cov(sys[-bad_index,1:3] - bf[-bad_index, ])
triv_cov <- cov(sys[20:2001,1:3] - triv[20:2001, ])
print(paste('UKF MSE', tr(ukf_cov), '| Trivial Estimator MSE', tr(triv_cov)))

write.csv(comp_rmse, 'comp_rmse.csv', row.names = FALSE)
write.csv(ukf_rmse, 'ukf_rmse.csv', row.names = FALSE)
write.csv(bf_rmse, 'bf_rmse.csv', row.names = FALSE)


