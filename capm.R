# Bayesian Capital Asset Pricing Model (CAPM)

# Libraries ---------------------------------------------------------------

library(rjags)
library(corrplot)

# Data Collection ---------------------------------------------------------

dat = read.csv("./datasets/Stock_Bond.csv")
cols = c("GM_AC","F_AC","UTX_AC","CAT_AC","MRK_AC","PFE_AC",
         "IBM_AC","MSFT_AC","C_AC","XOM_AC","S.P_AC")
price_ac = dat[4543:4794, cols]
tickers = c("GM","F","UTX","CAT","MRK","PFE",
            "IBM","MSFT","C","XOM","SP500")
colnames(price_ac) = tickers
str(price_ac)

# Data Exploration --------------------------------------------------------

n = dim(price_ac)[1]; m = dim(price_ac)[2]-1
r = price_ac[-1,]/price_ac[-n,] - 1
Cor = cor(r)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d", tl.cex=0.7)
corrplot(Cor, type="lower", method="number", col="black",
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

# Model Postulation -------------------------------------------------------

# OLS Regression
betas_ols = numeric(m)
for(i in 1:m)
{
  betas_ols[i] = coefficients(lm(r[,i] ~ r[,m+1]))[2]
}
names(betas_ols) = tickers[1:m]
round(betas_ols,4)

# Bayesian Hierarchical Linear Regression
bayesian_capm = "model
{
  for(t in 1:N)
  {
    for(j in 1:m)
    {
      R[t,j] ~ dnorm(beta[j]*mkt[t], tau_e[j])
    }
  }
  for(j in 1:m)
  {
    beta[j] ~ dnorm(mean_b, tau_b)
    tau_e[j] ~ dgamma(0.1, 0.001)
  }
  mean_b ~ dnorm(1, 1e-6)
  tau_b ~ dunif(1, 100)
}"

set.seed(42)

data = list(R=r[,1:m], N=dim(r)[1], mkt=r$SP500, m=m)

params = c("beta")
inits_capm = function() { list(beta=rep(1,m)) }

jags_capm = jags.model(
  textConnection(bayesian_capm),
  inits=inits_capm, data=data, n.chains=3
)
update(jags_capm, 1e3)

capm_sim = coda.samples(
  model=jags_capm,
  variable.names=params,
  n.iter=5e3
)

capm_csim = as.mcmc(do.call(rbind, capm_sim))

# Model Checking

plot(capm_sim)
gelman.diag(capm_sim)
autocorr.diag(capm_sim)
effectiveSize(capm_sim)

betas = colMeans(capm_csim)
rhat = matrix(r$SP500, ncol=1) %*% betas
resid = r - rhat
plot(rhat[,10], resid[,10], main=colnames(r)[10])

# Model Usage

res_betas = t(round(summary(capm_sim)[[1]],4)[,1:2])
colnames(res_betas) = tickers[1:m]
res_betas

res_ret = colMeans(rhat)
names(res_ret) = tickers[1:m]
round(res_ret * 252, 5)
