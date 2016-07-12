# load seroprevalence data and metadata
x = read.csv('../../data/seroprev_metadata.csv')
x$gdp_pcppp2005 = log(x$gdp_pcppp2005)
x$aedes = pmax(x$aegypti,x$albopictus)
x$avg_temp = rowMeans(as.matrix(x[,8:19]))
x$seroprev[x$seroprev==0] = 1e-3
x$seroprobit = qnorm(x$seroprev)


# compute an average model with average mosquito data
x$aegypti = pmax(x$aegypti,x$albopictus)
lm.est = lm(
  seroprobit ~
    (aegypti + avg_temp + gdp_pcppp2005) ^ 2 +
    I(aegypti^2) + I(avg_temp^2) + I(gdp_pcppp2005^2), data = x)
lm.step = step(lm.est)
save(lm.est,lm.step,file='../../generated/functions_stat_mean.RData')


# load replicate aegypti occurrence probabilities
aegypti.maps.pts.data = as.matrix(read.csv('../../data/aegypti_pt_vals_all.csv')[,-(1:4)])


# fit statistical models for 1,000 different samples with replacement from seroprevalence data
reps = 1000
lm.rand = list()
rep.master = cbind(
  1:reps,
  sample(ncol(aegypti.maps.pts.data),reps,replace=T))
for(rr in 1 : reps){
  x$aegypti = pmax(aegypti.maps.pts.data[,rep.master[rr,2]],x$albopictus)
  lm.est = lm(
    seroprobit ~
      (aegypti + avg_temp + gdp_pcppp2005) ^ 2 +
      I(aegypti^2) + I(avg_temp^2) + I(gdp_pcppp2005^2), data = x)
  
  lm.step = step(lm.est)
  lm.step$coefficients = rmvn(1,coef(lm.step),vcov(lm.step))
  
  lm.rand[[rr]] = lm.step
}

# save data
save(lm.rand,rep.master,reps,file='../../generated/functions_stat_random_draws.RData')
