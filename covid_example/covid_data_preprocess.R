library(mvtnorm)
library(tmvtnorm)
library(sde)
library(fChange)
library(invgamma)
library(reshape2)
library(ggplot2)
library(scpt)
library(splines)
library(matrixcalc)
library(stringr)
library(fda)

setwd("covid_data") # Replace with your directory containing the data
# the county center for all counties
load("countyCenter.RData") 
# A vector of county names
load("CountyName.RData")
# The daily newly increased covid-19 cases after filtering
# Each column represents a county
# "dates" and "ages" saves the corresponding date and age group for each row.
load("covidData.RData")
load("dates.RData")
load("ages.RData")

ns = 102 # total number of counties in IL
Nt = 95 # total number of time points
factor = 20 # the number to scale the functional data in order to get similar SNR as in the simulation

# Scale the spatial domain to 10X10
location = countyCenter[, c(2, 3)]
location[,1] = location[,1] - min(location[,1])
location[,2] = location[,2] - min(location[,2])
location[,1] = location[,1]/max(location[,1])*10
location[,2] = location[,2]/max(location[,2])*10
ds = as.matrix(dist(location))

################################################################################################
D = 7 # Number of Fourier basis functions
x0 = (0:Nt)/Nt # scale the time points into [0, 1]
nt = Nt - 1
n = ns * nt

setwd("covidresult") # Replace the directory where you would like to save the results
ovibasis = create.fourier.basis(rangeval = c(0, 1), nbasis = D)
psummary = data.frame(matrix(ncol = 8, nrow = 0))
colnames(psummary) = c("index", "County", "change_est", "CI_l",
                       "CI_u","day_bef", "day_aft", "pvalue")

first = TRUE
pic = TRUE
if(first){
  pvalue = change_est = trs = SNR = rep(0, ns)
  a_expect = b_expect = k_expect = a0_expect = b0_expect = beta_expect = beta_old = rep(0, ns)
  S = matrix(0, ncol = ns, nrow = Nt+1)
  confs = matrix(0, nrow = ns, ncol =2)
  set.seed(100)
  for(i in 1:ns){
    s = i
    mat.S = matrix(0, D, Nt)
    t = 0
    for(date in unique(dates)){
      index = which(dates==date)
      case = covidData[index, s]
      dage = ages[index]
      if(sum(case) > 0) {case = case/sum(case)}
      ovifourier.fd = smooth.basis(argvals = (1:length(case))/length(case), y = factor*case, fdParobj = ovibasis)$fd
      t = t + 1
      mat.S[, t] = ovifourier.fd$coefs
    }
    fdata_S = fd(mat.S, ovibasis)
    h = round(opt_bandwidth(fdata_S, "BT", "PR")$hat_h_opt)
    resultiid = change_FF(fdata_S,h=h)
    confint = Conf_int(fdata_S, h=h) #level 0.05
    confs[i, 1] = confint[1]/Nt
    confs[i, 2] = confint[3]/Nt
    pvalue[i] = resultiid$pvalue
    change_est[i] = resultiid$change
    fdata = center.fd(fdata_S) ; basis = fdata$basis; samp = fdata$coefs;N = ncol(samp)
    S[2,i] = sum((samp[, 1]- (1/N) * rowSums(samp[, 1:N]))^2)/N
    for (j in (2:N)) {
      S[j+1, i] = sum((rowSums(samp[, 1:j]) - (j/N) * rowSums(samp[, 1:N]))^2)/N
    }
    cat("County", s, "\n")
    result = change_FF(fdata_S, h=h, plot = F)
    if(pic){
      png(paste("CUSUM", s,".png", sep=""))
      plot(0:N, S[,i], xlab = "time", ylab = "CUSUM", type = "l", main = paste("County", s, "with p-value", result$pvalue))
      abline(v = result$change, col = "red")
      dev.off()
    }
    CI = try(Conf_int(fdata_S, h=h), silent=TRUE)
    tryres = try(CI[[1]])
    if(tryres == "Error in if (upper > N) { : missing value where TRUE/FALSE needed\n" ||
       tryres == "Error in if (lower < 0) { : missing value where TRUE/FALSE needed\n"){
      psummary = rbind(psummary, data.frame(index = s, County = CountyName[s]
                                            ,change_est = result$change, CI_l = NA, 
                                            CI_u = NA, day_bef = unique(dates)[result$change], 
                                            day_aft = unique(dates)[result$change + 1], 
                                            pvalue = result$pvalue))
      next
    } 
    psummary = rbind(psummary, data.frame(index = s, County = CountyName[s]
                                          ,change_est = result$change, CI_l = CI[[1]], 
                                          CI_u = CI[[3]], day_bef = unique(dates)[result$change], 
                                          day_aft = unique(dates)[result$change + 1], 
                                          pvalue = result$pvalue))

    k_star_guess = min(which(S[,i] == max(S[,i])))-1
    k_expect[i] = k_star_guess/Nt
    beta_y = max(S[,i])
    beta_x = (which.max(S[,i])-1)/Nt
    beta_old[i] = beta_y/(beta_x*(beta_x-1))
    
    LongRunC = LongRun(fdobj = fdata, h = h)
    lambda = LongRunC$e_val
    phi = LongRunC$e_fun$coefs
    dat.b = fdata_S[1:k_star_guess]
    dat.a = fdata_S[(k_star_guess + 1):N]
    mean.b = rowMeans(dat.b$coefs)
    mean.a = rowMeans(dat.a$coefs)
    delta_basis = mean.a - mean.b
    
    a_expect[i] = 2 * sum(lambda^2)
    for(m in 1:D){
      b_expect[i] = b_expect[i]+4*lambda[m]*(sum(phi[,m]*delta_basis))^2
    }
    a0_expect[i] = sum(lambda)
    b0_expect[i] = sum(delta_basis^2)
    
    trs[i] = sum(diag(LongRunC$covm))
    SNR[i] = k_expect[i]*(1-k_expect[i])*b0_expect[i]/trs[i]
    m = beta_old[i]*((k_expect[i]-1)*x0+(x0-k_expect[i])*(ifelse(x0>k_expect[i],1,0)))
    v = a_expect[i]*x0^2*(1-x0)^2+b_expect[i]*Nt*(k_expect[i])^2*x0*(1-x0)^3*(ifelse(x0>k_expect[i],1,0))+
      b_expect[i]*Nt*(1-k_expect[i])^2*x0^3*(1-x0)*(ifelse(x0<=k_expect[i],1,0))
    betas = seq(-1, floor(beta_old[i]), length = (abs(floor(beta_old[i]))-1)*10+1)
    res = rep(0, length(betas))
    for(j in 1:length(betas)){
      means = betas[j]*((k_expect[i]-1)*x0+(x0-k_expect[i])*(ifelse(x0>k_expect[i],1,0)))
      res[j] = sum(((S[,i]-means)[-c(1, length(means))]/sqrt(v[-c(1, length(means))]))^2)
    }
    beta_expect[i] = betas[which.min(res)]
  }
  
  cat(pvalue, "\n")
  save(S, file = "S.RData")
  save(beta_expect, file = "beta_expect.RData")
  save(beta_old, file = "beta_old.RData")
  save(a_expect, file = "a_expect.RData")
  save(b_expect, file = "b_expect.RData")
  save(k_expect, file = "k_expect.RData")
  save(a0_expect, file = "a0_expect.RData")
  save(b0_expect, file = "b0_expect.RData")
  save(pvalue, file = "pvalue.RData")
  save(trs, file = "trs.RData")
  save(SNR, file = "SNR.RData")
  write.csv(confs, file = "confs.csv")
}else{
  load("S.RData")
  load("a_expect.RData")
  load("b_expect.RData")
  load("beta_expect.RData")
  load("beta_old.RData")
  load("k_expect.RData")
  load("a0_expect.RData")
  load("b0_expect.RData")
  load("SNR.RData")
  load("trs.RData")
  load("pvalue.RData")
  confs = read.csv("confs.csv")[,-1]
}

#####################################################################
####################### Rejection Region ############################
#####################################################################
alt_index = which(round(p.adjust(pvalue, "BH"), 2)<=0.1)
save("alt_index", file = "alt_index.Rdata")