## Comparison of covid-19 vaccine efficacy of Pfizer_BioNTech,
## J&J and Novavax vaccines by weekly or biweekly data 
## in initial trial reports in FDA briefing reports.

## Pfizer-BioNTech Vaccine Data ##
## Data Source: Pfizer-BioNTech Vaccine document Figure 13 ##
## Figure 13 (P.58)	in FDA Pfizer-BioNTech Covid-19 Briefing report ##
## https://www.fda.gov/media/144246/download December 10, 2020.##

days = seq(0, 112, 7)
casecum.v = c(0, 21, 37, 39, 41, 42, 42, 43, 44, 47, 48, 48, 49, 49, 50, 50, 50)
casecum.u = c(0, 25, 55, 73, 97, 123, 143, 166, 192, 212, 235, 249, 257, 267, 274, 275, 275)

risk.v = c(21314, 21230, 21054, 20481, 19314, 18377, 17702, 17186, 15464, 14038, 12169, 9591, 6403, 3374, 1463, 398, 0)
risk.u = c(21258, 21170, 20970, 20366, 19209, 18218, 17578, 17025, 15290, 13876, 11994, 9471, 6294, 3301, 1449, 398, 0)

## exposure is calculated of 1) 7days * each subject at the end of current week + 3.5 days * each
## subject at the end of previous week and lost in the current week (including positives in this week).
## then divided by 365 to become person-year (unit for exposure). 
expos.v= (-diff(risk.v)/2 + risk.v[-1])*7/365
expos.u= (-diff(risk.u)/2 + risk.u[-1])*7/365

rate.v = diff(casecum.v) / expos.v
rate.u = diff(casecum.u) / expos.u

effica = 1-rate.v/rate.u
effica.loess = loess(effica[1:15]~days[2:16])

eff.pf_b = effica[1:15]
eff.loess.pf_b=effica.loess$fitted[1:15]
eff.loess.pf_b.x = effica.loess$x[1:15]

### Plot change parameters for plots ###
par(cex=1.5, lwd=2)

### plot of efficacy of  Pf-B  ###
plot(days[2:16], effica[1:15], col=4, xlim=c(0, 120),  ylim=c(-.3, 1.1),
     xlab="days after first dose", ylab = "efficacy")
points(days[2:4], effica[1:3], col=2)
lines(effica.loess$x[-1:-2],effica.loess$fitted[-1:-2],col=4)
lines(effica.loess$x[1:3],effica.loess$fitted[1:3],col=2)
title(main="Pfizer-BioNTech: efficacy, loess smoothing, 95% CI")
##legend("bottomright",pch=c(1,1), col=c(2,4), legend=c("Before 21 days","After 21 days"))
legend("bottomright",pch=c(1,1), col=c(2,4), legend=c("Before 21 days","After 21 days"))
abline(h=.95, col=1, lty=2)
abline(h=.903, col=1, lty=4)
abline(h=.976, col=1, lty=4)

### Empirical mean efficacy and standard deviation
mean.effi = mean(effica[4:15])
sd.effi = sqrt(var(effica[4:15]))
c(mean.effi, sd.effi)
## sd.effi = sqrt(var(effica[3:15])/length(effica[3:15]))

### Empirical mean efficacy and standard deviation in stable time-window 28 - 105 days for reliable estimation
mean.effi.window = mean(effica[4:15])
sd.effi.window = sqrt(var(effica[4:15]))
c(mean.effi.window, sd.effi.window)

c(mean.effi, sd.effi)
abline(h=mean.effi, col="orange", lty=2)
abline(h=mean.effi-1.96*sd.effi, col="orange", lty=4)
abline(h=min(1, mean.effi+1.96*sd.effi), col="orange", lty=4)
#legend("bottom",lty=4,col=c("blue", "orange"), legend=c("Reported 95% C.I.", "Empirical 1.96*SD"))
#legend("bottom", lty=c(1,1,4), col=c(2,4, 5), legend=c("loess smooth before 21 days", "reported efficacy"))
legend("bottom",lty=c(2,4,2,4),col=c("blue", "blue", "orange", "orange"), legend=c("Reported efficacy",  "Reported 95% C.I.", "Empirical efficacy", "Empirical 1.96*SD"))

## Examining normality of weekly efficacy ##
qqnorm((effica[4:15]-mean(effica[4:15]))/sqrt(var(effica[4:15])), ylim=c(-2.2, 2.2), xlim=c(-2.2, 2.2))
abline(a=0, b=1)

############ End of Pfizer-BioNTech data ################


### J&J document figure 1 
### Data Source: Figure 1 (P.31) in FDA (2021) FDA Briefing Document Janssen Ad26.COV2.S Vaccine.
### website https://www.fda.gov/media/146217/download , February 26, 2021. 
days = seq(0,126,7)
risk.v = c(19744, 19725, 19669, 19642, 19612, 19578, 18541, 14909, 10930, 7831, 3998, 1468, 713, 484, 483, 482, 142, 31, 0)
risk.u = c(19822, 19804, 19745, 19652, 19579, 19488, 18411, 14814, 10823, 7740, 3876, 1439, 708, 485, 482, 480, 133, 27, 0)
casecum.v = c(0, 27, 76, 96, 126, 151, 168, 178, 184, 188, 189, 191, 191, 192, 193, 193, 193, 193, 193)
casecum.u = c(0, 22, 81, 168, 237, 299, 351, 387, 407, 416, 423, 425, 430, 432, 432, 432, 432, 432, 432)

expos.v = (-diff(risk.v)/2*7+risk.v[-1]*7)/365
expos.u = (-diff(risk.u)/2*7+risk.u[-1]*7)/365
case.v = diff(casecum.v)
case.u = diff(casecum.u)

rate.v = case.v/expos.v
rate.u = case.u/expos.u

effi = 1-rate.v/rate.u
eff.loess = loess(effi[1:13]~days[2:14])

eff.jj = effi[1:13]
eff.loess.jj=eff.loess$fitted
eff.loess.jj.x = eff.loess$x

### plot of efficacy of J&J ###
plot(days[2:14],effi[1:13],xlim=c(0,100), ylim=c(-.3, 1.1), 
     pch=1, xlab="days after first dose", ylab="efficacy",col=4)
lines(eff.loess$x[1:3], eff.loess$fitted[1:3],col=2)
lines(eff.loess$x[-1:-2], eff.loess$fitted[-1:-2],col=4)
points(days[2:4],effi[1:3],col=2)
points(days[5:14],effi[4:13],col=4)
abline(h=.661, col=1, lty=2)
abline(h=.55, col=1, lty=4)
abline(h=.748, col=1, lty=4)
title(main="J&J: efficacy, loess smoothing, 95% CI")
legend("bottomright",pch=c(1,1), col=c(2,4), legend=c("Before 21 days","After 21 days"))

### Empirical mean efficacy and standard deviation
mean.effi = mean(effi[4:13])
sd.effi = sqrt(var(effi[4:13]))
c(mean.effi, sd.effi)

### Empirical mean efficacy and standard deviation in time-window 28-63 days
mean.effi.window = mean(effi[4:9])
sd.effi.window = sqrt(var(effi[4:9]))
c(mean.effi.window, sd.effi.window)

abline(h=mean.effi, col="orange", lty=2)
abline(h=mean.effi-1.96*sd.effi, col="orange", lty=4)
abline(h=min(1, mean.effi+1.96*sd.effi), col="orange", lty=4)
#legend("bottom",lty=4,col=c("blue", "orange"), legend=c("Reported 95% C.I.", "Empirical 1.96*SD"))
legend("bottom",lty=c(2,4,2,4),col=c("blue", "blue", "orange", "orange"), legend=c("Reported efficacy",  "Reported 95% C.I.", "Empirical efficacy", "Empirical 1.96*SD"))


qqnorm((effi[4:13]-mean(effi[4:13]))/sqrt(var(effi[4:13])), ylim=c(-2.2, 2.2), xlim=c(-2.2, 2.2))
abline(a=0, b=1)

############ End of J & J data ################

### Novavax Vaccine Data ### 
### Data source: Figure 1 (P.26) in	FDA briefing Document Novavax COVID-19 Vaccine
### https://www.fda.gov/media/158912/download
# Placebo cumcases  0,18,41,58,73,86,91,101,107,113,121,129,141,145,152,152,153,156,156
# vaccine cumcases  0,46,88,100,109,116,117,120,122,125,126,128,128,129,132,133,133,133,133
# Placebo at risk  9868,9800,9720,9653,9410,9180,8773,8421,8082,7790,7482,6730,5654,4513,3445,2544,1500,739,524,
# vaccine at risk  19714,19609,19483,19372,19000,18755,18311,17931,17568,17239,16863,15389,13015,10532,7925,5808,3181,1335,753

days = seq(0,126,7)
risk.v = c(19714,19609,19483,19372,19000,18755,18311,17931,17568,17239,16863,15389,13015,10532,7925,5808,3181,1335,753)
risk.u = c(9868,9800,9720,9653,9410,9180,8773,8421,8082,7790,7482,6730,5654,4513,3445,2544,1500,739,524)
casecum.v = c(0,46,88,100,109,116,117,120,122,125,126,128,128,129,132,133,133,133,133)
casecum.u = c(0,18,41,58,73,86,91,101,107,113,121,129,141,145,152,152,153,156,156)

  expos.v = (-diff(risk.v)/2*7+risk.v[-1]*7)/365
  expos.u = (-diff(risk.u)/2*7+risk.u[-1]*7)/365
  case.v = diff(casecum.v)
  case.u = diff(casecum.u)
  
  ### adjustment of last 4 weeks due to 0 cases in both vaccinated and unvaccinated groups
  length.case = length(case.v)  ### in case it is needed.
  
  ### change the weekly data to biweekly data of the last 4 weeks.
  case.u_4 = case.u [-3:0+length.case]
  case.v_4 = case.v [-3:0+length.case]
  
  case.u_biw2 = case.u_4[c(1,3)] + case.u_4[c(2,4)]
  case.v_biw2 = case.v_4[c(1,3)] + case.v_4[c(2,4)]
  
  
  expos.u_4 = expos.u [-3:0+length.case]
  expos.v_4 = expos.v [-3:0+length.case]
  
  expos.u_biw2 = expos.u_4[c(1,3)] + expos.u_4[c(2,4)]
  expos.v_biw2 = expos.v_4[c(1,3)] + expos.v_4[c(2,4)]
  
  case.v = case.v[-(-3:0+length.case)]
  case.u = case.u[-(-3:0+length.case)]
  
  expos.v = expos.v[-(-3:0+length.case)]
  expos.u = expos.u[-(-3:0+length.case)]
  
  rate.u = case.u/expos.u
  rate.v = case.v/expos.v
  
  rate.u_2 = case.u_biw2/expos.u_biw2
  rate.v_2 = case.v_biw2/expos.v_biw2


effi = 1-rate.v/rate.u
effi_biw2 = 1-rate.v_2 /rate.u_2

## Handle the missing data in risk issue due to 0 cases in control group.
effi.plot = c(effi[1:(length.case-3)], effi_biw2)
days.plot = c(days[1+1:(length.case-3)], days[1+length.case-c(2,0)])
eff.loess = loess(effi.plot~days.plot)

eff.novax = effi.plot
eff.loess.novax = eff.loess$fitted
eff.loess.novax.x = eff.loess$x
days.plot.novax =days.plot

### plot of efficacy of Novavax ###
plot (days.plot, effi.plot, xlim=c(0, 130), ylim=c(-.3, 1.1), 
      xlab="days after first dose", ylab = "efficacy", col=4)
lines(eff.loess.novax.x[1:3], eff.loess.novax[1:3], col=2)
lines(eff.loess.novax.x[-1:-2], eff.loess.novax[-1:-2],col=4)
points(days.plot[1:3], effi.plot[1:3], col=2)
abline(h=.9041, lty=2)  #reported efficacy
abline(h=.9432, lty=4)  #reported upper limit of 95% CI
abline(h=.8331, lty=4)  #reported lower limit of 95% CI
title(main="Novavax: efficacy, loess smoothing, 95% CI")
legend("bottomright",pch=c(1,1), col=c(2,4), legend=c("Before 21 days","After 21 days"))

### Empirical mean efficacy and standard deviation
mean.effi = mean(effi.plot[-1:-3], na.rm = T)
sd.effi = sqrt(var(effi.plot[-1:-3], na.rm = T))
c(mean.effi, sd.effi)

### Empirical mean efficacy and standard deviation within time window 28 - 98 days
mean.effi.window = mean(effi.plot[4:14], na.rm = T)
sd.effi.window = sqrt(var(effi.plot[4:14], na.rm = T))
c(mean.effi.window, sd.effi.window)

c(mean.effi, sd.effi)
abline(h=mean.effi, col="orange", lty=2)
abline(h=mean.effi-1.96*sd.effi, col="orange", lty=4)
abline(h=min(1, mean.effi+1.96*sd.effi), col="orange", lty=4)
legend("bottom",lty=c(2,4,2,4),col=c("blue", "blue", "orange", "orange"), legend=c("Reported efficacy",  "Reported 95% C.I.", "Empirical efficacy", "Empirical 1.96*SD"))

## Checking normality assumption of weekly efficacies with qqnorm plot. 
qqnorm((effi.plot[-c(1:3, 15)]-mean(effi.plot[-c(1:3, 15)]))/sqrt(var(effi.plot[-c(1:3, 15)])), xlim=c(-2.2, 2.2), ylim=c(-2.2, 2.2))
abline(a=0, b=1)

####### End of Novavax vaccine data analysis ######

#############################################
### Comparison of vaccine efficacies ###
plot(eff.loess.pf_b.x, eff.loess.pf_b, xlab="days after first dose", ylab="efficacy", 
     xlim=c(0, 130), ylim=c(-.3, 1.1), type='l', col=1)
lines(eff.loess.jj.x, eff.loess.jj, col=2)
lines(eff.loess.novax.x, eff.loess.novax, col=4)
title(main="Comparison of vaccine efficacy")
legend("bottomright", lty=c(1,1,1), col=c(1,2,4), pch=c(1,2,4), 
       legend=c("Pf-B", "J&J", "Novavax") )
points(eff.loess.pf_b.x[1:15], effica[1:15], pch=1)
points(eff.loess.jj.x[1:13], eff.jj,pch=2,col=2)
points(days.plot, eff.novax, pch=4, col=4)

############ End of vaccine efficacy comparison ############