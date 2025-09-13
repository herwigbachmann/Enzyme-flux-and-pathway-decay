library(sgolay)

time<-data$Time
Fl<-data[,-1:-2]

#convert fluorescence signal to pH  - basd on a standard curve
pH<- ((2.6765*(log10(Fl)))-5.6744)
#convert the pH the the proton concentration
proton<-10^(-pH)
#calculate the cumulative amount of protos produced in each well
protoncum<-apply(proton, MARGIN = 2, FUN = function(x) x-x[1])
#calculate the amount of lactic acid produced in each well based on a standard curve that was made in assay medium with 2e7 cells/ml (they buffer as well) (lactic acid (mM) = 0.0041ln(proton) + 0.0804)
lacadd<-apply(protoncum, MARGIN = c(1,2),FUN = function (x) 0.0041*log(x)+0.0804)
lacaddc<-lacadd
lacaddc[lacaddc=="NaN"]<-0
lacaddc[lacaddc=="-Inf"]<-0
lacaddc[lacaddc=="Inf"]<-0

#smoothing filter without derivation
smooth<-apply(lacaddc, MARGIN = 2, FUN = function (x) sgolayfilt(x = x,p = 1, m=0, n = 51,ts = .5))

#production rate
lacaddrate = log10(lacaddc/time)
smoothrate = log10(smooth/time)

lacaddrate[lacaddrate=="NaN"]<-NA
lacaddrate[lacaddrate=="-Inf"]<-NA
lacaddrate[lacaddrate=="Inf"]<-NA

smoothrate[smoothrate=="NaN"]<-NA
smoothrate[smoothrate=="-Inf"]<-NA
smoothrate[smoothrate=="Inf"]<-NA


#make dataframe for plotting
repi=length(time)
a=data.frame(
             time=rep(data$Time, each = 384),
             strain=rep(f$strain,times=repi), 
             carbon_pc=rep(f$csource_percent,times=repi),
             initialpH=rep(f$initial_ph,times=repi),
             pH=matrix(t(pH),ncol=1),
             protoncum=matrix(t(protoncum),ncol=1),
             proton=matrix(t(proton),ncol=1),
             lacadd=matrix(t(lacaddc),ncol=1),
             inoc=rep(f$inoc_density,times=repi),
             antibiotics=rep(f$antibiotics,times=repi),
             medium=rep(f$medium,times=repi),
             addition=rep(f$addition1,times=repi), 
             preculture=rep(f$preculture,times=repi), 
             manganese=rep(f$manganese,times=repi), 
             well=rep(f$well,times=repi),
             contamination=rep(f$contaminated,times=repi), #during the weeks of incubation in some wells growth started occuring (we do not know if this is contamination or spontaneous antibiotic resistance) ->remove
             smooth=matrix(t(smooth),ncol=1),
             smoothrate=matrix(t(smoothrate),ncol=1),
             lacaddrate=matrix(t(lacaddrate),ncol=1), 
             trans=rep(f$transition,times=repi))

a=na.omit(a[a$contamination=="no" & a$inoc=="1.00E+07" & a$trans=="1",])

  