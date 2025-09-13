library(sgolay)

time<-data$Time
Fl<-data[,-1:-2]

#convert fluorescence signal to pH  - basd on a standard curve
pH<- ((2.8974*(log10(Fl)))-7.0201)
corpH<-apply(pH, MARGIN = 2, FUN = function(x) x-x[1])
#convert the pH the the proton concentration
proton<-10^(-pH)
#calculate the cumulative amount of protos produced in each well by subracting the first value from each column from all values in the column
protoncum<-apply(proton, MARGIN = 2, FUN = function(x) x-x[1]) 
#calculate the amount of lactic acid produced in each well based on a standard curve that was made in assay medium with 2e7 cells/ml (they buffer as well)  (lactic acid (mM) = 0.0041ln(proton) + 0.0804)
lacadd<-apply(protoncum, MARGIN = c(1,2),FUN = function (x) 0.0041*log(x)+0.0804) 
lacaddc<-lacadd
lacaddc[lacaddc=="NaN"]<-0
lacaddc[lacaddc=="-Inf"]<-0
lacaddc[lacaddc=="Inf"]<-0


#smoothing filter without derivation
smooth<-apply(lacaddc, MARGIN = 2, FUN = function (x) sgolayfilt(x = x,p = 1, m = 1, n = 21,ts = .5))


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
             sugar_pc=rep(f$sugars_percent,times=repi),
             pH=matrix(t(pH),ncol=1),
             corpH=matrix(t(corpH),ncol=1),
             protoncum=matrix(t(protoncum),ncol=1),
             proton=matrix(t(proton),ncol=1),
             lacadd=matrix(t(lacaddc),ncol=1),
             inoc=rep(f$inoc_density,times=repi),
             antibiotics=rep(f$antibiotics,times=repi),
             medium=rep(f$medium,times=repi),
             addition=rep(f$addition,times=repi), 
             preculture=rep(f$preculture,times=repi),
             Bgal=rep(f$gal_activity,times=repi), 
             manganese=rep(f$mn_conc,times=repi), 
             well=rep(f$well,times=repi),
             smooth=matrix(t(smooth),ncol=1),
             smoothrate=matrix(t(smoothrate),ncol=1),
             lacaddrate=matrix(t(lacaddrate),ncol=1)
             )



#test plot individual wells
ggplot(a[a$well=="A3",], aes(x=time, y=pH, group = "well")) +geom_point()
  
