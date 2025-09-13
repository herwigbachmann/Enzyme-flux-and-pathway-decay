
##############COMBINED Mn and ATPase
#Plotting all curves in lactate production in time


a<-a_ori



a=na.omit(a[a$Bgal!="na" & a$Bgal!="173" & a$manganese!=20 & a$well!="F20" & a$well!="H19" & a$well!="M20" & a$well!="N19",])

a_mean <- a %>%
  group_by(time, Bgal,manganese) %>%
  summarise(n = n(),
            rate = mean(lacadd),
            lacY = mean(lacaddrate),
            sd = sd(lacadd)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = rate + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = rate - qt((1-0.95)/2, n - 1) * sem)







library(wesanderson)
pal <- wes_palette("Zissou1", 9, type = "continuous")

ggplot(a_mean, aes(x=time/24, y=rate, color = as.factor(Bgal))) +
  geom_line(aes(x=time/24, y=rate, color=as.factor(Bgal))) +
 # geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=as.factor(Bgal)),color=NA,alpha=0.1)+
  facet_grid(.~manganese)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=16),
        plot.title = element_text(hjust = 0.5,size=16))+
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(name="Time (days)") +
  scale_y_continuous(name="Lactic acid (mM)", limits=c(0, 0.02))+
  labs(title = expression("LA accumulation ATPase overexpression strains \n at different Mn concentrations"),)


#below is just to check replicates make sense
ggplot(na.omit(a), aes(x=time, y=lacadd, color = as.factor(Bgal))) +
  geom_line(aes(x=time, y=lacadd, color=as.factor(Bgal))) +
  geom_smooth(method="lm", linetype="dotted",se=FALSE,fullrange=TRUE)+
  facet_grid(manganese~Bgal)



#plot LA production rate vs yield
rate_mean <- a %>%
  group_by(time, manganese,Bgal) %>%
  summarise(n = n(),
            rate = mean(lacaddrate),
            lacY = mean(lacadd),
            sd = sd(lacaddrate))%>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = rate + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = rate - qt((1-0.95)/2, n - 1) * sem)

rate_mean = rate_mean[rate_mean$lacY>0.006 & rate_mean$rate>-4,]

library(broom) # to change model results into data frame

value = rate_mean %>% 
  group_by(manganese,Bgal) %>%
  do(glance(lm(rate ~ lacY, data = .)))

write.csv(value, file="Figures/Fits_linear_model_Mn-ATPase_data.csv")




#Supplementary Figure S1----
pdf("Figures/Fig.S1.pdf",width=11,height=6)

FigS1<-ggplot(na.omit(rate_mean[rate_mean$lacY>0.006 & rate_mean$rate>-4,]), aes(x=lacY, y=rate, color = as.factor(manganese))) +
 # geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=as.factor(manganese)),color=NA,alpha=0.4)+
  geom_point()+
  geom_smooth(method="lm", linetype="dashed",linewidth = 0.5,se=FALSE,fullrange=TRUE)+
  facet_wrap(.~Bgal,ncol=5)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=16),
        plot.title = element_text(hjust = 0.5,size=16))+
  scale_color_manual(values = pal, name="Manganese (mM)") +
  scale_fill_manual(values = pal, name="Manganese (mM)") +
  scale_x_continuous(name="Lactic acid (M)",limits=c(0,0.025),breaks = seq(0,0.025,by=0.01)) +
  scale_y_continuous(name=LA_rate,limits=c(-5.5,0),breaks = seq(-5,0,by=2))+
  labs(title = expression("F1-ATPase activity (Miller Units)"),)

# geom_line(aes(x=median, y=mean, color=as.factor(manganese))) +

print(FigS1)
FigS1

dev.off()

pdf("Figures/Fig.S1_grid.pdf",width=11,height=6)
FigS1_grid<-ggplot(na.omit(rate_mean[rate_mean$lacY>0.005 & rate_mean$rate>-4,]), aes(x=lacY, y=rate, color = as.factor(manganese))) +
  # geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=as.factor(manganese)),color=NA,alpha=0.4)+
  geom_point()+
  geom_smooth(method="lm", linetype="dashed",linewidth = 0.5,se=FALSE,fullrange=TRUE)+
  facet_grid(manganese~Bgal)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=16),
        plot.title = element_text(hjust = 0.5,size=16))+
  scale_color_manual(values = pal, name="Manganese (mM)") +
  scale_fill_manual(values = pal, name="Manganese (mM)") +
  scale_x_continuous(name="Lactic acid (M)",limits=c(0,0.016),breaks = seq(0,0.025,by=0.01)) +
  scale_y_continuous(name=LA_rate,limits=c(-5.5,0),breaks = seq(-5,0,by=2))+
  labs(title = expression("F1-ATPase activity (Miller Units)"),)

# geom_line(aes(x=median, y=mean, color=as.factor(manganese))) +

print(FigS1_grid)
dev.off()
FigS1_grid



###############################Correlation plot
library(broom) # to change model results into data frame

fit = a[a$lacadd>0.006 & a$lacaddrate>-4,]


fit = fit %>% 
  group_by(well,manganese,Bgal) %>%
  do(tidy(lm(lacaddrate ~ lacadd, data = .)))


fit_short=fit[1:5]
fit_short=dcast(fit_short,well+manganese+Bgal~term)
names(fit_short)[4] = "Intercept"
fit_short$lacadd=-fit_short$lacadd



###########Plot the combination off all conditions combining F1-ATPase overexpression strains and Mn variations

library(wesanderson)
pal <- wes_palette("Zissou1", 9, type = "continuous")


###define color palette for each 

library(scales)
library(ggnewscale)
library(ggpmisc)

#Figure 3----
pdf("Figures/Fig.3.pdf",width=5,height=5)
Fig3 = ggplot(fit_short,aes(x=lacadd,y=Intercept))+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=20),
        plot.title = element_text(hjust = 0.5,size=20),axis.title.x = element_text(margin=margin(10,0,0,0)),axis.title.y = element_text(margin=margin(0,10,0,0)), panel.spacing = unit(1, "lines"))+
  scale_shape_manual(values = 0:9, name="Manganese (mM)")+
  scale_color_manual(values = pal, name="F1-ATPase activity \n (Miller Unit)") +
  scale_fill_manual(values = pal) +
  scale_x_continuous(name=DC, limits=c(0,400),breaks=seq(0,400,by=100)) +
  scale_y_continuous(Max_LA_rate,limits=c(-4,0),breaks = seq(-4,0,by=2))+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE,size=7)+
  geom_jitter(data=fit_short[fit_short$manganese=="0",], size=5,alpha=0.6,aes(col=(Bgal)))+ scale_colour_steps(high = "#f5e9ea", low = "#f50014",n.breaks=9)+
  new_scale_color() +
  geom_jitter(data=fit_short[fit_short$manganese=="0.25",], size=5,alpha=0.6,aes(col=(Bgal)))+ scale_colour_steps(high = "#faf4ed", low = "#fa8901",n.breaks=9)+
  new_scale_color() +
  geom_jitter(data=fit_short[fit_short$manganese=="0.5",], size=5,alpha=0.6,aes(col=(Bgal)))+ scale_colour_steps(high = "#faf8ed", low = "#fad400",n.breaks=9)+
  new_scale_color() +
  geom_jitter(data=fit_short[fit_short$manganese=="0.75",], size=5,alpha=0.6,aes(col=(Bgal)))+ scale_colour_steps(high = "#f2fffa", low = "#00ba70",n.breaks=9)+
  new_scale_color() +
  geom_jitter(data=fit_short[fit_short$manganese=="1",], size=5,alpha=0.6,aes(col=(Bgal)))+ scale_colour_steps(high = "#d3dcde", low = "#00c0de",n.breaks=9)+
  new_scale_color() +
  geom_jitter(data=fit_short[fit_short$manganese=="1.5",], size=5,alpha=0.6,aes(col=(Bgal)))+ scale_colour_steps(high = "#f2f8ff", low = "#00428d",n.breaks=9)+
  new_scale_color() +
  geom_jitter(data=fit_short[fit_short$manganese=="2",], size=5,alpha=0.6,aes(col=(Bgal)))+ scale_colour_steps(high = "#e2dae6", low = "#5f2879",n.breaks=9)+
  geom_smooth(method = lm,se = FALSE)

print(Fig3)
Fig3
dev.off()

fit_s<-fit_short
fit_s <- rename(fit_s, Decay_constant = Intercept)
summary(lm(Decay_constant~lacadd, data = fit_s))

