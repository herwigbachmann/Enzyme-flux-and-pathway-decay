library(wesanderson)
library(ggsci)
library(ggpmisc)
library(maditr)
library(gridExtra)
library(broom) 
library(ggpmisc)

a$Bgal<-as.numeric(a$Bgal)
a_ori<-a


###############MANGANESE EFFECT select relevant wells########

a=a[(a$strain=="MG1363" | a$strain=="NCDO712") & a$well!="J24" & a$well!="A21" & a$well!="A22"& a$well!="B21"& a$well!="B22"& a$well!="C21"& a$well!="C22" & a$well!="D21"& a$well!="D22" & a$well!="E21"& a$well!="E22" & a$well!="F21" & a$well!="F22" &  a$manganese!="20" ,]


####calculate mean values per welll
rate_mean <- a %>%
  group_by(time, strain,manganese) %>%
  summarise(n = n(),
            rate = mean(lacaddrate),
            lacY = mean(lacadd),
            sd = sd(lacaddrate))%>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = rate + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = rate - qt((1-0.95)/2, n - 1) * sem)

rate_mean_Mn<-rate_mean

###########Plot LA production over time


a_mean <- a %>%
  group_by(time, strain,manganese) %>%
  summarise(n = n(),
            rate = mean(lacaddrate),
            lacY = mean(lacadd),
            sd = sd(lacadd)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = lacY + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = lacY - qt((1-0.95)/2, n - 1) * sem)

value = rate_mean %>% 
  group_by(manganese,strain) %>%
  do(glance(lm(rate ~ lacY, data = .)))

write.csv(value, file="Figures/Fits_linear_model_Mn_data.csv")

pal <- wes_palette("Zissou1", 7, type = "continuous")

ggplot(a_mean, aes(x=time/24, y=lacY, color = as.factor(manganese))) +
  geom_point(na.rm = TRUE) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=as.factor(manganese)),color=NA,alpha=0.1)+
  facet_grid(.~strain)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=16),
        plot.title = element_text(hjust = 0.5,size=16))+
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(name="Time (days)") +
  scale_y_continuous(name="Lactic acid (mM)", limits=c(0, 0.02))+
  labs(title = expression("Sugar preculture"),)



####Plot rate and yield graph 


#Fig2_PanelA - the selection of rate_mean>0.006 selects for values above the noise leves (see method paper)----
pdf("Figures/Fig.2A.pdf",width=9,height=5)
Fig2a<-ggplot(na.omit(rate_mean[rate_mean$lacY>0.006 & rate_mean$rate>-4,]), aes(x=lacY, y=rate, color = as.factor(manganese))) +
  geom_point()+
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=as.factor(manganese)),color=NA,alpha=0.5)+
  geom_smooth(method="lm", linetype="dashed",linewidth = 0.5,se=FALSE,fullrange=TRUE)+
  facet_grid(.~strain)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=14),
        plot.title = element_text(hjust = 0.5,size=14), axis.text = element_text(size=12),axis.title.x = element_text(margin=margin(10,0,0,0)),axis.title.y = element_text(margin=margin(0,10,0,0)), panel.spacing = unit(1, "lines"))+
  scale_color_manual(values = pal,name="Manganese (mM)") +
  scale_fill_manual(values = pal,name="Manganese (mM)") +
  scale_x_continuous(name="Lactic acid (M)",limits=c(0.005,0.02),breaks = seq(0.005,0.025,by=0.0075),labels = scales::number_format(accuracy = 0.001)) +
  scale_y_continuous(
    name=LA_rate,
    limits=c(-5,-1.5),
    breaks = seq(-5,-1.5,by=1.75),
    labels = scales::number_format(accuracy = 0.01))+
  labs(title = expression("Strain"))
print(Fig2a)
dev.off()
Fig2a 


###############################Correlation plot

fit = a[a$lacadd>0.006 & a$lacaddrate>-4,]


fit = fit %>% 
  group_by(well,strain, manganese) %>%
  do(tidy(lm(lacaddrate ~ lacadd, data = .)))

fit_short=fit[1:5]
fit_short=dcast(fit_short,well+strain+manganese~term)
names(fit_short)[4] = "Intercept"
fit_short$lacadd=-fit_short$lacadd


#Fig2_PanelB----
pdf("Figures/Fig.2B.pdf",width=8,height=5)
Fig2b<-ggplot(fit_short,aes(x=lacadd,y=Intercept,col='#3C5488FF',fill='#3C5488FF'))+
  geom_point(size=2,alpha=0.75,shape=16)+
  geom_smooth(method="lm",se=FALSE)+
  facet_grid(.~strain)+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=12),
        plot.title = element_text(hjust = 0.5,size=16),legend.position = "none",axis.title.x = element_text(margin=margin(10,0,0,0)),axis.title.y = element_text(margin=margin(0,10,0,0)), panel.spacing = unit(1, "lines"))+
  scale_x_continuous(name=DC, limits=c(0,275),breaks=seq(0,300,by=100)) +
  scale_y_continuous(Max_LA_rate,limits=c(-4,0),breaks = seq(-4,0,by=2))

print(Fig2b)
dev.off()
Fig2b
###################    DONE for manganese variations    




##########   Plot strains with F1-ATPase overexpression      

a<-a_ori

#calculate mean
a_mean <- a %>%
  group_by(time, strain, Bgal,manganese) %>%
  summarise(n = n(),
            lacY = mean(lacadd),
            rate = mean(lacaddrate),
            sd = sd(smooth)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = lacY + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = lacY - qt((1-0.95)/2, n - 1) * sem)



pal <- wes_palette("Zissou1", 11, type = "continuous")

#identify wells wells that contain outliers to remove

o<-ggplot(a_mean, aes(x=time/24, y=lacY, color = as.factor(Bgal))) +
  geom_point(aes(x=time/24, y=lacY, color=as.factor(Bgal))) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=as.factor(Bgal)),color=NA,alpha=0.1)+
  facet_grid(manganese~Bgal)+
  theme_light(base_size = 6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=16),
        plot.title = element_text(hjust = 0.5,size=16))+
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(name="Time (days)") +
  scale_y_continuous(name="Lactic acid (M)", limits=c(0, 0.02))+
  labs(title = expression("Sugar preculture"))
#print(o)



###remove outliers and on samples without additional manganese
a=na.omit(a[a$Bgal!="na" & a$Bgal!="173" & a$manganese==0 & a$well!="F20" & a$well!="H19" & a$well!="A19" &a$well!="M20" & a$well!="N19",])


####rate and yield data mean values and CI

rate_mean <- a %>%
  group_by(time,strain,Bgal) %>%
  summarise(n = n(),
            rate = mean(lacaddrate),
            lacY = mean(lacadd),
            sd = sd(lacaddrate))%>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = rate + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = rate - qt((1-0.95)/2, n - 1) * sem)

rate_mean_ATPase<-rate_mean
pal <- wes_palette("Zissou1", 10, type = "continuous")


#Fig2_PanelC----
pdf("Figures/Fig.2C.pdf",width=5,height=5)
Fig2c<-ggplot(na.omit(rate_mean[rate_mean$lacY>0.008,]), aes(x=lacY, y=rate, color = as.factor(Bgal))) +
  geom_point()+
 # geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=as.factor(Bgal)),color=NA,alpha=0.5)+
  geom_smooth(method="lm", linetype="dotted",se=FALSE,fullrange=TRUE)+
  facet_grid(.~.)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=12),
        plot.title = element_text(hjust = 0.5,size=12), axis.text = element_text(size=12),axis.title.x = element_text(margin=margin(10,0,0,0)),axis.title.y = element_text(margin=margin(0,10,0,0)), panel.spacing = unit(1, "lines"))+
  scale_color_manual(name="F1-ATPase\n(Miller Units)", values = pal) +
  scale_fill_manual(name="F1-ATPase\n(Miller Units)", values = pal) +
  scale_x_continuous(name="Lactic acid (M)",limits=c(0.005,0.02),breaks = seq(0.005,0.025,by=0.0075),labels = scales::number_format(accuracy = 0.001)) +
  scale_y_continuous(
    name=LA_rate,
    limits=c(-5,-2.2),
    breaks = seq(-5,-1.5,by=1.75),
    labels = scales::number_format(accuracy = 0.01))

print(Fig2c)
dev.off()
Fig2c



###############################Correlation plot
library(broom) 

fit = a[a$lacadd>0.007 ,]
fit = fit[fit$Bgal!="na" & fit$Bgal!="173",]

fit = fit %>% 
  group_by(well,Bgal) %>%
  do(tidy(lm(lacaddrate ~ lacadd, data = .)))



fit_short=fit[1:4]
fit_short=dcast(fit_short,well+Bgal~term)
names(fit_short)[3] = "Intercept"
fit_short$lacadd=-fit_short$lacadd


####Fig2_PanelD----
pdf("Figures/Fig.2D.pdf",width=5,height=5)
Fig2d<-ggplot(fit_short, aes(x = lacadd, y = Intercept, col = '#3C5488FF', fill = '#3C5488FF')) +
  geom_point(size = 2, alpha = 0.75, shape = 16) +
 # geom_text(aes(label = well), hjust = 0, vjust = -1, size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = y ~ x,
               aes(label = paste(..rr.label.., sep = "~~~")),
               parse = TRUE) +
  theme_light(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(colour = 'black', size = 12),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none",
    axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    panel.spacing = unit(1, "lines")
  ) +
  scale_x_continuous(
    name = DC, 
    limits = c(0, 250), 
    breaks = seq(0, 300, by = 100)) +
  scale_y_continuous(
    Max_LA_rate, 
    limits = c(-3, 0), 
    breaks = seq(-4, 0, by = 2))

print(Fig2d)
dev.off()
Fig2d




########DONE FOR ATPASE


