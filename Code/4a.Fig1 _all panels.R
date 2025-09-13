library(ggsci)
library(ggpmisc)
library(maditr)
library(gridExtra)
library(broom) 
library(ggpmisc)
library(wesanderson)
library(scales)

#calculate mean values
a_mean <- a %>%
  group_by(time, preculture,addition) %>%
  summarise(n = n(),
            lacY = mean(lacadd),
            rate = mean(lacaddrate),
            sd = sd(lacadd)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = lacY + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = lacY - qt((1-0.95)/2, n - 1) * sem)


rate_mean <- a %>%
  group_by(time, preculture,addition) %>%
  summarise(n = n(),
            rate = mean(lacaddrate),
            lacY = mean(lacadd),
            sd = sd(lacaddrate))%>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = rate + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = rate - qt((1-0.95)/2, n - 1) * sem)



#individual panels LA production per condition
ggplot(a_mean, aes(x=time/24, y=lacY,  color = addition)) +
  geom_line() +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=addition),color=NA,alpha=0.4,na.rm = TRUE)+
  facet_grid(addition~preculture)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=12),
        plot.title = element_text(hjust = 0.5,size=12), axis.text = element_text(size=12),axis.title.x = element_text(margin=margin(10,0,0,0)),axis.title.y = element_text(margin=margin(0,10,0,0)), panel.spacing = unit(1, "lines"))+
  scale_color_npg(name="Transitioned sugar", labels=c("galactose","glucose","lactose"))+
  scale_fill_npg(name="Transitioned sugar", labels=c("galactose","glucose","lactose"))+
  scale_x_continuous(name="Time (days)",breaks = seq(0,25,by=10)) +
  scale_y_continuous(name="Lactic acid (M)", limits=c(0, 0.025),breaks = seq(0,0.025,by=0.01))+
  labs(title = expression("Sugar preculture"),)

#Fig. 1 panel A----
pdf("Figures/Fig.1A.pdf",width=7,height=5)
pal <- wes_palette("Zissou1", 4, type = "continuous")

#LA production per condition - precultures combined per panel
panelA <- ggplot(a_mean, aes(x = time / 24, y = lacY, color = addition)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = addition), color = NA, alpha = 0.2) +
  geom_line() +
  facet_grid(. ~ preculture) +
  theme_light(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(colour = 'black', size = 12),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    panel.spacing = unit(1, "lines"),
    legend.position = "right"  
  ) +
 # scale_color_npg(name = "Transitioned sugar", labels = c("galactose", "glucose", "lactose")) +
 # scale_fill_npg(name = "Transitioned sugar", labels = c("galactose", "glucose", "lactose")) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(name = "Time (days)", breaks = seq(0, 25, by = 10)) +
  scale_y_continuous(name = "Lactic acid (M)  \n                                 ", limits = c(0, 0.025), breaks = seq(0, 0.025, by = 0.01)) +
  labs(title = expression("Sugar preculture"))
print(panelA)
dev.off()
panelA 


#define a lables for plots

LA_rate <- bquote("LA production rate" ~ (Log[10] ~ M ~ h^-1))
Max_LA_rate<-bquote("Max. LA production rate" ~ (Log[10] ~ M ~ h^-1))
DC<- bquote("Decay constant" ~( h^-1))

#LA producton rate vs LA yield
#Fig. 1 panel B----
pdf("Figures/Fig.1B.pdf",width=7,height= 5)

panelB<-ggplot(na.omit(rate_mean[rate_mean$lacY>0.006,]), aes(x=lacY, y=rate, color = addition)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=addition),color=NA,alpha=0.2)+ 
   geom_point()+
  geom_smooth(method="lm", linetype="dashed", linewidth=0.5, se=FALSE,fullrange=TRUE)+
  facet_grid(.~preculture)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_blank(),legend.position = "none", axis.text = element_text(size=12),axis.title.x = element_text(margin=margin(10,0,0,0)),axis.title.y = element_text(margin=margin(0,10,0,0)), panel.spacing = unit(1, "lines"))+
 # scale_color_npg(name="Transitioned sugar", labels=c("galactose","glucose","lactose"))+
#scale_fill_npg(name="Transitioned sugar", labels=c("galactose","glucose","lactose"))+
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(name="Lactic acid (M)",limits=c(0.005,0.025),breaks = seq(0.005,0.025,by=0.0075),labels = scales::number_format(accuracy = 0.001)) +
  scale_y_continuous(
    name = LA_rate,
    limits = c(-5, -1.5),
    breaks = seq(-5, -1.5, by = 1.75),
    labels = number_format(accuracy = 0.01)
  )
print(panelB)
dev.off()
panelB


###############################Correlation plot - Fig. 1 panel C

fit = a[a$lacadd>0.006,]

fit = fit %>% 
  group_by(well) %>%
  do(tidy(lm(lacaddrate ~ lacadd, data = .)))

fit_short=fit[1:3]
fit_short=dcast(fit_short,well~term)
names(fit_short)[2] = "Intercept"
fit_short$lacadd=-fit_short$lacadd

#Fig.1C----
pdf("Figures/Fig.1C.pdf",width=5,height=5)

panelC = ggplot(fit_short,aes(x=lacadd,y=Intercept,col='#3C5488FF',fill='#3C5488FF'))+
  geom_point(size=2,alpha=0.75,shape=16)+
  geom_smooth(method="lm",se=FALSE)+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_light(base_size = 16)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_rect(
          fill=NA),strip.text = element_text(colour = 'black',size=12),
        plot.title = element_text(hjust = 0.5,size=16),legend.position = "none",axis.title.x = element_text(margin=margin(10,0,0,0)),axis.title.y = element_text(margin=margin(0,10,0,0)), panel.spacing = unit(1, "lines"))+
  scale_x_continuous(
    name=DC, 
    limits=c(0,275),
    breaks=seq(0,300,by=100)) +
  scale_y_continuous(
    name= Max_LA_rate,limits=c(-4,0),
    breaks = seq(-4,0,by=2))
print(panelC)
dev.off()
panelC

#grid.arrange(panelA,panelB, panelC)





#Write fits to file
#Fitting a linear model to the rate data
value = rate_mean %>% 
  dplyr::filter(lacY > 0.006) %>% #remove noise 
  group_by( preculture,addition) %>%
  do(glance(lm(rate ~ lacY, data = .)))



write.csv(value, file="Figures/Fits_linear_model_Fig.1.csv")
 
