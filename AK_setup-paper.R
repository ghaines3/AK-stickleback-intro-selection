library(dplyr)
library(vegan)
library(ggplot2)
library(MASS)
library(tibble)
library(RColorBrewer)
library(caret)
library(jmv)
library(broom)
library(patchwork)


AKSL<-(read.csv("../AKStickleData - Sheet1_5-19-19.csv")[
  ,c("FishID","Lake","Mass_g","SL_mm","BodyDepth",
     "BuccalCavityLength","GapeWidth","EpWidth","RakerNum")]%>%
    full_join(read.csv("../AKStickleData _7-30-19.csv")[,c(
      "FishID","Caudal.peduncle","Jaw.Length","Snout.Length",
      "Eye.Diameter","Head.Length")], join_by(FishID==FishID)))%>%
  mutate(Caudal.peduncle=Caudal.peduncle*10)%>%
  mutate(Jaw.Length=Jaw.Length*10)%>%
  mutate(Snout.Length=Snout.Length*10)%>%
  mutate(Eye.Diameter=Eye.Diameter*10)%>%
  mutate(Head.Length=Head.Length*10) # multiplications by 10 here are convert cm to mm

# replaces erroneous value w/NA
AKSL[which(AKSL$FishID == "S180159"),"Head.Length"]<-NA
#fixes misplaced decimal
AKSL[which(AKSL$FishID=="S180648"),"Head.Length"]<-
  AKSL[which(AKSL$FishID=="S180648"),"Head.Length"]*10


#### Size correction 
# get regression slopes for log(trait)~log(SL)+Lake
# some measurements multiplied by 10 to change units to mm from cm
ancovaMass.AKSL <- lm(data=AKSL, log(Mass_g)~log(SL_mm)+Lake)
ancovaBCLength.AKSL <- lm(data=AKSL, log(BuccalCavityLength)~log(SL_mm)+Lake)
ancovaBD.AKSL <- lm(data=AKSL, log(BodyDepth)~log(SL_mm)+Lake)
ancovaGapeWidth.AKSL <- lm(data=AKSL, log(GapeWidth)~log(SL_mm)+Lake)
ancovaEpWidth.AKSL <- lm(data=AKSL, log(EpWidth)~log(SL_mm)+Lake)
ancovaCPed.AKSL <- lm(data=AKSL, log(Caudal.peduncle)~log(SL_mm)+Lake)
ancovaJL.AKSL <- lm(data=AKSL, log(Jaw.Length)~log(SL_mm)+Lake)
ancovaSnL.AKSL <- lm(data=AKSL, log(Snout.Length)~log(SL_mm)+Lake)
ancovaED.AKSL <- lm(data=AKSL, log(Eye.Diameter)~log(SL_mm)+Lake)
ancovaHL.AKSL <- lm(data=AKSL, log(Head.Length)~log(SL_mm)+Lake)


#makes trait.coefs df including vector w/ lake
#name and vectors with coefficients from trait ancovas
Coefs.df<-data.frame(Trait=c(
  "Mass","BuccalCavityL","BodyDepth","GapeWidth","EpWidth","RakerNum",
  "CPed","JawL","SnL","EyeD","HeadL"),
  Coef=c(ancovaMass.AKSL$coefficients["log(SL_mm)"],
         ancovaBCLength.AKSL$coefficients["log(SL_mm)"],
         ancovaBD.AKSL$coefficients["log(SL_mm)"],
         ancovaGapeWidth.AKSL$coefficients["log(SL_mm)"],
         ancovaEpWidth.AKSL$coefficients["log(SL_mm)"],NA,
         ancovaCPed.AKSL$coefficients["log(SL_mm)"],
         ancovaJL.AKSL$coefficients["log(SL_mm)"],
         ancovaSnL.AKSL$coefficients["log(SL_mm)"],
         ancovaED.AKSL$coefficients["log(SL_mm)"],
         ancovaHL.AKSL$coefficients["log(SL_mm)"]))


#creates vectors in AKSL of trait values that have been adjusted for allometry
AKSL<-AKSL%>%
  mutate(adj.Mass=Mass_g*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="Mass"),"Coef"])%>%
  mutate(adj.BC=BuccalCavityLength*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="BuccalCavityL"),"Coef"])%>%
  mutate(adj.BD=BodyDepth*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="BodyDepth"),"Coef"])%>%
  mutate(adj.GapeWidth=GapeWidth*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="GapeWidth"),"Coef"])%>%
  mutate(adj.EpWidth=EpWidth*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="EpWidth"),"Coef"])%>%
  mutate(adj.Cped=Caudal.peduncle*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="CPed"),"Coef"])%>%
  mutate(adj.JawL=Jaw.Length*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="JawL"),"Coef"])%>%
  mutate(adj.SnLength=Snout.Length*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="SnL"),"Coef"])%>%
  mutate(adj.EyeD=Eye.Diameter*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="EyeD"),"Coef"])%>%
  mutate(adj.HeadL=Head.Length*(mean(SL_mm, na.rm = T)/SL_mm)^
           Coefs.df[which(Coefs.df$Trait=="HeadL"),"Coef"])

AKSL%>%group_by(Lake)%>%summarise_at(c(
  "SL_mm","adj.Mass","adj.BC","adj.BD",
  "adj.GapeWidth","adj.EpWidth","RakerNum"), 
  list(~mean(.,na.rm = TRUE)%>%round(2) ,~sd(., na.rm=T)%>%round(2)  ))%>%
  kable(col.names = c(
    "Lake","SL mean", "Adj. Mass mean", 
    "Adj. BC mean", "Adj. BD mean","Adj. GW mean", 
    "Adj. PW mean", "Raker Num mean","SL SD", 
    "Adj. Mass SD" , "Adj. BC SD", "Adj. BD SD",
    "Adj. GW SD", "Adj. PW SD", "Raker Num SD"),
        label = "Population means and SDs, Part I",
        caption = "All summary statistics are presented in mm, except for mass, which is preseented in g.
Measurements are adjusted to the overall mean, not the population mean.")

AKSL%>%group_by(Lake)%>%summarise_at(c(
  "adj.Cped","adj.JawL","adj.SnLength","adj.EyeD","adj.HeadL"), 
  list(~mean(.,na.rm = TRUE)%>%round(2),~sd(., na.rm=T)%>%round(2) ))%>%
  kable(col.names = c(
    "Lake","Adj. CPed mean", "Adj. JL mean", "Adj. SnL mean", 
    "Adj. EyeD mean", "Adj. HeadL mean", "Adj. CPed SD", 
    "Adj. JL SD", "Adj. SnL SD", "Adj. EyeD SD", "Adj. HeadL SD"),
        label = "Population means and SDs, Part II",
        caption = "All summary statistics are presented in mm.
Measurements are adjusted to the overall mean, not the population mean.")


#####
theme_set(theme_classic())
AKSL<-(AKSL%>%mutate(Ecomorph=as.factor(if_else(
  Lake=="Tern" | Lake=="Corcoran" | Lake=="Watson" | Lake=="Walby",
  "Benthic", if_else(
    Lake=="South_Rolly"|Lake=="Long", "Limnetic",NA)))))

lda.lake <- lda(data = AKSL%>%mutate_if(is.numeric, scale), 
                Lake ~ adj.Mass+adj.BC + adj.BD + SL_mm + 
               adj.Cped + adj.JawL + adj.SnLength + 
               adj.EyeD + adj.HeadL)
predictions.morph.lake <- na.exclude(predict(
  lda.lake, AKSL%>%mutate_if(is.numeric, scale)))
AKSL.LDlake<-cbind(AKSL,predictions.morph.lake$x)

#prediction_table gives predicted ecomorph counts for each lake
prediction_table.lake<-table(Lake=AKSL$Lake, predictions.morph.lake$class)

#accuracy rate (correct assignment probability)#
assignment_accuracy.lake<-mean(na.omit(
  AKSL$Lake == predictions.morph.lake$class))
assignment_accuracy.lake%>%round(3)
kable(lda.lake$prior%>%round(3), caption="LDA Priors", align = "c")
kable(lda.lake$means%>%round(2), caption="LDA Group Means", align = "c")
kable(lda.lake$scaling%>%round(2), caption="LDA Trait Scaling", align = "c")
kable(prediction_table.lake, caption = "Individual Assignments", align = "c")

plot.lake<-ggplot(AKSL.LDlake[which(!is.na(AKSL.LDlake$LD1)),]%>%
                 mutate(Ecotype=as.factor(if_else(
                   is.na(Ecomorph),"NA",Ecomorph)))%>%
                 mutate(Region = as.factor(ifelse(
                   Lake=="Long" | Lake=="Walby" | 
                     Lake=="South_Rolly" | Lake=="Finger" |
                     Lake == "Corcoran", "Mat-Su", "Kenai"))), aes(x = reorder(
                       Lake, LD1, FUN=mean), y = LD1, fill = Ecotype))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("seagreen4", "skyblue2","coral2"))+
  geom_beeswarm(inherit.aes=FALSE,aes(x = Lake, y = LD1),
                method = 'center', shape = 16,
                size = .8, alpha = .8)+
  scale_colour_manual(values = c("coral4","seagreen4"))+
  #geom_point(alpha = .3, position=position_jitterdodge())+ 
  labs(x = NULL, 
       y = "LD 1 Score")+
  theme(axis.text.x=element_text(size = 4.5),
        axis.text.y=element_text(size = 6, angle = 0),
        axis.title.y =element_text(size = 8 ),
        legend.title = element_text(size = 10))


##########
##########
lda.a <- lda(data = AKSL, Ecomorph ~ adj.Mass
                         + adj.BC + adj.BD + SL_mm + 
                           adj.Cped + adj.JawL + adj.SnLength + 
                           adj.EyeD + adj.HeadL)
predictions.morph.a <- na.exclude(predict(lda.a, AKSL))
AKSL.LDa<-cbind(AKSL,predictions.morph.a$x)

#subsampling and LDA with subsampling
AKSL.sub <- AKSL %>% group_by(Lake) %>% slice_sample(n=23)
lda.sub <- lda(data = AKSL.sub, Ecomorph ~ adj.Mass
             + adj.BC + adj.BD + SL_mm + 
               adj.Cped + adj.JawL + adj.SnLength + 
               adj.EyeD + adj.HeadL)
predictions.morph.sub <- na.exclude(predict(lda.sub, AKSL))
AKSL.LDsub<-cbind(AKSL,predictions.morph.sub$x)

#prediction_table gives predicted ecomorph counts for each lake
prediction_table.a<-table(Lake=AKSL$Lake, predictions.morph.a$class)

#accuracy rate (correct assignment probability)#
assignment_accuracy.a<-mean(na.omit(
  AKSL$Ecomorph == predictions.morph.a$class))
cat("Assignment accuracy when ecomorph known:" )
assignment_accuracy.a%>%round(3)
kable(lda.a$prior%>%round(3), caption="LDA Priors", align = "c")
kable(lda.a$means%>%round(2), caption="LDA Group Means", align = "c")
kable(lda.a$scaling%>%round(2), caption="LDA Trait Scaling", align = "c")
kable(prediction_table.a, caption = "Individual Assignments", align = "c")

AKSL.LDa%>%describe()
AKSL.LDa%>%describeBy(AKSL.LDa$Lake)


plot.a<-ggplot(AKSL.LDa[which(!is.na(AKSL.LDa$LD1)),]%>%
         mutate(Ecotype=as.factor(if_else(
           is.na(Ecomorph),"NA",Ecomorph)))%>%
           mutate(Region = as.factor(ifelse(
             Lake=="Long" | Lake=="Walby" | 
               Lake=="South_Rolly" | Lake=="Finger" |
               Lake == "Corcoran", "Mat-Su", "Kenai"))), aes(x = reorder(
             Lake, LD1, FUN=mean), y = LD1, fill = Ecotype))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("seagreen4", "skyblue2","coral2"))+
  geom_beeswarm(inherit.aes=FALSE,aes(x = Lake, y = LD1),
                method = 'center', shape = 16,
                size = .8, alpha = .8)+
  scale_colour_manual(values = c("coral4","seagreen4"))+
  #geom_point(alpha = .3, position=position_jitterdodge())+ 
  labs(x = NULL, 
       y = "LD Score")+
  theme_update(axis.text.x=element_text(size = 6),
               axis.text.y=element_text(size = 6, angle = 0),
               axis.title.y =element_text(size = 8 ),
               legend.title = element_text(size = 10))





################
################
AKSL.b<-AKSL%>%
  mutate(Ecomorph=if_else(
    Lake=="Tern" | Lake=="Corcoran" |
      Lake=="Watson" | Lake=="Walby",
    "Benthic", if_else(
      Lake=="South_Rolly"|Lake=="Long"|
        Lake=="Wik"|Lake=="Spirit", "Limnetic",NA)))%>%
  filter(Lake == "Jean" | Lake == "Tern" | 
              Lake == "Watson" | Lake == "Engineer" | 
              Lake == "Engineer" | Lake == "Long" | 
              Lake == "Spirit" | Lake == "Wik" | Lake == "Walby" |
              Lake == "South_Rolly" | Lake == "Finger" | 
              Lake == "Corcoran")

lda.b <- lda(data = AKSL.b, Ecomorph ~ adj.Mass
             + adj.BC + adj.BD + SL_mm + 
               adj.Cped + adj.JawL + adj.SnLength + 
               adj.EyeD + adj.HeadL)
predictions.morph.b <- na.exclude(predict(lda.b, AKSL.b))
AKSL.LDb<-cbind(AKSL.b,predictions.morph.b$x)

#prediction_table gives predicted ecomorph counts for each lake
prediction_table.b<-table(Lake=AKSL.b$Lake, predictions.morph.b$class)

#accuracy rate (correct assignment probability)#
assignment_accuracy.b<-mean(na.omit(
  AKSL.b$Ecomorph == predictions.morph.b$class))
cat("Assignment accuracy when ecomorph known:" )
assignment_accuracy.b%>%round(3)
kable(lda.b$prior%>%round(3), caption="LDA Priors", align = "c")
kable(lda.b$means%>%round(2), caption="LDA Group Means", align = "c")
kable(lda.b$scaling%>%round(2), caption="LDA Trait Scaling", align = "c")
kable(prediction_table.b, caption = "Individual Assignments", align = "c")


plot.b<-ggplot(AKSL.LDb[which(!is.na(AKSL.LDb$LD1)),]%>%
                 mutate(Region = as.factor(ifelse(
                    Lake=="Long" | Lake=="Walby" | 
                      Lake=="South_Rolly" | Lake=="Finger" |
                      Lake == "Corcoran", "Mat-Su", "Kenai")))%>%
         mutate(Ecotype=as.factor(if_else(
           is.na(Ecomorph),"NA",Ecomorph))), aes(x = reorder(
             Lake, LD1, FUN=mean), y = LD1 ))+
  geom_boxplot(aes(fill = Ecotype), outlier.shape = NA)+
  scale_fill_manual(values = c("seagreen4", "skyblue2","coral2"))+
  geom_beeswarm(inherit.aes=FALSE,aes(x = Lake, y = LD1),
               method = 'center', shape = 16,
               size = .8, alpha = .8)+
  scale_colour_manual(values = c("coral4","seagreen4"))+
  labs(x = NULL, y = "LD Score")+
  theme(axis.text.x=element_text(size = 6, angle = 0),
        axis.text.y=element_text(size = 6, angle = 0),
        axis.title.y =element_text(size = 8 ),
        legend.title = element_text(size = 10 ))




################
################
AKSL.c<-AKSL%>%
  mutate(Ecomorph=if_else(
    Lake=="Tern" | Lake=="Finger" |
      Lake=="Watson" | Lake=="Walby",
    "Benthic", if_else(
      Lake=="South_Rolly"|Lake=="Long"|
        Lake=="Wik"|Lake=="Spirit", "Limnetic",NA)))%>%
  filter(Lake == "Tern" | Lake == "Watson" | Lake == "Long" | 
           Lake == "Spirit" | Lake == "Wik" | Lake == "Walby" |
           Lake == "South_Rolly" | Lake == "Finger")

lda.c <- lda(data = AKSL.c, Ecomorph ~ adj.Mass
             + adj.BC + adj.BD + SL_mm + 
               adj.Cped + adj.JawL + adj.SnLength + 
               adj.EyeD + adj.HeadL)
predictions.morph.c <- na.exclude(predict(lda.c, AKSL.c))
AKSL.LDc<-cbind(AKSL.c,predictions.morph.c$x)

#prediction_table gives predicted ecomorph counts for each lake
prediction_table.c<-table(Lake=AKSL.c$Lake, predictions.morph.c$class)

#accuracy rate (correct assignment probability)#
assignment_accuracy.c<-mean(na.omit(
  AKSL.LDc$Ecomorph == predictions.morph.c$class))
cat("Assignment accuracy when ecomorph known:" )
assignment_accuracy%>%round(3)
kable(lda.c$prior%>%round(3), caption="LDA Priors", align = "c")
kable(lda.c$means%>%round(2), caption="LDA Group Means", align = "c")
kable(lda.c$scaling%>%round(2), caption="LDA Trait Scaling", align = "c")
kable(prediction_table.c, caption = "Individual Assignments", align = "c")


plot.c<-ggplot(AKSL.LDc[which(!is.na(AKSL.LDc$LD1)),]%>%
         mutate(Region = as.factor(ifelse(
                    Lake=="Long" | Lake=="Walby" | 
                      Lake=="South_Rolly" | Lake=="Finger",
                      "Mat-Su", "Kenai")))%>%
         mutate(Ecotype=as.factor(if_else(
           is.na(Ecomorph),"NA",Ecomorph))), aes(x = reorder(
             Lake, LD1, FUN=mean), y = LD1, fill = Ecotype))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("seagreen4", "skyblue2","coral2"),
                    guide="none")+
  geom_beeswarm(inherit.aes=FALSE,aes(x = Lake, y = LD1),
                method = 'center', shape = 16,
                size = .8, alpha = .8)+
  scale_colour_manual(values = c("coral4","seagreen4"))+
  labs(x = NULL, 
       y = "LD Score")+
  theme(axis.text.x=element_text(size = 6, angle = 0),
        axis.text.y=element_text(size = 6, angle = 0),
        axis.title.y = element_text(size = 8 ),
        legend.title = element_text(size = 10 ))


plot.a/plot.b/plot.c+ plot_layout(guides = 'collect')


# This is the LDA figure for the paper
plot.lake/plot.b/plot.c+ plot_layout(guides = 'collect')



final.selection<-tibble(Lake=c("Tern","Watson","Finger","Walby",
                   "Spirit","Long","South_Rolly","Wik","G"),
            Ecotype=c("Benthic","Benthic","Benthic","Benthic",
                      "Limnetic","Limnetic","Limnetic","Limnetic",
                      "Transplant Lakes"))

#map
#loads map data

library(vegan)
library(corrplot)
library(Hmisc)
library(pracma)
# following packages required for map
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

eco.dat<-read.csv("/AK_ENV_zoo_invert_2018-2019_v2AmNat.csv")

map.dat<-merge(final.selection,eco.dat, by = "Lake")

world <- ne_countries(scale=10,returnclass = 'sf')
usa <- subset(world, admin == "United States of America")
can <- subset(world, admin == "Canada")
rus <- subset(world, admin == "Russia")
usacanrus<-rbind(usa,can,rus)

#inset map of entire state of AK
alaska<-ggplot(data = usacanrus) +
  geom_sf(fill = "honeydew2", color="grey60", size=.2) +
  panel_border(color = "grey50")+ theme_grey()+
  geom_rect(xmin = -152, xmax = -148.5, ymin = 59.5, ymax = 62, 
            fill = NA, colour = "black", linewidth = .35)+
  coord_sf(xlim = c(-173, -135), 
           ylim = c(55, 73), expand = FALSE, datum = NA)+
  theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))
# Kenai Peninsula local map
alaskalocal <- ggplot(data = usa) + theme_cowplot(font_size = 9)+
  geom_sf(fill = "honeydew2",color="grey50", size=.3) +
  coord_sf(xlim = c(-152, -148.5), ylim = c(59.5, 62), expand = F)
# combines inset and local Kenai map, and adds points and labels
AKlocal_inset<-alaskalocal + annotation_custom(
  grob = ggplotGrob(alaska),
  xmin = -150,
  xmax = -148.3,
  ymin = 59.5,
  ymax = 60.2) +
  geom_point(data=map.dat,aes(
    x=LONGITUDE,y=LATITUDE,color = Ecotype),size = 2)+
  scale_color_manual(values = c(
    "seagreen4", "skyblue2","coral2"),guide="none")+
  geom_text_repel(data = ~ mutate(map.dat, Lake= if_else(
    Lake=="G","Restoration\nArea",Lake)),aes(
    x=LONGITUDE,y=LATITUDE,label = Lake),
    box.padding=.4,point.padding = .1,size = 3,
    fontface = "bold")+labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "none",
        axis.text = element_text(size = 6.5))
ggsave("AK_intro-map.jpg", width = 5, height = 6, dpi = 600)
ggsave("AK_intro-map.pdf", width = 5, height = 6, dpi = 600)

## FIGURE 2
AKlocal_inset 
