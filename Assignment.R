library('corrr')
library("corrplot")
library(rquery)
library(stats)
library(car)
library(rgl)
library(MASS)
library(DescTools)
library(dplyr)
library(fuzzySim)
library(vegan)
library(gmodels)
library(tidyr)
library(ggplot2)
library(bestNormalize)
library(Hmisc)
library(caret)
library('betapart')
library(labdsv)
library(ggvegan)
library(ggrepel)
library(egg)
library(devtools)
library(pairwiseAdonis)


inselb_pred <- read.delim("~/Documents/Belgium Master/Third semester/Analysis of biological data/Analysis-of-Biological-Data---002905-2020-Oct-22_16-13-07/viewer/files/Inselberg_predictors.txt", header=TRUE)

inselb_sp <- read.delim("~/Documents/Belgium Master/Third semester/Analysis of biological data/Analysis-of-Biological-Data---002905-2020-Oct-22_16-13-07/viewer/files/Inselberg_species_composition.txt", header=TRUE)

#pres_abs = 1*(inselb_sp[,4:length(inselb_sp)] > 0)

pres_abs = decostand(inselb_sp[,4:length(inselb_sp)], "pa")

Speciesrichness = rowSums(pres_abs)

inselb_sp = cbind(inselb_sp, Speciesrichness)

data = cbind(inselb_pred, inselb_sp[,4:length(inselb_sp)])

arc_sp = cbind(Archipelago=data[,1], inselb_sp[,4:(length(inselb_sp)-1)])

#arc_sp_group = arc_sp %>% 
#  group_by(Archipelago) %>% 
#  summarise_all(sum)

#arc_sp_group_pa = cbind(arc_sp_group[,1],(decostand(arc_sp_group[,2:length(arc_sp_group)], 'pa')))

#arc_sp_group_ord = arc_sp_group %>% 
 #   pivot_longer(!Archipelago, names_to= 'Species', values_to = 'Abundance')

# separating env_variables

env_data = inselb_pred[,6:(length(inselb_pred)-4)]

env_data = env_data[-c(25)]

colnames(env_data) = c('Annual_Mean_Temperature','Mean_Diurnal_Range','Isothermality','Temperature_Seasonality','Max_Temperature_of_Warmest Month','Min_Temperature_of_Coldest_Month','Temperature_Annual_Range','Mean_Temperature_of_Wettest_Quarter','Mean_Temperature_of_Driest_Quarter','Mean_Temperature_of_Warmest_Quarter','Mean_Temperature_of_Coldest_Quarter','Annual_Precipitation','Precipitation_of_Wettest-Month','Precipitation_of_Driest_Month','Precipitation_Seasonality_(Coefficient_of_Variation)','Precipitation_of_Wettest_Quarter','Precipitation_of_Driest_Quarter','Precipitation_of_Warmest_Quarter','Precipitation_of_Coldest_Quarter',"pet","aridity","srad",'altitude',"area","isolation1","isolation5","isolation15","dist_coast")

climate_var = env_data[,1:22]
regiona_var = env_data[,23:length(env_data)]
regiona_var = cbind(regiona_var, log_area=log(regiona_var$area))
# remember log_area

cor_clim = cor(climate_var, method = 's')
corrplot(cor_clim, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.4)

# To know the correlation value bt each pair of variables 
clim_cor_var = climate_var %>% 
  correlate(method='s') %>%
  shave() %>%
  stretch() %>% 
  arrange(r) %>%
  na.exclude(r)

clim_cor_var

#write.csv(clim_cor_var,'clim_cor_var_new.csv')

#Removing variables - Our cutoff 0.75
rem_clim_index = findCorrelation(cor_clim, cutoff=0.75) # putt any value as a "cutoff" 
rem_clim_index = sort(rem_clim_index)

#reduced clim variables
red_clim_var = climate_var[,-c(rem_clim_index)]

# The same procedure for the regionalism variables
cor_regi = cor(regiona_var, method = 's')
corrplot(cor_regi, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.7)

# To know the correlation value bt each pair of variables 
region_cor_var = regiona_var %>% 
  correlate(method='s') %>%
  shave() %>%
  stretch() %>% 
  arrange(r) %>%
  na.exclude(r)

#write.csv(region_cor_var,'region_cor_var_new.csv')

#Removing variables - Our cutoff 0.75
rem_reg_index = findCorrelation(cor_regi, cutoff=0.75) # putt any value as a "cutoff" 
rem_reg_index = sort(rem_reg_index)

#reduced reg variables
red_reg_var = regiona_var[,-c(rem_reg_index)]

# 6 variables for climate, 4 for regionalism

### SPECIES MATRIX

# To test regionalism we will use the differences between species diversity in the different archipelagos 'beta diversity’ is applied in a broad sense to any measure of variation in species composition

# We need to separate the archipelagos and to assign P and A, also we will include the name of the inselbergs for each row

#A1 = arc_sp %>% 
#  filter(Archipelago == 'A1') 
#row.names(A1) = as.matrix(data %>% 
#  filter(Archipelago =='A1') %>% 
#    select(inselberg))
#A1 = A1[,2:length(A1)]

#A2 = arc_sp %>% 
#  filter(Archipelago == 'A2') 
#row.names(A2) = as.matrix(data %>% 
#         filter(Archipelago =='A2') %>% 
#         select(inselberg))
#A2 = A2[,2:length(A2)]

#A3 = arc_sp %>% 
#  filter(Archipelago == 'A3') 
#row.names(A3) = as.matrix(data %>% 
#        filter(Archipelago =='A3') %>% 
#        select(inselberg))
#A3 = A3[,2:length(A3)]

#A4 = arc_sp %>% 
#  filter(Archipelago == 'A4') 
#row.names(A4) = as.matrix(data %>% 
#  filter(Archipelago =='A4') %>% 
#  select(inselberg))
#A4 = A4[,2:length(A4)]

# get betapart object for beta diversity
# read this paper http://webspersoais.usc.es/export9/sites/persoais/persoais/andres.baselga/pdfs/Baselga-Orme2012.pdf
# Explains why we are applying this, we want to make a difference in the reasons of biodiversity changes 

#pa data
#A1.pa = decostand(A1, "pa")
#A2.pa = decostand(A2, "pa")
#A3.pa = decostand(A3, "pa")
#A4.pa = decostand(A4, "pa")

#A1.core <- betapart.core(A1.pa)
#A2.core <- betapart.core(A2.pa)
#A3.core <- betapart.core(A3.pa)
#A4.core <- betapart.core(A4.pa)

# multiple site measures - Sorensen index 
#A1.multi <- beta.multi(A1.core)
#A2.multi <- beta.multi(A2.core)
#A3.multi <- beta.multi(A3.core)
#A4.multi <- beta.multi(A4.core)

#A1.abund = beta.multi.abund(A1, index.family="bray")
#A2.abund = beta.multi.abund(A2, index.family="bray")
#A3.abund = beta.multi.abund(A3, index.family="bray")
#A4.abund = beta.multi.abund(A4, index.family="bray")

# All species
# You need to reduce the variation (transformation)
# Use square root or proportions to minimize influence of most abundant groups
# This part is the important

all_sp   = inselb_sp[,4:(length(inselb_sp)-1)]
row.names(all_sp) = inselb_sp$inselberg
all_sp_sqrt   = sqrt(all_sp) #square root transformation

all_sp_pa= decostand(all_sp, 'pa')

# all variation - Bray curtis with separated index 
# You obtained 3 matrix with different index
# - Bray-Curtis dissimilarity (abundance weighted)
all.abund.matrix = beta.pair.abund(all_sp_sqrt, index.family="bray")

# Abundance-based multiple-site dissimilarities - Bray curtis
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12693

#PERMANOVA Non-paramentric, based on dissimilarities. Allows for partitioning of variability, similar to ANOVA, allowing for complex design (multiple factors, nested design, interactions, covariates). Uses permutation to compute F-statistic (pseudo-F).

set.seed(1609) #reproducible ;)

#matrix total 
insel.div = adonis2(all.abund.matrix[["beta.bray"]]~Archipelago, data=data, permutations = 999, method="bray")
insel.div

pairwise.adonis2(all.abund.matrix[["beta.bray"]]~Archipelago, data=data, permutations = 999, method="bray")

#matrix gra # no significant
insel.div1 = adonis2(all.abund.matrix[["beta.bray.gra"]]~Archipelago, data=data, permutations = 999, method="bray")
insel.div1



#matrix bal (balanced variation)
insel.div2 = adonis2(all.abund.matrix[["beta.bray.bal"]]~Archipelago, data=data, permutations = 999, method="bray")
insel.div2

# significant p - value of 0.001

### Post hoc test permanova ??? 

# Now multivariate dispersion - The average distance to group centroid. Used as a measure of multivariate beta diversity.

dispersion1<-betadisper(all.abund.matrix[["beta.bray"]], group=data$Archipelago)
permutest(dispersion1)

# https://rstudio-pubs-static.s3.amazonaws.com/246172_1930ddfb5f064b2bab54b11016ab407e.html

# matrix total
source("myplotbetadisp.r")

insel = inselb_sp$inselberg
col.fill.rect <- apply(col2rgb(c("black", 'red','green','blue')), 2, addAlpha, alpha=0.3)
col.text.rect <- apply(col2rgb(c("black")), 2, addAlpha, alpha=0.9)
transp.centroids <- 1

myplotbetadisper(dispersion1, ellipse = TRUE, hull = FALSE, fillrect=col.fill.rect, coltextrect=col.text.rect, labPoints=insel,main= '(a) Total dissimilarity', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n',poslabPoints = 4, alphaPoints=0)

dispersion1.HSD = TukeyHSD(dispersion1)


# matrix grad

dispersion2<-betadisper(all.abund.matrix[["beta.bray.gra"]], group=data$Archipelago)
permutest(dispersion2)

plot(dispersion2, hull=FALSE, ellipse=TRUE, main = '(b) Dissimilarity due to abundance gradients', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n')##sd ellipse()

myplotbetadisper(dispersion2, ellipse = TRUE, hull = FALSE, fillrect=col.fill.rect, coltextrect=col.text.rect, labPoints=insel,main= '(b) Dissimilarity due to abundance gradients', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n',poslabPoints = 4, alphaPoints=0)

dispersion2.HSD = TukeyHSD(dispersion2)


# matrix balanced var
dispersion3<-betadisper(all.abund.matrix[["beta.bray.bal"]], group=data$Archipelago)
permutest(dispersion3)

plot(dispersion3, hull=FALSE, ellipse=TRUE, main = '(c) Dissimilarity due to balanced variation', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n')##sd ellipse()

myplotbetadisper(dispersion3, ellipse = TRUE, hull = FALSE, fillrect=col.fill.rect, coltextrect=col.text.rect, labPoints=insel,main= '(c) Dissimilarity due to balanced variation', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n',poslabPoints = 4, alphaPoints=0)

dispersion3.HSD = TukeyHSD(dispersion3)

par(mfrow=c(1,3))

myplotbetadisper(dispersion1, ellipse = TRUE, hull = FALSE, fillrect=col.fill.rect, coltextrect=col.text.rect, labPoints=insel,main= '(a) Total dissimilarity', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n',poslabPoints = 4, alphaPoints=0)

myplotbetadisper(dispersion2, ellipse = TRUE, hull = FALSE, fillrect=col.fill.rect, coltextrect=col.text.rect, labPoints=insel,main= '(b) Dissimilarity due to abundance gradients', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n',poslabPoints = 4, alphaPoints=0)

myplotbetadisper(dispersion3, ellipse = TRUE, hull = FALSE, fillrect=col.fill.rect, coltextrect=col.text.rect, labPoints=insel,main= '(c) Dissimilarity due to balanced variation', xlab = '', ylab = '', sub= '', labels=FALSE, xaxt= 'n', yaxt= 'n',poslabPoints = 4, alphaPoints=0)
#  RDA 

## ABUNDANCES 
# hellinger transformation

all_sp  = hellinger(all_sp)
row.names(red_clim_var) = inselb_sp$inselberg
row.names(red_reg_var)  = inselb_sp$inselberg

# Clim_var

red_clim_var = as.data.frame(scale(red_clim_var))

colnames(red_clim_var) = c("temp_seas", "mtemp_warmest_quarter","pp_wettest_quarter","pp_driest_quarter","pp_warmest_quarter", "srad")

rda_tree1 = rda(all_sp ~ . , data=red_clim_var)
rda_tree1

# The output above provides us with some useful information. Inertia is another name for variation or variance in this case. “Total” refers to total variance, “Constrained” refers to the amount of variance explained by the explanatory variables, “Unconstrained” refers to the residual variance.

RsquareAdj(rda_tree1)

aplot(rda_tree1, type='n', scaling=1)
orditorp(rda_tree1, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_tree1, display='cn', col='red')

autoplot(rda_tree1, geom = c("point",'text', 'colour') , data= inselb_pred, colour='Archipelago', legend='none', layers = c("species",'biplot'), arrows = FALSE, title = 'Ordination RDA - Climatic variables') +theme_bw()+ theme(legend.position='none') 

smry <- summary(rda_tree1)
df1  <- data.frame(smry$species[,1:2])  

numbers_sp = paste(rep("sp",110), seq(1,110,1), sep = "_")
row.names(df1) = numbers_sp

df2  <- data.frame(smry$biplot[,1:2])     # loadings for PC1 and PC2
rda.plot <- ggplot(df1, aes(x=RDA1, y=RDA2)) + 
  geom_text(aes(label=rownames(df1)),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()
rda.plot

rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2), color="darkred", arrow=arrow(length=unit(0.01,"npc"))) + geom_text_repel(data=df2, aes(x=RDA1,y=RDA2,label=rownames(df2), hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), color="darkred", size=4)+ xlab('RDA1 (12.09%)') + ylab('RDA2 (6.95%)') + annotate('label', label ='italic(Adj.R)^2 == 0.24',x = 0.75, y=0.45, parse=T) + theme_bw() + ggtitle('(a) Climatic variables')

clim = rda.biplot + theme(plot.title = element_text(color="black", size=12, face='bold', hjust = 0.5))


anova(rda_tree1, by = 'term')

# Regionalism 

red_reg_var = as.data.frame(scale(red_reg_var))
rda_tree2 = rda(all_sp ~ . , data=c(red_reg_var))
rda_tree2

RsquareAdj(rda_tree2)

plot(rda_tree2, type='n', scaling=1)
orditorp(rda_tree2, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_tree2, display='cn', col='red')

smry1 <- summary(rda_tree2)
df3  <- data.frame(smry1$species[,1:2])  

#numbers_sp = paste(rep("sp",110), seq(1,110,1), sep = "_")
row.names(df3) = numbers_sp

df4  <- data.frame(smry1$biplot[,1:2])     # loadings for PC1 and PC2
rda.plot1 <- ggplot(df3, aes(x=RDA1, y=RDA2)) + 
  geom_text(aes(label=rownames(df3)),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()
rda.plot1

rda.biplot1 <- rda.plot1 +
  geom_segment(data=df4, aes(x=0, xend=RDA1, y=0, yend=RDA2), color="darkred", arrow=arrow(length=unit(0.01,"npc"))) + geom_text_repel(data=df4, aes(x=RDA1,y=RDA2,label=rownames(df4), hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), color="darkred", size=4)+ xlab('RDA1 (9.21%)') + ylab('RDA2 (6.00%)') + annotate('label', label ='italic(Adj.R)^2 == 0.10',x = 0.2, y=-0.9, parse=T) + theme_bw() + ggtitle('(b) Topographic variables')


reg = rda.biplot1 + theme(plot.title = element_text(color="black", size=12, face='bold', hjust = 0.5))

ggarrange(clim, reg,
          ncol = 2, nrow = 1)

anova(rda_tree2, by = 'term')

# We interpret the plot above as we have interpreted the previously ordination plots with one important difference. The environmental variables are now displayed and their placement indicates their loading on the two displayed RDA axes. elev is loading heavily on RDA1 indicating that this variable explains a larger portion of the variance associated with axis 1. The location of the species relative to the environmental variables indicates how strongly a species is associated with a particular environmental variable

#  RDA 

## PRESENCE and ABSENCE  
# jaccard p/a

all_sp_pa  = vegdist(all_sp_pa, method='jaccard')
row.names(red_clim_var) = inselb_sp$inselberg
row.names(red_reg_var)  = inselb_sp$inselberg


# Clim_var

red_clim_var = as.data.frame(scale(red_clim_var))

rda_tree3 = rda(all_sp_pa ~ . , data=red_clim_var)
rda_tree3

# The output above provides us with some useful information. Inertia is another name for variation or variance in this case. “Total” refers to total variance, “Constrained” refers to the amount of variance explained by the explanatory variables, “Unconstrained” refers to the residual variance.

RsquareAdj(rda_tree1)

plot(rda_tree1, type='n', scaling=1)
orditorp(rda_tree1, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_tree1, display='cn', col='red')

anova(rda_tree1, by = 'term')

# Regionalism 

red_reg_var = as.data.frame(scale(red_reg_var))
rda_tree2 = rda(all_sp ~ . , data=red_reg_var)
rda_tree2

RsquareAdj(rda_tree2)

plot(rda_tree2, type='n', scaling=1)
orditorp(rda_tree2, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_tree2, display='cn', col='red')

anova(rda_tree2, by = 'term')

# We interpret the plot above as we have interpreted the previously ordination plots with one important difference. The environmental variables are now displayed and their placement indicates their loading on the two displayed RDA axes. elev is loading heavily on RDA1 indicating that this variable explains a larger portion of the variance associated with axis 1. The location of the species relative to the environmental variables indicates how strongly a species is associated with a particular environmental variable



