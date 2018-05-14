library(nlme)
library(multcomp)
library(MuMIn)
library(vegan)
library(ade4)
library(Hmisc)
library(lsmeans)
library(car)
library(plyr)
library(cowplot)
library(grid)

Field_data = read.csv("Field_Data.csv", sep = '\t')
Plot1m2_data = read.csv('Plot1m2_Data.csv', sep = '\t')
Plot1m2_community = read.csv('Plot1m2_Community.csv', sep = '\t')
Plot1m2_data$LU = factor(Plot1m2_data$LU, c('ULR', 'M', 'YRM', 'OR'))
Plot1m2_community$LU = factor(Plot1m2_community$LU, c('ULR', 'M', 'YRM', 'OR'))
## Topographical data
Slope_per_field = tapply(Plot1m2_data$Slope, Plot1m2_data$Plot, mean) ### check

#### Variations of soil properties with land use ####
# Bulk density
mean(Plot1m2_data$Bulk_density);sd(Plot1m2_data$Bulk_density)
anova(lme(Bulk_density~LU, data = Plot1m2_data, random = ~1|Field))
# Humidity
mean(Plot1m2_data[Plot1m2_data$LU == 'OR',]$Humidity); sd(Plot1m2_data[Plot1m2_data$LU == 'OR',]$Humidity)
mean(Plot1m2_data[Plot1m2_data$LU != 'OR',]$Humidity); sd(Plot1m2_data[Plot1m2_data$LU != 'OR',]$Humidity)
anova(lme(Humidity~LU, data = Plot1m2_data, random = ~1|Field))
# C N
mean(Plot1m2_data$C); sd(Plot1m2_data$C)
mean(Plot1m2_data$N); sd(Plot1m2_data$N)
anova(lme(C~LU, data = Plot1m2_data, random = ~1|Field))
anova(lme(N~LU, data = Plot1m2_data, random = ~1|Field))
# Interactions between soil characteristics
anova(lme(N~LU, data = Plot1m2_data, random = ~1|Field))
anova(lme(Bulk_density~C, data = Plot1m2_data, random = ~1|Field))


#### Weed richness, abundance and composition ####
# Rice fields richer ????
LU_rice_field = ifelse(Field_data$LU == 'ULR', 'ULR','notULR' )
anova(lm(Richness~LU_rice_field, data = Field_data))

Plot1m2_data$LU_rice = ifelse(Plot1m2_data$LU == 'ULR', 'ULR','notULR' )
# Living biomass
summary(lme(Fresh~LU_rice, data = Plot1m2_data, random = ~1|Field))
# Living soilcover
summary(lme(Living_cover~LU_rice, data = Plot1m2_data, random = ~1|Field, na.action = na.omit))
# Litter biomass
mod = lme(Litter~LU, data = Plot1m2_data, random = ~1|Field, na.action = na.omit)
summary(glht(mod, linfct = mcp(LU = "Tukey")))
# Interactions between plant characteristics
mod = lme(Fresh~Living_cover, data = Plot1m2_data, random = ~1|Field, na.action = na.omit)
mod = lme(Living_cover~Abundance, data = Plot1m2_data, random = ~1|Field, na.action = na.omit)
r.squaredGLMM(mod)

## Between Class Analysis
Community =  Plot1m2_community[, tolower(colnames(Plot1m2_community)) %in% tolower(list_herb)]

Community_hellinger = decostand(Community, 'hellinger')
PCA = dudi.pca(Community_hellinger, scannf = FALSE, nf = 4)
BCA = bca(PCA, as.factor(Plot1m2_community$LU), scannf = FALSE, nf = 3)
randtest(BCA, 1000, "two-sided")
# Test significance of variation on 2 first axes
BCA_data = data.frame('CS1' = BCA$ls[,1], 'CS2' = BCA$ls[,2], 'LU'=as.factor(Plot1m2_community$LU), 'Field' = as.factor(Plot1m2_community$Field))
rownames(BCA_data)
mod_cs1 = lme(CS1~LU, data = BCA_data, random = ~1|Field)
mod_cs2 = lme(CS2~LU, data = BCA_data, random = ~1|Field)
summary(glht(mod_cs1, linfct = mcp(LU = "Tukey")))
summary(glht(mod_cs2, linfct = mcp(LU = "Tukey")))


#### Interactions between weeds and soil properties ####
# Pearson correlation coefficients
rcorr(as.matrix(Plot1m2_data[,2:14]))
for (i in c( 'Bulk_density','Humidity', 'C', 'N')){
  for (j in c('Fresh','Richness', 'Litter', 'Living_cover')){
    print('----------------------------------------------')
    print(i)
    print(j)
    Y = Plot1m2_data[,i]
    X = Plot1m2_data[,j]
    LU = factor(Plot1m2_data[,'LU'], levels = c('ULR', 'M', 'YRM', 'OR'))
    Field = Plot1m2_data[,'Field']
    DATA = data.frame(Y,X,LU, Field)
    
    ## Model A: removes effect of land-use
    mod0 = lm(Y~LU, data = DATA, na.action = na.omit)
    res = mod0$residuals
    modA = lme(res ~ X, data = DATA, random = ~1|Field, na.action = na.omit)
    print(summary(modA))
    print(r.squaredGLMM(modA))
    
    ## Model B: interaction with land-use
   #  modB1 = lme(Y~LU*X, data = DATA, random = ~1|Field, na.action = na.omit)
   #  modB2 = lme(Y~LU/X, data = DATA, random = ~1|Field, na.action = na.omit)

   # print(summary(modB))

   #  print(cld(lsmeans(modB1, ~LU), alpha = 0.1))
   #   print(cld(lstrends(modB2, specs = 'LU', var = 'X', alpha = 0.05)))
    #print(r.squaredGLMM(modB1))
    
    
    # plot(residuals(modA)~Plot1m2_data$Plot, main = i, xlab = j)
    # plot(fitted(modA), residuals(modA), xlab = 'Fitted Values', ylab = 'Residuals')
    # abline(h=0, lty=2)
    # lines(smooth.spline(fitted(modA), residuals(modA)))
    
  }}




##### Figures #####

#### Figure 3,4 ####
niceBP = function(dat, var, fact, group, sig, Ylab, main, ypos, ntch = TRUE, symb = TRUE, coul){
  #  if (fact == cluster){
  #    coul =  alpha(coul3, 0.6)
  #    names(coul) = c('A', 'B', 'C')}
  #  else {  coul = c('orange2','steelblue3','palegreen3', 'mediumorchid3')
  #  names(coul) = c('M', 'OR', 'ULR', 'YR')}
  #    print(levels(dat[,fact]))
  #  print(coul)
  # dat[,fact] = droplevels(dat[,fact])
  if (symb == TRUE){
    D = data.frame(fact = dat[, fact], var = dat[, var], group = dat[, fact], supp = gsub('[A-z]+', '',dat[, 'Plot']))
    p<-ggplot(D, aes(x=fact, y=var, fill= fact)) +
      geom_boxplot(outlier.colour = NA, aes(x=fact, y=var), notch = ntch, width=1.2,fill = coul[levels(dat[, fact])]) 
    p = p+   geom_jitter(position=position_dodge(0.6), data = D, aes(x=fact, y=var, shape=supp), cex = 3.5, show_guide = FALSE)
    print(D)
  }
  else {  
    
    D = data.frame(fact = dat[, fact], var = dat[, var], group = dat[, fact])
    p<-ggplot(D, aes(x=fact, y=var, fill= fact)) +
      geom_boxplot(outlier.colour = NA, aes(x=fact, y=var), notch = ntch, width=0.8,fill = coul[levels(dat[, fact])])
    p = p+ geom_jitter(position=position_jitter(width = 0.2), data = D, aes(x=fact, y=var), cex = 3.5, show.legend = FALSE)
  }
  # p <- ggplot(D, aes(fact, var))+ 
  #  geom_boxplot(outlier.shape = 3) +
  #  geom_boxplot(outlier.colour = NA, fill = coul[levels(Com_agg$Crop)])+
  #  geom_point(position = position_jitter(width = 0.2))+
  p = p + 
    ylab(Ylab) +
    xlab(group) +
    theme_bw() +
    theme(axis.title=element_text(size=15,face="bold"), axis.text=element_text(size=13)) 
  # print(dat[,fact])
  p = p + annotate("text", x= 1:length(levels(dat[,fact])), y = rep(ypos, length(levels(dat[,fact]))), label=sig, colour = 'black', cex = 12)
  # plot(p)
  return(p)
}

lm.BP = function(DATA, var, fact, rand = '', colors = couleurs, ymax, YLAB = '', random_effect = TRUE, plim = 0.05, plot_hypo = TRUE, lambda = 1){
  library(DHARMa)
  DATA[,fact] = droplevels(DATA[,fact])
  if (random_effect == TRUE){
    D = data.frame('var' =  DATA[,var], 
                   'var_trans' =  yjPower(DATA[,var], lambda),
                   'fact' = DATA[,fact],
                   'rand' = DATA[,rand])
    print('a')
    mod = lme(var_trans~ fact,data = D,random = ~1|rand, method = "ML", na.action = na.omit)
    print('b')}
  if (random_effect == FALSE){
    D = data.frame('var' =  DATA[,var],
                   'var_trans' =  yjPower(DATA[,var], lambda),
                   'fact' = DATA[,fact])
    mod = lm(var~ fact,data = D, method = "ML") }
  PW  = summary(glht(mod, mcp(fact ="Tukey")))
  print(cld(PW, level = plim))
  let = cld(PW, level = plim)$mcletters$Letters
  par(xpd=TRUE)
  # if (plot_hypo == TRUE){
  #simulationOutput <- simulateResiduals(fittedModel = mod, refit = T)
  #plotSimulatedResiduals(simulationOutput)
  #testUniformity(simulationOutput)
  #testOverdispersion(simulationOutput)
  # }
  BP= niceBP(DATA, var, fact, fact, let, YLAB, '', ymax, ntch = F, coul = colors,symb = F)
  return(list("BP" = BP, "mod" = list('m' = mod, 'v' = var, 'f' = fact)))
}

cbbPalette = c('gray', 'gray','gray', 'gray' )
names(cbbPalette) = c("ULR", "M", "YRM", 'OR')

ric = lm.BP(Plot1m2_data, 'Richness', 'LU', 'Field', colors = cbbPalette, ymax = 15,expression(bold(paste("Number of species / ", m^{2}))), 'TRUE', lambda = 0.3)
pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/richness_crop.pdf", width = 7.5, height = 5)
ric$BP = ric$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                                   labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14)) + xlab('Land use')
plot(ric$BP)
dev.off()

Lit = lm.BP(Plot1m2_data, 'Litter', 'LU', 'Field', colors = cbbPalette, ymax = 920,expression(bold(paste("Litter biomass (g/", m^{2}, ")"))), 'TRUE', lambda = 0.3)
pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/Lit_crop.pdf", width = 7.5, height = 5)

Lit$BP = Lit$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                                   labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14)) + xlab('Land use')
plot(Lit$BP)
dev.off()

Fr= lm.BP(Plot1m2_data, 'Fresh', 'LU', 'Field', colors = cbbPalette, ymax = 300,expression(bold(paste("Living biomass (g/", m^{2}, ")"))), 'TRUE', lambda = 0.3)
pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/Fr_crop.pdf", width = 7.5, height = 5)

Fr$BP = Fr$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                                 labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14))+ xlab('Land use')
plot(Fr$BP)
dev.off()

li= lm.BP(Plot1m2_data[!is.na(Plot1m2_data$Living_cover),], 'Living_cover', 'LU', 'Field', colors = cbbPalette, ymax = 80,expression(bold(paste('Living soil cover (%)'))), 'TRUE', lambda = 0.25)
pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/li_crop.pdf", width = 7.5, height = 5)
li$BP = li$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                                 labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14))+ xlab('Land use')
plot(li$BP)
dev.off()



hum= lm.BP(Plot1m2_data, 'Humidity', 'LU', 'Field', colors = cbbPalette, ymax = 30,'Soil water content (%)', 'TRUE', lambda = 0.25)
pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/hum_crop.pdf", width = 7.5, height = 5)

hum$BP = hum$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                                   labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14))+ xlab('Land use')
plot(hum$BP)
dev.off()


bd= lm.BP(Plot1m2_data, 'Bulk_density', 'LU', 'Field', colors = cbbPalette, ymax = 1.5,expression(bold(paste("Bulk density (g/", cm^{3}, ")"))), 'TRUE', lambda = 0.25)
pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/bd_crop.pdf", width = 7.5, height = 5)

bd$BP = bd$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                                 labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14))+ xlab('Land use')
plot(bd$BP)
dev.off()

c= lm.BP(Plot1m2_data, 'C', 'LU', 'Field', colors = cbbPalette, ymax = 5,'Carbon content (%)', 'TRUE', lambda = 0.25)
pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/c_crop.pdf", width = 7.5, height = 5)
c$BP = c$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                               labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14))+ xlab('Land use')
plot(c$BP)
dev.off()

pdf("/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/n_crop.pdf", width = 7.5, height = 5)
n= lm.BP(Plot1m2_data, 'N', 'LU', 'Field', colors = cbbPalette, ymax = 0.5,'Nitrogen content (%)', 'TRUE', lambda = 0.25)
n$BP = n$BP + scale_x_discrete(breaks=c("ULR", "M", "YRM", 'OR'),
                               labels=c("Rice", "Maize", "Young RT + Maize", 'Mature RT'))+
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=14))+ xlab('Land use')
plot(n$BP)
dev.off()

#### Figure 5 ####
BCA_data$CS1 = -BCA_data$CS1
BCA_data$CS2 = -BCA_data$CS2
BCA_data$LU = factor(BCA_data$LU, levels = c('ULR', 'M', 'YRM','OR'))

ggbca = ggplot(BCA_data, aes(CS1, CS2, shape = LU, linetype = LU)) +  
  scale_shape_manual(name="Land use", values = c(1,2,19,17), breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT'), guide = FALSE) +
  scale_linetype_manual(name="Land use", values = c("solid",'dashed','dotted','twodash'), breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT')) +
  scale_fill_manual(values=cbbPalette) +
  geom_point(aes(x=CS1, y=CS2)) +
  stat_ellipse(level=0.95)+
  theme(#legend.position="none", 
    panel.background=element_blank(),panel.border=element_blank(),
    plot.background=element_blank()) +
  xlab('CS1 (7.9%)') + ylab('CS2 (7.1%)') 

save_plot('/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/BCA_leg2.pdf',ggbca, base_aspect_ratio = 1.85 )


bca_pairtest = cbind(-BCA$ls[,1:3],Plot1m2_community$LU)
colnames(bca_pairtest)[4] = 'LU'
mod1 = lm(CS1~LU, data = bca_pairtest)
mod2 = lm(CS2~LU, data = bca_pairtest)
#cld(lsmeans(mod1, 'Crop'))
cs1_means <- ddply(bca_pairtest, "LU", summarise, cs1.mean=mean(CS1))
cs1 = cs1_means[cs1_means$LU %in% c('ULR', 'M', 'OR'),'cs1.mean']#+ c(0, -0.1,0)

cs2_means <- ddply(bca_pairtest, "LU", summarise, cs2.mean=mean(CS2))
cs2 = cs2_means[cs2_means$LU %in% c('ULR', 'M', 'OR'),'cs2.mean']#+ c(0, -0.1,0)
cs2_means

#plot.cs1 =
ggplot(bca_pairtest, aes(x=CS1,  linetype = LU)) + geom_density(alpha = 0.1) +
  #geom_segment(data=cs1_means, aes(x=cs1.mean, y=0, xend=cs1.mean, yend=0.75), linetype = 'solid')+
#  geom_histogram(aes(x=CS1, fill = LU))+
 # geom_segment(data=cs1_means, aes(x=cs1.mean[3]-0.06, y=0.8, xend=cs1.mean[4]+0.06, yend=0.8), color = 'black' ,linetype = 'solid', alpha=0.95)+
#  annotate('text',label = c('a', 'b'), x = c(-2.7,1), y = 1, cex = 7)+ ylab('Density (CS1)') +
#  ylim(c(0, 1.1))+
#  scale_linetype_manual(name="Land use", values = c("solid",'dashed','dotted','twodash'), breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT')) +
  theme(axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),
        plot.background=element_blank())+
  theme(legend.position="none") + xlab('') 

plot.cs2 = ggplot(bca_pairtest, aes(x=CS2, linetype=LU)) + geom_density(alpha = 0.1) +
  ylim(c(0, 1.1))+
  scale_linetype_manual(name="Land use", values = c("solid",'dashed','dotted','twodash'), breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT')) +
  
  geom_segment(data=cs2_means, aes(x=cs2.mean, y=0, xend=cs2.mean, yend=0.75), linetype = 'solid')+
#  geom_segment(data=cs2_means, aes(x=cs2.mean[1]+0.06, y=0.8, xend=cs2.mean[2]-0.06, yend=0.8,linetype = 'solid')
#             , color = 'black', alpha=0.95) +
 # geom_segment(data=cs2_means, aes(x=cs2.mean[2]+0.06, y=0.8, xend=cs2.mean[3]-0.06, yend=0.8,linetype = 'solid')
   #            , color = 'black', alpha=0.95) +
  theme(axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),
        plot.background=element_blank(), legend.position = c(0.4, 1.2))+
  coord_flip()  + xlab('') + ylab('Density (CS2)') +
  annotate('text',label = c('a', 'b'), x = c(1,-0.3), y = 0.95, cex = 7)

a = ggdraw() +
  draw_plot(ggbca,     0,    0, .75, .75) +
  draw_plot(plot.cs1 , 0.05, 0.75,  .70, .25) +
  draw_plot(plot.cs2, .75,   0.05,  .25, .70)
a
save_plot('/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/BCA_plot.pdf',a, base_aspect_ratio = 1.85 )



circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}



tot_names = c('Con.sum', 'Lyg.fle', 'The.sub' ,'Pen.pol', 'Pas.con', 
              'Chr.odo', 'Blu.lac', 'Sel.hel' ,'Cra.cre' ,'Thy.lat',
              'Dig.rad' ,'Eup.hir' ,'Man.uti' ,'Ele.ind', 'Age.con',
              'Cyn.dac' ,'Spi.pan', 'Mel.rep' ,'Cen.asi' )
XLAB = (BCA$co[tot_names,]$Comp1)
YLAB = (BCA$co[tot_names,]$Comp2)    
l_seg = sqrt(XLAB^2 + YLAB^2) + 0.05
XLAB_prim = XLAB/l_seg 
YLAB_prim = YLAB/l_seg        

XLAB_lab =  XLAB_prim   # + c(0,      -0.0,       -0.21,     0.0,     0,
                            -0,       0.1 ,       0.15,      +0.1,     0.23,
                            0.15,     -0.02,       -0.13,     0.1,         0, 
                            -0.12,    0.0,       -0.1,    - 0.12)
YLAB_lab =  YLAB_prim  # + c(0,      - 0.05,       0.13,      -0.05,   -0.1,
                          -0,        -0.04 ,      0.12,      0.02,    0.12, 
                          -0.1,      -0.1,       0.21,     0.02,        0,    
                          0.15,     -0.03,    0.1,     0.1)
circle <- circleFun(c(0,0),2,npoints = 100)

      #  ggcircle =
          ggplot(BCA$co[tot_names,], aes(x=Comp1, y=-Comp2)) +
  #annotation_custom(grob=circleGrob(r=unit(1,"npc")), xmin=-1, xmax=1, ymin=-1, ymax=1) +
  geom_segment(data=BCA$co[tot_names,], aes(x=0, y=0, xend=Comp1, yend=-Comp2),
               colour = 'black', arrow=arrow(length=unit(0.2,"cm")), alpha=1)+
  geom_segment(data=BCA$co[tot_names,], aes(x=0, y=0, xend=XLAB_prim, yend=-YLAB_prim),
               colour = 'black', alpha=0.25)+
  ylim(c(-1,1)) + xlim(c(-1,1)) + theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  annotate('text', label = rownames(BCA$co[tot_names,]), x = XLAB_lab,y = -YLAB_lab, cex = 7)+
  geom_path(data = circle, aes(x=x,y=y), color = 'lightgray') + xlab('CS1 (7.8%)') + ylab('CS2 (6.9%)') +
            xlim(c(-1.5,1)) + ylim(c(-1, 1.5))

pdf('/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/circle.pdf', width = 4, height = 4)
plot(ggcircle)
dev.off()


#### Figure 6 ####
trend_plot = function(data, Y, X, Crop_to_draw, LP, xtext, ytext, angl){
  dat = data[, c(Y, X, 'LU')]
  colnames(dat) = c('Y', 'X', 'LU')
  gg = ggplot(data = dat, aes(y=Y, x=X,colour=LU, fill = LU)) + 
    theme_bw() +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5))+
    xlab('') + ylab('') + 
    scale_color_manual(name="Land use", values=c('black', 'black', 'black', 'black'), breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT'), guide = FALSE) + 
    scale_fill_manual(name="Land use", values=c('black', 'black', 'black', 'black'), breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT'), guide = FALSE) +
    #scale_linetype_manual(name="Land use", values=LP, breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT'), guide = FALSE) +
    geom_point(size = 0.65, aes(shape = LU), alpha = 0.8) +
    scale_shape_manual(name="Land use", values = c(1,2,19,17), breaks=c("ULR", "M", "YRM", 'OR'), labels = c("Rice", "Maize", "Young RT + Maize", 'Mature RT'), guide = FALSE)
    
  Crop_to_write = recode(Crop_to_draw,"'ULR'='Rice';'M'='Maize';'OR'='Mature RT';'YR'='Young RT'")
  
  if (!is.na(Crop_to_draw)){
    for (i in 1:length(Crop_to_draw)){
      gg = gg  + geom_smooth(data=subset(dat, LU == Crop_to_draw[i]),method='lm', linetype = LP[i], alpha = 0.2) +
        # geom_smooth(data = dat2, aes(y=Y2, x=X2), )#, linetype = LP[i])  
        annotate('text', label = Crop_to_write[i], y = ytext[i], x = xtext[i], cex = 3, angle = angl[i])  
    }}
  return(gg)
}
hum_li = trend_plot(Plot1m2_data[!is.na(Plot1m2_data$Living_cover),], 'Humidity', 'Living_cover', c('OR', 'ULR'), c('dashed', 'solid'),  c(40, 64),c(17,10), c(0,22))
hum_ric = trend_plot(Plot1m2_data, 'Humidity', 'Richness', c('M', 'OR'), c('dashed', 'solid'),  c(9, 6),c(5.8,15.8), c(0,-10))
hum_lit = trend_plot(Plot1m2_data, 'Humidity', 'Litter', c('OR'), c( 'solid'),  c(500),c(18), c(12))
hum_fr = trend_plot(Plot1m2_data, 'Humidity', 'Fresh', c('ULR'), c( 'solid'),  c(220),c(10), c(15))
bd_fr = trend_plot(Plot1m2_data, 'Bulk_density', 'Fresh', NA, NA,  NA,NA,NA)
bd_ric = trend_plot(Plot1m2_data, 'Bulk_density', 'Richness', NA, NA,  NA,NA,NA)
bd_lit = trend_plot(Plot1m2_data, 'Bulk_density', 'Litter', NA, NA,  NA,NA,NA)
bd_li = trend_plot(Plot1m2_data[!is.na(Plot1m2_data$Living_cover),], 'Bulk_density', 'Living_cover', NA,  NA,NA,NA)
n_fr = trend_plot(Plot1m2_data, 'N', 'Fresh', 'ULR', 'dashed',  245,0.21,8)
n_ric = trend_plot(Plot1m2_data, 'N', 'Richness', 'OR', 'solid',  10,0.105,0)
n_lit = trend_plot(Plot1m2_data, 'N', 'Litter', NA, NA,  NA,NA,NA)
n_li = trend_plot(Plot1m2_data[!is.na(Plot1m2_data$Living_cover),], 'N', 'Living_cover', NA,  NA,NA,NA)
c_fr = trend_plot(Plot1m2_data, 'N', 'Fresh', NA, NA,  NA,NA,NA)
c_ric = trend_plot(Plot1m2_data, 'N', 'Richness', NA, NA,  NA,NA,NA)
c_lit = trend_plot(Plot1m2_data, 'N', 'Litter', NA, NA,  NA,NA,NA)
c_li = trend_plot(Plot1m2_data[!is.na(Plot1m2_data$Living_cover),], 'N', 'Living_cover', NA,  NA,NA,NA)

library(cowplot)

tot_plot = ggdraw() + 
  draw_plot(bd_ric, 0.1, 0.7, 0.2, .2) +
  draw_plot(bd_lit, 0.3, 0.7, 0.2, .2) +
  draw_plot(bd_fr, 0.5, 0.7, 0.2, .2) +
  draw_plot(bd_li, 0.7, 0.7, 0.2, .2) +
  draw_plot(hum_ric, 0.1, 0.5, 0.2, .2) +
  draw_plot(hum_lit, 0.3, 0.5, 0.2, .2) +
  draw_plot(hum_fr, 0.5, 0.5, 0.2, .2) +
  draw_plot(hum_li, 0.7, 0.5, 0.2, .2) +
  draw_plot(n_ric, 0.1, 0.3, 0.2, .2) +
  draw_plot(n_lit, 0.3, 0.3, 0.2, .2) +
  draw_plot(n_fr, 0.5, 0.3, 0.2, .2) +
  draw_plot(n_li, 0.7, 0.3, 0.2, .2) +         
  draw_plot(c_ric, 0.1, 0.1, 0.2, .2) +
  draw_plot(c_lit, 0.3, 0.1, 0.2, .2) +
  draw_plot(c_fr, 0.5, 0.1, 0.2, .2) +
  draw_plot(c_li, 0.7, 0.1, 0.2, .2) +
  annotate('text', label = c('C content (%)','N content (%)','Water content (%)', 'Bulk density (g.cm  )',  
                             'Species richness', 'Litter biomass (g.m  )', "Living biomass (g.m  )", 'Living soil cover (%)'),
           x = c(0.1,0.1,0.1,0.1,0.22,0.42,0.62,0.82), y = c(0.22,0.42,0.62,0.82, 0.1,0.1,0.1,0.1),
           angle = c(90,90,90,90, 0,0,0,0), cex = 5)  +
  draw_plot_label(tolower(LETTERS)[1:16], x= 0.02+ rep(c(0.1, 0.3, 0.5, 0.7), 4), y= 1.02- rep(c(0.1, 0.3, 0.5, 0.7), each= 4), size = 15)
pdf('/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/tot_plot.pdf', height = 10, width = 10)
plot(tot_plot)
dev.off()





##### Fig S1 #####

meteo = read.table('monthly_meteo2015-2016.csv', sep = '\t', head = TRUE)
meteo$month = 1:22
meteo$Season = c('Dry season', rep('Rainy season', 8), rep('Dry season',4), rep('sainy Season',8), 'Dry season')

rf = ggplot(meteo, aes(y = Rainfall, x = month))  + theme_bw()+
  #  color = c(rep('black', 12), 'red', rep('black', 9)))+
  labs(x = "", y = 'Rainfall (mm)') +
  #scale_fill_manual(labels = c('Dry season', 'Rainy season'), values = c('gray50', 'gray10'))+
  # color = c(rep('black', 12), 'red', rep('black', 9)))+
  scale_x_continuous(breaks = meteo$month[seq(1,22, 2)] , labels = meteo$X[seq(1,22, 2)]) +
  annotate("rect", xmin = 2, xmax = 9, ymin = -5, ymax = 420,  alpha = .8,fill = 'gray90') +
  annotate("rect", xmin = 14, xmax = 21, ymin = -5, ymax = 420,  alpha = .8,fill = 'gray90') +
  geom_vline(aes(xintercept=meteo$month),linetype="dotted", 
             size = c(rep(0.3, 12), 1, rep(0.3, 9)))+ 
  geom_bar(stat="identity") +
  annotate("rect", xmin = 4, xmax = 8, ymin = -35, ymax = -5, alpha = 1,fill = 'white') +
  annotate("rect", xmin = 16, xmax = 20, ymin = -35, ymax = -5, alpha = 1,fill = 'white') #+
  #annotate("text", x = c(6,18),y = c(-25, -25),  label = c('2015', '2016'))



rh = ggplot(meteo, aes(y = RH, x = month)) +  
  annotate("rect", xmin = 2, xmax = 9, ymin = 25, ymax = 100,  alpha = .8, fill = 'gray90') +
  annotate("rect", xmin = 14, xmax = 21, ymin = 25, ymax = 100,  alpha = .8,fill = 'gray90') +
  geom_line() +
  geom_vline(aes(xintercept=meteo$month),linetype="dotted", 
             size = c(rep(0.3, 12), 1, rep(0.3, 9))) +
  scale_x_continuous(breaks = meteo$month[seq(1,22, 2)] , labels = meteo$X[seq(1,22, 2)])  +
  labs(x = "", y = 'Relative air humidity (%)') + theme_bw()  

rh

rdens = ggplot(meteo, aes(y = Rdens, x = month))  +
  annotate("rect", xmin = 2, xmax = 9, ymin = 0, ymax = 0.4,  alpha = .8,fill = 'gray90') +
  annotate("rect", xmin = 14, xmax = 21, ymin = 0, ymax = 0.4, alpha = .8,fill = 'gray90') +
  geom_vline(aes(xintercept=meteo$month),linetype="dotted", 
             size = c(rep(0.3, 12), 1, rep(0.3, 9)))+
  scale_x_continuous(breaks = meteo$month[seq(1,22, 2)] , labels = meteo$X[seq(1,22, 2)])  +
  labs(x = "", y = 'Rainfall erosivity') + theme_bw() +
  geom_bar(stat="identity") 

library(ggplot2)
library(grid)
library(gridExtra)

pdf('/Users/Margot/Desktop/Projet_M2/Documents/Papers/Paper1_SUM/Figures/meteo.pdf', width = 9, height = )
grid.arrange(rf, rh, rdens, ncol = 1)
dev.off()

##### Fig S2 #####

colnames(Plot1m2_community)
names_spe = names(sort(colSums(Plot1m2_community[,1:44]), decreasing = TRUE))[1:9]
names_spe = names_spe[names_spe !=  "Not_identified"]
Field_com = aggregate(Plot1m2_community[,names_spe ], list(Plot1m2_community$Field), mean )
Field_com$LU = factor(rep(c('M', 'OR', 'ULR', 'YR'), each = 5), levels = c('ULR', 'M', 'YR', 'OR'))
levels(Field_com$LU) <- c("Upland rice", "Maize", "Young RT + Maize", 'Mature RT')
Field_com_melt = melt(Field_com, id.vars = c('Group.1', 'LU'))
ggplot(Field_com_melt, aes(x = variable, y = value+1)) + geom_boxplot()+theme_bw()+
#  geom_jitter(position=position_dodge(0.6)) +
  facet_wrap(~LU, ncol = 2) + scale_y_log10() + xlab('Species') +
  ylab(bquote(Plants~per~m^2))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

