##   R Core Team (2015). R: A language and environment for
##   statistical computing. R Foundation for Statistical Computing,
##   Vienna, Austria. URL http://www.R-project.org/.

## load packages ##
require(readxl)
require(ggplot2)
require(ggpmisc)
require(MCMCglmm)
require(ape)
require(phytools)
require(cowplot)

## Import data
#dataframe <- read_excel("~/Dropbox/Skolenotater/Uni/5. År/Q3+Q4/Niels Chr Speciale/dataframe.xlsx" ,  sheet = "Data")
#setwd("Dropbox/Skolenotater/Uni/5. År/Q3+Q4/Niels Chr Speciale")
setwd(dir = "/Users/au231308/Desktop/")
dataframe <- read_excel(path = "./dataframe.xlsx" ,  sheet = "Data")
dataframe$airbreathing = as.factor(dataframe$airbreathing)
dataframe$species = as.factor(dataframe$species)
dataframe$fishid  = as.factor(dataframe$fishid)
dataframe$po2<-
  (dataframe$O2/100)*(dataframe$daily_pressure-26.70) 
# multiply dailypres and correct for H2O vapor pressure (26.70) at 27°C
dataframe$log_po2<-log10(dataframe$po2)
dataframe$log_Y <- log10((dataframe$sat / (100-dataframe$sat))) 




## Analyze each oxygen equilibrium curve

# Setup data frame for data storage
df.sum<-
  data.frame(
    species = 
      aggregate(
        species~curve,
        data = dataframe,
        FUN = unique)[2],
      
    fishid = 
      aggregate(
        fishid~curve,
        data = dataframe,
        FUN = unique)[2],
      
    pH = 
      aggregate(
        pH~curve,
        data = dataframe,
        FUN = unique)[2],
      
    P50 = rep(x = NA,max(dataframe$curve)),
    
    n50 = rep(x = NA,max(dataframe$curve)),
    
    r2adj = rep(x = NA,max(dataframe$curve))
  )

for (k in 1:length(unique(df.sum$fishid))){
  print(k)
  pH7 <- df.sum[df.sum$pH==7,]
  pH7.6 <- df.sum[df.sum$pH==7.6,]
  print(pH7)
}

for (k in 1:max(dataframe$curve)){
  # subset dataframe
  df<-
    dataframe[dataframe$curve==k,]
  
  i <- ggplot(df, aes(log_po2, log_Y))
  fit <- lm(df$log_Y ~ df$log_po2)
  P50 <- 10^((-fit$coef[1]) / (signif(fit$coef[2])))
  # calculates P50 from Hill Plot
  n50<-signif(fit$coef[2]) 
  # calculates n50 from Hill plot
  
  i + geom_point(colour="red")+
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
    theme_classic()+
    xlab("log(PO2 [mmHg])") +
    ylab("log(Y/1-Y)")+
    labs(title = paste(df$fishid[1], "pH",df$pH[1]),
         subtitle = 
           paste(
             "P50 =",signif(P50,5), 
             " n50 =",signif(n50,5), 
             " y =", signif(fit$coef[2], 3), "x", signif(fit$coef[1], 3),
             " Adj R2 =", signif(summary(fit)$adj.r.squared, 3)))
  
  ggsave(
    filename = paste("./Figures/",df.sum$fishid[k],"_pH",df.sum$pH[k],"_Hill.pdf",sep = ""),
    width = 4,
    height = 4
  )
  #plotting Hill plot with n50, P50, intercept,slope and adjR2 shown
  
  po2<-seq(1,100, by=1)
  # create PO2 from 1 to 100mmHg for hemoglobin binding curve
  new_sat<- 100*(po2^n50) / ((po2^n50)+(P50^n50)) 
  # calculated Y from the hill eq by using P50 and n50
  p <- ggplot(mapping = aes(po2,new_sat))
  p + geom_line(colour="black")+
    theme_classic()+
    geom_point(data = df, aes(po2,sat),colour="red",size=2)+
    geom_text(data = df, aes(po2,sat, label = paste(po2,",", sat)), hjust = -0.3, vjust = 0)+
    labs(title = paste(df$fishid[1], "pH",df$pH[1]),
         subtitle = paste("P50 = ", signif(P50,5),"     ","n50 = ", signif(n50,3)))+
    xlab("PO2[mmHg]")+
    ylab("Saturation (%)")
  # plot of Hemoglobin oxygen binding curve withindividual data, P50 and n50 are shown.
  
  ggsave(
    filename = paste("./Figures/",df.sum$fishid[k],"_pH",df.sum$pH[k],"_OEC.pdf",sep = ""),
    width = 4,
    height = 4
  )
  
  # export derived parameters to df.sum
  df.sum$r2adj[k]<-summary(fit)$adj.r.squared
  df.sum$P50[k]<-10^((-fit$coef[1]) / (fit$coef[2]))
  df.sum$n50[k]<-fit$coef[2]
}

df.sum

#importing P50 and n50 values to dataframe
#define airbreathing in df.sum
df.sum$airbreathing<-
  ifelse(df.sum$species=="Corydoras_paleatus", 
         0,
         ifelse(df.sum$species=="Auchenoglanis_occidentalis",
                0,
                ifelse(df.sum$species=="Ancistrus_sp.", 
                       1,
                       ifelse(df.sum$species=="Synodontis_petricula",
                              1,
                              1))))#else == Pangasius
#convert to factor
df.sum$airbreathing<-
  as.factor(df.sum$airbreathing)


# Import tree
NC_tree<-read.tree("/Users/au231308/Dropbox/Projects/Air-breathing review_full/mcc.nexus") # MCC tree from Rabosky, download from github page
df.sum$sp <- as.character(df.sum$species)
dataframe$sp <- as.character(dataframe$species)
df.sum$sp[which(df.sum$sp=="Ancistrus_sp.")]<-"Ancistrus_cirrhosus"
dataframe$sp[which(dataframe$sp=="Ancistrus_sp.")]<-"Ancistrus_cirrhosus"
NC_tree<-keep.tip(NC_tree,df.sum$sp)
NC_tree<-ladderize(NC_tree)
plot(NC_tree)

prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
#lav nodelabels
NC_tree$node.label<-seq(1,length(NC_tree$node.label),1)
#inserted tol=1
inv.phylo<-MCMCglmm::inverseA(NC_tree,nodes="TIPS",scale=TRUE, tol = 1, reduced = FALSE)

#creates bohr data set
bohr_data <- df.sum[df.sum$pH==7.6,]
bohr_data["bohreffect"] <- (log10(df.sum$P50[df.sum$pH==7.6])-log10(df.sum$P50[df.sum$pH==7.0]))/(7.6-7.0)



df.sum$airbreathing

fit.p50 <- 
  MCMCglmm::MCMCglmm(
    log10(P50)~airbreathing,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df.sum[which(df.sum$pH==7.6),],
    nitt=500000,
    burnin=1000,
    thin=500)
summary(fit.p50)


fit.n50 <- 
  MCMCglmm::MCMCglmm(
    log10(n50)~airbreathing,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df.sum[which(df.sum$pH==7.6),],
    nitt=500000,
    burnin=1000,
    thin=500)
summary(fit.n50)

bohr_data$bohreffect
fit.bohr <- 
  MCMCglmm::MCMCglmm(
    bohreffect~airbreathing,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=bohr_data,
    nitt=500000,
    burnin=1000,
    thin=500)
summary(fit.bohr)


# Phylogenetic signal
p50<-aggregate(P50~sp, df.sum[which(df.sum$pH==7.6),],mean)
n50<-aggregate(n50~sp, df.sum[which(df.sum$pH==7.6),],mean)
bohr<-aggregate(bohreffect~sp, bohr_data, mean)
range(bohr$bohreffect)
phylosig(tree = NC_tree,x = setNames(p50$P50,p50$sp),method = "lambda",test = T)
phylosig(tree = NC_tree,x = setNames(n50$n50,n50$sp),method = "lambda",test = T)
phylosig(tree = NC_tree,x = setNames(bohr$bohreffect,bohr$sp),method = "lambda",test = T)
ace <- signif(fastAnc(NC_tree, setNames(p50$P50,p50$sp)),3)


# Fig. 1A
NC_tree$root.edge <- 25
plot.phylo(NC_tree,root.edge = T,show.tip.label = F)
nodelabels(ace,frame = "none",adj = -0.1)
edgelabels(frame = "circle",edge = c(8,5,4), col = "#e41a1c", bg = "#e41a1c")


# Fig. 1BCD
#species order to fit phylogeny
species_order <- c("Auchenoglanis_occidentalis",
                   "Synodontis_petricola",
                   "Pangasianodon_hypophthalmus",
                   "Corydoras_paleatus", 
                   "Ancistrus_cirrhosus")


#P50 pointplot
ggplot(bohr_data, aes(y=factor(sp,level=species_order), x=P50, colour=airbreathing)) + 
  geom_point(shape = 1) +
  scale_x_continuous(limits = c(0,6),breaks = c(1,3,5))+
  scale_color_manual(breaks = c(0, 1),
                     values=c("#377eb8", "#e41a1c"),
                     labels=c("Water breathing", "Air breathing")) +
  xlab(label = expression("P"["50"]*" at pH 7.6"*" (mmHg)")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y =element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.line.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())->p1;p1

#n50 pointplot
ggplot(bohr_data, aes(y=factor(sp,level=species_order), x=n50, colour=airbreathing)) + 
  geom_point(shape = 1) +
  scale_x_continuous(limits = c(1,3),breaks = c(1,2,3))+
  scale_color_manual(breaks = c(0, 1),
                     values=c("#377eb8", "#e41a1c"),
                     labels=c("Water breathing", "Air breathing")) +
  xlab(label = expression("n"["50"]*" at pH 7.6")) +
  theme_classic()+
  theme(legend.position = "none",
        panel.grid.major.y = element_line(colour = "grey90"),
        axis.title.y =element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.line.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())->p2;p2

#bohreffekt pointplot
ggplot(bohr_data, aes(y=factor(sp,level=species_order), x=bohreffect, colour=airbreathing)) + 
  geom_point(shape=1) +
  scale_color_manual(breaks = c(0, 1),
                     values=c("#377eb8", "#e41a1c"),
                     labels=c("Water breathing", "Air breathing")) +
  scale_x_reverse(limits = c(0,-1),breaks = c(0,-0.5,-1)) +
  xlab(expression(phi*" ("*Delta*"log"["10"]*"(P"["50"]*") "*Delta*"pH"^"-1"*")")) +
  theme_classic()+
  theme(axis.title.y =element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())->p3;p3


pdf("./PanelB.pdf",width = 3.75,height = 2,useDingbats = F)
cowplot::plot_grid(p1,p2,p3,nrow = 1,ncol =3,align = "h",labels = c("B","C","D"),label_size = 12)
dev.off()

