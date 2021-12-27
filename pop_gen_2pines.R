
###################### MAP of sampled P. pungens and P. rigida populations ####################
install.packages("ggplot2")
library(ggplot2)
install.packages("sf")
library(sf)
theme_set(theme_bw())

install.packages("rnaturalearth")
library(rnaturalearth)
install.packages("rnaturalearthdata")
library(rnaturalearthdata)
install.packages("ggthemes")
library(ggthemes)
install.packages("ggspatial")
library(ggspatial)
install.packages("ggrepel")
library(ggrepel)
install.packages("ggsn")
library(ggsn)


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

setwd("~/path/to/dir")
sites <- read.csv("rigida_mean_coordinates_perPOP.csv")
sites_sf <- st_as_sf(sites, coords = c("mean_long", "mean_lat"), 
                     crs = 4326, agr = "constant")
View(sites)
poly <- st_read("~/path/to/pinurigi.shp", crs=4326)

rigida_map <- ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  geom_sf(data=poly, size=0.5, color="darkorange", fill= "orange") +
  geom_sf(data = sites_sf, size = 6,shape = 21, color="white", fill = "black") +
  coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") +
  coord_sf(xlim = c(-90, -65), ylim = c(25, 48), expand = FALSE)+
  ggsn::scalebar(x.min = -73, x.max = -68, y.min= 27, y.max=27.5, dist = 200, transform=TRUE,model='WGS84', height=0.5,st.dist = 2, dist_unit ="km")

north2(rigida_map, 0.8, 0.25, symbol = 3)

ggsave("~/path/to/Rigida_pops_sampled_map.pdf", width = 6, height = 7)



sites_pung <- read.csv("pungens_mean_coordinates_perPOP.csv")
sites_sf_pung <- st_as_sf(sites_pung, coords = c("pung_mean_long", "pung_mean_lat"), 
                          crs = 4326, agr = "constant")
View(sites_pung)
poly_pung <- st_read("~/path/to/pinupung.shp", crs=4326)

pungens_map <- ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  geom_sf(data=poly_pung, size=0.5, color="blue", fill= "royalblue") +
  geom_sf(data = sites_sf_pung, size = 6,shape = 21, color="white", fill = "black") +
  coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") +
  coord_sf(xlim = c(-90, -65), ylim = c(25, 48), expand = FALSE)+
  ggsn::scalebar(x.min = -73, x.max = -68, y.min= 27, y.max=27.5, dist = 200, transform=TRUE,model='WGS84', height=0.5,st.dist = 2, dist_unit ="km")

north2(pungens_map, 0.8, 0.25, symbol = 3)

ggsave("~/path/to/Pungens_pops_sampled_map.pdf", width = 6, height = 7)

################################## STRUCTURE ##########################################
##(1) fastSTRUCTURE 
## assessment for selecting the best K from log likelihood results files across 10 runs. 

####You may want to also run STRUCTURE separately for each species ####
# This is optional, but could be useful to describe genetic structure within species. 


remotes::install_github('royfrancis/pophelper')
library(pophelper)
install.packages("corrplot")
library(corrplot)

setwd("~/Documents/Ch1_files/SNP_data/K2_runs")
Run1 <- readQ("Run1_fastSTRUCTK2.2.meanQ")
Run2 <- readQ("Run2_fastSTRUCTK2.2.meanQ")
Run3 <- readQ("Run3_fastSTRUCTK2.2.meanQ")
Run4 <- readQ("Run4_fastSTRUCTK2.2.meanQ")
Run5 <- readQ("Run5_fastSTRUCTK2.2.meanQ")
Run6 <- readQ("Run6_fastSTRUCTK2.2.meanQ")
Run7 <- readQ("Run7_fastSTRUCTK2.2.meanQ")
Run8 <- readQ("Run8_fastSTRUCTK2.2.meanQ")
Run9 <- readQ("Run9_fastSTRUCTK2.2.meanQ")
Run10 <- readQ("Run10_fastSTRUCTK2.2.meanQ")

sfiles <- list.files(path = "~/path/to/K2_runs", full.names = T)
slist <- readQ(files=sfiles, indlabfromfile = T)
head(slist)
attributes(slist)
align <- alignK(slist)
write.csv(align, "aligned_K2_pophelper.csv")
clusters <- read.csv("aligned_K2_pophelper.csv")
mislabeled <- clusters[104,] # sample 104 was incorrectly identified in the field
reassign_mislabel <- rbind(clusters[1:103,], clusters[105:301,])
cluster1 <- reassign_mislabel[,c(1,2,4,6,8,10,12,14,16,18,20)]
cluster2 <- reassign_mislabel[, c(1,3,5,7,9,11,13,15,17,19,21)]
avg_cluster1 <- rowMeans(cluster1[,-1])
View(avg_cluster1)
avg_cluster2 <- rowMeans(cluster2[,-1])
mean_clusters <- cbind(reassign_mislabel[,1], avg_cluster1, avg_cluster2)
View(mean_clusters)
write.csv(mean_clusters, "averaged_K2_assignments.csv")

p1 <- plotQ(align, imgoutput = "join", returnplot = T, exportplot = F, basesize = 11)
p1

K_results <- read.csv("~/Documents/SNPdata_Ch1/SNP_data/Pop_ID_with_location_data_AND_K_assignments_minus104.csv")
head(K_results)
barplot(t(as.matrix(K_results[,8:9])), col= c("blue", "orange"), border=NA, space=FALSE, ylab="Ancestry (Q)", ylim=c(0.0,1.0), legend.text=c(expression(italic("P. pungens")), expression(italic("P. rigida"))))
dev.copy2pdf(file="~/path/to/K2_structurePLOT_minus104.pdf", useDingbats=FALSE, family="sans")


## K plot by species by latitude

K_byLat <- read.csv("~/path/to/Pop_ID_K_by_latitude.csv")
head(K_byLat)
barplot(t(as.matrix(K_byLat[,8:9])), col= c("blue","orange"), border=NA, space=FALSE, ylab="Ancestry (Q)", ylim=c(0.0,1.0), legend.text=c(expression(italic("P. pungens")), expression(italic("P. rigida"))))
dev.copy2pdf(file="~/path/to/K2_structurePLOT_minus104_BYlatitude.pdf", useDingbats=FALSE, family="sans")

avg_assignments <- read.csv("averaged_K2_assignments_minus104.csv")
head(avg_assignments)
pung_avg_ass <- avg_assignments[1:125,]
dim(pung_avg_ass)
head
ADMIX_2 <- which(pung_avg_ass$avg_cluster1 <(0.98))
ADMIX_2 #28
pung_avg_ass[pung_avg_ass$avg_cluster1 <= 0.98 & pung_avg_ass$avg_cluster1 >= 0.90,]
#16
rig_avg_ass <- avg_assignments[126:300,]
ADMIX_2r <- which(rig_avg_ass$avg_cluster2 <(0.98))
ADMIX_2r #29
rig_avg_ass[rig_avg_ass$avg_cluster2 <= 0.98 & rig_avg_ass$avg_cluster2 >= 0.90,]
#25

## correlation with lat, long, elev
Pop_ID_Sum <- read.csv("~/path/to/Pop_ID_Sum.csv")
env_data <- read.csv("~/path/to/env_final_5var.csv")
head(env_data)
elev <- env_data[,6]
Pop_coord_k2 <- Pop_ID_Sum[,6:8]
cor(Pop_coord_k2$Latitude, Pop_coord_k2$anc_k2 ) #-0.4140138 (CI= -0.50364 - -0.3155679)
cor(Pop_coord_k2$Longitude, Pop_coord_k2$anc_k2) #-0.2910068 (CI= -0.3913511 - -0.1838236 )
cor(elev, Pop_coord_k2$anc_k2) #0.5091966 (CI= 0.42018 - 0.58850)

data_cor <- cbind(Pop_coord_k2, elev)
head(data_cor)

signif_cor <- cor.mtest(data_cor, conf.level = .95)
signif_cor # all significant
write.csv(signif_cor, "~/path/to/Admixture_correlation_to_ENV_significance_and_CI.csv")
#$p
#Longitude     Latitude       anc_k2         elev
#Longitude 0.000000e+00 4.532503e-94 2.884947e-07 3.074140e-34
#Latitude  4.532503e-94 0.000000e+00 7.473537e-14 9.981061e-44
#anc_k2    2.884947e-07 7.473537e-14 0.000000e+00 3.403278e-21
#elev      3.074140e-34 9.981061e-44 3.403278e-21 0.000000e+00




##################################### PCA ######################################## 


###PCA for 012 coded  vcf files
##Following method in Patterson et al 2006

library(data.table)
library(tidyverse)
library(ggsci)
library(ggforce)
library(devtools)
library(viridis)
install.packages("LEA")
library(LEA)
library(lme4)
library(lmerTest)
library(readr)
library(plyr)
library(readxl)
library("wesanderson")

source('~/path/to/PCA/Imports.R')

########### PCA_gen ##############

##### input files #####

### df_gen: genotypic data with individuals as rows and snps as coloumns.
###         Can include missing data. Either genotype probabilities or 012 format

### indv:  dataframe with rows corresponding to individuals in df_gen file 
###        Must have Pop and ID coloumn 

##### output files ####

### pca_out: 

### $ 'pca_df': dataframe with rows as individuals and coloumns as PC1-30, Pop, ID
### $ 'pve': list of proportion of variance explained for each PC

####### function ###########

PCA_gen <- function(df_gen,indv,num=5){ #add ggplot, add tw, add # of output
  pkgTest("ggplot2")
  
  df_gen <- apply(df_gen, 2, function(df) gsub(-1,NA,df,fixed=TRUE))
  df_gen <- apply(df_gen, 2, function(df) as.numeric(df))
  
  colmean<-apply(df_gen,2,mean,na.rm=TRUE)
  
  normalize<-matrix(nrow = nrow(df_gen),ncol=ncol(df_gen))
  af<-colmean/2
  
  for (m in 1:length(af)){
    nr<-df_gen[ ,m]-colmean[m]
    dn<-sqrt(af[m]*(1-af[m]))
    normalize[ ,m]<-nr/dn
  }
  
  normalize[is.na(normalize)]<-0
  
  method1<-prcomp(normalize, scale. = FALSE,center = FALSE)
  pve <- summary(method1)$importance[2,]
  print(pve[1:300])
  
  if(nrow(df_gen) < num){
    num <- nrow(df_gen)
  }
  
  pca_X<- method1$x[,1:num]
  
  pca_X <- as.data.frame(pca_X)
  pca_X$Pop <- indv$Sp_Pop
  pca_X$ID <- indv$ID
  
  pca_out <- list("pca_df"=pca_X,"pve"=pve)
  
  print(PCA_fig(pca_out))
  
  return(pca_out)
  return(method1)
}


############ Read in files and make PCA #########
df012<-fread("~/path/to/single_dadi_minus104.recode.vcf.recode.vcf.gz.012",sep="\t", data.table=F) 
df012 <- df012[,-1]
dim(df012)
#[1]  300 2168
df012_na <- apply(df012, 2, function(df) gsub('-1','NA',df,fixed=TRUE))
View(df012_na[1:20,1:20])
formated<-apply(df012_na,2,function(df) as.numeric(df))
View(formated[1:20, 1:20])
Pop_ID_Sum <- read.csv("~/path/to/Pop_ID_Sum.csv")


############# PCA results ##############
pca_out <- PCA_gen(formated,Pop_ID_Sum)
head(pca_out)
pve <- pca_out$pve[1:5]
pve
#    PC1     PC2     PC3     PC4     PC5 
#0.04232 0.00849 0.00707 0.00658 0.00623

pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(Pop_ID_Sum,pca_df)

write.csv(pca_df,'pca_df.csv',row.names = FALSE)




########## PCA PLOT ###########

pca_df <- read.csv('pca_df.csv')
pve <- c(0.04232, 0.00849, 0.00707, 0.00658, 0.00623)

#Sp

colSp <- c('#376BD1','#FAA722')
pca_sp <- ggplot(data = pca_df, aes(x=-1*PC1,y=PC2,fill=as.character(Sp))) + 
  geom_point(colour='black',size = 5,pch=21) + 
  xlim(-30,30) +
  ylim(-20,20) +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
  scale_fill_manual(values=colSp) +
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pca_sp

ggsave('PCA_Sp.pdf',pca_sp,height=5,width=6,units='in')




#### PCA color coded by population assignment

#colors chosen from iwanthue.com hard(force vector)
col33 <- c("#b4a4ff",
           "#237f00",
           "#a843c9",
           "#00d78c",
           "#930994",
           "#bad154",
           "#6b4bd1",
           "#7b7d00",
           "#b57eff",
           "#01702c",
           "#be0081",
           "#02d3c0",
           "#d00728",
           "#0194c1",
           "#ff762e",
           "#0058b5",
           "#c2cd6c",
           "#ff97e8",
           "#425b15",
           "#e7acff",
           "#8b6900",
           "#90b9ff",
           "#9b3f00",
           "#436c42",
           "#ff5a55",
           "#d7c58c",
           "#8c315c",
           "#ffb570",
           "#ff77a3",
           "#844600",
           "#df94a2",
           "#823f39",
           "#ff8894")

pca_sp_pop <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Sp_Pop),shape=as.character(Sp))) + 
  geom_point(colour='black',size = 4) + 
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
  scale_fill_manual(values=col33) + 
  scale_shape_manual(values=c(21,24)) + 
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pca_sp_pop

ggsave('~/path/to/PCA_Sp_Pop.pdf',pca_sp_pop,height=5,width=6,units='in')


sp_pop_df <- unique(pca_df[,c(1,2,5)])
shape_list <- rep(NA,times=nrow(sp_pop_df))
shape_list[which(sp_pop_df$Sp == 'PU')] <- 21
shape_list[which(sp_pop_df$Sp == 'RI')] <- 24

#sp_pop_ with legend
pca_sp_pop_leg <- ggplot(data = pca_df, aes(x=PC1,y=PC2,fill=as.character(Sp_Pop),shape=as.character(Sp))) + 
  geom_point(colour='black',size = 4) + 
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
  scale_fill_manual(name="Species_Pop:",values=col33,
                    guide = guide_legend(override.aes = list(shape = shape_list))) + 
  scale_shape_manual(name="Species_Pop:",values=c(21,24)) + 
  #guide = guide_legend(override.aes = list(fill=col33))) + 
  theme_bw() + guides(shape = 'none')  +
  theme(#legend.position = 'none',
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    legend.text = element_text(size=13), 
    legend.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
pca_sp_pop_leg

ggsave('Figures/PCA_Sp_Pop_leg.pdf',pca_sp_pop_leg,height=5,width=7,units='in')






############################# Hierarchical F-statistics ###########################
##(3) Hierarchical F-statistics as estimated using hierfstat, 
#where grouping levels are species, populations nested in species, and individuals 
#nested into populations.
install.packages('hierfstat')
library(hierfstat)


setwd("~/path/to/SNP_data")
data_vcf <- read.table("single_dadi.recode.vcf", header=F)
View(data_vcf[299:350,1:20])
SNP_info <- data_vcf[,1:8]
head(SNP_info)
write.table(SNP_info, "SNP_info_from_vcf_file.txt", sep="\t")
data_012 <- read.table("single_dadi.recode.vcf.012")
data_012 <- data_012[,-1]
data_012[1:20,1:20]
dim(data_012)
View(data_012[100:142,1:20])


MinorToHierf<-function(df,n){
  
  #this script will use a dataframe with snp and ind and pop information
  # the snps need to be coded as counts of the minor allele
  #n indicates the column at which the loci information begins
  
  df<-df[ ,n:ncol(df)]
  
  formated <- apply(df, 2, function(df) gsub('-1','NA',df,fixed=TRUE))
  formated<- apply(formated, 2, function(df) gsub('1','het',df))
  formated <- apply(formated, 2, function(df) gsub('2','min',df))
  formated <- apply(formated, 2, function(df) gsub('0','maj',df))
  
  formated <- apply(formated, 2, function(df) gsub('NA', 'NA', df))
  formated<- apply(formated, 2, function(df) gsub('het','12',df))
  formated <- apply(formated, 2, function(df) gsub('maj','22',df))
  formated <- apply(formated, 2, function(df) gsub('min','11',df))
  
  formated<-apply(formated,2,function(df) as.numeric(df))
  return(data.frame(formated))
}

formatedGdf<-MinorToHierf(data_012,n = 1)


## format "loci" to have sample IDs as column names

data_part <- data_vcf[,10:310]
indv <- read.table("single_dadi.recode.vcf.012.indv")
View(indv)
indv_names <- as.character(indv$V1)
colnames(data_part) <- indv_names
View(data_part[299:350,100:120])



setwd("~/Documents/SNPdata_Ch1")
contigs <- read.table("single_dadi.recode.vcf.012.pos")
View(contigs)
write.table(contigs, "SNP_ID_and_position.txt", sep="\t")
contig_names <- as.character(contigs$V1)
colnames(formatedGdf) <- contig_names
View(formatedGdf[280:302,1:20])
dim(formatedGdf)


loci <- formatedGdf
dim(loci)
View(loci[100:140,1:20])
df_minus104 <- loci[-104,] #misidentified sample
dim(df_minus104)

write.table(df_minus104, "loci_for_hierfstats_minus104.txt", sep="\t")
View(df_minus104[102:120,1:20])


#factors
factors <- read.table("~/path/to/levels_for_hierFstats.txt", sep="\t")
dim(factors)
factors_edit <- factors[-104,] ##  deleting PU_SH_7 because it was misidentified. It was a pitch pine, but no other pitch pines
## were sampled from this location.
View(factors_edit)

colnames(factors_edit) <- c("Sp", "Pop", "Indv")
factors_edit1 <- factors_edit[,1:2]
write.table(factors_edit1, "levels_for_hierFstats_minus104.txt", sep="\t")

factors_edit$Sample_ID <- paste(factors_edit$Pop,factors_edit$Indv, sep="_")
write.csv(factors_edit, "Sample_IDs_Sp_pop_indv_minus104.csv")



Fstats <- varcomp.glob(levels = factors_edit1, loci= df_minus104, diploid = TRUE)
Fstats$F
#$F
#           Sp              Pop              Ind
#Total      0.1171143   0.123325587   0.16079417
#Sp         0.0000000   0.007035155   0.04947392
#Pop        0.0000000   0.000000000   0.04273945(Fis)
### above values not exactly like AMOVA because it was done in a more linear mixed model way (Yang, 1998)
# see Arlequin program for nominclature (Fct, etc)


Fstats$overall
write.csv(Fstats$overall, "Fstats_overall_minus104.csv")
#  Sp           Pop        Ind       Error 
#35.526527   1.884173  11.366058 254.572294 


write.csv(Fstats$F, "Fstats_minus104.csv")

write.csv(Fstats$loc, "Fstats_byLOCUS_minus104.csv")


#### calculate using the variance components
### Fstats_by_Locus.csv was modified within Excel to calculate stats per locus
setwd("~/path/to/dir")
F_loci <- read.csv("Fstats_byLOCUS_minus104_v2.csv")

F_loci$Fraction.of.Structure.due.to.species[which(F_loci$Fraction.of.Structure.due.to.species < 0)] <- 0

hist(F_loci$Fraction.of.Structure.due.to.species)

F_loci$F.Sp.total.[which(F_loci$F.Sp.total.<0)] <- 0
hist(F_loci$F.Sp.total.)
max(F_loci$F.Sp.total.)

F_loci$F.pop.Sp.[which(F_loci$F.pop.Sp.<0)] <- 0 
hist(F_loci$F.pop.Sp.)

F_loci$F.pop.total.[which(F_loci$F.pop.total.<0)] <- 0
hist(F_loci$F.pop.total.)

write.csv(F_loci, "Fstats_byLOCUS_minus104_noNegatives.csv")




## estimate confidence intervals for multilocus F-statistics
boot <- boot.vc(levels=factors_edit1,loci=df_minus104,diploid=TRUE,nboot=1000,quant=c(0.025,0.5,0.975))
df_boot <- boot$boot
write.csv(df_boot, "df_boot_Fstats_minus104.csv")
bott_stats <- boot$res
write.csv(bott_stats, "Fstats_bootstrap_minus104.csv")
boot_ci <- boot$ci
boot$ci
#                    Fct          Fst                          Fsc
#       H-Total F-Sp/Total F-Pop/Total F-Ind/Total   H-Sp F-Pop/Sp F-Ind/Sp  H-Pop F-Ind/Pop
# 2.5%   0.1345     0.0993      0.1057      0.1382 0.1188   0.0055   0.0358 0.1180    0.0285
# 50%    0.1399     0.1167      0.1228      0.1602 0.1235   0.0070   0.0492 0.1226    0.0424
# 97.5%  0.1459     0.1364      0.1431      0.1851 0.1282   0.0088   0.0637 0.1273    0.0570
#Hobs
# 2.5%  0.1122
# 50%   0.1174
# 97.5% 0.1230
write.csv(boot_ci, "Fstats_confidence_intervals_minus104.csv")



View(factors_edit1)

############## Summaries of genetic diversity ###############
#(4) Summaries of genetic diversity by population for each species, including at least observed 
#and expected heterozygosities for each population. 

install.packages('hierfstat')
library(hierfstat)

# all the files in the hierfstat_outputs_ALLpopsTWOspecies folder are likely.. 
#irrelevant/not useful
 #setwd("~/Documents/SNPdata_Ch1")
#View(factors_edit)
data_pop <- cbind(factors_edit1$Pop, df_minus104)
colnames(data_pop)[which(names(data_pop)== "factors_edit1$Pop")] <- "Pop_ID"
View(data_pop[100:130,1:20])




#### genetic diversity for P pungens populations ####
data_pung <- data_pop[1:125,]
View(data_pung[103:125,1:20])
View(df_minus104[1:20,1:20])


data_pung[] <- lapply(data_pung, type.convert)
Pung_diversity <- basic.stats(data_pung, diploid = TRUE, digits = 4)
Pung_diversity$overall
# Ho     Hs      Ht     Dst    Htp    Dstp    Fst   Fstp    Fis   Dest 
# 0.1267 0.1306 0.1313 0.0007 0.1314 0.0008 0.0055 0.0059 0.0296 0.0009 

write.csv(Pung_diversity$overall, "pungens_genetic_diveristy_estimates.csv")
write.csv(Pung_diversity$Ho, "pungens_obs_heterozygosity_perPOP_perLOCUS.csv")
write.csv(Pung_diversity$Hs, "pungens_obs_gene_diversity_perPOP_perLOCUS.csv")
write.csv(Pung_diversity$Fis, "pungens_obs_Fis_perPOP_perLOCUS.csv")
write.csv(Pung_diversity$perloc, "pungens_basic_stats_perLOCUS.csv")

#### Find CI for species Fst ####
## estimate confidence intervals for multilocus F-statistics
data_P <- data_pung[,-1]
View(data_P[1:20,1:20])
factors_edit_pung <- factors_edit1[1:125,2]
View(factors_edit_pung)

Fstats_pung <- varcomp.glob(levels = factors_edit_pung, loci= data_P, diploid = TRUE)
Fstats_pung$F
#         Pop        Ind
#Total 0.005750188 0.03435611
#Pop   0.000000000 0.02877137
Fstats_pung$overall
#Pop        Ind      Error 
#1.635907   8.138278 274.722074 
Fstats_pung$loc

boot <- boot.vc(levels=factors_edit_pung,loci=data_P,diploid=TRUE,nboot=1000,quant=c(0.025,0.5,0.975))
df_boot <- boot$boot
write.csv(df_boot, "df_boot_Pungens_Fstats_minus104.csv")
bott_stats <- boot$res
write.csv(bott_stats, "Pungens_Fstats_bootstrap_minus104.csv")
boot_ci <- boot$ci
boot$ci

#       H-Total F-Pop/Total F-Ind/Total  H-Pop  F-Ind/Pop   Hobs
#2.5%   0.1525      0.0032      0.0110   0.1517    0.0054  0.1461
#50%    0.1591      0.0058      0.0347   0.1581    0.0291  0.1537
#97.5%  0.1657      0.0084      0.0586   0.1650    0.0527  0.1610



snpval <- numeric(2168)

for(i in 1:2168) {
  snpval[i] <- length(table(data_P[,i]))
}
b <- which(snpval <2)

data_P_onlySNPs <- data_P[,-b]
View(data_P_onlySNPs[1:20,1:20])
write.table(data_P_onlySNPs, "Pungens_SNPs_with_variation.txt", sep="\t")

Fstats_pung_onlySNPs <- varcomp.glob(levels = factors_edit_pung, loci= data_P_onlySNPs, diploid = TRUE)
Fstats_pung_onlySNPs$F


#### genetic diversity for P rigida populations ####
data_rig <- data_pop[126:300,]
View(data_rig[103:125,1:20])
View(df_minus104[1:20,1:20])
data_R <- data_rig[,-1]
View(data_R[1:20,1:20])
factors_edit_rig <- factors_edit1[126:300,2]
View(factors_edit_rig)

Fstats_rig <- varcomp.glob(levels = factors_edit_rig, loci= data_R, diploid = TRUE)
Fstats_rig$F
#         Pop        Ind
#Total 0.005658343 0.1116902
#Pop   0.000000000 0.1066352
Fstats_rig$overall
#Pop        Ind      Error 
#1.414395  26.504383 222.047502 
Fstats_rig$loc

boot <- boot.vc(levels=factors_edit_rig,loci=data_R,diploid=TRUE,nboot=1000,quant=c(0.025,0.5,0.975))
df_boot_rig <- boot$boot
write.csv(df_boot_rig, "df_boot_Rigida_Fstats_minus104.csv")
bott_stats_rig <- boot$res
write.csv(bott_stats_rig, "Rigida_Fstats_bootstrap_minus104.csv")
boot_ci_rig <- boot$ci
boot$ci
#       H-Total F-Pop/Total  F-Ind/Total  H-Pop   F-Ind/Pop   Hobs
#2.5%   0.1215      0.0032      0.0930    0.1208    0.0874   0.1069
#50%    0.1267      0.0056      0.1121    0.1260    0.1070   0.1125
#97.5%  0.1323      0.0082      0.1302    0.1316    0.1247   0.1181



data_rig[] <- lapply(data_rig, type.convert)
Rig_diversity <- basic.stats(data_rig, diploid = TRUE, digits = 4)
Rig_diversity$overall
#  Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.1018 0.1144 0.1150 0.0006 0.1150 0.0007 0.0054 0.0057 0.1102 0.0007 

write.csv(Rig_diversity$overall, "rigida_genetic_diveristy_estimates_final.csv")
write.csv(Rig_diversity$Ho, "rigida_obs_heterozygosity_perPOP_perLOCUS_final.csv")
write.csv(Rig_diversity$Hs, "rigida_obs_gene_diversity_perPOP_perLOCUS_final.csv")
write.csv(Rig_diversity$Fis, "rigida_obs_Fis_perPOP_perLOCUS_final.csv")
write.csv(Rig_diversity$perloc, "rigida_basic_stats_perLOCUS_final.csv")



snpval <- numeric(2168)

for(i in 1:2168) {
  snpval[i] <- length(table(data_R[,i]))
}
a <- which(snpval <2)

data_R_onlySNPs <- data_R[,-a]
View(data_R_onlySNPs[1:20,1:20])
write.table(data_R_onlySNPs, "Rigida_SNPs_with_variation.txt", sep="\t")

Fstats_rig_onlySNPs <- varcomp.glob(levels = factors_edit_rig, loci= data_R_onlySNPs, diploid = TRUE)
Fstats_rig_onlySNPs$F # output same as above with data_R used for loci





##################### Calculate exp and obs heterozygosity per SNP per population #####################


View(df012_pop[1:30,1:30])


# replace -1 with NA in SNP data
gdat <- df012_pop
df_gen <- apply(gdat, 2, function(df) gsub(-1,NA,df,fixed=TRUE))
dim(df_gen)
df_gen <- apply(df_gen[1:300,2:2169], 2, function(df) as.numeric(df))
View(df_gen[1:20,1:20])
df_gen2 <- cbind(Pop_ID$Sp_Pop, df_gen)
colnames(df_gen2) <- SNP_names
View(df_gen2[1:20,1:20])
dim(df_gen2)
View(df_gen2[,1])
write.table(df_gen2, "Two_SP_data_with_NAs_minus104.txt", sep="\t")

df_gen3 <- read.table("~/path/to/Two_SP_data_with_NAs_minus104.txt", sep="\t")
View(df_gen3[1:20,1:20])



##### Function to get Ho for a vector of 0,1,2 genotypes #####

obsH <- function(dat) {
  
  if(length(which(is.na(dat) == "TRUE")) == 0) {nsamps <- length(dat)} else {nsamps <- length(dat[-which(is.na(dat) == "TRUE")])}
  
  nhets <- length(dat[which(dat == 1)])
  
  if(nsamps > 0) {return(nhets/nsamps)} else {return(NA)} ## This should account for cases where all samples are missing for a given population.
  
}

expH <- function(dat) {
  
  if(length(which(is.na(dat) == "TRUE")) == 0) {nsamps <- length(dat)} else {nsamps <- length(dat[-which(is.na(dat) == "TRUE")])}
  
  acounts <- sum(dat, na.rm = T)
  
  minaf <- acounts/(2*nsamps)
  
  majaf <- 1 - minaf
  
  if(nsamps > 0) {return(2*minaf*majaf)} else{return(NA)} ## This should account for cases where all samples are missing for a given population.
  
}

### Function to take input file and parse into populations using a list, where each element of the list is a population.

popParse <- function(dat) {
  
  popsplits <- vector("list", length(unique(dat[,1])))
  
  labs <- unique(dat[,1])
  
  for(i in 1:length(labs)) {
    
    popsplits[[i]] <- dat[which(dat[,1] == labs[i]),2:ncol(dat)]
    
  }
  
  setNames(popsplits, labs)
  
  return(popsplits)
  
}

### Now, use the functions above to cycle through the full dataset to get Ho and He by pop for each SNP.

### create empty output objects to fill... one for Ho and one for He
### The dims of each object are: nrows = no. pops., ncols = no. snps.

labs <- unique(df_gen3[,1])

ho_out <- matrix(0, nrow = length(labs), ncol = ncol(df_gen3)-1)

row.names(ho_out) <- labs

colnames(ho_out) <- colnames(df_gen3[2:ncol(df_gen3)])

he_out <- matrix(0, nrow = length(labs), ncol = ncol(df_gen3)-1)

row.names(he_out) <- labs

colnames(he_out) <- colnames(df_gen3[2:ncol(df_gen3)])

### Next, parse full dataset.

pops_parse <- popParse(dat = df_gen3)

### Next to last (or penultimately if you are exporting the resulting objects), use simple for() loop to cycle through pops to get each statistic for each SNP.

for(i in 1:length(labs)) {
  
  ho_out[i,] <- apply(pops_parse[[i]], 2, obsH) ### Applies obsH across all SNPs at once.
  
  he_out[i,] <- apply(pops_parse[[i]], 2, expH) ### Applies expH across all SNPs at once.
  
}


### Last, use write.table() to output the results as needed. 
write.table(he_out, "exp_heterozygosity_perSNPperPOP_df_gen3_minus104.txt", sep="\t")
write.table(ho_out, "obs_heterozygosity_perSNPperPOP_df_gen3_minus104.txt", sep="\t")




##### Ho mean and SD #####
H_means <- read.csv("~/path/to/mean_perPOP_heterozygosity.csv")
H_PU <- Ho_means[1:14,]
H_RI <- Ho_means[15:33,]

RI_mean <- mean(H_RI$mean_Ho) #0.101796
RI_sd <- sd(H_RI$mean_Ho) #0.008929
PU_mean <- mean(H_PU$mean_Ho) #0.1268225
PU_sd <- sd(H_PU$mean_Ho) #0.0150094

##### He mean and SD #####
RI_mean2 <- mean(H_RI$mean_He) #0.104102
RI_sd2 <- sd(H_RI$mean_He) #0.0047928
PU_mean2 <- mean(H_PU$mean_He) #0.11821
PU_sd2 <- sd(H_PU$mean_He) #0.0080857

#### correlations of Ho with lat, long, elevation ####
Ho <- read.csv("~/path/to/obs_heterozygosity_perSNPperPOP_df_gen3_minus104.csv")
Ho_pung <- Ho[1:14,]
Ho_rig <- Ho[15:33,]
cor_Ho_pung_lat <- cor(Ho_pung$Latitude, Ho_pung$mean_Ho) #0.08009
cor_Ho_pung_long <- cor(Ho_pung$Longitude, Ho_pung$mean_Ho) #0.1750
cor_Ho_rig_lat <- cor(Ho_rig$Latitude, Ho_rig$mean_Ho) #-0.008109
cor_Ho_rig_long <- cor(Ho_rig$Longitude, Ho_rig$mean_Ho) #0.1126
elev <- read.csv("~/path/to/Pungens_WorldClim_values_perPOP_minus104.csv")
cor_Ho_pung_elev <- cor(elev$elevation.asc, Ho_pung$mean_Ho) #-0.16846
rig_elev <- read.csv("~/path/to/Rigida_WorldClim_values_perPOP_minus104.csv")
cor_Ho_rig_elev <- cor(rig_elev$elevation.asc, Ho_rig$mean_Ho) #0.25304



################################ Missing data ##############################3#
## missing data per SNP
total_NA <- colSums(is.na(data_pop))
write.csv(total_NA, "missing_data_per_SNP_minus104.csv")

## missing data per individual calc
ind_NA <- rowSums(is.na(data_pop))
write.csv(ind_NA, "missing_data_per_INDV_minus104.csv")

## not missing
total_rep <- colSums(!is.na(data_pop))
write.csv(total_rep, "Has_data_per_SNP_minus104.csv")

## missing data for pungens 
pung_NA <- colSums(is.na(data_pung))
write.csv(pung_NA, "missing_SNPdata_for_pungens_minus104.csv")
## not missing
pung_rep <- colSums(!is.na(data_pung))
write.csv(pung_rep, "Has_SNPdata_for_pungens_minus104.csv")

## missing data for rigida
rig_NA <- colSums(is.na(data_rig))
write.csv(rig_NA, "missing_SNPdata_for_rigida_minus104.csv")
## not missing
rig_rep <- colSums(!is.na(data_rig))
write.csv(rig_rep, "Has_SNPdata_for_rigida_minus104.csv")



missingP_anc <- read.csv("~/path/to/pop_N_pungens.csv")
head(missingP_anc)
missing_data <- missingP_anc[,2]
anc <- missingP_anc[,4]
miss_anc <- cbind(missing_data,anc)
miss_anc <- as.data.frame(miss_anc)


cor_missingP <- cor(miss_anc)
cor_missingP
signif_missing_ancP <- cor.mtest(miss_anc, conf.level = .95)
signif_missing_ancP


missingR_anc <- read.csv("~/path/to/pop_N_rigida.csv")
head(missingR_anc)
missing_dataR <- missingR_anc[,2]
anc <- missingR_anc[,4]
miss_ancR <- cbind(missing_dataR,anc)
miss_ancR <- as.data.frame(miss_ancR)


subset90_p <- subset(miss_anc, anc < 0.90)
subset90plus_P <- subset(miss_anc, anc > 0.90)
subset95_p <- subset(miss_anc, anc < 0.95)
subset95plus_P <- subset(miss_anc, anc > 0.95)

boxplot(subset90_p$missing_data)
boxplot(miss_anc$missing_data)

subset10_r <- subset(miss_ancR, anc > 0.10)
subset5_r <- subset(miss_ancR, anc > 0.05)
subset10plus_R <- subset(miss_ancR, anc < 0.10)
subset5plus_R <- subset(miss_ancR, anc < 0.05)

boxplot(subset10_r$missing_data)
boxplot(miss_ancR$missing_data,subset10_r$missing_data, miss_anc$missing_data, subset90_p$missing_data, ylab= "Missing data counts", ylim=c(0,2000))

boxplot(subset5plus_R$missing_data,subset5_r$missing_data, subset95plus_P$missing_data, subset95_p$missing_data, ylab= "Missing data counts", ylim=c(0,2000))

together <- rbind(missingP_anc, missingR_anc)

Pop_ID_Sum <- read.csv("~/path/to/Pop_ID_Sum.csv")
SP <- Pop_ID_Sum$Sp
together_sp <- cbind(SP, together)
plot(together_sp$ancestry,together_sp$missing_data, colour=together_sp$SP, data= together_sp)
library(ggplot2)

ggplot(data = together_sp, aes(ancestry, missing_data, color = SP)) +
  geom_point() +
  xlim(0,1) +
  ylim(0,1300) +
  scale_color_manual(values = c("PU" = "blue", "RI" = "orange"))


cor_missingR <- cor(miss_ancR)
cor_missingR
signif_missing_ancR <- cor.mtest(miss_ancR, conf.level = .95)
signif_missing_ancR

total_data <- rbind(miss_anc, miss_ancR)
cor_total <- cor(total_data)
cor_total
sign_total <- cor.mtest(total_data, conf.level = .95)
sign_total




#### minor allele counts per pop ####

## pungens pops
pop1 <- subset(data_pop, Pop_ID == "PU_BB")
View(data_pop[1:10,1:10])
View(pop1[1:5,1:5])
pop2 <- subset(data_pop, Pop_ID == "PU_BN")
pop3 <- subset(data_pop, Pop_ID == "PU_BV")
pop4 <- subset(data_pop, Pop_ID == "PU_DT")
pop5 <- subset(data_pop, Pop_ID == "PU_EG")
pop6 <- subset(data_pop, Pop_ID == "PU_EK")
pop7 <- subset(data_pop, Pop_ID == "PU_GA")
pop8 <- subset(data_pop, Pop_ID == "PU_LG")
pop9 <- subset(data_pop, Pop_ID == "PU_NM")
pop10 <- subset(data_pop, Pop_ID == "PU_PM")
pop11 <- subset(data_pop, Pop_ID == "PU_SC")
pop12 <- subset(data_pop, Pop_ID == "PU_SH")
pop13 <- subset(data_pop, Pop_ID == "PU_SV")
pop14 <- subset(data_pop, Pop_ID == "PU_TR")

## rigida pops
pop1 <- subset(data_pop, Pop_ID == "RI_BR")
View(data_pop[1:10,1:10])
View(pop1[1:5,1:5])
pop2 <- subset(data_pop, Pop_ID == "RI_CT")
pop3 <- subset(data_pop, Pop_ID == "RI_DT")
pop4 <- subset(data_pop, Pop_ID == "RI_GA")
pop5 <- subset(data_pop, Pop_ID == "RI_GW")
pop6 <- subset(data_pop, Pop_ID == "RI_HH")
pop7 <- subset(data_pop, Pop_ID == "RI_JF")
pop8 <- subset(data_pop, Pop_ID == "RI_KY")
pop9 <- subset(data_pop, Pop_ID == "RI_ME")
pop10 <- subset(data_pop, Pop_ID == "RI_MI")
pop11 <- subset(data_pop, Pop_ID == "RI_NJ")
pop12 <- subset(data_pop, Pop_ID == "RI_NY")
pop13 <- subset(data_pop, Pop_ID == "RI_OH")
pop14 <- subset(data_pop, Pop_ID == "RI_RS")
pop15 <- subset(data_pop, Pop_ID == "RI_SH")
pop16 <- subset(data_pop, Pop_ID == "RI_SP")
pop17 <- subset(data_pop, Pop_ID == "RI_TN")
pop18 <- subset(data_pop, Pop_ID == "RI_TR")
pop19 <- subset(data_pop, Pop_ID == "RI_VT")

# counting homozygous calls for minor allele
pop1_minors <- colSums(pop1=="11", na.rm=TRUE)
pop2_minors <- colSums(pop2== "11", na.rm=TRUE)
pop3_minors <- colSums(pop3=="11", na.rm=TRUE)
pop4_minors <- colSums(pop4=="11", na.rm=TRUE)
pop5_minors <- colSums(pop5=="11", na.rm=TRUE)
pop6_minors <- colSums(pop6=="11", na.rm=TRUE)
pop7_minors <- colSums(pop7=="11", na.rm=TRUE)
pop8_minors <- colSums(pop8=="11", na.rm=TRUE)
pop9_minors <- colSums(pop9=="11", na.rm=TRUE)
pop10_minors <- colSums(pop10=="11", na.rm=TRUE)
pop11_minors <- colSums(pop11=="11", na.rm=TRUE)
pop12_minors <- colSums(pop12=="11", na.rm=TRUE)
pop13_minors <- colSums(pop13=="11", na.rm=TRUE)
pop14_minors <- colSums(pop14=="11", na.rm=TRUE)
### above for just pungens analysis... add the five others below for rigida analysis
pop15_minors <- colSums(pop15=="11", na.rm=TRUE)
pop16_minors <- colSums(pop16=="11", na.rm=TRUE)
pop17_minors <- colSums(pop17=="11", na.rm=TRUE)
pop18_minors <- colSums(pop18=="11", na.rm=TRUE)
pop19_minors <- colSums(pop19=="11", na.rm=TRUE)


# counting heterozygous calls
pop1_het <- colSums(pop1=="12", na.rm=TRUE)
pop2_het <- colSums(pop2== "12", na.rm=TRUE)
pop3_het <- colSums(pop3=="12", na.rm=TRUE)
pop4_het  <- colSums(pop4=="12", na.rm=TRUE)
pop5_het  <- colSums(pop5=="12", na.rm=TRUE)
pop6_het  <- colSums(pop6=="12", na.rm=TRUE)
pop7_het  <- colSums(pop7=="12", na.rm=TRUE)
pop8_het <- colSums(pop8=="12", na.rm=TRUE)
pop9_het <- colSums(pop9=="12", na.rm=TRUE)
pop10_het  <- colSums(pop10=="12", na.rm=TRUE)
pop11_het  <- colSums(pop11=="12", na.rm=TRUE)
pop12_het  <- colSums(pop12=="12", na.rm=TRUE)
pop13_het  <- colSums(pop13=="12", na.rm=TRUE)
pop14_het <- colSums(pop14=="12", na.rm=TRUE)
### above for just pungens analysis... add the five others below for rigida analysis
pop15_het  <- colSums(pop15=="12", na.rm=TRUE)
pop16_het  <- colSums(pop16=="12", na.rm=TRUE)
pop17_het  <- colSums(pop17=="12", na.rm=TRUE)
pop18_het  <- colSums(pop18=="12", na.rm=TRUE)
pop19_het <- colSums(pop19=="12", na.rm=TRUE)

View(pop1_het)
View(pop1_minors)
see <- pop1_minors*2 + pop1_het
see[30:50]

### make data frame where homozygous calls for minor allele are multiplied by 2 for each population
### and the heterozygous calls remain worth the value of 1 minor allele count
PUNG_pops_minorAlleles <- cbind(pop1_minors*2 + pop1_het, pop2_minors*2 + pop2_het, pop3_minors + pop3_het, 
                                pop4_minors*2 + pop4_het, pop5_minors*2 + pop5_het, pop6_minors*2 + pop6_het, 
                                pop7_minors*2 + pop7_het, pop8_minors*2 + pop8_het,pop9_minors*2 + pop9_het, pop10_minors*2 + pop10_het,
                                pop11_minors*2 + pop11_het, pop12_minors*2 + pop12_het, pop13_minors*2 + pop13_het, pop14_minors*2 + pop14_het)

PUNG_pops_minorAlleles <- PUNG_pops_minorAlleles[-1,]
colnames(PUNG_pops_minorAlleles) <- c( "PU_BB", "PU_BN", "PU_BV", "PU_DT", "PU_EG", "PU_EK", "PU_GA", "PU_LG", "PU_NM", "PU_PM", "PU_SC", "PU_SH", "PU_SV", "PU_TR")
View(PUNG_pops_minorAlleles[1:10,1:14])
write.csv(PUNG_pops_minorAlleles, "minor_allele_counts_for_pungens_pops.csv")


RIG_pops_minorAlleles <- cbind(pop1_minors*2 + pop1_het, pop2_minors*2 + pop2_het, pop3_minors + pop3_het, 
                               pop4_minors*2 + pop4_het, pop5_minors*2 + pop5_het, pop6_minors*2 + pop6_het, 
                               pop7_minors*2 + pop7_het, pop8_minors*2 + pop8_het,pop9_minors*2 + pop9_het, pop10_minors*2 + pop10_het,
                               pop11_minors*2 + pop11_het, pop12_minors*2 + pop12_het, pop13_minors*2 + pop13_het, pop14_minors*2 + pop14_het,
                               pop15_minors*2 + pop15_het, pop16_minors*2 + pop16_het, pop17_minors*2 + pop17_het, pop18_minors*2 + pop18_het, pop19_minors*2 + pop19_het)

RIG_pops_minorAlleles <- RIG_pops_minorAlleles[-1,]
colnames(RIG_pops_minorAlleles) <- c( "RI_BR", "RI_CT", "RI_DT", "RI_GA", "RI_GW", "RI_HH", "RI_JF", "RI_KY", "RI_ME", "RI_MI", "RI_NJ", "RI_NY", "RI_OH", "RI_RS", "RI_SH", "RI_SP", "RI_TN", "RI_TR", "RI_VT")
View(RIG_pops_minorAlleles[1:10,1:19])
write.csv(RIG_pops_minorAlleles, "minor_allele_counts_for_rigida_pops_minus104.csv")


############################### Pairwise Fst ##############################

# https://rdrr.io/cran/hierfstat/src/R/pairwise.fst.R
pung_pairFst <- pairwise.WCfst(dat= data_pung, diploid=TRUE)
View(pung_pairFst)
write.csv(pung_pairFst, "pungens_pairwise_WC_Fst.csv")


rig_pairFst <- pairwise.WCfst(dat= data_rig, diploid=TRUE)
View(rig_pairFst)
write.csv(rig_pairFst, "rigida_pairwise_WC_Fst_minus104.csv")

rigidaFstpair <- read.csv("rigida_pairwise_WC_Fst_minus104.csv")
r_Fpw <- rigidaFstpair[,2:19]
r_Fpw <- as.data.frame(r_Fpw)
r_Fpw[is.na(r_Fpw)] <- 0
as.data.frame(r_Fpw)
rf <- as.numeric(r_Fpw)
r_Fpw[which(r_Fpw < 0)] <- 0
min(rigidaFstpair[,2:19], na.rm=T) #-0.008867819 RI_NY and RI_CT
max(rigidaFstpair[,2:19], na.rm=T) # 0.02566537 RI_SH and RI_HH

pungensFstpair <- read.csv("~/Downloads/pungens_pairwise_WC_Fst.csv")
pungensFstpair[,2:15]
min(pungensFstpair, na.rm=T) #-0.01057421 
max(pungensFstpair[,2:15], na.rm=T) # 0.04573286 PU_DT and PU_BB... although all
# PU_DT comparisons were high range = 0.0146-0.0457
pung_num <- pungensFstpair[,2:15]
pung_num <- as.data.frame(pung_num)
pung_num <- as.numeric(pung_num)
mean_Pung_Fstpair <- apply(pung_num, 2, mean, na.rm=TRUE)
View(mean_Pung_Fstpair)
avg_pung <- mean(mean_Pung_Fstpair) 

## lineralize pairwise distances
pung_linear <- apply(pung_pairFst, 1:2, function(x) x/(1-x))
View(pung_linear)

rig_linear <- apply(rig_pairFst, 1:2, function(x) x/(1-x))
View(rig_linear)

write.csv(pung_linear, "pungens_linearized_Fst.csv")
write.csv(rig_linear, "rigida_linearized_Fst_minus104.csv")


## figuring out PU_DT
View(df_gen3[1:50,1:5])
PU_DT <- df_gen3[26:32,]
rig_data <- df_gen3[126:300,]
rig_PU_DT <- rbind(PU_DT,rig_data)
View(rig_PU_DT[1:20,1:20])

MinorToHierf<-function(df,n){
  
  #this script will use a dataframe with snp and ind and pop information
  # the snps need to be coded as counts of the minor allele
  #n indicates the column at which the loci information begins
  
  df<-df[ ,n:ncol(df)]
  
  formated <- apply(df, 2, function(df) gsub('-1','NA',df,fixed=TRUE))
  formated<- apply(formated, 2, function(df) gsub('1','het',df))
  formated <- apply(formated, 2, function(df) gsub('2','min',df))
  formated <- apply(formated, 2, function(df) gsub('0','maj',df))
  
  formated <- apply(formated, 2, function(df) gsub('NA', 'NA', df))
  formated<- apply(formated, 2, function(df) gsub('het','12',df))
  formated <- apply(formated, 2, function(df) gsub('maj','22',df))
  formated <- apply(formated, 2, function(df) gsub('min','11',df))
  
  formated<-apply(formated,2,function(df) as.numeric(df))
  return(data.frame(formated))
}

hierf_rig_PU_DT<-MinorToHierf(rig_PU_DT[,2:2169],n = 1)
View(hierf_rig_PU_DT[1:20,1:20])
hierf_rig_PU_DT <- cbind(rig_PU_DT$Pop, hierf_rig_PU_DT)
rig_PU_DT_pairFst <- pairwise.WCfst(dat= hierf_rig_PU_DT, diploid=TRUE)
write.csv(rig_PU_DT_pairFst, "~/Documents/Ch1_files/PU_DT_FstPairwise_with_RigidaPops.csv")




############################ Mantel test ##########################################                  

############# IBD ##############

install.packages('geosphere')
library(geosphere)
install.packages('vegan')
library(vegan)

# calculate mean locations for each population 
data_loc <- read.csv("~/path/to/Pop_ID_with_location_data.csv")
View(data_loc)
coords <- data_loc[-104,]
View(coords)
mean_lat <- tapply(coords$Latitude, coords$Sp_Pop, mean)
View(mean_lat)
mean_long <- tapply(coords$Longitude, coords$Sp_Pop, mean)
View(mean_long)
pop_mean_coords <- cbind(mean_long, mean_lat)
View(pop_mean_coords)


## P pungens Mantel test
## mean location data 
pung_geo <- pop_mean_coords[1:14,]
View(pung_geo)
write.csv(pung_geo, "pungends_mean_coordinates_perPOP.csv", row.names = TRUE)

#### change to kilometers ####
# geographic distances
pung_dist_matrix <- distm(pung_geo, fun= distVincentyEllipsoid)
pung_in_km <- pung_dist_matrix/1000
write.csv(pung_in_km, "pungens_geo_distances.csv")
View(pung_in_km)

# genetic distances
pung_linear <- read.csv("pungens_linearized_Fst.csv")
View(pung_linear)
pung_linear <- pung_linear[,-1]

pung_mantel <- mantel(pung_in_km, pung_linear, method = "pearson", permutations = 999, strata=NULL)
View(pung_mantel)
pung_mantel$statistic # The Mantel statistic
# -0.07886864
pung_mantel$signif #Empirical significance level from permutations.
# 0.638



## P rigida Mantel test
# mean location data
rig_geo <- pop_mean_coords[15:33,]
View(rig_geo)
write.csv(rig_geo, "~/path/to/rigida_mean_coordinates_perPOP_minus104.csv")

# geographic distances
rig_dist_matrix <- distm(rig_geo, fun=distVincentyEllipsoid)
rig_in_km <- rig_dist_matrix/1000
write.csv(rig_in_km, "rigida_geo_distances_minus104.csv")
View(rig_in_km)

# genetic distances
rig_linear <- read.csv("rigida_linearized_Fst_minus104.csv")
View(rig_linear)
rig_linear <- rig_linear[,-1]

rig_mantel <- mantel(rig_in_km, rig_linear, method = "pearson", permutations = 999, strata=NULL)
View(rig_mantel)
rig_mantel$statistic # The Mantel statistic
# 0.1757708
rig_mantel$signif # Empirical significance level from permutations
# 0.055


############ IBE #############

#calc environmental distance for pungens
library(raster)
library(rgdal)



#Calc IBE from worldClim data

##rigida
setwd("~/Documents/TGG_SDM/WorldClim_data/wc2.1_30s_bio/")
pts_rig <- read.csv("~/Documents/SNPdata_Ch1/rigida_mean_coordinates_perPOP_minus104.csv")
R_loc <- pts_rig[,2:3] 
Data_R <- R_loc
View(Data_R)
tifFiles <- list.files (pattern=".tif")

stacked_tifs <- stack(tifFiles)
Rigida_clim <- raster::extract(stacked_tifs, Data_R)

Rigida_clim <- cbind(Data_R, Rigida_clim)


## calc IBE for Pinus rigida occurrence data

write.csv(Rigida_clim, "~/path/to/Rigida_WorldClim_values_perPOP_minus104.csv")


rigida_env_data <- read.csv("~/path/to/Rigida_WorldClim_values_perPOP_minus104.csv")
rigida_env_data <- rigida_env_data[,4:22]
rigida_env_data
rigENV_scaled <- scale(rigida_env_data, center = TRUE, scale = TRUE)
rigENV_dist <- dist(rigENV_scaled, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

# genetic distances
rig_linear <- read.csv("~/path/to/rigida_linearized_Fst_minus104.csv")
View(rig_linear)
rig_linear <- rig_linear[,-1]

rigENV_mantel <- mantel(rigENV_dist, rig_linear, method = "pearson", permutations = 999, strata=NULL)
View(rigENV_mantel)
rigENV_mantel$statistic
# -0.06694538
rigENV_mantel$signif
#  0.633



## pungens
pts_pung <- read.csv("~/path/to/pungens_mean_coordinates_perPOP.csv")
View(pts_pung)
P_loc <- pts_pung[,2:3]
Data_P <- P_loc

crs(stacked_tifs) <-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

Pungens_clim <- raster::extract(stacked_tifs, Data_P)

Pungens_clim <- cbind(Data_P, Pungens_clim)
write.csv(Pungens_clim, "~/path/to/Pungens_WorldClim_values_perPOP_minus104.csv")

pung_env_data <- read.csv("~/path/to/Pungens_WorldClim_values_perPOP_minus104.csv")
pung_env_data <- pung_env_data[,4:22]
pungENV_scaled <- scale(pung_env_data, center = TRUE, scale = TRUE)
pungENV_dist <- dist(pungENV_scaled, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

# genetic distances
pung_linear <- read.csv("~/path/to/pungens_linearized_Fst.csv")
View(pung_linear)
pung_linear <- pung_linear[,-1]

pungENV_mantel <- mantel(pungENV_dist, pung_linear, method = "pearson", permutations = 999, strata=NULL)
View(pungENV_mantel)
pungENV_mantel$statistic
#  0.013068
pungENV_mantel$signif
# 0.411




#### vcf recode_ removed PU_SH_7 for dadi runs ####
# performed in vcftools
fixed_data <- read.table("single_dadi_minus104.recode.vcf.recode.vcf")
View(fixed_data[1:20,1:20])



############################### popgraph #######################################

install.packages("devtools")
require(devtools)
install.packages("popgraph")
install_github("dyerlab/popgraph")
library(popgraph)


##The genotype data have to be run through to_mv() from gstudio first 
##and then the resulting object can be used with the popgraph() function in popgraph to infer the graph. 


install.packages( c("RgoogleMaps",
                    "geosphere",
                    "proto",
                    "sampling",
                    "seqinr",
                    "spacetime",
                    "spdep"), 
                  dependencies=TRUE )

install_github("dyerlab/gstudio")
library(gstudio)



txt_data2 <- read_population("~/path/to/Two_SP_data_with_NAs_minus104.txt", sep="\t", type = "snp", locus.columns = c(2:2169), header = TRUE)
View(txt_data2[1:20,1:20])
off_mv <- to_mv(txt_data2)

graph <- popgraph( off_mv, groups=txt_data2$Pop)
#Warning message:
#  In popgraph(off_mv, groups = txt_data2$Pop) :
#  32 variables are collinear and being dropped from the discriminant rotation.
plot(graph)







################################## RDA ########################################

########################## Prep files for RDA ##################################


#### extract WorldClim data for all individuals for RDA
library(raster)
setwd("~/Downloads")
elev <- raster("wc2.1_30s_elev.tif")
setwd("~/Downloads/wc2.1_30s_bio/")
bioclim <- list.files(pattern = ".tif")
stack_bioclim <- stack(bioclim)
crs (stack_bioclim) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


sample_coords <- read.csv("~/path/to/Pop_ID_with_location_data.csv")
Longitude <- as.numeric(sample_coords[,6])
Latitude <- as.numeric(sample_coords[,7])
xy <- cbind(Longitude, Latitude)
elev_coords <- cbind(xy, raster::extract(elev, xy))
colnames(elev_coords) <- c("Longitude", "Latitude", "Elevation")

env_coords <- cbind(elev_coords, raster::extract(stack_bioclim, xy))

colnames(env_coords) <- c("Longitude", "Latitude", "Elevation", "Bio1", "Bio10", "Bio11","Bio12", "Bio13", "Bio14", "Bio15", "Bio16",
                          "Bio17", "Bio18", "Bio19", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7", "Bio8", "Bio9")
ID_coords_env <- cbind(sample_coords[,1:4], env_coords)
write.csv(ID_coords_env, "~/Documents/SNPdata_Ch1/TwoSP_WorldClim_30s_values_perINDV.csv")

ID_coords_5env <- cbind(sample_coords[,1:4], ID_coords_env$Bio1, ID_coords_env$Bio12, ID_coords_env$Bio17, ID_coords_env$Bio4, ID_coords_env$Bio17)
colnames(ID_coords_5env) <- c("Sp", "Pop", "ID", "All", "Bio1", "Bio12", "Bio17","Bio4", "Bio7")
write.csv(ID_coords_5env, "~/Documents/SNPdata_Ch1/Pop_ID_env_5vars_30s.csv")

all_RDA_data <- cbind(sample_coords[,1:4], env_coords[,1:3], ID_coords_5env[,5:9])
write.csv(ID_coords_5env, "~/Documents/SNPdata_Ch1/Pop_ID_latlong_env_5vars_30s.csv")


ID_coords_env <- read.csv("~/Documents/SNPdata_Ch1/TwoSP_WorldClim_30s_values_perINDV.csv")


library(corrplot)

env <- env_coords[,4:22]
corr_all <- cor(env, use= "complete.obs", method= "pearson")
corr_all_NA <- round(corr_all, 2)
corr_all_NA[which(abs(corr_all) < 0.75)] <- NA
abs(corr_all_NA)


## 2.5 minute world clim data correlations
setwd("~/Documents/TGG_SDM/WorldClim_data/wc2.1_2.5m_bio/")
bioclim_2 <- list.files(pattern = ".tif")
stack_bioclim_2 <- stack(bioclim_2)
crs (stack_bioclim_2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

env_2 <- raster::extract(stack_bioclim_2, xy)
colnames(env_2) <- c( "Bio1", "Bio10", "Bio11","Bio12", "Bio13", "Bio14", "Bio15", "Bio16",
                          "Bio17", "Bio18", "Bio19", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7", "Bio8", "Bio9")

corr_all_2 <- cor(env_2, use= "complete.obs", method= "pearson")
corr_all_NA2 <- round(corr_all_2, 2)
corr_all_NA2[which(abs(corr_all_2) < 0.75)] <- NA
abs(corr_all_NA2)


#### RDA_connie_TGG.R ####


#straight from forester paper
#https://popgen.nescent.org/2018-03-27_RDA_GEA.html
library(data.table)
install.packages("tidyverse")
library(tidyverse)
library(vegan)
install.packages("psych")
library(psych)
library(viridis)
library(ggsci)
install.packages("ggcorrplot")
library(ggcorrplot)
require(ggrepel)
library(lme4)
library(readr)

source('~/path/to/PCA/Imports.R')

# genetic data 
pos <- read_table2("~/Documents/Ch1_files/SNP_data/single_dadi.recode.vcf.012.pos",col_names = FALSE)
snp_ID <- paste(pos$X1,pos$X2,sep='_')
df012<-fread("~/Desktop/TrevorAnalyses/data/single_dadi_minus104.recode.vcf.recode.vcf.gz.012",sep="\t", data.table=F) 
df012 <- df012[,-1]
dim(df012)
#[1]  300 2168
View(df012[1:20,1:20])
df012_na <- apply(df012, 2, function(df) gsub('-1','NA',df,fixed=TRUE))
View(df012_na[1:20,1:20])
formated<-apply(df012_na,2,function(df) as.numeric(df))


sum(is.na(formated)) #128238
sum(!is.na(formated)) #522162

128238/650528 #0.197... 19.7% missing data

colmean<-apply(formated,2,mean,na.rm=TRUE)
  
normalize<-matrix(nrow = nrow(formated),ncol=ncol(formated))

af<-colmean/2
  
for (m in 1:length(af)){
    nr<-formated[ ,m]-colmean[m]
    dn<-sqrt(af[m]*(1-af[m]))
    normalize[ ,m]<-nr/dn
  }
  
normalize[is.na(normalize)]<-0
scaled_df012 <- as.data.frame(normalize)
View(scaled_df012[1:20,1:20])





#environmental variables 

all_env <- read.csv('~/path/to/TwoSP_WorldClim_30s_values_perINDV.csv')
names(all_env)
all_env <- all_env[,9:27]

corr_all <- cor(all_env, use="complete.obs", method="pearson")

corr_all_NA <- round(corr_all,2)
corr_all_NA[which(abs(corr_all) < .75)] <- NA
abs(corr_all_NA)
write.csv(corr_all, "~/Desktop/TrevorAnalyses/data/ENV_correlations.csv")
write.csv(corr_all_NA, "~/Desktop/TrevorAnalyses/data/ENV_correlations_NAs_r75.csv")


#env correlations
corr <- cor(env, use="complete.obs", method="pearson")
ggcorrplot(corr)

corr_NA <- round(corr,2)
corr_NA[which(abs(corr) < .75)] <- NA
abs(corr_NA)

Bio2 <- all_env[,12]
Bio10 <- all_env[,2]
Bio11 <- all_env[,3]
Bio15 <- all_env[,7]
Bio17 <- all_env[,9]

climate_5var <- cbind(Bio2,Bio10, Bio11, Bio15, Bio17)
clim5_scaled <- apply(climate_5var, 2, scale)
clim_df <-as.data.frame(clim5_scaled)



#### partitioning of variance... varpart() genotype, climate and geo####
geo_env <- read.csv("~/path/to/indv_summary_WorldClim.csv")
View(geo_env)

lat <- geo_env$Latitude
long <- geo_env$Longitude
elev <- geo_env$Elevation
geo_df <- cbind(long, lat, elev)
geo_scaled <- apply(geo_df,2, scale)
geo_scaled <- as.data.frame(geo_scaled)


vp<- varpart(scaled_df012, ~ as.matrix(clim_df), ~ as.matrix(geo_scaled))
vp
#Partition of variance in RDA 

#Call: varpart(Y = scaled_df012, X = ~as.matrix(clim_df),
#~as.matrix(geo_scaled))

#Explanatory tables:
#  X1:  ~as.matrix(clim_df)
#  X2:  ~as.matrix(geo_scaled) 

#No. of explanatory tables: 2 
#Total variation (SS): 1148706 
#Variance: 3841.8 
#No. of observations: 300 

#Partition table:
#                     Df R.squared Adj.R.squared Testable
#[a+b] = X1            5   0.02659       0.01003     TRUE
#[b+c] = X2            3   0.02029       0.01036     TRUE
#[a+b+c] = X1+X2       8   0.04157       0.01522     TRUE
#Individual fractions                                    
#[a] = X1|X2           5                 0.00486     TRUE
#[b]                   0                 0.00518    FALSE
#[c] = X2|X1           3                 0.00519     TRUE
#[d] = Residuals                         0.98478    FALSE
#---
#  Use function ‘rda’ to test significance of fractions of interest

climate_ind_accounts <- (0.00486/0.01522)*100 # % PVE
climate_ind_accounts #31.93167
geo_ind_accounts <- (0.00519/0.01522)*100 # % PVE
geo_ind_accounts # 34.09987
confounded_effect <- ((0.01522-0.01003-0.01036)/0.01522)*100
confounded_effect # -33.96846






####  run rda with bioclim and geo####

env_scaled <- cbind(clim_df,geo_scaled)

m <- rda(formula = scaled_df012 ~.,scale=FALSE, center=TRUE, data = env_scaled)
m
# Call: rda(formula = scaled_df012 ~ Bio2 + Bio10 + Bio11 + Bio15 + Bio17 +
#long + lat + elev, data = env_scaled, scale = FALSE, center = TRUE)

#Inertia Proportion Rank
#Total         3.842e+03  1.000e+00     
#Constrained   1.597e+02  4.157e-02    8
#Unconstrained 3.682e+03  9.584e-01  291
#Inertia is variance 

#Eigenvalues for constrained axes:
#  RDA1  RDA2  RDA3  RDA4  RDA5  RDA6  RDA7  RDA8 
#67.55 15.61 14.43 13.10 12.84 12.74 11.89 11.55 

#Eigenvalues for unconstrained axes:
#  PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
#108.19  32.07  26.93  24.84  23.74  23.48  23.12  22.60 
#(Showing 8 of 291 unconstrained eigenvalues)


anova(m)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = scaled_df012 ~ Bio2 + Bio10 + Bio11 + Bio15 + Bio17 + long + lat + elev, data = env_scaled, scale = FALSE, center = TRUE)
#Df Variance      F Pr(>F)    
#Model      8    159.7 1.5777  0.001 ***
#  Residual 291   3682.1                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m,by="axis")
#Permutation test for rda under reduced model
#Forward tests for axes
#Permutation: free
#Number of permutations: 999

#Model: rda(formula = scaled_df012 ~ Bio2 + Bio10 + Bio11 + Bio15 + Bio17 + long + lat + elev, data = env_scaled, scale = FALSE, center = TRUE)
#         Df Variance      F Pr(>F)    
#RDA1       1     67.6 5.3387  0.001 ***
#  RDA2       1     15.6 1.2335  0.216    
#RDA3       1     14.4 1.1404  0.487    
#RDA4       1     13.1 1.0356  0.970    
#RDA5       1     12.8 1.0146  0.972    
#RDA6       1     12.7 1.0067  0.921    
#RDA7       1     11.9 0.9394  0.997    
#RDA8       1     11.5 0.9125  0.958    
#Residual 291   3682.1                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


plot(m)

RsquareAdj(m)
#$r.squared
# 0.04156999

#$adj.r.squared
# 0.0152214


summary(eigenvals(m,model='constrained'))
#Importance of components:
#                        RDA1     RDA2     RDA3     RDA4     RDA5     RDA6     RDA7  RDA8
#Eigenvalue            67.552 15.60838 14.43046 13.10422 12.83757 12.73846 11.88649  11.5467
#Proportion Explained   0.423  0.09773  0.09036  0.08205  0.08038  0.07976  0.07443  0.0723
#Cumulative Proportion  0.423  0.52072  0.61107  0.69313  0.77351  0.85327  0.92770  1.0000

#summary(eigenvals(m,model='unconstrained'))

screeplot(m) #1

RDA1_gen <- sort(abs(summary(m)$biplot[,1]),decreasing=TRUE)
RDA1_gen
#          elev        lat      Bio15      Bio11       long       Bio2      Bio17 
#       0.76604557 0.63993638 0.61361573 0.54587654 0.44497951 0.38483591 0.32138114 
#  Bio10 
#  0.01682639 

Bio15_pung <- Bio15[1:125]
mean_15pung <- mean(Bio15_pung)#11.3276
median_15pung <- median(Bio15_pung) #10.9869
Bio15_sd_pung <- sd(Bio15_pung) # 1.82639

Bio15_rig <- Bio15[126:300]
mean_15_rig <- mean(Bio15_rig) #14.2307
median_15rig <- median(Bio15_rig) # 12.4749
Bio15_sd_rig <- sd(Bio15_rig) # 3.968987


bio11_pung <- Bio11[1:125]
mean_11_pung <- mean(bio11_pung)# 1.2497
median_11_pung <- median(bio11_pung) # 0.84999
bio11_sd_pung <- sd(bio11_pung) # 2.005

bio11_rig <- Bio11[126:300]
mean_11_rig <- mean(bio11_rig) # -0.78523
median_11_rig <- median(bio11_rig) # 0.28333
bio11_sd_rig <- sd(bio11_rig) # 3.08468

elev_pung <- elev[1:125]
mean_E_pung <- mean(elev_pung)# 724.648
median_E_pung <- median(elev_pung) # 711
elev_sd_pung <- sd(elev_pung) # 224.1727

elev_rig <- elev[126:300]
mean_E_rig <- mean(elev_rig) #399.6857
median_E_rig <- median(elev_rig) # 316
elev_sd_rig <- sd(elev_rig) # 292.2594

RDA2_gen <- sort(abs(summary(m)$biplot[,2]),decreasing=TRUE)
RDA2_gen
#long       Bio2        lat      Bio11      Bio10      Bio17       elev 
#0.79750970 0.76666974 0.65583062 0.62156004 0.52255556 0.47143457 0.23884969 
#Bio15 
#0.07187071 
> 

sort(abs(summary(m)$biplot[,1]) + abs(summary(m)$biplot[,2]),decreasing=TRUE)
#       lat      long     Bio11      Bio2      elev     Bio17     Bio15     Bio10 
#   1.2957670 1.2424892 1.1674366 1.1515056 1.0048953 0.7928157 0.6854864 0.5393820 

setwd("~/Desktop/TrevorAnalyses/Figures/")
saveRDS(m,'RDA_TGG_clim_geo_scaled.RDS')

##### make csv of RDA axis env loadings ####

RDA_df <- as.data.frame(summary(m)$biplot)
write.csv(RDA_df,'RDA_TGG_clim_geo_scaled.csv',row.names = F)

############# PLOT RDA #################

######## ggplot
library(ggsci)
require(ggrepel)



#raw
m <- readRDS('~/Desktop/TrevorAnalyses/Figures/RDA_TGG_clim_geo_scaled.RDS')
env <- read.csv('~/Desktop/TrevorAnalyses/data/TwoSP_WorldClim_30s_values_perINDV.csv')
Pop_ID_Sum <- read.csv('~/Desktop/TrevorAnalyses/data/Pop_ID_Sum.csv')

rda_df <- data.frame(Sp=as.character(Pop_ID_Sum$Sp),Pop = as.character(Pop_ID_Sum$Pop),
                     ID=as.character(Pop_ID_Sum$ID),anc=as.character(Pop_ID_Sum$anc_k2))

m_sum_pve <- summary(m)
RDA1_pve <- paste("RDA1 (",round((m_sum_pve$concont$importance[2,1]*100),digits=2),"%)", sep="")
RDA2_pve <- paste("RDA2 (",round((m_sum_pve$concont$importance[2,2]*100),digits=2),"%)", sep="")
m_sum <- summary(m)
m_sum <- scores(m, display=c("sp","sites", "bp")) 
rda_snp <- as.data.frame(m_sum$species)
rda_indv <- as.data.frame(m_sum$sites)
rda_indv <- cbind(rda_df,rda_indv)
write.csv(rda_indv, "~/path/to/RDA_indv_assignment_loadings_m_clim_geo_scaled.csv")
rda_biplot <- as.data.frame(m_sum$biplot)
rda_biplot$var <- row.names(rda_biplot)

basplot <- plot(m)
mult <- attributes(basplot$biplot)$arrow.mul

colSp <- c('#376BD1','#FAA722')

rda_plot <- ggplot(data=rda_indv, aes(x=RDA1, y=RDA2)) + 
  geom_point(data=rda_indv, aes(x=RDA1, y=RDA2,fill=Sp),pch=21,col='black',size=4) +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult * RDA1 * .8,y = 0, yend = mult * RDA2 * .8),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black",size=1) +
  theme_bw()  

rda_plot_nolab <- rda_plot  + 
  scale_fill_manual(values=colSp) +
  xlab(RDA1_pve) + ylab(RDA2_pve) +
  theme(legend.position = "none",
        axis.text = element_text(size=13), 
        axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

rda_plot_nolab

ggsave('RDA_nolab_5vars_corrected.pdf',rda_plot_nolab,height=5,width=6,units='in')



### With label repel  ###

rda_plot <- ggplot(data=rda_indv, aes(x=RDA1, y=RDA2)) + 
  geom_point(data=rda_indv, aes(x=RDA1, y=RDA2,fill=Sp),pch=21,col='black',size=4) +
  xlim(-6, 6) +
  ylim(-6,6) +
  scale_fill_manual(values=colSp) +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult * RDA1 * 0.7,y = 0, yend = mult * RDA2 * 0.7),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot,
                   aes(x= (mult + mult/8) * RDA1, y = (mult + mult/8) * RDA2, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 4,fontface = "bold") + 
  xlab(RDA1_pve) + ylab(RDA2_pve) +theme_bw()  + 
  theme(#legend.position = "none",
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
rda_plot  

ggsave('RDA_label_clim_geo_scaled.pdf',rda_plot,height=6,width=7,units='in')


#### plot RDA lat long elev####



#### correlate RDA and SDM loadings for climate variables ####

clim_results <- read.csv("~/path/to/RDA_SDM_loadings.csv")
rownames(clim_results) <- clim_results$Variable
clim_results <- clim_results[,-1]

signif_cor2 <- cor.mtest(clim_results, conf.level = .95)
signif_cor2 # only SDM of rigida to pungens was significant (p < 0.05)
#$p
#             RDA1      RDA2   SDM_rigida  SDM_pungens
#RDA1        0.0000000 0.5615040 0.30855785  0.50442858
#RDA2        0.5615040 0.0000000 0.11100532  0.25942625
#SDM_rigida  0.3085578 0.1110053 0.00000000  0.04941535
#SDM_pungens 0.5044286 0.2594263 0.04941535  0.00000000



#### rda with lat long transformations ####

lat_sq <- lat^2
lat_cube <- lat^3
long_sq <- long^2
long_cube <- long^3


geo_df2 <- cbind(long_sq, lat_sq)
geo_df3 <- cbind(long_cube, lat_cube)
geo_df2s <- scale(geo_df2)
geo_df3s <- scale(geo_df3)

geo_df2e <- cbind(long_sq, lat_sq, elev)
geo_df2es <- scale(geo_df2e)
geo_df2es <- as.data.frame(geo_df2es)

geo_df3e <- cbind(long_cube, lat_cube, elev)
geo_df3es <- scale(geo_df3e)
geo_df3es <- as.data.frame(geo_df3es)





geo_trans <- cbind(geo_df, geo_df2, geo_df3)
geo_trans <- scale(geo_trans)

geo_log <- cbind(long_log, lat_log, elev)
geo_log <- scale(geo_log)
geo_log <- as.data.frame(geo_log)

geo_log2 <- cbind(long_log2, lat_log2, elev)
geo_log2 <- scale(geo_log2)
geo_log2 <- as.data.frame(geo_log2)

geo_log10 <- cbind(long_log10, lat_log10, elev)
geo_log10 <- scale(geo_log10)
geo_log10 <- as.data.frame(geo_log10)

geo_logb <- cbind(long_logb, lat_logb, elev)
geo_logb <- scale(geo_logb)
geo_logb <- as.data.frame(geo_logb)

geo_logTrans <- cbind(geo_log[,1:2], geo_log2[,1:2], geo_logb[,1:2], geo_log10)


### when NA is the mean and lat long squared
vp_mean_sq <- varpart(scaled, ~ as.matrix(clim_df), ~ as.matrix(geo_df2es))
vp_mean_sq
#Partition of variance in RDA 



### when NA is the mean and lat long cubed
vp_mean_cube <- varpart(df012_scaled, ~ as.matrix(clim_df), ~ as.matrix(geo_df3es))
vp_mean_cube
#Partition of variance in RDA 

#Call: varpart(Y = df012_scaled, X = ~as.matrix(clim_df),
#              ~as.matrix(geo_df3es))

#Explanatory tables:
#  X1:  ~as.matrix(clim_df)
#X2:  ~as.matrix(geo_df3es) 

#No. of explanatory tables: 2 
#Total variation (SS): 519994 
#Variance: 1739.1 
#No. of observations: 300 

#Partition table:
#  Df R.squared Adj.R.squared Testable
#[a+b] = X1            5   0.02608       0.00952     TRUE
#[b+c] = X2            3   0.02013       0.01020     TRUE
#[a+b+c] = X1+X2       8   0.03954       0.01313     TRUE
#Individual fractions                                    
#[a] = X1|X2           5                 0.00293     TRUE
#[b]                   0                 0.00659    FALSE
#[c] = X2|X1           3                 0.00362     TRUE
#[d] = Residuals                         0.98687    FALSE



geo_trans <- as.data.frame(geo_trans)
n2<- rda(formula = df012_scaled ~ lat + long + elev + lat_sq + lat_cube + long_sq + long_cube ,scale=FALSE, center=TRUE, data = geo_trans)
n2 
plot(n2)



