library(reshape)
library(ggplot2)
library(diptest)
library(mousetrap)

###############################################################################
###############################################################################

#Archean carbon data analysis
#Written by: ssz
#Last revised: 1/4/2021

###############################################################################
###############################################################################

setwd("~/Documents/Caltech/Research/Kerogen/LiteratureReview/data analysis")

# Load the data into R
tocData <- read.csv("TOCdata.csv", header = T, fileEncoding = "cp932", stringsAsFactors=FALSE)
kerogenData <- read.csv("kerogenData.csv", header = T, fileEncoding = "cp932", stringsAsFactors=FALSE)
phanerozoicData <-read.csv("phanerozoic_d13C.csv", header = T, stringsAsFactors=FALSE)
carbonateData <-read.csv("carbonates.csv", header = T, stringsAsFactors=FALSE)
aliphaticData <-read.csv("aliphatics.csv", header = T, stringsAsFactors=FALSE)
aromaticData <-read.csv("aromatics.csv", header = T, stringsAsFactors=FALSE)
murchAromaticData <- read.csv("murchison_aromatics.csv", header = T, stringsAsFactors=FALSE)

#split carbonate data into phanerozoic and archean
archeanCarbonateData <- carbonateData[carbonateData$Age..gya. >= 2.25, ]
phanerozoicCarbonateData <- carbonateData[carbonateData$Age..gya. < 2.25, ]

###############################################################################
###############################################################################

# Create arrays for plotting and data analysis
TOC_age <- as.numeric(tocData$Age..Ga.)
TOC_d13C <- as.numeric(tocData$d13C)
TOC_lithology <- tocData$Lithology
TOC_metaGrade <- tocData$Metamorphic.grade
kerogen_age <- as.numeric(kerogenData$Age.age.estimate)
kerogen_d13C <- as.numeric(kerogenData$d13C)
kerogen_lithology <- kerogenData$Lithology
kerogen_metaGrade <- kerogenData$Metamorphic.grade
kerogen_HtoC <- kerogenData$H.C

###############################################################################
###############################################################################

# Compute initial statistics

#mean and standard deviation
mean(TOC_d13C, na.rm=TRUE)
mean(kerogen_d13C, na.rm=TRUE)
sd(TOC_d13C,na.rm=TRUE)
sd(kerogen_d13C,na.rm=TRUE)

#median
median(TOC_d13C, na.rm=TRUE)
median(kerogen_d13C, na.rm=TRUE)

#mean by lithology (TOC, all data)
aggregate(x = tocData$d13C, by = list(tocData$Lithology),FUN = mean)                          

#mean by meta grade (kerogen)
df_meta <- with(kerogenData, kerogenData[!(kerogenData$Metamorphic.grade == "" | is.na(kerogenData$Metamorphic.grade)), ]) #remove rows without meta grade specified
aggregate(x = df_meta$d13C, by = list(df_meta$Metamorphic.grade),function(x) c(mean = mean(x), sd = sd(x))) 

mean(carbonateData$d13C.carb)
sd(carbonateData$d13C.carb)

#comparison to phanerozoic data
phanerozoic_d13C = as.numeric(phanerozoicData$d13C.org)
mean(phanerozoic_d13C)
sd(phanerozoic_d13C)
median(phanerozoic_d13C)

# Use t.test to quantify different assertions of difference
#Are the kerogen and TOC values different from each other and from phanerozoic data?
ks.test(kerogen_d13C, TOC_d13C)
ks.test(kerogen_d13C, phanerozoic_d13C)
ks.test(TOC_d13C, phanerozoic_d13C)


###############################################################################
###############################################################################

# Initial plots of d13C versus age

#plot d13C versus age for TOC and kerogen separately
plot(TOC_age, TOC_d13C,  xlab="Age (Gyr)", ylab="d13C (PDB, per mille)", main="d13C versus age, TOC", xlim = rev(range(TOC_age)))
plot(kerogen_age, kerogen_d13C, xlab="Age (Gyr)", ylab="d13C (PDB, per mille)", main="d13C versus age, kerogen",xlim = rev(range(kerogen_age)))

#comparison to carbonate data (toc+kerogen)
ggplot() + theme_bw() +
  geom_point(data=archeanCarbonateData, aes(x=Age..gya., y=as.numeric(d13C.carb)), color='gray') + 
  geom_point(data=tocData, aes(x=Age..Ga., y=d13C), color='black') + 
  geom_point(data=kerogenData, aes(x=Age.age.estimate, y=d13C), shape=2, color="black")+ 
  labs(x="Age (Gya)", y="d13C (PDB, per mille)", title="Archean organics (TOC + kerogen) versus carbonates")+
  scale_x_reverse()

#plot as a box and whisker plot

#subsample the data from df
#TODO: finish this
austrailiaToc2pt7 <- austrailiaToc[austrailiaToc$Age..Ga.>2.65 & austrailiaToc$Age..Ga.<2.8,]
sAfricaToc2pt7 <- sAfricaToc[sAfricaToc$Age..Ga.>2.65 & sAfricaToc$Age..Ga.<2.8,]
indiaToc2pt7 <- indiaToc[indiaToc$Age..Ga.>2.65 & indiaToc$Age..Ga.<2.8,]
zimbabweToc2pt7 <- zimbabweTOC[zimbabweTOC$Age..Ga.>2.65 & zimbabweTOC$Age..Ga.<2.8,]

ggplot(df_meta, aes(x=factor(Metamorphic.grade), y=d13C)) + geom_boxplot() +
  theme(text=element_text(size=20))+
  scale_x_discrete(limits = c("zeolite", "sub greenschist", "lower greenschist", "greenschist", "upper greenschist", "lower amphibolite", "amphibolite")) +
  labs(x="Metamorphic grade", y="d13C (PDB, per mille)", title="d13C versus metamorphic grade")+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()

###############################################################################
###############################################################################

# Kernel density estimation for kerogen versus TOC measurements
ggplot() + theme_bw()+
  geom_density(mapping=aes(x=as.numeric(phanerozoicData$d13C.org),y=..scaled..), color=rgb(0,0,0), size=1.5, linetype="dotted")+
  geom_density(mapping=aes(x=as.numeric(tocData$d13C),y=..scaled..), color=rgb(0,0,0), size=1.5)+
  geom_density(mapping=aes(x=as.numeric(kerogenData$d13C),y=..scaled..), color=rgb(0,0,0), size=1.5, linetype="dashed")+
  labs(x="d13C", y="density", title="d13C kernel density plot")

#TODO: Tests of bimodality
# Use diptest to evaluate multimodality

dip.test(kerogen_d13C)
dip.test(TOC_d13C)
dip.test(phanerozoic_d13C)

bimodality_coefficient(kerogen_d13C)
bimodality_coefficient(TOC_d13C)

###############################################################################
###############################################################################

# Look at effects of lithology

# plot versus lithology
# for all TOC data
df_lith_TOC = with(tocData, tocData[!(tocData$Lithology == "" | is.na(tocData$Lithology)), ])
ggplot(df_lith_TOC, aes(x=Age..Ga., y=as.numeric(d13C), color=Lithology)) + geom_point()+
  theme_bw()+ theme(text=element_text(size=14))+scale_x_reverse()+
  labs(x="Age (Gya)", y="d13C (PDB, per mille)", title="d13C versus age, by lithology (TOC)") +
  theme(panel.background = element_blank())

# for kerogen measurements
df_lith_kerogen = with(kerogenData, kerogenData[!(kerogenData$Lithology == "" | is.na(kerogenData$Lithology)), ])
ggplot(df_lith_kerogen, aes(x=Age.age.estimate, y=as.numeric(d13C), color=Lithology)) + geom_point()+
  theme_bw()+theme(text=element_text(size=14))+scale_x_reverse()+
  labs(x="Age (Gya)", y="d13C (PDB, per mille)", title="d13C versus age, by lithology (Kerogen)")+
  theme(panel.background = element_blank())
#TODO: change color to match the lithology chart 
###############################################################################
###############################################################################
# Is the kerogen density plot just because certain units/groups have been sampled more?

#(1) First, test by subsampling
#subsample kerogen dataset and plot kernel density functions
N = 50
number_units = 100
M <- matrix(ncol=N, nrow=number_units)
for (i in 1:N) {
  M[,i]<- sample(kerogen_d13C, number_units, replace = FALSE)
}
subset_df = as.data.frame(M)
melted_subset_df = melt(subset_df)

kerogen_density <- ggplot() + theme_bw() +theme(text=element_text(size=14))
kerogen_density <- kerogen_density +  geom_density(melted_subset_df, mapping=aes(x=as.numeric(value), y=..scaled..,color=variable), alpha=0.1)+
  theme(legend.position = "none")+
  labs(x="d13C (PDB, per mille)", y="density", title="kerogen density subsampling")
kerogen_density <- kerogen_density + geom_density(mapping=aes(x=kerogen_d13C,y=..scaled..), color=rgb(0,0,0), size=1.5, linetype="dashed") 
kerogen_density <- kerogen_density + geom_density(mapping=aes(x=TOC_d13C,y=..scaled..), color=rgb(0,0,0), size=1.5)
kerogen_density

#create subsample for TOC to compare to kerogen density subsampling
M2 <- matrix(ncol=N, nrow=number_units)
for (j in 1:N) {
  M2[,j]<- sample(TOC_d13C, number_units, replace=FALSE)
}
subset_df_TOC = as.data.frame(M2)
melted_subset_df_TOC = melt(subset_df_TOC)

#(2) Second, remove some formations at 2.7 to see how that affects the distribution
# look at distributions without tumbiana

sansTumbianaDF = kerogenData[- grep("Tumbiana", kerogenData$Formation),]
sansFortescueDF = sansTumbianaDF[- grep("Fortescue", sansTumbianaDF$Group),]
sansHamersleyDF = sansFortescueDF[- grep("Hamersley", sansFortescueDF$Group),]

ggplot() + theme_bw() +theme(text=element_text(size=20))+
  geom_density(mapping=aes(x=as.numeric(sansTumbianaDF$d13C),y=..scaled..), size=1.5,color="red", linetype="dotted")+
  geom_density(mapping=aes(x=as.numeric(sansFortescueDF$d13C),y=..scaled..), size=1.5,color="blue", linetype="dotted")+
  geom_density(mapping=aes(x=as.numeric(sansHamersleyDF$d13C),y=..scaled..),size=1.5,color="orange", linetype="dotted")+
  geom_density(mapping=aes(x=kerogen_d13C,y=..scaled..), color=rgb(0,0,0), size=1.5, linetype="dashed")+
  geom_density(mapping=aes(x=TOC_d13C,y=..scaled..), color=rgb(0,0,0), size=1.5)+
  labs(x="d13C (PDB, per mille)", y="density", title="density subsampling")
  # = c("ker without Tumbiana", "ker without Tumb+Fortescue", "ker without Tumb+Fort+Hamersley", "all kerogen", "TOC"))
#TODO: add a legend here

# Now is this different from TOC data? do a t-test!
ks.test(as.numeric(sansHamersleyDF$d13C), TOC_d13C)

# (3) Other questions: is this 2.7 gya, tumbiana or shallow environment? 
# are there deep water units that could be tested

###############################################################################
###############################################################################

# Look at Tumbiana specifically  -- red
# TODO: change the colors to match the lithology figure

# are the kerogens versus toc different?
kerogenTumbianaDF = kerogenData[grep("Tumbiana", kerogenData$Formation),]
tocTumbianaDF = tocData[grep("Tumbiana", tocData$Group.Formation),]

tumbianaShale = tocTumbianaDF[grep("shale", tocTumbianaDF$Lithology),]
tumbianaCarbonate = tocTumbianaDF[grep("carbonate", tocTumbianaDF$Lithology),]
tumbianaStrom = tocTumbianaDF[grep("stromatolite", tocTumbianaDF$Lithology),]

ggplot() + theme_bw() +theme(text=element_text(size=20))+
  geom_density(mapping=aes(x=as.numeric(kerogenTumbianaDF$d13C),y=..scaled..), color="black", size=1, linetype="dashed", alpha=0.5)+
  geom_density(mapping=aes(x=as.numeric(tocTumbianaDF$d13C),y=..scaled..), color="black", size=1, alpha=0.5)+
  geom_density(mapping=aes(x=as.numeric(tumbianaShale$d13C),y=..scaled..), color=rgb(52/255,52/255,52/255), linetype="dotted", size=1.5)+
  geom_density(mapping=aes(x=as.numeric(tumbianaCarbonate$d13C),y=..scaled..), color=rgb(156/255,156/255,156/255), linetype="dotted", size=1.5)+
  geom_density(mapping=aes(x=as.numeric(tumbianaStrom$d13C),y=..scaled..), color=rgb(47/255,193/255,153/255), linetype="dotted", size=1.5)+
  labs(x="d13C (PDB, per mille)", y="density", title="tumbiana density")
ks.test(as.numeric(kerogenTumbianaDF$d13C), as.numeric(tocTumbianaDF$d13C))

min(as.numeric(tocTumbianaDF$d13C))

# Fortescue -- blue
kerogenFortescueDF = kerogenData[grep("Fortescue", kerogenData$Group),]
tocFortescueDF = tocData[grep("Fortescue", tocData$Group.Formation),]

fortescueShale = tocFortescueDF[grep("shale", tocFortescueDF$Lithology),]
fortescueCarbonate = tocFortescueDF[grep("carbonate", tocFortescueDF$Lithology),]
fortescueChert = tocFortescueDF[grep("chert", tocFortescueDF$Lithology),]

ggplot() + theme_bw() +theme(text=element_text(size=20))+
  geom_density(mapping=aes(x=as.numeric(kerogenFortescueDF$d13C),y=..scaled..), color="black", size=1, linetype="dashed")+
  geom_density(mapping=aes(x=as.numeric(tocFortescueDF$d13C),y=..scaled..), color="black", size=1)+
  geom_density(mapping=aes(x=as.numeric(fortescueShale$d13C),y=..scaled..), color=rgb(52/255,52/255,52/255), linetype="dotted", size=1.5)+
  geom_density(mapping=aes(x=as.numeric(fortescueCarbonate$d13C),y=..scaled..), color=rgb(156/255,156/255,156/255), linetype="dotted", size=1.5)+
  geom_density(mapping=aes(x=as.numeric(fortescueChert$d13C),y=..scaled..), color=rgb(75/255,0,0), linetype="dotted", size=1.5)+
  labs(x="d13C (PDB, per mille)", y="density", title="fortescue density")

ks.test(as.numeric(kerogenFortescueDF$d13C), as.numeric(tocFortescueDF$d13C))

# Hamersley - orange

kerogenHamersleyDF = kerogenData[grep("Hamersley", kerogenData$Group),]
tocHamersleyDF = tocData[grep("Hamersley", tocData$Group.Formation),]

hamShale = tocHamersleyDF[grep("shale", tocHamersleyDF$Lithology),]
hamCarbonate = tocHamersleyDF[grep("carbonate", tocHamersleyDF$Lithology),]
hamChert = tocHamersleyDF[grep("chert", tocHamersleyDF$Lithology),]

ggplot() + theme_bw() +theme(text=element_text(size=20))+
  geom_density(mapping=aes(x=as.numeric(kerogenFortescueDF$d13C),y=..scaled..),color="black", size=1, linetype="dashed")+
  geom_density(mapping=aes(x=as.numeric(tocHamersleyDF$d13C),y=..scaled..), color="black", size=1)+
  geom_density(mapping=aes(x=as.numeric(hamShale$d13C),y=..scaled..), color=rgb(52/255,52/255,52/255), linetype="dotted", size=1.5)+
  geom_density(mapping=aes(x=as.numeric(hamCarbonate$d13C),y=..scaled..), color=rgb(156/255,156/255,156/255), linetype="dotted", size=1.5)+
  geom_density(mapping=aes(x=as.numeric(hamChert$d13C),y=..scaled..), color=rgb(75/255,0,0), linetype="dotted", size=1.5)+
  labs(x="d13C (PDB, per mille)", y="density", title="hamersley density")

ks.test(as.numeric(kerogenHamersleyDF$d13C), as.numeric(tocHamersleyDF$d13C))

# compare Australia, South Africa and Indian rock TOC at ~2.7Gya

austrailiaToc = tocData[grep("Austrailia", tocData$Location),]
sAfricaToc = tocData[grep("South Africa", tocData$Location),]
indiaToc =tocData[grep("India", tocData$Location),]
zimbabweTOC= tocData[grep("Zimbabwe", tocData$Location),]

austrailiaToc2pt7 <- austrailiaToc[austrailiaToc$Age..Ga.>2.65 & austrailiaToc$Age..Ga.<2.8,]
sAfricaToc2pt7 <- sAfricaToc[sAfricaToc$Age..Ga.>2.65 & sAfricaToc$Age..Ga.<2.8,]
indiaToc2pt7 <- indiaToc[indiaToc$Age..Ga.>2.65 & indiaToc$Age..Ga.<2.8,]
zimbabweToc2pt7 <- zimbabweTOC[zimbabweTOC$Age..Ga.>2.65 & zimbabweTOC$Age..Ga.<2.8,]

ggplot() + theme_bw() +theme(text=element_text(size=20))+
  geom_density(mapping=aes(x=as.numeric(austrailiaToc2pt7$d13C),y=..scaled..), color=rgb(153/255,197/255,83/255), size=1.5)+
  geom_density(mapping=aes(x=as.numeric(sAfricaToc2pt7$d13C),y=..scaled..), color=rgb(38/255,173/255,155/255), size=1.5)+
  geom_density(mapping=aes(x=as.numeric(indiaToc2pt7$d13C),y=..scaled..), color=rgb(44/255,20/255,255/255), size=1.5)+
  geom_density(mapping=aes(x=as.numeric(zimbabweToc2pt7$d13C),y=..scaled..), color=rgb(0/255,18/255,112/255), size=1.5)+
  labs(x="d13C (PDB, per mille)", y="density", title="~2.7ga rocks, austrailia (red) versus s africa (black), india(blue), and zimbabwe (purple)")

###############################################################################
###############################################################################

# Look at effects of metamorphic grades

# Plot versus H/C ratio
plot(kerogen_HtoC,kerogen_d13C, xlab="H/C ratio", ylab="d13C (PDB, per mille)", main="kerogen d13C versus H/C")

# Plot versus metamorphic grade
ggplot(df_meta, aes(x=factor(Metamorphic.grade), y=d13C)) + geom_boxplot() +
  theme(text=element_text(size=20))+
  scale_x_discrete(limits = c("zeolite", "sub greenschist", "lower greenschist", "greenschist", "upper greenschist", "lower amphibolite", "amphibolite")) +
  labs(x="Metamorphic grade", y="d13C (PDB, per mille)", title="d13C versus metamorphic grade")+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()



###############################################################################
###############################################################################

# Plot aromatics and aliphatics separately

#aromatics and aliphatics
ggplot() + theme_bw()+ scale_color_grey()+theme(text=element_text(size=20))+
  geom_point(aliphaticData, mapping=aes(x=aliphaticData$Exact.age..or.average.of.range.if.given..Ga., y=as.numeric(aliphaticData$d13C..VPBD.), color=aliphaticData$Soluble.Insoluble..HyPy.), shape=15, size=6) + 
  geom_point(aromaticData, mapping=aes(x=aromaticData$Specific.age..or.avg.if.range.given..Ga., y=as.numeric(aromaticData$d13C..VPDB.),  color=aromaticData$Soluble.Insoluble..HyPy.),shape=18, size=6)+
  labs(x="Age (Gya)", y="d13C (PDB, per mille)", title="d13C versus age, aliphatics and aromatics")+
  scale_x_reverse() + theme(legend.position="none")

#plot aromatic by type
#ggplot() +  theme_bw()+scale_color_grey()+theme(text=element_text(size=20))+
#  geom_point(aromaticData,mapping=aes(x=aromaticData$PAH, y=aromaticData$d13C..VPDB., color=aromaticData$Soluble.Insoluble..HyPy.), shape=18, size=3)+
#  labs(x="Aromatic molecule", y="d13C (PDB, per mille)", title="d13C aromatics by molecule, soluble and hy-py produced")+
#  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="none")

#plot aliphatic by type
ggplot() +  theme_bw()+scale_color_ordinal()+theme(text=element_text(size=20))+
  geom_point(aliphaticData,mapping=aes(x=aliphaticData$X..of.carbons, y=aliphaticData$d13C..VPBD., color=aliphaticData$Formation),  shape=15, size=6)+
  labs(x="Carbon number", y="d13C (PDB, per mille)", title="d13C aliphatics, Carawine Dolomite and Dresser Fm")+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="none")

# Comparison of aromatics with Murchison measurements
ggplot() +  theme_bw()+scale_color_grey()+theme(text=element_text(size=16))+
  geom_point(aromaticData,mapping=aes(x=aromaticData$PAH, y=aromaticData$d13C..VPDB., color=aromaticData$Soluble.Insoluble..HyPy.), shape=18, size=4)+
  geom_point(murchAromaticData,mapping=aes(x=murchAromaticData$PAH, y=murchAromaticData$d13C, color=murchAromaticData$type), shape=23,size=3)+
  labs(x="Aromatic molecule", y="d13C (PDB, per mille)", title="d13C aromatics by molecule, comparison with murch")+
  theme(axis.text.x = element_text(angle = 90))+ theme(legend.position="none")


