################# Deriving a model for quantitative song variation in a hybrid zone chickadee population #################

# Get acoustic measurement data from master's thesis
list.files()
data0 <- read.csv("SPPUA_song_level_measurements_final.csv")
# need to sample to 100 songs per bird
data <- data0[sample(which(data0$ind_ID==unique(data0$ind_ID)[1]), size = 100, replace = FALSE),]

for (i in 2:length(unique(data0$ind_ID))) {
  a <- data0[sample(which(data0$ind_ID==unique(data0$ind_ID)[i]), size = 100, replace = FALSE),]
  data <- rbind(data, a)
}
dim(data)

#### PCA with all measurements ####
colnames(data)
pca.all <- prcomp(data[,c(3:8,10:14)],
                  scale. = T,
                  center = T)
pca.all
summary(pca.all)
pcscores <- as.data.frame(pca.all$x)


#### Get distance data ####
library(kde1d)
library(vegan)

## make a function that will get mst distances using n songs for each individual x
gettrees2 <- function(x, n) {
  score <- pcscores[which(data$ind_ID==x),]
  rare <- score[sample(1:nrow(score), n), ]
  distmat <- dist(rare[,1:2])
  tree <- spantree(distmat)
  return(tree$dist)
}
gettrees2(x=unique(data$ind_ID)[1], 
          n=100) # it works



#### KDE for all individuals and sample sizes ####
# Make an empty 2-layered list. Outer layer is individuals. Inner layer is rarefaction levels.
# NOTE: rarefaction analysis not completed; will only use the full dataset moving forward
inds <- unique(data$ind_ID)
final.trees <- vector("list",10)
names(final.trees) <- inds

rare.levels <- c(100,90,80,70,60,50,40,30)
empty.rare.list <- vector("list",8)
names(empty.rare.list) <- as.character(rare.levels)

# now fill the list with values for the corresponding individual and rarefaction level
for (i in 1:length(final.trees)) {
  final.trees[[i]] <- empty.rare.list
  for (j in 1:length(rare.levels)) {
    final.trees[[i]][[j]] <- gettrees2(x=inds[i],
                                       n=rare.levels[j])
  }
}

final.trees


################################### Model Selection ###################################
#### Null model: one distribution for all individuals ####
alltrees <- c()
for (i in 1:length(inds)) {
  current <- final.trees[[inds[i]]][[1]]
  alltrees <- c(alltrees, current)
}
alltrees

H0kde <- kde1d(alltrees
)
# separate the steps
H0LL <- logLik(H0kde)
H0p <- H0kde$edf
H0SS <- H0kde$nobs
BICH0_0 <- -2*H0LL + H0p*log(H0SS)
# use BIC function
BICH0 <- BIC(H0kde)

BICH0_0
BICH0 # same result

#### Model 1: Individual distributions ####
all.kdes <- vector("list",10)
names(all.kdes) <- inds

for (i in 1:length(inds)) {
  all.kdes[[i]] <- kde1d(final.trees[[i]][[1]]
  )
}

# overall likelihood score
H1LL <- 0
for(i in 1:10){
  current <- logLik(all.kdes[[i]])
  H1LL <- H1LL + current
}
H1LL

# number of parameters estimated
H1p <- 0
for(i in 1:10){
  current <- all.kdes[[i]]$edf
  H1p <- H1p + current
}
H1p

# number of observations
H1SS <- 0
for(i in 1:10){
  current <- all.kdes[[i]]$nobs
  H1SS <- H1SS + current
}
H1SS

# H1 BIC
BICH1 <- -2*H1LL + H1p*log(H1SS)
BICH1

#### Model 2: Distributions based on genetic ancestry ####
gendata <- read.csv("All_Data_1.csv")
colnames(gendata)
gendata <- gendata[,c(2,3,5)]

# we will use P(assignment to CACH) < 0.9 to classify hybrids
inds.hy <- gendata$group[which(gendata$prob_CA<0.9)]
inds.cc <- gendata$group[which(gendata$prob_CA>0.9)]

# group the trees by genetic grouping
hytrees <- c()
for (i in 1:length(inds.hy)) {
  current <- final.trees[[which(names(all.kdes) %in% inds.hy)[i]]][[1]]
  hytrees <- c(hytrees, current)
}
hytrees

cctrees <- c()
for (i in 1:length(inds.cc)) {
  current <- final.trees[[which(names(all.kdes) %in% inds.cc)[i]]][[1]]
  cctrees <- c(cctrees, current)
}
cctrees

# get kdes
hykde <- kde1d(hytrees
)
cckde <- kde1d(cctrees
)

# log-likelihood
H2LL <- logLik(hykde) + logLik(cckde)

# number of parameters estimated
H2p <- hykde$edf + cckde$edf

# number of observations
H2SS <- hykde$nobs + cckde$nobs

# H2 BIC
BICH2 <- -2*H2LL + H2p*log(H2SS)
BICH2

