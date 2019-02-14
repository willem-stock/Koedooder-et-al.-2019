#Koedooder et al. diversity partitioning
#original code provided by J. Baert
data=read.table('growth_rates.txt', sep='\t', header=T, dec='.')
data[,4:6][is.na(data[,4:6])]<-0

# convert denisties to biovolumes
vol               <- c(226.9,501.28,189.56)                                     # cell volumes cyl; sem; nav in µm²
data[,4:6]       <- as.matrix(data[,4:6])%*%diag(vol)

# calculate diversity levels
diversity         <- sapply(1:nrow(data),function(x) 3-length(which(data[x,4:6]==0)) )
data              <- cbind(data,diversity)

#plot all datapoints
plot(data$diversity,log10(rowSums(data[,4:6])),col=data$treatment,xlab="Species richness",ylab="Log productivity",main='Diversity partitioning')
reg       <- lm(log10(rowSums(data[,4:6]))~data$diversity*data$treatment)            
summary(reg)

#axenic####
data_ax      <- data[which(data$treatment=='axenic'),]  

# plot biodiveristy-productivity relationship
plot(data_ax$diversity,log10(rowSums(data_ax[,4:6])),xlab="Species richness",ylab="Log productivity",main='Diversity partitioning',ylim=c(1.8,3.2))
reg       <- lm(log10(rowSums(data_ax[,4:6]))~data_ax$diversity)            
summary(reg)
abline(coef(reg)[1],coef(reg)[2],lty=2)

# calculate the average monoculture function
data.mon              <- data_ax[which(data_ax$diversity==1),]                                        #  extract monoculture data
data.mon$combination  <- factor(data.mon$combination)  
mon                   <- tapply(rowSums(data.mon[,4:6]),data.mon$combination,mean)             #  calculate monoculture functions
mon<-mon[c(1,3,2)]                                             #  reorder monoculture functions according to dataframe

# calculate biodiversity effects
data.com              <- data_ax[which(data_ax$diversity>1),] 

bio.eff               <- array(dim=c(nrow(data.com),5))

for(i in 1:nrow(data.com)){
  RY.e           <- rep(1/data.com$diversity[i],data.com$diversity[i])
  RY.o           <- data.com[i,4:6]/mon
  M              <- mon[which(RY.o>0)]
  RY.o           <- unlist(RY.o[which(RY.o>0)])
 
  delta.Y        <- sum(data.com[i,4:6])-sum(RY.e*M)
  delta.RY       <- RY.o-RY.e
  n.spec         <- data.com$diversity[i]
  
  
  # Loreau & Hector (2001)
  sel            <- n.spec*cov(delta.RY,M)
  comp           <- n.spec*mean(delta.RY)*mean(M)

  
  # store biodiversity effects
  bio.eff[i,1]   <- toString(data.com$combination[i])
  bio.eff[i,2:5] <- c(data.com$diversity[i],delta.Y,sel,comp)
}

bio.eff<-as.data.frame(bio.eff)
names(bio.eff) <- c("combination","diversity","dY","sel","comp")

#with bacterial inoculum (non-axenic)####
data_in      <- data[which(data$treatment=='bacteria'),]  

# plot biodiveristy-productivity relationship for combinations with bacteria (non-axenic)
points(data_in$diversity,log10(rowSums(data_in[,4:6])),col='red')
reg2       <- lm(log10(rowSums(data_in[,4:6]))~data_in$diversity)            
summary(reg2)
abline(coef(reg2)[1],coef(reg2)[2],lty=2,col='red')

# calculate the average monoculture function
data.mon              <- data_in[which(data_in$diversity==1),]
data.mon$combination  <- factor(data.mon$combination)  
mon                   <- tapply(rowSums(data.mon[,4:6]),data.mon$combination,mean)            
mon<-mon[c(1,3,2)]

# calculate biodiversity effects
data.com              <- data_in[which(data_in$diversity>1),] 

bio.eff2             <- array(dim=c(nrow(data.com),5))

for(i in 1:nrow(data.com)){
  RY.e           <- rep(1/data.com$diversity[i],data.com$diversity[i])
  RY.o           <- data.com[i,4:6]/mon
  M              <- mon[which(RY.o>0)]
  RY.o           <- unlist(RY.o[which(RY.o>0)])
  
  delta.Y        <- sum(data.com[i,4:6])-sum(RY.e*M)
  delta.RY       <- RY.o-RY.e
  n.spec         <- data.com$diversity[i]
  
  
  # Loreau & Hector (2001)
  sel            <- n.spec*cov(delta.RY,M)
  comp           <- n.spec*mean(delta.RY)*mean(M)
  
  
  # store biodiversity effects
  bio.eff2[i,1]   <- toString(data.com$combination[i])
  bio.eff2[i,2:5] <- c(data.com$diversity[i],delta.Y,sel,comp)
}

bio.eff2<-as.data.frame(bio.eff2)
names(bio.eff2) <- c("combination","diversity","dY","sel","comp")

#combine the results for the axenic and non-axenic data
bio.eff$composition<-'axenic'
bio.eff2$composition<-'non-axenic'
bio.eff.tot<-rbind(bio.eff,bio.eff2)

