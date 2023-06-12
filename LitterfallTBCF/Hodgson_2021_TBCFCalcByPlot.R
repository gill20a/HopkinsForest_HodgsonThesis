urlfile = "https://raw.githubusercontent.com/gill20a/HopkinsForest_HodgsonThesis/main/LitterfallTBCF/Hodgson_2021_Litterfall.csv"
litter <- read.csv(url(urlfile))
litter_sum<-aggregate(LitterMass ~ Plot, data = litter, FUN = sum)
litter_sum$LitterMass_gC_m2<-(litter_sum$LitterMass/0.145)/2
litter_sum$Site<-c(rep("Lower", 24), rep("Upper", 24))


urlfile = "https://raw.githubusercontent.com/gill20a/HopkinsForest_HodgsonThesis/main/LitterfallTBCF/Hodgson_SeasonalSoilRespiration.csv"
resp <- read.csv(url(urlfile))
head(resp)
colnames(resp)<-c("Plot", "SeasonalResp")
head(resp)
resp$Site<-c(rep("Lower", 24), rep("Upper", 24))

resp$LitterMass_gC_m2<-litter_sum$LitterMass_gC_m2
hist(resp$SeasonalResp)
resp$SeasonalResp[resp$SeasonalResp > 3000] <- NA

resp$TBCF<-resp$SeasonalResp-litter_sum$LitterMass_gC_m2
boxplot(resp$TBCF~resp$Site)
summary(lm(resp$TBCF~resp$Site))

library(Matrix)
meanresp<-tapply(resp$SeasonalResp, resp$Site, mean, na.rm=T)
upper<-tapply(resp$SeasonalResp, resp$Site, quantile, 0.975, na.rm=T)
lower<-tapply(resp$SeasonalResp, resp$Site, quantile, 0.025, na.rm=T)
meanresp
upper
lower

fit1<-lm(resp$SeasonalResp~resp$Site)
summary(fit1)

meanlitter<-tapply(litter_sum$LitterMass_gC_m2, litter_sum$Site, mean, na.rm=T)
upper<-tapply(litter_sum$LitterMass_gC_m2, litter_sum$Site, quantile, 0.975, na.rm=T)
lower<-tapply(litter_sum$LitterMass_gC_m2, litter_sum$Site, quantile, 0.025, na.rm=T)
meanlitter
upper
lower

meanTBCF<-tapply(resp$TBCF, resp$Site, mean, na.rm=T)
upper<-tapply(resp$TBCF, resp$Site, quantile, 0.975, na.rm=T)
lower<-tapply(resp$TBCF, resp$Site, quantile, 0.025, na.rm=T)
meanTBCF
upper
lower

