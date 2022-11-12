
##setting the working directory and loading the data spreadhseet, as well as calling the required packages
urlfile = "https://raw.githubusercontent.com/gill20a/HopkinsForest_HodgsonThesis/main/2021_11_01_Timesandfluxes.csv"
Flux<-read.csv(url(urlfile))


head(Flux)
# install.packages("Hmisc")
# install.packages("Matrix")
library(Hmisc)
library(Matrix)



##Combine the separate date and time columns into one which R can easily work with
Flux$Time_1 <- as.POSIXct(paste(Flux$date1, Flux$time1), format="%Y-%m-%d %H:%M")
Flux$Time_2 <- as.POSIXct(paste(Flux$date2, Flux$time2), format="%Y-%m-%d %H:%M")
Flux$Time_3 <- as.POSIXct(paste(Flux$date3, Flux$time3), format="%Y-%m-%d %H:%M")


#Calculate amount of time elapsed during measurements
Flux$TimeInt_0<-as.numeric(as.character(difftime(Flux$Time_1, Flux$Time_1), format="%Y-%m-%d %H:%M")) #Should be a vector full of zeros
Flux$TimeInt_1<-as.numeric(as.character(difftime(Flux$Time_2, Flux$Time_1), format="%Y-%m-%d %H:%M"))
Flux$TimeInt_2<-as.numeric(as.character(difftime(Flux$Time_3, Flux$Time_1), format="%Y-%m-%d %H:%M"))

##Calculate change in C in between measurements
Flux$CInt_0<-(Flux$c1-Flux$c1)
Flux$CInt_1<-(Flux$c2-Flux$c1)
Flux$CInt_2<-(Flux$c3-Flux$c1)


##Assigning a site variable

Flux$Site <-c(rep("Lower", 36), rep("Upper", 36))


##making sub-parts of the larger matrix, first by time relative to addition
BeforeFlux<-subset(Flux, Flux[, "Set"]=="Before")
DuringFlux<-subset(Flux, Flux$Set=="During")
AfterFlux<-subset(Flux, Flux$Set=="After")


###Calculate slope and R2 for all jars - This is calculating the flux rate and a measure of error on that
## for each individual jar used in this experiment.
JarSlope <- numeric(length = length(Flux$CInt_1))
for (i in seq_along(JarSlope)) {
        JarSlope[i] <- coef(lm(unlist(Flux[i,c("CInt_0","CInt_1","CInt_2")])~
                                       unlist(Flux[i,c("TimeInt_0","TimeInt_1","TimeInt_2")]), na.action=na.omit))[2] 
        
}
JarSlope
Flux$JarSlope<-JarSlope

JarR2 <- numeric(length = length(Flux$CInt_1))
for (i in seq_along(JarR2)) {
        JarR2[i] <- summary(lm(unlist(Flux[i,c("CInt_0","CInt_1","CInt_2")])~
                                       unlist(Flux[i,c("TimeInt_0","TimeInt_1","TimeInt_2")]), na.action=na.omit))$r.squared 
        
        
}
JarR2
min(JarR2)
Flux$JarR2<-JarR2

##Find slopes per g of dry soil, and propagate the R2 as error assuming no error on the mass
Flux$Slope_gsoil <- (Flux$JarSlope/Flux$DryMass_g)

Flux$R2_gsoil <- (Flux$JarR2/Flux$DryMass_g) 

##Find slopes in nmol co2 per g of dry soil (per hour)
Flux$molslope <- ((Flux$Slope_gsoil*(4.0874*10^(-4)))*1000)
Flux$R2_mol <- ((Flux$R2_gsoil*(4.0874*10^(-4)))*1000)



##Now find response ratios for the during and after measurements relative to the before measurement for each jar

Flux$Rep <- rep(c(rep(c(1,2),2),1,1),36)

Flux_mean <- aggregate(Flux$molslope, by = data.frame(Flux$Set,Flux$Addition, Flux$Plot), mean, na.rm=T)
head(Flux_mean)
Flux_R2 <- aggregate(Flux$R2_mol, by = data.frame(Flux$Set,Flux$Addition, Flux$Plot), mean, na.rm=T)
head(Flux_R2)
Flux_mean$R2<-Flux_R2$x
Flux_mean$molslope<-Flux_mean$x

Flux_mean_G <- subset(Flux_mean, Flux_mean$Flux.Addition=="G")
Flux_mean_V <- subset(Flux_mean, Flux_mean$Flux.Addition=="V")
Flux_mean_L <- subset(Flux_mean, Flux_mean$Flux.Addition=="L")
Flux_mean_X <- subset(Flux_mean, Flux_mean$Flux.Addition=="X")

Flux_mean_G$SlopeRR <- (Flux_mean_G$molslope/Flux_mean_X$molslope)
Flux_mean_V$SlopeRR <- (Flux_mean_V$molslope/Flux_mean_X$molslope)
Flux_mean_L$SlopeRR <- (Flux_mean_L$molslope/Flux_mean_X$molslope)
Flux_mean_X$SlopeRR <- (Flux_mean_X$molslope/Flux_mean_X$molslope)

##corresponds to supplemental figure X?
plot(Flux_mean_G$Flux.Plot[Flux_mean_G$Flux.Set=="During"], Flux_mean_G$SlopeRR[Flux_mean_G$Flux.Set=="During"], col = "forest green", pch=16,xlab = "Plot Number",ylab = "Flux Rate Response Ratio", ylim=c(0,max(Flux_mean_G$SlopeRR)))
points(Flux_mean_V$Flux.Plot[Flux_mean_V$Flux.Set=="During"], Flux_mean_V$SlopeRR[Flux_mean_V$Flux.Set=="During"], col = "orange", pch=16)
points(Flux_mean_L$Flux.Plot[Flux_mean_L$Flux.Set=="During"], Flux_mean_L$SlopeRR[Flux_mean_L$Flux.Set=="During"], col = "purple", pch=16)
abline(h=1, lty=2, col = "gray", lwd=3)

plot(Flux_mean_G$Flux.Plot[Flux_mean_G$Flux.Set=="After"], Flux_mean_G$SlopeRR[Flux_mean_G$Flux.Set=="After"], col = "forest green", pch=16,xlab = "Plot Number",ylab = "Flux Rate Response Ratio", ylim=c(0,max(Flux_mean_G$SlopeRR)))
points(Flux_mean_V$Flux.Plot[Flux_mean_V$Flux.Set=="After"], Flux_mean_V$SlopeRR[Flux_mean_V$Flux.Set=="After"], col = "orange", pch=16)
points(Flux_mean_L$Flux.Plot[Flux_mean_L$Flux.Set=="After"], Flux_mean_L$SlopeRR[Flux_mean_L$Flux.Set=="After"], col = "purple", pch=16)
abline(h=1, lty=2, col = "gray", lwd=3)



FluxResponseB <- subset(Flux, Flux$Set == "Before")
FluxResponseD <- subset(Flux, Flux$Set == "During")
FluxResponseA <- subset(Flux, Flux$Set == "After")

Fluxduringx <- subset(FluxResponseD, FluxResponseD$Addition=="X")
ControlFlux<-c(rep(Fluxduringx$Slope_gsoil[1],6), rep(Fluxduringx$Slope_gsoil[2],6), rep(Fluxduringx$Slope_gsoil[3],6), rep(Fluxduringx$Slope_gsoil[4],6),
  rep(Fluxduringx$Slope_gsoil[5],6),rep(Fluxduringx$Slope_gsoil[6],6), rep(Fluxduringx$Slope_gsoil[7],6), rep(Fluxduringx$Slope_gsoil[8],6), 
  rep(Fluxduringx$Slope_gsoil[9],6), rep(Fluxduringx$Slope_gsoil[10],6), rep(Fluxduringx$Slope_gsoil[11],6), rep(Fluxduringx$Slope_gsoil[12],6))

FluxResponseB$SlopeRR <- (FluxResponseB$Slope_gsoil/FluxResponseB$Slope_gsoil)

FluxResponseD$SlopeRR <- (FluxResponseD$Slope_gsoil/ControlFlux)

FluxResponseA$SlopeRR <- (FluxResponseA$Slope_gsoil/FluxResponseB$Slope_gsoil)

## This code generates Manuscript Figure 4

boxplot(FluxResponseD$SlopeRR~FluxResponseD$Site*FluxResponseD$Addition, xlim=c(0.7 ,6.3), xaxt="n", xlab="Substrate", ylab="Flux Response Ratio", col=c(rgb(0.65,0.15,0.2), rgb(0,0.447,0.698)), cex.axis=1.5, cex.lab=2)
axis(side=1, at=c(1.65,3.65,5.45), labels=c("Glucose", "Lignin", "Vanillin"), cex.axis=1.5)
legend("topright", c("Upper", "Lower"), col=c(rgb(0, 0.447, 0.698, 1), rgb(0.65, 0.15, 0.2, 1)), pch=16, bty="n", cex=1.5)

##Average duplicates into a new matrix


colnames(Flux)
Flux$Site<-rep((c(rep("Lower", 36), rep("Upper", 36))),3)
Lower = subset(Flux, Flux$Site=="Lower")
Upper = subset(Flux, Flux$Site=="Upper")

pick1<-Lower$Set=="Before"
BeforeLower<-tapply(Lower[pick1,'molslope'], Lower[pick1,"Addition"], mean, na.rm=T)
pick1<-Lower$Set=="During"
DuringLower<-tapply(Lower[pick1,'molslope'], Lower[pick1,"Addition"], mean, na.rm=T)
pick1<-Lower$Set=="After"
AfterLower<-tapply(Lower[pick1,'molslope'], Lower[pick1,"Addition"], mean, na.rm=T)


pick1<-Lower$Set=="Before"
BeforeLower_se<-tapply(Lower[pick1,'molslope'], Lower[pick1,"Addition"], sd, na.rm=T)/sqrt(36)
pick1<-Lower$Set=="During"
DuringLower_se<-tapply(Lower[pick1,'molslope'], Lower[pick1,"Addition"], sd, na.rm=T)/sqrt(36)
pick1<-Lower$Set=="After"
AfterLower_se<-tapply(Lower[pick1,'molslope'], Lower[pick1,"Addition"], sd, na.rm=T)/sqrt(36)

BeforeLower_Upper<-BeforeLower+BeforeLower_se
DuringLower_Upper<-DuringLower+DuringLower_se
AfterLower_Upper<-AfterLower+AfterLower_se

BeforeLower_Lower<-BeforeLower-BeforeLower_se
DuringLower_Lower<-DuringLower-DuringLower_se
AfterLower_Lower<-AfterLower-AfterLower_se

summary<-as.data.frame(rbind(BeforeLower, DuringLower, AfterLower,
      BeforeLower_Upper, DuringLower_Upper, AfterLower_Upper,
      BeforeLower_Lower, DuringLower_Lower, AfterLower_Lower))

summary_means<-summary[1:3,]
summary_UB<-summary[4:6,]
summary_LB<-summary[7:9,]
position=c(1,2,3)

## This code generates two figures - supplementals??

par(mar=c(4,12,4,12))
plot(summary_means$G, ylim=c(0,5), col="forest green", pch=16, xlab = ".." , ylab = expression(paste("", n, "mol ",CO[2]," ",hr^-1," ",per, " ", g, " ", dry, " ", soil, "", sep="")), xaxt="n", cex.axis=1.5, cex.lab=2)
arrows(position, summary_means$G, position, summary_UB$G, angle=90,length=0.05, lwd=2, col="forest green")
arrows(position,summary_means$G, position, summary_LB$G, angle=90,length=0.05, lwd=2, col="forest green")
points(summary_means$G,  col="forest green", pch=18, cex=2)


points(summary_means$V, ylim=c(min(AfterLower_Lower), max(DuringLower_Upper)), col="Orange", pch=16)
arrows(position, summary_means$V, position, summary_UB$V, angle=90,length=0.05, lwd=2, col="orange")
arrows(position,summary_means$V, position, summary_LB$V, angle=90,length=0.05, lwd=2, col="orange")
points(summary_means$V,  col="Orange", pch=18, cex=2)

points(summary_means$L, ylim=c(min(AfterLower_Lower), max(DuringLower_Upper)), col="purple", pch=16)
arrows(position, summary_means$L, position, summary_UB$L, angle=90,length=0.05, lwd=2, col="purple")
arrows(position,summary_means$L, position, summary_LB$L, angle=90,length=0.05, lwd=2, col="purple")
points(summary_means$L,  col="purple", pch=18, cex=2)

points(summary_means$X, ylim=c(min(AfterLower_Lower), max(DuringLower_Upper)), col="gray40", pch=16)
arrows(position, summary_means$X, position, summary_UB$X, angle=90,length=0.05, lwd=2, col="gray40")
arrows(position,summary_means$X, position, summary_LB$X, angle=90,length=0.05, lwd=2, col="gray40")
points(summary_means$X,  col="Gray40", pch=18, cex=2)

legend("topright", c("Glucose", "Vanillin", "Lignin", "Control"), col=c("forest green", "Orange", "purple", "gray40"), pch=18, bty="n", cex=1.4)
axis(1, c("Before", "0 Hours", "72 Hours"), at=c(1,2,3), cex.axis=1.5)

pick1<-Upper$Set=="Before"
BeforeUpper<-tapply(Upper[pick1,'molslope'], Upper[pick1,"Addition"], mean, na.rm=T)
pick1<-Upper$Set=="During"
DuringUpper<-tapply(Upper[pick1,'molslope'], Upper[pick1,"Addition"], mean, na.rm=T)
pick1<-Upper$Set=="After"
AfterUpper<-tapply(Upper[pick1,'molslope'], Upper[pick1,"Addition"], mean, na.rm=T)


pick1<-Upper$Set=="Before"
BeforeUpper_se<-tapply(Upper[pick1,'molslope'], Upper[pick1,"Addition"], sd, na.rm=T)/sqrt(36)
pick1<-Upper$Set=="During"
DuringUpper_se<-tapply(Upper[pick1,'molslope'], Upper[pick1,"Addition"], sd, na.rm=T)/sqrt(36)
pick1<-Upper$Set=="After"
AfterUpper_se<-tapply(Upper[pick1,'molslope'], Upper[pick1,"Addition"], sd, na.rm=T)/sqrt(36)

BeforeUpper_Upper<-BeforeUpper+BeforeUpper_se
DuringUpper_Upper<-DuringUpper+DuringUpper_se
AfterUpper_Upper<-AfterUpper+AfterUpper_se

BeforeUpper_Lower<-BeforeUpper-BeforeUpper_se
DuringUpper_Lower<-DuringUpper-DuringUpper_se
AfterUpper_Lower<-AfterUpper-AfterUpper_se

summary<-as.data.frame(rbind(BeforeUpper, DuringUpper, AfterUpper,
                             BeforeUpper_Upper, DuringUpper_Upper, AfterUpper_Upper,
                             BeforeUpper_Lower, DuringUpper_Lower, AfterUpper_Lower))

summary_means<-summary[1:3,]
summary_UB<-summary[4:6,]
summary_LB<-summary[7:9,]
position=c(1,2,3)


par(mar=c(4,12,4,12))
plot(summary_means$G, ylim=c(0,5), col="forest green", pch=16, xaxt="n", xlab="", ylab = expression(paste("", n, "mol ",CO[2]," ",hr^-1," ",per, " ", g, " ", dry, " ", soil, "", sep="")), cex.axis=1.5, cex.lab=2)
arrows(position, summary_means$G, position, summary_UB$G, angle=90,length=0.05, lwd=2, col="forest green")
arrows(position,summary_means$G, position, summary_LB$G, angle=90,length=0.05, lwd=2, col="forest green")
points(summary_means$G,  col="forest green", pch=18, cex=2)
axis(1, c("Before", "0 Hours", "72 Hours"), at=c(1,2,3), cex.axis=1.5)
legend("topright", c("Glucose", "Vanillin", "Lignin", "Control"), col=c("forest green", "Orange", "purple", "Gray40"), pch=18, bty="n", cex=1.4)

points(summary_means$V, ylim=c(min(AfterUpper_Lower), max(DuringUpper_Upper)), col="Orange", pch=16)
arrows(position, summary_means$V, position, summary_UB$V, angle=90,length=0.05, lwd=2, col="orange")
arrows(position,summary_means$V, position, summary_LB$V, angle=90,length=0.05, lwd=2, col="orange")
points(summary_means$V,  col="Orange", pch=18, cex=2)

points(summary_means$L, ylim=c(min(AfterUpper_Lower), max(DuringUpper_Upper)), col="purple", pch=16)
arrows(position, summary_means$L, position, summary_UB$L, angle=90,length=0.05, lwd=2, col="purple")
arrows(position,summary_means$L, position, summary_LB$L, angle=90,length=0.05, lwd=2, col="purple")
points(summary_means$L,  col="purple", pch=18, cex=2)

points(summary_means$X, ylim=c(min(AfterUpper_Lower), max(DuringUpper_Upper)), col="Gray40", pch=16)
arrows(position, summary_means$X, position, summary_UB$X, angle=90,length=0.05, lwd=2, col="gray40")
arrows(position,summary_means$X, position, summary_LB$X, angle=90,length=0.05, lwd=2, col="gray40")
points(summary_means$X,  col="Gray40", pch=18, cex=2)


###########################################################################
##making sub-parts of the big matrix, first by time relative to addition, then by type of addition 
BeforeFlux<-subset(Flux, Flux[, "Set"]=="Before")
DuringFlux<-subset(Flux, Flux$Set=="During")
AfterFlux<-subset(Flux, Flux$Set=="After")


##Now testing different model structures to best describe the data
## we have variables time, addition/substrate, and site
#model1(b)(d)(a) is separating by time and then doing F~s*a
##This info goes to make manuscript Table 3

modelcrazy <- lm(Slope_gsoil~Site*Addition*Set, data=Flux)
anova(modelcrazy)
summary(modelcrazy)

model1b <- lm(Slope_gsoil~Site*Addition, data=BeforeFlux)
anova(model1b)
summary(model1b)
model1d <- lm(Slope_gsoil~Site*Addition, data=DuringFlux)
anova(model1d)
summary(model1d)

model1a <- lm(Slope_gsoil~Site*Addition, data=AfterFlux)
summary(model1a)
anova(model1a)

