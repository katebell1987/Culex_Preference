##### Set working directory ####
setwd("~/Dropbox/UMD_Research/Preference/HostCompExp/")

### Load Libraries ####
library(rjags)
library(coda)
#install.packages("kableExtra")
library(kableExtra)

#### Read in the data ####
data<-read.csv("RawData/KateVsAnna_RDforR.csv",header=T)


########### Get data ready to run the model for Evanston ####################
countsE<-data$NumberHuman[data$Population=="Evanston" ]
trialsE<-data$NumberChick[data$Population=="Evanston"] + data$NumberHuman[data$Population=="Evanston"]
chickE<-data$Chick[data$Population=="Evanston"]
humanE<-as.numeric(as.factor(data$HumanHost[data$Population=="Evanston"]))
datamodE<-as.data.frame(cbind(countsE,trialsE,chickE,humanE))
colnames(datamodE)<-c("counts","trials","chick","human")
datamodE$chick<-as.numeric(as.factor(datamodE$chick))
datamodE$counts<-as.numeric(datamodE$counts)
datamodE$trials<-as.numeric(datamodE$trials)

basicDataE<-list(Counts=datamodE$counts[],Responders=datamodE$trials[],chick=datamodE$chick[],human=datamodE$human[],N=length(countsE))

######################## Specify the model for Evanston #######################
model<-"model{
for(i in 1:N){
Counts[i] ~ dbin(p[i], Responders[i])
p[i] <- 1 / (1 + exp(-1 * alpha[i]))
alpha[i] <- int + beta1[human[i]]}

#prior on intercept 
int ~ dnorm(0,0.001)

beta1[1]<-0 ## Anna
beta1[2]~ dnorm(0,0.001) ## Kate

}"

basicModelE<-jags.model(textConnection(model),data=basicDataE,n.chains=4)
samplesE<-jags.samples(model=basicModelE,variable.names=c("beta1","int","p"),n.iter=100000,n.burnin=500000)

#save.image("test2.RData")

#################### Check for stable sampling distribution ###################

beta1EES<-effectiveSize(samplesE$beta1[2,,])
intEES<-effectiveSize(samplesE$int[,,])
#pEES<-effectiveSize(samplesE$p[1,,])
gelman.diag(mcmc.list(as.mcmc(samplesE$beta1[2,,1]),as.mcmc(samplesE$beta1[2,,2]),as.mcmc(samplesE$beta1[2,,3]),as.mcmc(samplesE$beta1[2,,3])))

plot(samplesE$beta1[2,,1],col="grey")
points(samplesE$beta1[2,,2],col="lightblue")
points(samplesE$beta1[2,,3],col="lightpink")
points(samplesE$beta1[2,,4],col="lightgreen")

plot(samplesE$int[1,,1],col="grey")
points(samplesE$int[1,,2],col="lightblue")
points(samplesE$int[1,,3],col="lightpink")
points(samplesE$int[1,,4],col="lightgreen")

######### Plots of P and beta's ########
pdf("Evanston_beta.pdf")
plot(density(samplesE$beta1[2,,]),col="#2166ac",xlim=c(-5,5),ylim=c(0,1.5),xlab="Posterior Probability",main="Posterior Probability Estimates\n of Effect of Human Host for Evanston") 
lines(density(samplesE$int[1,,]),col="#b2182b")
legend(x=2.75,y=1.5,legend=c("Human Host 1",
                             "Human Host 2")
       ,cex=0.9,col=c("#b2182b","#2166ac"),pch=19)
dev.off()
combSamplesEbeta<-c(samplesE$beta1[2,,1],samplesE$beta1[2,,2],samplesE$beta1[2,,3],samplesE$beta1[2,,4])
plogis(summary(combSamplesEbeta))
summary(combSamplesEbeta)
combSamplesEint<-c(samplesE$int[1,,1],samplesE$int[1,,2],samplesE$int[1,,3],samplesE$int[1,,4])
summary(combSamplesEint)
plogis(summary(combSamplesEint))

pdf("Evanston_p.pdf")
plot(density(samplesE$p[,,]),col="black",xlim=c(0,1),ylim=c(0,15),xlab="Posterior Probability of Landing on a Human Host",main="Posterior Probability Estimates\n of Preference for Evanston") 
dev.off()

######################## Model for Cal1 ###################
countsC<-data$NumberHuman[data$Population=="Cal1"]
trialsC<-data$NumberChick[data$Population=="Cal1"] + data$NumberHuman[data$Population=="Cal1"]
chickC<-data$Chick[data$Population=="Cal1"]
humanC<-as.numeric(as.factor(data$HumanHost[data$Population=="Cal1"]))
datamodC<-as.data.frame(cbind(countsC,trialsC,chickC,humanC))
colnames(datamodC)<-c("counts","trials","chick","human")
datamodC$chick<-as.numeric(as.factor(datamodC$chick))
datamodC$counts<-as.numeric(datamodC$counts)
datamodC$trials<-as.numeric(datamodC$trials)
#datamod$dayPair<-data$DayPair[data$Population=="Evanston"]

basicDataC<-list(Counts=datamodC$counts[],Responders=datamodC$trials[],human=datamodC$human[],N=length(countsC))
############### specify the model #################
model<-"model{
for(i in 1:N){
Counts[i] ~ dbin(p[i], Responders[i])
p[i] <- 1 / (1 + exp(-1 * alpha[i]))
alpha[i] <- int + beta1[human[i]]}

#prior on intercept 
int ~ dnorm(0,0.001)

beta1[1]<-0 ## Anna
beta1[2]~ dnorm(0,0.001) ## Kate

}"

basicModelC<-jags.model(textConnection(model),data=basicDataC,n.chains=4)
samplesC<-jags.samples(model=basicModelC,variable.names=c("beta1","int","p"),n.iter=100000,n.burnin=500000)


#################### Check for stable sampling distribution ###################

beta1CES<-effectiveSize(samplesC$beta1[2,,])
intCES<-effectiveSize(samplesC$int[,,])
#pEES<-effectiveSize(samplesE$p[1,,])
gelman.diag(mcmc.list(as.mcmc(samplesC$beta1[2,,1]),as.mcmc(samplesC$beta1[2,,2]),as.mcmc(samplesC$beta1[2,,3]),as.mcmc(samplesC$beta1[2,,3])))

plot(samplesC$beta1[2,,1],col="grey")
points(samplesC$beta1[2,,2],col="lightblue")
points(samplesC$beta1[2,,3],col="lightpink")
points(samplesC$beta1[2,,4],col="lightgreen")

plot(samplesC$int[1,,1],col="grey")
points(samplesC$int[1,,2],col="lightblue")
points(samplesC$int[1,,3],col="lightpink")
points(samplesC$int[1,,4],col="lightgreen")

######### Plots of P and beta's ########
pdf("Cal1_beta.pdf")
plot(density(samplesC$beta1[2,,]),col="#2166ac",xlim=c(-6,7),ylim=c(0,1.5),xlab="Posterior Probability",main="Posterior Probability Estimates\n of Effect of Human Host for Cal1") 
lines(density(samplesC$int[1,,]),col="#b2182b")
legend(x=2.75,y=1.5,legend=c("Human Host 1",
                             "Human Host 2")
       ,cex=0.9,col=c("#b2182b","#2166ac"),pch=19)
dev.off()

pdf("Cal1_p.pdf")
plot(density(samplesC$p[,,]),col="black",xlim=c(0,1),ylim=c(0,15),xlab="Posterior Probability of Landing on a Human Host",main="Posterior Probability Estimates\n of Preference for Cal1") 
dev.off()

combSamplesCbeta<-c(samplesC$beta1[2,,1],samplesC$beta1[2,,2],samplesC$beta1[2,,3],samplesC$beta1[2,,4])
plogis(summary(combSamplesCbeta))
summary(combSamplesCbeta)
combSamplesCint<-c(samplesC$int[1,,1],samplesC$int[1,,2],samplesC$int[1,,3],samplesC$int[1,,4])
summary(combSamplesCint)
plogis(summary(combSamplesCint))

######################## Model for Northfield ###################
countsN<-data$NumberHuman[data$Population=="Northfield"]
trialsN<-data$NumberChick[data$Population=="Northfield"] + data$NumberHuman[data$Population=="Northfield"]
chickN<-data$Chick[data$Population=="Northfield"]
humanN<-as.numeric(as.factor(data$HumanHost[data$Population=="Northfield"]))
datamodN<-as.data.frame(cbind(countsN,trialsN,chickN,humanN))
colnames(datamodN)<-c("counts","trials","chick","human")
datamodN$chick<-as.numeric(as.factor(datamodN$chick))
datamodN$counts<-as.numeric(datamodN$counts)
datamodN$trials<-as.numeric(datamodN$trials)
#datamod$dayPair<-data$DayPair[data$Population=="Evanston"]

basicDataN<-list(Counts=datamodN$counts[],Responders=datamodN$trials[],human=datamodN$human[],N=length(countsN))
############### specify the model #################
model<-"model{
for(i in 1:N){
Counts[i] ~ dbin(p[i], Responders[i])
p[i] <- 1 / (1 + exp(-1 * alpha[i]))
alpha[i] <- int + beta1[human[i]]}

#prior on intercept 
int ~ dnorm(0,0.001)

beta1[1]<-0 ## Anna
beta1[2]~ dnorm(0,0.001) ## Kate

}"

basicModelN<-jags.model(textConnection(model),data=basicDataN,n.chains=4)
samplesN<-jags.samples(model=basicModelN,variable.names=c("beta1","int","p"),n.iter=100000,n.burnin=500000)


#################### Check for stable sampling distribution ###################

beta1NES<-effectiveSize(samplesN$beta1[2,,])
intNES<-effectiveSize(samplesN$int[,,])
#pEES<-effectiveSize(samplesE$p[1,,])
gelman.diag(mcmc.list(as.mcmc(samplesN$beta1[2,,1]),as.mcmc(samplesN$beta1[2,,2]),as.mcmc(samplesN$beta1[2,,3]),as.mcmc(samplesN$beta1[2,,3])))

plot(samplesC$beta1[2,,1],col="grey")
points(samplesC$beta1[2,,2],col="lightblue")
points(samplesC$beta1[2,,3],col="lightpink")
points(samplesC$beta1[2,,4],col="lightgreen")

plot(samplesC$int[1,,1],col="grey")
points(samplesC$int[1,,2],col="lightblue")
points(samplesC$int[1,,3],col="lightpink")
points(samplesC$int[1,,4],col="lightgreen")

######### Plots of P and beta's ########
pdf("North_beta.pdf")
plot(density(samplesN$beta1[2,,]),col="#2166ac",xlim=c(-6,7),ylim=c(0,1.5),xlab="Posterior Probability",main="Posterior Probability Estimates\n of Effect of Human Host for Northfield") 
lines(density(samplesN$int[1,,]),col="#b2182b")
legend(x=2.75,y=1.5,legend=c("Human Host 1",
                             "Human Host 2")
       ,cex=0.9,col=c("#b2182b","#2166ac"),pch=19)
dev.off()

pdf("Northfield_p.pdf")
plot(density(samplesN$p[,,]),col="black",xlim=c(0,1),ylim=c(0,6),xlab="Posterior Probability of Landing on a Human Host",main="Posterior Probability Estimates\n of Preference for Northfield") 
dev.off()

combSamplesNbeta<-c(samplesN$beta1[2,,1],samplesN$beta1[2,,2],samplesN$beta1[2,,3],samplesN$beta1[2,,4])
plogis(summary(combSamplesNbeta))
summary(combSamplesNbeta)
combSamplesNint<-c(samplesN$int[1,,1],samplesN$int[1,,2],samplesN$int[1,,3],samplesN$int[1,,4])
summary(combSamplesNint)
plogis(summary(combSamplesNint))


########################### Combined Plots #####################
par(mfrow=c(3,2))
plot(density(samplesE$beta1[2,,]),col="#2166ac",xlim=c(-5,5),ylim=c(0,1.5),xlab="Posterior Probability",main="Posterior Probability Estimates\n of Effect of Human Host for Evanston") 
lines(density(samplesE$int[1,,]),col="#b2182b")

plot(density(samplesE$p[,,]),col="black",xlim=c(0,1),ylim=c(0,15),xlab="Posterior Probability of Landing on a Human Host",main="Posterior Probability Estimates\n of Preference for Evanston") 

plot(density(samplesC$beta1[2,,]),col="#2166ac",xlim=c(-6,7),ylim=c(0,1.5),xlab="Posterior Probability",main="Posterior Probability Estimates\n of Effect of Human Host for Cal1") 
lines(density(samplesC$int[1,,]),col="#b2182b")

plot(density(samplesC$p[,,]),col="black",xlim=c(0,1),ylim=c(0,15),xlab="Posterior Probability of Landing on a Human Host",main="Posterior Probability Estimates\n of Preference for Cal1") 

plot(density(samplesN$beta1[2,,]),col="#2166ac",xlim=c(-6,7),ylim=c(0,1.5),xlab="Posterior Probability",main="Posterior Probability Estimates\n of Effect of Human Host for Northfield") 
lines(density(samplesN$int[1,,]),col="#b2182b")

plot(density(samplesN$p[,,]),col="black",xlim=c(0,1),ylim=c(0,6),xlab="Posterior Probability of Landing on a Human Host",main="Posterior Probability Estimates\n of Preference for Northfield") 






####################### Model for chick nested within human ##################
## I think to do this we need to index each human/chick combination ###

basicData<-list(Counts=datamod$counts[],Responders=datamod$trials[],chick=datamod$chick[],human=datamod$human[],N=length(counts),chickNum=length(unique(chick)),humanhostNum=length(unique(human)))

#### Specify the model ####
model<-"model{
for(i in 1:N){
Counts[i] ~ dbin(p[i], Responders[i])
p[i] <- 1 / (1 + exp(-1 * alpha[i]))
alpha[i] <- int[human[i]] + beta1[human[i] & chick[i]]}

#prior on intercept
int ~ dnorm(0,0.0001)

for(j in 1:chickNum){
beta1[j] ~ dnorm(beta1mu, beta1tau) #### Hierarchical prior on chick
}

for(k in 1:humanhostNum){
beta2[k] ~ dnorm(beta2mu, beta2tau)
}

#hyperprior on precision 
beta1tau ~ dgamma(0.1,0.001)
beta2tau ~ dgamma(0.1,0.001)

#hyperprior on mean 
beta1mu ~ dnorm(0,0.00001)
beta2mu ~ dnorm(0,0.00001)

}"

basicModel<-jags.model(textConnection(model),data=basicData,n.chains=4)
samples<-jags.samples(model=basicModel,variable.names=c("beta1","beta2","int","beta1mu","beta2mu","beta1tau","beta2tau"),n.iter=500000,n.burnin=100000)

#gelman.diag(samples$beta2)


#### Model for just one human host  ####
data<-read.csv("RawData/KateVsAnna_RDforR.csv",header=T)

counts<-data$NumberHuman[data$Population=="Evanston" & data$HumanHost=="Kate"]
trials<-data$NumberChick[data$Population=="Evanston"& data$HumanHost=="Kate"] + data$NumberHuman[data$Population=="Evanston"& data$HumanHost=="Kate"]
chick<-data$Chick[data$Population=="Evanston"& data$HumanHost=="Kate"]
human<-as.numeric(as.factor(data$HumanHost[data$Population=="Evanston"& data$HumanHost=="Kate"]))
datamod<-as.data.frame(cbind(counts,trials,chick,human))
colnames(datamod)<-c("counts","trials","chick","human")
datamod$chick<-as.numeric(as.factor(datamod$chick))
datamod$counts<-as.numeric(datamod$counts)
datamod$trials<-as.numeric(datamod$trials)
#datamod$dayPair<-data$DayPair[data$Population=="Evanston"]

basicData<-list(Counts=datamod$counts[],Responders=datamod$trials[],chick=datamod$chick[],N=length(counts),chickNum=length(unique(chick)))
#### Specify the model ####
model<-"model{
for(i in 1:N){
Counts[i] ~ dbin(p[i], Responders[i])
p[i] <- 1 / (1 + exp(-1 * alpha[i]))
alpha[i] <- int + beta1[chick[i]]}

#prior on intercept
int ~ dnorm(0,0.001)

for(j in 1:chickNum){
beta1[j] ~ dnorm(beta1mu, beta1tau) #### Hierarchical prior on chick
}

#for(k in 1:humanhostNum){
#beta2[k] ~ dnorm(beta2mu, beta2tau)
#}

#hyperprior on precision 
beta1tau ~ dgamma(0.1,0.001)
#beta2tau ~ dgamma(0.1,0.001)

#hyperprior on mean 
beta1mu ~ dnorm(0,0.1)
#beta2mu ~ dnorm(0,0.001)

}"

basicModel<-jags.model(textConnection(model),data=basicData,n.chains=4)
samples<-jags.samples(model=basicModel,variable.names=c("beta1","int","beta1mu","beta1tau","p"),n.iter=500000,n.burnin=100000)

plot(density(samples$beta1[1,,]))
lines(density(samples$beta1[2,,]))
lines(density(samples$beta1[3,,]))
lines(density(samples$beta1[4,,]))
lines(density(samples$beta1[5,,]))
lines(density(samples$beta1[6,,]))
lines(density(samples$beta1[7,,]))
lines(density(samples$beta1[8,,]))
lines(density(samples$beta1[9,,]))
lines(density(samples$beta1[10,,]))
lines(density(samples$beta1[11,,]))
lines(density(samples$beta1[12,,]))

################################### OLD CODE ##################################

#for(i in 1:length(human)){
#  if(human[i]==1){
#    human[i]<-0}
#    else(human[i]<-1)
#}
#datamod$dayPair<-data$DayPair[data$Population=="Evanston"]


#for(j in 1:chickNum){
#beta1[j] ~ dnorm(beta1mu, beta1tau) #### Hierarchical prior on chick
#}

#for(k in 1:2){
#beta1[k] ~ dnorm(beta1mu, beta1tau)
#}

#hyperprior on precision 
#beta1tau ~ dgamma(0.1,0.001)
#beta2tau ~ dgamma(0.1,0.001)

#hyperprior on mean 
#beta1mu ~ dnorm(0,0.001)
#beta2mu ~ dnorm(0,0.001)
