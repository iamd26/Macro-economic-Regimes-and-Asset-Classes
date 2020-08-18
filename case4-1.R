################################
##############Case 4############
################################
library(dplyr)
#load the data
load("data_4.RData")
head(dt)

################################
#Direction 1
#five samples
whole <- dt
i1g1 <- dt %>% filter(INF1GRW1 == 1)
i2g1 <- dt %>% filter(INF2GRW1 == 1)
i1g2 <- dt %>% filter(INF1GRW2 == 1)
i2g2 <- dt %>% filter(INF2GRW2 == 1)

#1a - whole
mean_whole <- apply(whole[,6:10],2,mean)
sd_whole <- apply(whole[,6:10],2,sd)
sr_whole <- vector(length=4)
for (i in 1:4){sr_whole[i] <- (mean_whole[i] - mean_whole[5])/sd_whole[i]}
whole_a <- rbind(mean_whole[1:4],sd_whole[1:4],sr_whole)
row.names(whole_a) <- c("mean","sd","sr")
whole_a

#1a - i1g1
mean_i1g1 <- apply(i1g1[,6:10],2,mean)
sd_i1g1 <- apply(i1g1[,6:10],2,sd)
sr_i1g1 <- vector(length=4)
for (i in 1:4){sr_i1g1[i] <- (mean_i1g1[i] - mean_i1g1[5])/sd_i1g1[i]}
i1g1_a <- rbind(mean_i1g1[1:4],sd_i1g1[1:4],sr_i1g1)
row.names(i1g1_a) <- c("mean","sd","sr")
i1g1_a

#1a - i2g1
mean_i2g1 <- apply(i2g1[,6:10],2,mean)
sd_i2g1 <- apply(i2g1[,6:10],2,sd)
sr_i2g1 <- vector(length=4)
for (i in 1:4){sr_i2g1[i] <- (mean_i2g1[i] - mean_i2g1[5])/sd_i2g1[i]}
i2g1_a <- rbind(mean_i2g1[1:4],sd_i2g1[1:4],sr_i2g1)
row.names(i2g1_a) <- c("mean","sd","sr")
i2g1_a

#1a - i1g2
mean_i1g2 <- apply(i1g2[,6:10],2,mean)
sd_i1g2 <- apply(i1g2[,6:10],2,sd)
sr_i1g2 <- vector(length=4)
for (i in 1:4){sr_i1g2[i] <- (mean_i1g2[i] - mean_i1g2[5])/sd_i1g2[i]}
i1g2_a <- rbind(mean_i1g2[1:4],sd_i1g2[1:4],sr_i1g2)
row.names(i1g2_a) <- c("mean","sd","sr")
i1g2_a

#1a - i2g2
mean_i2g2 <- apply(i2g2[,6:10],2,mean)
sd_i2g2 <- apply(i2g2[,6:10],2,sd)
sr_i2g2 <- vector(length=4)
for (i in 1:4){sr_i2g2[i] <- (mean_i2g2[i] - mean_i2g2[5])/sd_i2g2[i]}
i2g2_a <- rbind(mean_i2g2[1:4],sd_i2g2[1:4],sr_i2g2)
row.names(i2g2_a) <- c("mean","sd","sr")
i2g2_a

#1b - vcv
head(whole)
vcv_whole <- cov(whole[,6:9])
vcv_i1g1 <- cov(i1g1[,6:9])
vcv_i2g1 <- cov(i2g1[,6:9])
vcv_i1g2 <- cov(i1g2[,6:9])
vcv_i2g2 <- cov(i2g2[,6:9])

#1c - weights
one <- c(1,1,1,1)
excess_whole <- whole_a[1,] - mean_whole[5]
excess_i1g1 <- i1g1_a[1,] - mean_i1g1[5]
excess_i1g2 <- i1g2_a[1,] - mean_i1g2[5]
excess_i2g1 <- i2g1_a[1,] - mean_i2g1[5]
excess_i2g2 <- i2g2_a[1,] - mean_i2g2[5]
w_whole <- solve(vcv_whole)%*%(excess_whole)%*%solve(one%*%solve(vcv_whole)%*%excess_whole)
w_i1g1 <- solve(vcv_i1g1)%*%(excess_i1g1)%*%solve(one%*%solve(vcv_i1g1)%*%excess_i1g1)
w_i1g2 <- solve(vcv_i1g2)%*%(excess_i1g2)%*%solve(one%*%solve(vcv_i1g2)%*%excess_i1g2)
w_i2g1 <- solve(vcv_i2g1)%*%(excess_i2g1)%*%solve(one%*%solve(vcv_i2g1)%*%excess_i2g1)
w_i2g2 <- solve(vcv_i2g2)%*%(excess_i2g2)%*%solve(one%*%solve(vcv_i2g2)%*%excess_i2g2)

#1c - expected return
er_whole <- t(as.matrix(whole_a[1,]))%*%w_whole
er_i1g1 <- t(as.matrix(i1g1_a[1,]))%*%w_i1g1
er_i1g2 <- t(as.matrix(i1g2_a[1,]))%*%w_i1g2
er_i2g1 <- t(as.matrix(i2g1_a[1,]))%*%w_i2g1
er_i2g2 <- t(as.matrix(i2g2_a[1,]))%*%w_i2g2

#1c - sd
std_whole <- (t(w_whole)%*%vcv_whole%*%w_whole)^0.5
std_i1g1 <- (t(w_i1g1)%*%vcv_i1g1%*%w_i1g1)^0.5
std_i1g2 <- (t(w_i1g2)%*%vcv_i1g2%*%w_i1g2)^0.5
std_i2g1 <- (t(w_i2g1)%*%vcv_i2g1%*%w_i2g1)^0.5
std_i2g2 <- (t(w_i2g2)%*%vcv_i2g2%*%w_i2g2)^0.5

#1c - sr
shr_whole <- (er_whole - mean_whole[5])/std_whole
shr_i1g1 <- (er_i1g1 - mean_i1g1[5])/std_i1g1
shr_i1g2 <- (er_i1g2 - mean_i1g2[5])/std_i1g2
shr_i2g1 <- (er_i2g1 - mean_i2g1[5])/std_i2g1
shr_i2g2 <- (er_i2g2 - mean_i2g2[5])/std_i2g2

#1d - gmv weights
wgmv_whole <- solve(vcv_whole)%*%one%*%solve(t(one)%*%solve(vcv_whole)%*%one)
wgmv_i1g1 <- solve(vcv_i1g1)%*%one%*%solve(t(one)%*%solve(vcv_i1g1)%*%one)
wgmv_i1g2 <- solve(vcv_i1g2)%*%one%*%solve(t(one)%*%solve(vcv_i1g2)%*%one)
wgmv_i2g1 <- solve(vcv_i2g1)%*%one%*%solve(t(one)%*%solve(vcv_i2g1)%*%one)
wgmv_i2g2 <- solve(vcv_i2g2)%*%one%*%solve(t(one)%*%solve(vcv_i2g2)%*%one)

#1d - gmv expected return
ergmv_whole <- t(as.matrix(whole_a[1,]))%*%wgmv_whole
ergmv_i1g1 <- t(as.matrix(i1g1_a[1,]))%*%wgmv_i1g1
ergmv_i1g2 <- t(as.matrix(i1g2_a[1,]))%*%wgmv_i1g2
ergmv_i2g1 <- t(as.matrix(i2g1_a[1,]))%*%wgmv_i2g1
ergmv_i2g2 <- t(as.matrix(i2g2_a[1,]))%*%wgmv_i2g2

#1d - gmv sd
sdgmv_whole <- (t(wgmv_whole)%*%vcv_whole%*%wgmv_whole)^0.5
sdgmv_i1g1 <- (t(wgmv_i1g1)%*%vcv_i1g1%*%wgmv_i1g1)^0.5
sdgmv_i1g2 <- (t(wgmv_i1g2)%*%vcv_i1g2%*%wgmv_i1g2)^0.5
sdgmv_i2g1 <- (t(wgmv_i2g1)%*%vcv_i2g1%*%wgmv_i2g1)^0.5
sdgmv_i2g2 <- (t(wgmv_i2g2)%*%vcv_i2g2%*%wgmv_i2g2)^0.5

#1d - gmv sr
srgmv_whole <- (ergmv_whole - mean_whole[5])/sdgmv_whole
srgmv_i1g1 <- (ergmv_i1g1 - mean_i1g1[5])/sdgmv_i1g1
srgmv_i1g2 <- (ergmv_i1g2 - mean_i1g2[5])/sdgmv_i1g2
srgmv_i2g1 <- (ergmv_i2g1 - mean_i2g1[5])/sdgmv_i2g1
srgmv_i2g2 <- (ergmv_i2g2 - mean_i2g2[5])/sdgmv_i2g2

#1e - wa
A <- c(1/1.3,1/2.8,1/6.5,1/10.5,1/16.9)
wa_whole <- solve(vcv_whole)%*%excess_whole%*%A
Tbills_whole <- 1 - apply(wa_whole,2,sum)
wa_whole <- rbind(wa_whole,Tbills_whole)
colnames(wa_whole) <- c("A=1.3","A=2.8","A=6.5","A=10.5","A=16.9")

wa_i1g1 <- solve(vcv_i1g1)%*%excess_i1g1%*%A
Tbills_i1g1 <- 1 - apply(wa_i1g1,2,sum)
wa_i1g1 <- rbind(wa_i1g1,Tbills_i1g1)
colnames(wa_i1g1) <- c("A=1.3","A=2.8","A=6.5","A=10.5","A=16.9")

wa_i1g2 <- solve(vcv_i1g2)%*%excess_i1g2%*%A
Tbills_i1g2 <- 1 - apply(wa_i1g2,2,sum)
wa_i1g2 <- rbind(wa_i1g2,Tbills_i1g2)
colnames(wa_i1g2) <- c("A=1.3","A=2.8","A=6.5","A=10.5","A=16.9")

wa_i2g1 <- solve(vcv_i2g1)%*%excess_i2g1%*%A
Tbills_i2g1 <- 1 - apply(wa_i2g1,2,sum)
wa_i2g1 <- rbind(wa_i2g1,Tbills_i2g1)
colnames(wa_i2g1) <- c("A=1.3","A=2.8","A=6.5","A=10.5","A=16.9")

wa_i2g2 <- solve(vcv_i2g2)%*%excess_i2g2%*%A
Tbills_i2g2 <- 1 - apply(wa_i2g2,2,sum)
wa_i2g2 <- rbind(wa_i2g2,Tbills_i2g2)
colnames(wa_i2g2) <- c("A=1.3","A=2.8","A=6.5","A=10.5","A=16.9")

wa_whole
wa_i1g1
wa_i1g2
wa_i2g1
wa_i2g2

#1e - expected return
er_allocation_whole <- mean_whole%*%wa_whole
er_allocation_i1g1 <- mean_i1g1%*%wa_i1g1
er_allocation_i1g2 <- mean_i1g2%*%wa_i1g2
er_allocation_i2g1 <- mean_i2g1%*%wa_i2g1
er_allocation_i2g2 <- mean_i2g2%*%wa_i2g2

#1e - sd
risk_whole <- apply(wa_whole[1:4,],2,sum)
sd_allocation_whole <- risk_whole*as.numeric(std_whole)
risk_i1g1 <- apply(wa_i1g1[1:4,],2,sum)
sd_allocation_i1g1 <- risk_i1g1*as.numeric(std_i1g1)
risk_i1g2 <- apply(wa_i1g2[1:4,],2,sum)
sd_allocation_i1g2 <- risk_i1g2*as.numeric(std_i1g2)
risk_i2g1 <- apply(wa_i2g1[1:4,],2,sum)
sd_allocation_i2g1 <- risk_i2g1*as.numeric(std_i2g1)
risk_i2g2 <- apply(wa_i2g2[1:4,],2,sum)
sd_allocation_i2g2 <- risk_i2g2*as.numeric(std_i2g2)

#1e - sr
sr_allocation_whole <- vector(length = 5)
for (i in 1:5){sr_allocation_whole[i] <- (er_allocation_whole[1,i]-mean_whole[5])/sd_allocation_whole[i]}
sr_allocation_i1g1 <- vector(length = 5)
for (i in 1:5){sr_allocation_i1g1[i] <- (er_allocation_i1g1[1,i]-mean_i1g1[5])/sd_allocation_i1g1[i]}
sr_allocation_i1g2 <- vector(length = 5)
for (i in 1:5){sr_allocation_i1g2[i] <- (er_allocation_i1g2[1,i]-mean_i1g2[5])/sd_allocation_i1g2[i]}
sr_allocation_i2g1 <- vector(length = 5)
for (i in 1:5){sr_allocation_i2g1[i] <- (er_allocation_i2g1[1,i]-mean_i2g1[5])/sd_allocation_i2g1[i]}
sr_allocation_i2g2 <- vector(length = 5)
for (i in 1:5){sr_allocation_i2g2[i] <- (er_allocation_i2g2[1,i]-mean_i2g2[5])/sd_allocation_i2g2[i]}

################################
#Direction 2
#summary table - whole
sr_whole <- (mean_whole - mean_whole[5])/sd_whole
sum_whole <- rbind(mean_whole,sd_whole,sr_whole,wa_whole[,3])
rownames(sum_whole) <- c("mean","sd","sr","A=6.5")
sum_whole
#summary table - i1g1
sr_i1g1 <- (mean_i1g1 - mean_i1g1[5])/sd_i1g1
sum_i1g1 <- rbind(mean_i1g1,sd_i1g1,sr_i1g1,wa_i1g1[,3])
rownames(sum_i1g1) <- c("mean","sd","sr","A=6.5")
sum_i1g1
#summary table - i1g2
sr_i1g2 <- (mean_i1g2 - mean_i1g2[5])/sd_i1g2
sum_i1g2 <- rbind(mean_i1g2,sd_i1g2,sr_i1g2,wa_i1g2[,3])
rownames(sum_i1g2) <- c("mean","sd","sr","A=6.5")
sum_i1g2
#summary table - i2g1
sr_i2g1 <- (mean_i2g1 - mean_i2g1[5])/sd_i2g1
sum_i2g1 <- rbind(mean_i2g1,sd_i2g1,sr_i2g1,wa_i2g1[,3])
rownames(sum_i2g1) <- c("mean","sd","sr","A=6.5")
sum_i2g1
#summary table - i2g2
sr_i2g2 <- (mean_i2g2 - mean_i2g2[5])/sd_i2g2
sum_i2g2 <- rbind(mean_i2g2,sd_i2g2,sr_i2g2,wa_i2g2[,3])
rownames(sum_i2g2) <- c("mean","sd","sr","A=6.5")
sum_i2g2

################################
#Direction 3
#Static - EW Portfolio
static_ew <- 0.25*sum_i1g1[4,] + 0.25*sum_i1g2[4,] + 0.25*sum_i2g1[4,] + 0.25*sum_i2g2[4,]
#Tilt - i1g1
tilt_i1g1 <- 0.5*sum_i1g1[4,] + (1/6)*sum_i1g2[4,] + (1/6)*sum_i2g1[4,] + (1/6)*sum_i2g2[4,]
#Tilt - i1g2
tilt_i1g2 <- (1/6)*sum_i1g1[4,] + (0.5)*sum_i1g2[4,] + (1/6)*sum_i2g1[4,] + (1/6)*sum_i2g2[4,]
#Tilt - i2g1
tilt_i2g1 <- (1/6)*sum_i1g1[4,] + (1/6)*sum_i1g2[4,] + (0.5)*sum_i2g1[4,] + (1/6)*sum_i2g2[4,]
#Tilt - i2g2
tilt_i2g2 <- (1/6)*sum_i1g1[4,] + (1/6)*sum_i1g2[4,] + (1/6)*sum_i2g1[4,] + (0.5)*sum_i2g2[4,]
#
static_ew
tilt_i1g1
tilt_i1g2
tilt_i2g1
tilt_i2g2

################################
#Direction 4
whole_6.5 <- sum_whole[4,]
i1g1_6.5 <- sum_i1g1[4,]
i1g2_6.5 <- sum_i1g2[4,]
i2g1_6.5 <- sum_i2g1[4,]
i2g2_6.5 <- sum_i2g2[4,]

#all the weight for A = 6.5
whole_6.5
i1g1_6.5
i1g2_6.5
i2g1_6.5
i2g2_6.5
static_ew
tilt_i1g1
tilt_i1g2
tilt_i2g1
tilt_i2g2
weight <- rbind(whole_6.5,i1g1_6.5,i1g2_6.5,i2g1_6.5,i2g2_6.5,static_ew,tilt_i1g1,tilt_i1g2,tilt_i2g1,tilt_i2g2)
weight

#expected return
unconditional_return <- mean_whole%*%t(weight)
i1g1_return <- mean_i1g1%*%t(weight)
i1g2_return <- mean_i1g2%*%t(weight)
i2g1_return <- mean_i2g1%*%t(weight)
i2g2_return <- mean_i2g2%*%t(weight)
################
unconditional_return
i1g1_return
i1g2_return
i2g1_return
i2g2_return



#sd
unconditional_sd <- matrix(nrow=1,ncol=10)
for (i in 1:10) {unconditional_sd[1,i] <- (t(weight[i,1:4])%*%vcv_whole%*%weight[i,1:4])^0.5}
colnames(unconditional_sd) <- c("whole_6.5","i1g1_6.5","i1g2_6.5","i2g1_6.5","i2g2_6.5",
                               "static_ew","tilt_i1g1","tilt_i1g2","tilt_i2g1","tilt_i2g2")
i1g1_sd <- matrix(nrow=1,ncol=10)
for (i in 1:10) {i1g1_sd[1,i] <- (t(weight[i,1:4])%*%vcv_i1g1%*%weight[i,1:4])^0.5}
colnames(i1g1_sd) <- c("whole_6.5","i1g1_6.5","i1g2_6.5","i2g1_6.5","i2g2_6.5",
                       "static_ew","tilt_i1g1","tilt_i1g2","tilt_i2g1","tilt_i2g2")
i1g2_sd <- matrix(nrow=1,ncol=10)
for (i in 1:10) {i1g2_sd[1,i] <- (t(weight[i,1:4])%*%vcv_i1g2%*%weight[i,1:4])^0.5}
colnames(i1g2_sd) <- c("whole_6.5","i1g1_6.5","i1g2_6.5","i2g1_6.5","i2g2_6.5",
                       "static_ew","tilt_i1g1","tilt_i1g2","tilt_i2g1","tilt_i2g2")
i2g1_sd <- matrix(nrow=1,ncol=10)
for (i in 1:10) {i2g1_sd[1,i] <- (t(weight[i,1:4])%*%vcv_i2g1%*%weight[i,1:4])^0.5}
colnames(i2g1_sd) <- c("whole_6.5","i1g1_6.5","i1g2_6.5","i2g1_6.5","i2g2_6.5",
                       "static_ew","tilt_i1g1","tilt_i1g2","tilt_i2g1","tilt_i2g2")
i2g2_sd <- matrix(nrow=1,ncol=10)
for (i in 1:10) {i2g2_sd[1,i] <- (t(weight[i,1:4])%*%vcv_i2g2%*%weight[i,1:4])^0.5}
colnames(i2g2_sd) <- c("whole_6.5","i1g1_6.5","i1g2_6.5","i2g1_6.5","i2g2_6.5",
                       "static_ew","tilt_i1g1","tilt_i1g2","tilt_i2g1","tilt_i2g2")
################
unconditional_sd
i1g1_sd
i1g2_sd
i2g1_sd
i2g2_sd

#rf
unconditional_rf <- as.numeric(mean_whole[5])
i1g1_rf <- as.numeric(mean_i1g1[5])
i1g2_rf <- as.numeric(mean_i1g2[5])
i2g1_rf <- as.numeric(mean_i2g1[5])
i2g2_rf <- as.numeric(mean_i2g2[5])

#sr
unconditional_sr <- (unconditional_return - unconditional_rf)/unconditional_sd
i1g1_sr <- (i1g1_return - i1g1_rf)/i1g1_sd
i1g2_sr <- (i1g2_return - i1g2_rf)/i1g2_sd
i2g1_sr <- (i2g1_return - i2g1_rf)/i2g1_sd
i2g2_sr <- (i2g2_return - i2g2_rf)/i2g2_sd
###########
unconditional_sr
i1g1_sr
i1g2_sr
i2g1_sr
i2g2_sr
Sharpe_Ratios <- rbind(unconditional_sr,i1g1_sr,i1g2_sr,i2g1_sr,i2g2_sr)
Sharpe_Ratios <- t(Sharpe_Ratios)
colnames(Sharpe_Ratios) <- c("Unconditional","i1g1","i1g2","i2g1","i2g2")

############
Sharpe_Ratios

############################Tables###################################
sum_whole
sum_i1g1
sum_i1g2
sum_i2g1
sum_i2g2 #direction 2

static_ew
tilt_i1g1
tilt_i1g2
tilt_i2g1
tilt_i2g2 #direction 3

Sharpe_Ratios #direction 4