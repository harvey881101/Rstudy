#정준상관분석 x1,x2,x3,...,xp 와 y1,y2... yq 까지의 관련성이 있어야만 정준상관분석 의미가 있음. 
#공분산행렬이 0이면 정준상관계수 r1, r2, r3... r_s는 통계적인 의미가 없음
install.packages("CCA")
library(CCA)

chem <- read.table("chemistry.txt", head=T, skip=3, sep="\t")
chem
Y<-chem[,1:2]
X<-chem[,3:5]

library(CCA)
print(matcor(X,Y), digits=3)

cc1<-cc(X,Y)
cc1[1] #정준상관계수 

cc1[3:4] #정준변수계수 
#u1= -0.136X1+ -0.121X2 -0.159X3 ##계수가 -이므로 아래에서도 -로 나올 것
#v1= 0.~~

cc2<-comput(X,Y,cc1)
cc2[3:6] #변수와 정준변수 간의 상관계수 

cov.X <- cov(X) # X의 공분산행렬 x1,x2,x3
sd.X <- sqrt(diag(cov.X)) #x1,x2,x3의 표준편차
S.X <- diag(sd.X) #X1, X2, X3의 표준편차 대각행렬 
S.X %*% cc1$xcoef #X의 표준화정준변수계수 

U1<-cc1$scores$xscores[,1]
V1<-cc1$scores$yscores[,1]
windows()
plot(U1,V1 ,pch=18,main="First Canonical Plot" )
plt.cc(cc1,type="v", var.label = T)


################################################################################
## 판별분석 ### 
ldahist(data, g, breaks, width, xlim=range(breaks),
        ymax=0, type=c("histogram", "density", "both"),
        )

lda prior = proportions # 사전확률 c(0.5, 0.3, 0.3) 알고있으면 쓸 것, 없으면 표본비율이 


#등분산검정(biotools 패키지) 
boxM(data, grouping) #data : 공분산행렬 데이터 // grouping : 그룹변수 
install.packages("bitools")

#판별그래프 klaR 패키지 
partimat(x, grouping, data, method="lda" or "qda", plot.matrix=FALSE, imageplot=TRUE)
# x= 독립변수, grouping = 집단변수 
install.packages("biotools")
library(klaR)
#판별분석 변수 선택 stepclass 결과가 다르게 나올수도있음. 

#method : 판별분석함수

skull <- read.table("skull.txt", header=TRUE, skip=6, sep="\t")
attach(skull)
n <- dim(skull)[1]
X <- skull[,1:4]
boxM(X,Year) # H0 : 귀무가설(집단의 분산이 같다.) -> 선형 판별분석 실시 

ldahist(X1, g=Year, type="histogram", col="black", border="white")
#x1는 히스토그램만그룹별로 히스토그램을 그림.
ldahist(X2, g=Year, type="density")
# x2는 밀도그림만
ldahist(X3, g=Year, type="both", col="black", border="white")
ldahist(X4, g=Year, type="both", col="black", border="white")
# 히스토그램을 그렸을 때 명확한 차이가 있어야 좋음 

##패키지 설치
library(MASS)
library(klaR)
library(biotools)
##

ld <- lda(Year~X1+X2+X3+X4, data=skull)
# or ld <- lda(X, Year, data=skull)
ld[1] # 각각 갯수가 30개씩있음. 
ld[3] # 
ld[4] # 0.1187*X1 + -0.084*X2 + -0.10842*X3 + 0.026*X4 

pred.ld <- predict(ld, X)
attributes(pred.ld)
pred.LG <- pred.ld$class
pred.LG <- as.numeric(pred.LG)
pred.LG[pred.LG==1]<-(-4000)
pred.LG[pred.LG==2]<-150

Lc.Table <- table(Year, pred.LG)
ld.correct.rate <- sum(diag(Lc.Table))/n

ld.error.rate <- 1-ld.correct.rate
Lc.Table;ld.correct.rate;ld.error.rate

windows()
partimat(X, factor(Year, labels=c("B","A")), data=skull, method="lda")
partimat(X, factor(Year, labels=c("B","A")), data=skull, method="lda",
         plot.matrix=TRUE, imageplot=FALSE)


###만약 2차 판별함수를 사용한다면,
qd <- qda(X, Year)
pred.qd <- predict(qd)
pred.qd

#####
ldc <- lda(Year~X1+X2+X3+X4, data=skull, CV=TRUE, prior=c(0.5,0.5))
names(ldc)
result <- data.frame(Year, ldc$class, ldc$posterior)
head(result, n=4) # 사후 확률이 큰쪽으로 분류를 한다. 
c.table <- table(Year, ldc$class)
prop.table(c.table, margin=1)

credit <- read.table("credit.txt", head=T, skip=13, sep=",")
credit$G <- factor(credit$G, labels=c("양호", "모호", "불량"))
head(credit, n=3)

attach(credit)
n<-dim(credit)[1]
X<-credit[,2:11] 

#동질성검정은 생략 



############### 군집분석 #####################

# 계층적 군집방법 : n개의 군집으로부터 점차 군집의 개수를 줄여 나가는 방법 (n개부터 1개까지)
# 최단연결법 / 최장연결법 / 평균연결법 / Ward 방법 
# C1, C2간의 거리를 C1와 C2거리중에서 가장 낮은 것으로 정의.. 

# 계층적 군집 방법으로 적절한 군집의 개수를 정하고 
# K-means 방법으로 같이 쓰일 수 있음. 

# 군집 개수의 결정 
# 1. 일단해보고... 사후에 결정
# 2. Ward의 방법에서는 ESS 증분값이 크게 변화하는 곳 
# 3. 통계적 방법은 아직 없음 (다변량 분산분석을 통해 차이가 가장 크게..)
# 4. 계층적 군집 방법 : Mojena(1977) 
# 5. Cross valuation 
# 6. 모형기반 군집 방법 : 확률 분포적 모형에 대한 정보를 이용 

dist(x, method="euclidean", diag=FALSE, upper=FALSE, p=2)
# diag (자기 자신의 거리는 0) 
# upper는 하한쪽값만 / (어차피 똑같기 때문에..) 
# p: 민코우스키 거리의 지수 (p>0)

hclust(d, method="complete")
# d : 거리 행렬 개체 method = "complete" 

plot(x, labels=NULL, hang=0.1 )
cutree(tree, k=NULL ) 

crime <- read.table("crime.txt",skip=4, header=T, sep="\t")
head(crime,n=3)
attach(crime)
X <- crime[,2:3]


#kmeans(x, centers)
# x : 군집분석을 위한 원자료 centers : 군집의 수 
n<-dim(crime)[1] #표본의 수 
cl.k <- n-1 #군집 최대수
wss <- rep(NA, cl.k) #군집의 수에 따른 군집내 제곱합 
wss[1] <- (n-1)*sum(sapply(X, FUN=var))


for(i in 2:cl.k){
  wss[i] <- sum(kmeans(X, centers=i)$withinss)
}

wss
plot(1:cl.k, wss, type="b", xlab="군집수",ylab="군집내 제곱합")
km <- kmeans(X, center=5)
attributes(km)
clus <- data.frame(City, X, G=km$cluster)
plot(X, pch=km$cluster, main="K-means")
text(X,labels=City, pos=0)
points(km$centers, pch=15)
LB <- paste("(", round(km$centers[,1],1),",",round(km$centers[,2],1),")",sep="")
text(km$centers, labels=LB, pos=3)
km$centers

library(mclust)
# Mclust(data, G=최대 군집수, modelNames=NUll) #BIC가 가장 큰 값으로 모형을 자동으로 함
mc <- Mclust(X, 2:5) #군집 2개에서 ~5개 적절히..
mc.clus <- data.frame(crime, G=mc$classification)
mc.clus

plot(mc, what="BIC")
plot(mc, what="classification")
plot(mc, what="uncertainty")
plot(mc, what="density")


## 다차원 척도법 : 객체간 얼마나 가까이에 있는지 
 # 데이터와 유사한 지도를 산
 # 비유사성 행렬 / 거리 행렬 
 # 근접성 행렬을 구함 
 # 유사성 판단을 설명하는 이론의 개발 : 
 # 다차원 척도법은 공간모형 
 

X<- read.table("MDS.txt", head=TRUE, sep="\t")
X
n<-nrow(X)
D<-dist(X)
round(D,1)
print(cmdscale(D, k=9, eig=TRUE),2)
#1부터 5까지의 값이 양수값. 
#따라서 points의 처음 5개 열만 유클리드 거리행렬을 나타냄 

max(abs(dist(X)-dist(cmdscale(D, k=5)))) # 거의 같다. 
max(abs(prcomp(X)$x)-abs(cmdscale(D,k=5)))


D.m <- dist(X, method="manhattan")
MDS.m <- cmdscale(D.m, k=n-1, eig=TRUE)
#유클리디안이 아니면 eig값이 음수가 나올 수 있음. 
#Pm1 대략3차원이 적당함. 
cumsum(abs(MDS.m$eig))/sum(abs(MDS.m$eig))
cumsum(abs(MDS.m$eig^2))/sum(abs(MDS.m$eig^2))
