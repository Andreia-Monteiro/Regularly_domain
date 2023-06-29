rm(list=ls())
library(INLA)
library(spatstat)
library(ggplot2)
library(inlabru)
#library(PStestR)
#library(maptools)
#library(sp)
library(sf)


#Caso beta=2 com 
#5% k vizinhos mais pr√≥ximos
#n=100

seed <- 217
set.seed(seed=seed)

# sample size
n=250

# mesh to simulate spatially structured random effect
domain_boundary <- matrix(c(0,0,1,0,1,1,0,1,0,0), byrow=TRUE, ncol=2) # the study region is set as a square domain: (0,1)x(0,1) \subset R^2
#mesh <- inla.mesh.2d(loc.domain=domain_boundary, max.edge=c(0.025, 0.15))
#mesh <- inla.mesh.2d(loc.domain=domain_boundary, max.edge=c(0.05, 0.15))
#actual mesh
#mesh <- inla.mesh.2d(loc.domain=domain_boundary, max.edge=c(0.1, 0.15))
#new mesh1
mesh <- inla.mesh.2d(loc.domain=domain_boundary, max.edge=c(0.3, 0.15))
ggplot() + gg(mesh) + theme_bw() 


# choose preferential parameter (beta=0 is non-preferential)
beta=-2

## SIMULATION STUDY FOR SEVERAL REPLICAS

total_rep<-100
rep<-0
while(rep < total_rep){
  rep=rep+1


### simulation of the spatial random effect that will be weighted
sigma0 <- sqrt(1.5)    #marginal sd of the field
range0 <- 0.15         #scale
kappa0 <- sqrt(8)/range0     
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
spde_matern <- inla.spde2.matern(mesh=mesh, B.tau=cbind(log(tau0), -1, +1),
                                 B.kappa=cbind(log(kappa0), 0, -1), theta.prior.mean=c(0,0),
                                 theta.prior.prec=c(0.1,0.1))

  
Q <- inla.spde.precision(spde=spde_matern, theta=c(0,0))
#u <- inla.qsample(n=1, Q=Q, mu=rep(4,nrow(Q)), seed=seed) #m√©dia 4; poderia ter sido acrescentada mais tarde 
u <- inla.qsample(n=1, Q=Q, mu=rep(4,nrow(Q))) #m√©dia 4; poderia ter sido acrescentada mais tarde 
  
### simulation grid
  
limlattice <- c(0,1) # limits of the lattice simulation 
grid <- as.matrix(expand.grid(
  seq(limlattice[1], limlattice[2], length.out=100),
  seq(limlattice[1], limlattice[2], length.out=100)))

#Campo simulado
u_ver <- as.vector(inla.spde.make.A(mesh=mesh, loc=grid)%*%as.vector(u))

  
### simulating a gaussian response variable

#intercept <- 4   #caso quisesse acrescentar a m√©dia mu=4 aqui
#predictor <- intercept + u_ver 
predictor <- u_ver 
prec.gaussian <- (0.1)^(-1)  #nugget precision

ysim <- rnorm(nrow(grid), mean=predictor, sd=prec.gaussian**(-1/2))
DFSim <- data.frame(x=grid[,1], y=grid[,2], ysim=ysim)

### Sample the data according to PP with intensity exp(beta*S(x))

sampData <- DFSim[sample(1:nrow(DFSim), size=n, prob = exp(beta*u_ver)),]  


#Vou usar a mesh que j√° tenho, √© parecida com a que usavamos nas versoes anteriores deste trabalho com o RandomFields
# ##Mesh
# fronteira<- cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
# mesh <- inla.mesh.2d(loc.domain = fronteira, max.edge = c(0.092, 0.2))
# nv <- mesh$n
# #plot(mesh,main="Mesh")
# #points(sampData$coords, pch=19, cex=.5)

##########################################################
#######Inclusion of Construct Covariate 5% neighbors######
##########################################################


###############################
###Luci test1 without covariate
###############################
fronteira<- cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
x <- fronteira
x.aux <- st_sf(st_sfc(st_polygon(list(x))))

#tamanho do lado da grelha
#library(distances)
dados<-sampData

#Decis„o Final: considerar a densidade dos pontos
area<-st_area(x.aux)
area
lado.grid<-sqrt(area/n)
lado.grid


library(mgcv)

res.1<-data.frame(lado.grid=lado.grid,corr=NA,p.value=NA)

for (lg in 1:length(lado.grid)){
  
  #grid
  grid <- st_make_grid(x.aux, cellsize = c(lado.grid[lg],lado.grid[lg]))
  
  #DeterminaÁ„o de quais os quadradinhos dentro do poligono: s„o os que tÍm todos os pontos de
  #grid[[i]] dentro da fronteira.pt
  #library(mgcv)
  grid.in<-NULL
  for (i in 1:length(grid)){
    aux<-in.out(as.matrix(x),as.matrix(grid[[i]])) #se algum ponto de grid[[i]] estiver dentro da fronteira vir· true(=1)
    if(sum(aux)==5) #5-1 vÈrtices do quadrado inÌcio para fechar
    {grid.in[i]<-"dentro"}
    else{
      if(sum(aux)==0) grid.in[i]<-"fora" else grid.in[i]<-"borda"}
  } #end for i
  
  #CriaÁ„o da data.frame onde colocar a info para o teste
  base.teste<-data.frame(ID=1:length(grid.in),pos.grelha=grid.in,n.pontos=NA,media.y=NA)
  
  #Areas da grid que est„o totalmente dentro da fronteira
  
  
  for (i in 1:length(grid)){
    if (base.teste$pos.grelha[i]=="dentro" | base.teste$pos.grelha[i]=="borda"){
      aux<-in.out(as.matrix(grid[[i]]),as.matrix(data.frame(dados$x,dados$y))) #ver que dados ent„o dentro do quadradinho da grelha
      base.teste$n.pontos[i]<-sum(aux)
      if (sum(aux)!=0) base.teste$media.y[i]<-mean(dados$ysim[aux],na.rm=T)
      else base.teste$n.pontos[i]<-NA
      # cat(paste("IteraÁ„o:",i,"\n"))
    }
  }
  base.teste<-base.teste[base.teste$pos.grelha=="dentro"|base.teste$pos.grelha=="borda",]
  
  teste<-cor.test(base.teste$n.pontos, base.teste$media.y,  method = "spearman", use="pairwise.complete.obs",exact=FALSE)
  res.1$corr[lg]<-teste$estimate
  res.1$p.value[lg]<-teste$p.value
  
  
} #end for lg



RES<-res.1
RES



####An·lise da grelha####

#Determinar n˙mero de quadrados da grelha
length(base.teste$n.pontos)
#Determinar n˙mero de obs
sum(base.teste$n.pontos, na.rm=T)
#Determinar n˙mero mÈdio de observaÁıes por quadrado
mean(base.teste$n.pontos,na.rm=T)
hist(base.teste$n.pontos)



#RES_NOcovariate<-c(l=RES$lado.grid, cor=RES$corr, pvalue=RES$p.value)
write.table(matrix(RES,ncol=3), file="RES_nocovariate_betaM2_250obs_5p_mesh1.txt", append=T,row.names=F, col.names=F)


#base.teste

#Constructed covariate nas observaÁıes
#Dist‚ncia aos k vizinhos mais prÛximos 5% dos vizinhos

x<-sampData[,1]
y<-sampData[,2]
marks<-sampData[,3]

viz<-0.05*n
nndistance1<-nndist(x,y, k=1:viz)
nndistance2<-as.data.frame(nndistance1)
nndistance<-apply(nndistance2,1,mean)

# # Model Fitting without Preferential Sampling with const. covariate

range0<-0.05
pr<-0.01
sigma0<-1
ps<-0.01


spde<-inla.spde2.pcmatern(mesh,alpha=2,prior.range=c(range0,pr),prior.sigma=c(sigma0,ps))

# # Define projector matrix from rough mesh for fast computation
proj_dados1 <- inla.spde.make.A(mesh = mesh, loc =cbind(x,y))
#proj_pred_dados1 <- inla.spde.make.A(mesh, loc = mesh$loc[,1:2]) # Identity matrix

#Create data matrix for inla

s.index<-inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)

stack_smooth_dados2 <- inla.stack(data=data.frame(y=marks),
                                  A=list(2,proj_dados1),
                                  effects=list(data.frame(Intercept=rep(1,length(x)),
                                                          cov1=nndistance),
                                               spatial.field=s.index),
                                  tag='obs')



stack_dados2 <- inla.stack(stack_smooth_dados2)
formula_smooth_dados2 <- y ~ -1 + Intercept + cov1 + f(spatial.field, model=spde)

result_smooth_dados2 <- inla(formula_smooth_dados2,
                             data=inla.stack.data(stack_dados2, spde = spde),
                             family="gaussian",
                             control.predictor = list(A = inla.stack.A(stack_dados2),
                                                      compute = TRUE),
                             num.threads = 2,
                             #control.mode=list(theta=c(2.3622, 3.7063, -0.9590),restart=T),
                             verbose = T)



summary(result_smooth_dados2)


aux1 <-result_smooth_dados2$summary.fitted.values[inla.stack.index(stack_dados2,"obs")$data, "mean"]

aux2 <-result_smooth_dados2$summary.fitted.values[inla.stack.index(stack_dados2,"obs")$data, "sd"]

resid1 <- (marks - aux1) /aux2



############################################
###Luci test2 with constructed covariate
############################################


#> Linking to GEOS 3.5.1, GDAL 2.1.3, proj.4 4.9.2
x <- fronteira
x.aux <- st_sf(st_sfc(st_polygon(list(x))))

#tamanho do lado da grelha
#library(distances)
dados2<-sampData
dados2$data<-resid1

#lado.grid continua a ser o definido pela densidade dos pontos

library(mgcv)

res.2<-data.frame(lado.grid=lado.grid,corr=NA,p.value=NA)

for (lg in 1:length(lado.grid)){
  
  #grid
  grid <- st_make_grid(x.aux, cellsize = c(lado.grid[lg],lado.grid[lg]))
  
  #DeterminaÁ„o de quais os quadradinhos dentro do poligono: s„o os que tÍm todos os pontos de
  #grid[[i]] dentro da fronteira.pt
  #library(mgcv)
  grid.in<-NULL
  for (i in 1:length(grid)){
    aux<-in.out(as.matrix(x),as.matrix(grid[[i]])) #se algum ponto de grid[[i]] estiver dentro da fronteira vir· true(=1)
    if(sum(aux)==5) #5-1 vÈrtices do quadrado inÌcio para fechar
    {grid.in[i]<-"dentro"}
    else{
      if(sum(aux)==0) grid.in[i]<-"fora" else grid.in[i]<-"borda"}
  } #end for i
  
  #CriaÁ„o da data.frame onde colocar a info para o teste
  base.teste<-data.frame(ID=1:length(grid.in),pos.grelha=grid.in,n.pontos=NA,media.y=NA)
  
  #Areas da grid que est„o totalmente dentro da fronteira
  
  
  for (i in 1:length(grid)){
    if (base.teste$pos.grelha[i]=="dentro" | base.teste$pos.grelha[i]=="borda"){
      aux<-in.out(as.matrix(grid[[i]]),as.matrix(data.frame(dados2$x,dados2$y))) #ver que dados ent„o dentro do quadradinho da grelha
      base.teste$n.pontos[i]<-sum(aux)
      if (sum(aux)!=0) base.teste$media.y[i]<-mean(dados2$data[aux],na.rm=T)
      else base.teste$n.pontos[i]<-NA
      # cat(paste("IteraÁ„o:",i,"\n"))
    }
  }
  
  base.teste<-base.teste[base.teste$pos.grelha=="dentro"|base.teste$pos.grelha=="borda",]
  teste<-cor.test(base.teste$n.pontos, base.teste$media.y,  method = "spearman", use="pairwise.complete.obs",exact=FALSE)
  res.2$corr[lg]<-teste$estimate
  res.2$p.value[lg]<-teste$p.value
  
  
} #end for lg



RES2<-res.2
RES2
#base.teste


#RES_covariate<-c(l=RES2$lado.grid, cor=RES2$corr, pvalue=RES2$p.value)
write.table(matrix(RES2,ncol=3), file="RES_covariate_betaM2_250obs_5p_mesh1.txt", append=T,row.names=F, col.names=F)


} # end for while

#To see results
#n=100
RES_nocovariate_betaM2_100obs_5p<-read.table(file="RES_nocovariate_betaM2_100obs_5p.txt")
RES_nocovariate_betaM2_100obs_5p
aa<-round(RES_nocovariate_betaM2_100obs_5p,5)
#library(xtable)
#xtable(aa, type = "latex", file = "filename2.tex")
RES1<-RES_nocovariate_betaM2_100obs_5p$V3<0.1
RES1
sum(RES1)


RES_covariate_betaM2_100obs_5p<-read.table(file="RES_covariate_betaM2_100obs_5p.txt")
RES_covariate_betaM2_100obs_5p
bb<-round(RES_covariate_betaM2_100obs_5p,5)

#xtable(bb, type = "latex", file = "filename3.tex")
RES2<-RES_covariate_betaM2_100obs_5p$V3<0.1
RES2
sum(RES2)

#n=50
RES_nocovariate_betaM2_50obs_5p<-read.table(file="RES_nocovariate_betaM2_50obs_5p.txt")
RES_nocovariate_betaM2_50obs_5p
aa<-round(RES_nocovariate_betaM2_50obs_5p,5)
#library(xtable)
#xtable(aa, type = "latex", file = "filename2.tex")
RES1<-RES_nocovariate_betaM2_50obs_5p$V3<0.1
RES1
sum(RES1)


RES_covariate_betaM2_50obs_5p<-read.table(file="RES_covariate_betaM2_50obs_5p.txt")
RES_covariate_betaM2_50obs_5p
bb<-round(RES_covariate_betaM2_50obs_5p,5)

#xtable(bb, type = "latex", file = "filename3.tex")
RES2<-RES_covariate_betaM2_50obs_5p$V3<0.1
RES2
sum(RES2)

#n=250
RES_nocovariate_betaM2_250obs_5p<-read.table(file="RES_nocovariate_betaM2_250obs_5p.txt")
RES1<-RES_nocovariate_betaM2_250obs_5p$V3<0.1
RES1
sum(RES1)


RES_covariate_betaM2_250obs_5p<-read.table(file="RES_covariate_betaM2_250obs_5p.txt")
RES2<-RES_covariate_betaM2_250obs_5p$V3<0.1
RES2
sum(RES2)


#To see results
#mesh1
#n=100
RES_nocovariate_betaM2_100obs_5p<-read.table(file="RES_nocovariate_betaM2_100obs_5p_mesh1.txt")
RES_nocovariate_betaM2_100obs_5p
aa<-round(RES_nocovariate_betaM2_100obs_5p,5)
#library(xtable)
#xtable(aa, type = "latex", file = "filename2.tex")
RES1<-RES_nocovariate_betaM2_100obs_5p$V3<0.05
RES1
sum(RES1)


RES_covariate_betaM2_100obs_5p<-read.table(file="RES_covariate_betaM2_100obs_5p_mesh1.txt")
RES_covariate_betaM2_100obs_5p
bb<-round(RES_covariate_betaM2_100obs_5p,5)

#xtable(bb, type = "latex", file = "filename3.tex")
RES2<-RES_covariate_betaM2_100obs_5p$V3<0.05
RES2
sum(RES2)


#n=50
RES_nocovariate_betaM2_50obs_5p<-read.table(file="RES_nocovariate_betaM2_50obs_5p_mesh1.txt")
RES_nocovariate_betaM2_50obs_5p
aa<-round(RES_nocovariate_betaM2_50obs_5p,5)
#library(xtable)
#xtable(aa, type = "latex", file = "filename2.tex")
RES1<-RES_nocovariate_betaM2_50obs_5p$V3<0.05
RES1
sum(RES1)


RES_covariate_betaM2_50obs_5p<-read.table(file="RES_covariate_betaM2_50obs_5p_mesh1.txt")
RES_covariate_betaM2_50obs_5p
bb<-round(RES_covariate_betaM2_50obs_5p,5)

#xtable(bb, type = "latex", file = "filename3.tex")
RES2<-RES_covariate_betaM2_50obs_5p$V3<0.05
RES2
sum(RES2)


#n=250
RES_nocovariate_betaM2_250obs_5p<-read.table(file="RES_nocovariate_betaM2_250obs_5p_mesh1.txt")
RES1<-RES_nocovariate_betaM2_250obs_5p$V3<0.01
RES1
sum(RES1)


RES_covariate_betaM2_250obs_5p<-read.table(file="RES_covariate_betaM2_250obs_5p_mesh1.txt")
RES2<-RES_covariate_betaM2_250obs_5p$V3<0.05
RES2
sum(RES2)



