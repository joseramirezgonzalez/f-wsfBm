
install.packages("pracma")
install.packages("MASS")
install.packages("mvtnorm")
install.packages("optimr")
install.packages("bbmle")
install.packages("mle.tools")
install.packages("plot3D")
install.packages("rworldmap")
install.packages("geosphere")
install.packages("maxLik")


#----Library Packages
library(maxLik)
library(pracma)
library(MASS)
library(mvtnorm)
library(optimr)
library(bbmle)
library(mle.tools)
library(plot3D)
library(rworldmap)
library(geosphere)


#-----------Proceso de trayectoria funcion polinomial
#-------------Process one---------------#
D_1<-function(u,s,t,sigma,a) #----Delta function
{
   return( (sigma^2) *(u^a) *
            (  -(t-u)*log(t-u)-(s-u)*log(s-u)+(t+s-2*u)*log(t+s-2*u) ) 
) 
} 
cov_1<-function(s,t,sigma,a)
{
  sapply(min(s,t),function(y){g<-function(z)
              {
                return(D_1(z,s=s,t=t,sigma=sigma,a=a))#---t_1 parametro de h
              }
         return(integrate(g,0,y)$val)
          }

) } 
cov_wfbm_1<-function(t,sigma,a){
n<-length(t)
K<-matrix(rep(0,n*n),n,n)
   for(i in 1:n) {
     for(j in 1:i) {
       K[i,j]<-cov_1(t[j],t[i],sigma,a)
       if(j<i) {
          K[j,i]<-K[i,j]
       }
} } 
return(K) } 




#-----------Proceso de trayectoria funcion exponencial------
#-------------Process one---------------#
D_2<-function(u,s,t,b) #----Delta function
{
   return(  (exp(-b*u)) *
            (  -(t-u)*log(t-u)-(s-u)*log(s-u)+(t+s-2*u)*log(t+s-2*u) ) 
) 
} 
cov_2<-function(s,t,b)
{
  sapply(min(s,t),function(y){g<-function(z)
              {
                return(D_2(z,s=s,t=t,b=b))#---t_1 parametro de h
              }
         return(integrate(g,0,y,rel.tol = 1e-6)$val)
          }

) } 


cov_wfbm_2<-function(t,b){
n<-length(t)
K<-matrix(rep(0,n*n),n,n)
   for(i in 1:n) {
     for(j in 1:i) {
       K[i,j]<-cov_2(t[j],t[i],b=b)
       if(j<i) {
          K[j,i]<-K[i,j]
       }
} } 
return(K) } 





log_verosimilitude_practical<-function(T,b,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n      
t<-seq(D,T,length=n)               
K<-cov_wfbm_2(t,b=b)
M<-K[index,index]
g<-((t(as.vector(mu[2:(k+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(k+1)]-mu[1]) )[1,1])/k
S<-g*M   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(c(dmvnorm(mu[2:(k+1)],rep(mu[1],k),S,log=TRUE),sqrt(b*g)))
}



log_verosimilitude_practical_2_2<-function(T,sigma,b,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n      
t<-seq(D,T,length=n)               
K<-cov_wfbm_2(t,b=b)
M<-((sigma^2)/b)*K[index,index]
   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(dmvnorm(mu[2:(k+1)],rep(mu[1],k),M,log=TRUE))
}





#-----------Proceso de trayectoria con logaritmo--------
D_3<-function(u,s,t) #----Delta function
{
   return(  
            (  -(t-u)*log(t-u)-(s-u)*log(s-u)+(t+s-2*u)*log(t+s-2*u) ) 
) 
} 
cov_3<-function(s,t)
{
  sapply(min(s,t),function(y){g<-function(z)
              {
                return(D_3(z,s=s,t=t))#---t_1 parametro de h
              }
         return(integrate(g,0,y,rel.tol = 1e-6)$val)
          }

) } 
cov_wfbm_3<-function(t){
K<-matrix(rep(0,n*n),n,n)
   for(i in 1:n) {
     for(j in 1:i) {
       K[i,j]<-cov_3(t[j],t[i])
       if(j<i) {
          K[j,i]<-K[i,j]
       }
} } 
return(K) } 

log_verosimilitude_practical_log<-function(T,sigma,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n      
t<-seq(D,T,length=n)               
K<-cov_wfbm_3(t)
M<-K[index,index]
S<-(sigma^2)*M   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(dmvnorm(mu[2:(k+1)],rep(mu[1],k),S,log=TRUE))
}







#----------Proceso de Jhoson

c_H<-function(u,v,h){
  return(0.5*(u^(2*h)+v^(2*h)-abs(u-v)^(2*h)))
 }


#-----Covariance funtion y

cov_OUH_zeta_2<-function(T,n,b,h) {
    D<-T/n
    index<-0
    G<-function(u,v){return(exp(-b*(u+v))*abs(v-u+D*index)^(2*h))}
    H<-function(u){return(exp(-b*u)*(D*(index+1)-u)^(2*h))}
    int_2<-rep(0,n)
    int_1<-rep(0,n)
    for(i in 1:n){
    int_2[i]<-integral2(G,0,D,0,D,reltol = 1e-6)$Q
    int_1[i]<-integrate(H,0,D)$val
    index<-i 
    }
    K<-matrix(rep(0,n*n),n,n)
    for(i in 1:n) {
        for(j in 1:i){
              K[i,j]<-0.5*((1-exp(-b*D))/b)*(int_1[i]+int_1[j])-0.5*int_2[i-j+1]
              K[j,i]<-K[i,j]
        }
    }
    return(K)
 }

#------Log-likelihood trajectorias velocidad--------#
                    
log_verosimilitude_practical_OU<-function(T,b,h,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n                     
K<-cov_OUH_zeta_2(T,n,b,h)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*D*b)
  } 
}
M<-(S_B%*%K%*%t(S_B))[index,index]
g<-(t(as.vector(mu[2:(k+1)]-mu[1]))%*%solve(M)%*%as.vector(mu[2:(k+1)]-mu[1]) )[1,1]/k
S<-g*M   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(dmvnorm(mu[2:(k+1)],rep(mu[1],k),S,log=TRUE))
}



log_verosimilitude_practical_OU_2<-function(T,sigma,b,mu,index){
k<-length(mu)-1
n<-index[k]                
D<-T/n                     
K<-cov_OUH_zeta_2(T,n,b,0.5)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*D*b)
  } 
}
M<-(S_B%*%K%*%t(S_B))[index,index]
S<-(sigma^2)*M   #----No hay estabilidad numerico con y por lo tanto consideramos S_B
return(dmvnorm(mu[2:(k+1)],rep(mu[1],k),S,log=TRUE))
}



#----------




#---Practical case: Study of whales
ruta<-file.choose()
library(readxl)
datos<-data.frame(read_excel(ruta)) #--Data

number<-1 #---Number of whale in data

cont<-0
for(j in 1:length(datos[,3*(number-1)+2])){
  if(datos[j,3*(number-1)+2]!="NA"){
  cont<-cont+1
 }
  
}

k<-cont
lon<-datos[1:k,3*(number-1)+2]
lat<-datos[1:k,3*number]
index<-2*datos[2:k,3*(number-1)+1] #--Index of day informations

D<-0.5
T<-index[k-1]*D #----initial + k-1 data
n<-T/D
t<-seq(D,T,length=n)


data <- st_as_sf(datos, coords = c("lon", "lat"), crs = 4326)
box <- st_bbox(data)



nc_osm <- get_tiles(box, crop = TRUE, zoom = 12)

pdf("bat.pdf", width = 6*0.965, height =7 )
plot(st_as_sfc(box), graticule = TRUE, axes = TRUE, xlab = "Longitude", ylab = "Latitude")
plot_tiles(nc_osm, add = TRUE)
plot(st_geometry(data), add = TRUE, type = "l", lwd = 3)
plot(st_geometry(data)[1], col = "red", add = TRUE, pch = 19)
plot(st_geometry(data)[62], col = "blue", add = TRUE, pch = 19)
legend("bottomright", c("Begin", "End", "Trajectory"),
    pch = c(19, 19, NA), lty = c(NA, NA, 1), col = c("red", "blue", "black"),
    inset = 0.01, bg = "white")
dev.off()

#------------Change to lon data---------------
log_vero_pract_OU<-function(x){
-log_verosimilitude_practical_OU(T,x,0.5,lat,index)
}







#------Profile Likelihood (b,H)

time <- proc.time()
block<-100
reg_b<-seq(0.1,20,length=block)
profile_b<-rep(0,length(reg_b))
for(i in 1:block){
   profile_b[i]<-log_vero_pract_OU(reg_b[i])
 }
 plot(reg_b,-profile_b,xlab=expression(beta),ylab=" ",main=expression(paste("Profile log-likelihood ",beta)),type="l")
beta_max_lon_OU<-opm(1,fn=log_vero_pract_OU, lower=0.1, upper=20,method="L-BFGS-B")$p1
lines(c(beta_max_lon_OU,beta_max_lon_OU),c(0,1350),col="red")

proc.time()-time


K_1<-cov_OUH_zeta_2(T,n,beta_max_lon_OU,0.5)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*beta_max_lon_OU*D)
  } 
}
M<-(S_B%*%K_1%*%t(S_B))[index,index]
sigma_max_lon_OU<-sqrt((t(as.vector(lat[2:k]-lat[1]))%*%solve(M)%*%as.vector(lat[2:k]-lat[1]) )[1,1]/(k-1))




beta_max_lon_OU
sigma_max_lon_OU

max_OU<-log_verosimilitude_practical_OU_2(T,sigma_max_lon_OU,beta_max_lon_OU,lat,index)




#------------Change to lon data---------------
log_vero_pract_OU_beta<-function(x){
	log_verosimilitude_practical_OU_2(T,sigma_max_lon_OU,x,lat,index)
	}



#------------Change to lon data---------------
log_vero_pract_OU_sigma<-function(x){
	log_verosimilitude_practical_OU_2(T,x,beta_max_lon_OU,lat,index)
	}


4-2*max_OU


par(mfrow=c(2,1))



block<-50
reg_aux<-seq(0.3,1.6,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_OU_beta(reg_aux[i])
}






ml<-maxLik(log_vero_pract_OU_beta,start=beta_max_lon_OU)
norm_1<-1/(-hessian(ml))
normal_approx_1<-dnorm(reg_aux,beta_max_lon_OU,sqrt(norm_1))



valc1_1<-qnorm(0.005,beta_max_lon_OU,sqrt(norm_1))
valc2_1<-qnorm(0.995,beta_max_lon_OU,sqrt(norm_1))
val1<-dnorm(valc1_1,beta_max_lon_OU,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",beta)),xlab=expression(beta),ylab=expression(paste("Profile Likelihood Ratio ",beta)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(beta_max_lon_OU,beta_max_lon_OU),c(0,3000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,3000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,3000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_OU)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)



block<-50
reg_aux<-seq(0.0031,0.0051,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_OU_sigma(reg_aux[i])
}


ml2<-maxLik(log_vero_pract_OU_sigma,start=sigma_max_lon_OU)
norm_1<-1/(-hessian(ml2))
normal_approx_1<-dnorm(reg_aux,sigma_max_lon_OU,sqrt(norm_1))





valc1_1<-qnorm(0.005,sigma_max_lon_OU,sqrt(norm_1))
valc2_1<-qnorm(0.995,sigma_max_lon_OU,sqrt(norm_1))
val1<-dnorm(valc1_1,sigma_max_lon_OU,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",sigma)),xlab=expression(sigma),ylab=expression(paste("Profile Likelihood Ratio ",sigma)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(sigma_max_lon_OU,sigma_max_lon_OU),c(0,3000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,3000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,3000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_OU)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)






#-------------Start latitude------------------


#----------------


#------------Change to lon data---------------
log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x,lat,index)[1]
}



#----------------
#------------Change to lon data---------------
log_vero_pract_2<-function(x){
log_verosimilitude_practical(T,x,lat,index)
}

par(mfrow=c(2,1))


#------Profile Likelihood (b,H)

time <- proc.time()
block<-50
reg_b_log<-seq(0.001,0.1,length=block)
reg_c_log<-rep(0,block)
profile_b<-rep(0,block)
aux<-rep(0,2)
for(i in 1:block){
   aux<-log_vero_pract_2(reg_b_log[i])
   profile_b[i]<-aux[1]
   reg_c_log[i]<-(aux[2]^2)/reg_b_log[i]
 }
plot(reg_b_log,profile_b,xlab=expression(beta),ylab=" ",main=expression(paste("Profile log-likelihood ",beta)),type="l")
lines(c(beta_max_lon,beta_max_lon),c(0,1350),col="red")


plot(reg_b_log,reg_c_log,,main="cociente",type="l")

proc.time()-time

beta_max_lon<-opm(0.002,fn=log_vero_pract, lower=0.001, upper=0.1,method="L-BFGS-B")$p1
sigma_max_lon<-log_vero_pract_2(beta_max_lon)[2]


beta_max_lon 
sigma_max_lon

max_2_2<--log_vero_pract(beta_max_lon)

#------------Change to lon data---------------
log_vero_pract_2_beta<-function(x){
	log_verosimilitude_practical_2_2(T,sigma_max_lon,x,lat,index)
	}


#------------Change to lon data---------------
log_vero_pract_2_sigma<-function(x){
	log_verosimilitude_practical_2_2(T,x,beta_max_lon,lat,index)
	}



4-2*max_2_2


par(mfrow=c(2,1))



block<-50
reg_aux<-seq(0.008,0.020,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_2_beta(reg_aux[i])
}






ml<-maxLik(log_vero_pract_2_beta,start=beta_max_lon)
norm_1<-1/(-hessian(ml))
normal_approx_1<-dnorm(reg_aux,beta_max_lon,sqrt(norm_1))


valc1_1<-qnorm(0.005,beta_max_lon,sqrt(norm_1))
valc2_1<-qnorm(0.995,beta_max_lon,sqrt(norm_1))
val1<-dnorm(valc1_1,beta_max_lon,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",beta)),xlab=expression(beta),ylab=expression(paste("Profile Likelihood Ratio ",beta)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(beta_max_lon,beta_max_lon),c(0,3000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,3000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,3000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_2_2)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


block<-50
reg_aux<-seq(0.00025,0.00046,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_2_sigma(reg_aux[i])
}



ml2<-maxLik(log_vero_pract_2_sigma,start=sigma_max_lon)
norm_1<-1/(-hessian(ml2))
normal_approx_1<-dnorm(reg_aux,sigma_max_lon,sqrt(norm_1))



valc1_1<-qnorm(0.005,sigma_max_lon,sqrt(norm_1))
valc2_1<-qnorm(0.995,sigma_max_lon,sqrt(norm_1))
val1<-dnorm(valc1_1,sigma_max_lon,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",sigma)),xlab=expression(sigma),ylab=expression(paste("Profile Likelihood Ratio ",sigma)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(sigma_max_lon,sigma_max_lon),c(0,50000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,50000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,50000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_2_2)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)













#-----End latitude



#-------------Start longitude------------------


#----------------


#------------Change to lon data---------------
log_vero_pract<-function(x){
-log_verosimilitude_practical(T,x,lon,index)[1]
}



#----------------
#------------Change to lon data---------------
log_vero_pract_2<-function(x){
log_verosimilitude_practical(T,x,lon,index)
}

par(mfrow=c(2,1))


#------Profile Likelihood (b,H)

time <- proc.time()
block<-100
reg_b_log<-seq(0.02,0.08,length=block)
reg_c_log<-rep(0,block)
profile_b<-rep(0,block)
aux<-rep(0,2)
for(i in 1:block){
   aux<-log_vero_pract_2(reg_b_log[i])
   profile_b[i]<-aux[1]
   reg_c_log[i]<-(aux[2]^2)/reg_b_log[i]
 }
plot(reg_b_log,profile_b,main="ACI zeta process",type="l",xlab=expression(beta),ylab=" ")



plot(reg_b_log,reg_c_log,,main="cociente",type="l")

proc.time()-time


beta_max_lat<-opm(0.05,fn=log_vero_pract, lower=0.02, upper=0.08,method="L-BFGS-B")$p1
sigma_max_lat<-log_vero_pract_2(beta_max_lat)[2]

beta_max_lat 
sigma_max_lat

max_2_2<--log_vero_pract(beta_max_lat)

4-2*max_2_2



#------------Change to lon data---------------
log_vero_pract_2_beta<-function(x){
	log_verosimilitude_practical_2_2(T,sigma_max_lat,x,lon,index)
	}


#------------Change to lon data---------------
log_vero_pract_2_sigma<-function(x){
	log_verosimilitude_practical_2_2(T,x,beta_max_lat,lon,index)
	}



4-2*max_2_2


par(mfrow=c(2,1))



block<-50
reg_aux<-seq(0.03,0.042,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_2_beta(reg_aux[i])
}


#ml<-maxLik(log_vero_pract_2_beta,start=beta_max_lat)
norm_1<-1/(-hessian(log_vero_pract_2_beta,beta_max_lat)[1,1])
normal_approx_1<-dnorm(reg_aux,beta_max_lat,sqrt(norm_1))


valc1_1<-qnorm(0.005,beta_max_lat,sqrt(norm_1))
valc2_1<-qnorm(0.995,beta_max_lat,sqrt(norm_1))
val1<-dnorm(valc1_1,beta_max_lat,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",beta)),xlab=expression(beta),ylab=expression(paste("Profile Likelihood Ratio ",beta)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(beta_max_lat,beta_max_lat),c(0,3000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,3000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,3000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_2_2)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


block<-50
reg_aux<-seq(0.0006,0.002,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_2_sigma(reg_aux[i])
}



#ml2<-maxLik(log_vero_pract_2_sigma,start=sigma_max_lat)
norm_1<-1/(-hessian(log_vero_pract_2_sigma,sigma_max_lat)[1,])
normal_approx_1<-dnorm(reg_aux,sigma_max_lat,sqrt(norm_1))



valc1_1<-qnorm(0.005,sigma_max_lat,sqrt(norm_1))
valc2_1<-qnorm(0.995,sigma_max_lat,sqrt(norm_1))
val1<-dnorm(valc1_1,sigma_max_lat,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",sigma)),xlab=expression(sigma),ylab=expression(paste("Profile Likelihood Ratio ",sigma)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(sigma_max_lat,sigma_max_lat),c(0,50000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,50000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,50000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_2_2)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)




#----------Poner lo del logaritmo, la perfil de b, el cociente entre la perfil de b y el maximo con el logaritmo


#------------Change to lon data---------------
log_vero_pract_log<-function(x){
-log_verosimilitude_practical_log(T,x,lon,index)
}




par(mfrow=c(2,1))


#------Profile Likelihood (b,H)

time <- proc.time()
block<-50
reg_s_log<-seq(0.001,0.01,length=block)
profile_s_log<-rep(0,block)

for(i in 1:block){
   profile_s_log[i]<-log_vero_pract_log(reg_s_log[i])
 }
plot(reg_s_log,-profile_s_log,xlab=expression(alpha),ylab=" ",main=expression(paste("Profile log-likelihood ",alpha)),type="l")



K<-cov_wfbm_3(t)
M<-K[index,index]
sigma_max_lon_log<-sqrt((t(as.vector(lon[2:k]-lon[1]))%*%solve(M)%*%as.vector(lon[2:k]-lon[1]) )[1,1]/(k-1))

lines(c(sigma_max_lon_log,sigma_max_lon_log),c(0,1350),col="red")

proc.time()-time







sigma_max_lon_log
max_log<-log_verosimilitude_practical_log(T,sigma_max_lon_log,lon,index)



#------------Change to lon data---------------
log_vero_pract_log_2<-function(x){
	exp(log_verosimilitude_practical_log(T,x,lon,index)-max_log)
	}



block<-50
reg_aux<-seq(0.0025,0.0042,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_log_2(reg_aux[i])
}


norm_1<-1/(-hessian(log_vero_pract_log_2,sigma_max_lon_log))
norm_1<-norm_1[1,1]
normal_approx_1<-dnorm(reg_aux,sigma_max_lon_log,sqrt(norm_1))

valc1_1<-qnorm(0.005,sigma_max_lon_log,sqrt(norm_1))
valc2_1<-qnorm(0.995,sigma_max_lon_log,sqrt(norm_1))
val1<-dnorm(valc1_1,sigma_max_lon_log,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",alpha)),xlab=expression(alpha),ylab=expression(paste("Profile Likelihood Ratio ",alpha)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(sigma_max_lon_log,sigma_max_lon_log),c(0,3000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,3000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,3000),col="green",lwd=2)
lines(reg_aux,val_perfil_aux*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)







#-------------------------



#-------------Poner le modelo de Jhoson--------------


log_vero_pract_OU<-function(x){
-log_verosimilitude_practical_OU(T,x,0.5,lon,index)
}


#------Profile Likelihood (b,H)

time <- proc.time()
block<-100
reg_b<-seq(0.1,20,length=block)
profile_b<-rep(0,length(reg_b))
for(i in 1:block){
   profile_b[i]<-log_vero_pract_OU(reg_b[i])
 }
 plot(reg_b,-profile_b,xlab=expression(beta),ylab=" ",main=expression(paste("Profile log-likelihood ",beta)),type="l")
beta_max_lat_OU<-opm(1,fn=log_vero_pract_OU, lower=0.1, upper=20,method="L-BFGS-B")$p1
lines(c(beta_max_lat_OU,beta_max_lat_OU),c(0,1350),col="red")

proc.time()-time


K_1<-cov_OUH_zeta_2(T,n,beta_max_lat_OU,0.5)
S_B<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  for(j in 1:i){
     S_B[i,j]<-exp((j-i)*beta_max_lat_OU*D)
  } 
}
M<-(S_B%*%K_1%*%t(S_B))[index,index]
sigma_max_lat_OU<-sqrt((t(as.vector(lon[2:k]-lon[1]))%*%solve(M)%*%as.vector(lon[2:k]-lon[1]) )[1,1]/(k-1))




beta_max_lat_OU
sigma_max_lat_OU

max_OU_lat<-log_verosimilitude_practical_OU_2(T,sigma_max_lat_OU,beta_max_lat_OU,lon,index)



#------------Change to lon data---------------
log_vero_pract_OU_beta<-function(x){
	log_verosimilitude_practical_OU_2(T,sigma_max_lat_OU,x,lon,index)
	}



#------------Change to lon data---------------
log_vero_pract_OU_sigma<-function(x){
	log_verosimilitude_practical_OU_2(T,x,beta_max_lat_OU,lon,index)
	}


4-2*max_OU_lat


par(mfrow=c(2,1))



block<-50
reg_aux<-seq(0.1,1,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_OU_beta(reg_aux[i])
}


ml<-maxLik(log_vero_pract_OU_beta,start=beta_max_lat_OU)
norm_1<-1/(-hessian(ml))
normal_approx_1<-dnorm(reg_aux,beta_max_lat_OU,sqrt(norm_1))



valc1_1<-qnorm(0.005,beta_max_lat_OU,sqrt(norm_1))
valc2_1<-qnorm(0.995,beta_max_lat_OU,sqrt(norm_1))
val1<-dnorm(valc1_1,beta_max_lat_OU,sqrt(norm_1))



plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",beta)),xlab=expression(beta),ylab=expression(paste("Profile Likelihood Ratio ",beta)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(beta_max_lat_OU,beta_max_lat_OU),c(0,3000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,3000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,3000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_OU_lat)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)



block<-50
reg_aux<-seq(0.0042,0.0072,length=block)
val_perfil_aux<-rep(0,block)

for(i in 1:block){
	val_perfil_aux[i]<-log_vero_pract_OU_sigma(reg_aux[i])
}


ml<-maxLik(log_vero_pract_OU_sigma,start=sigma_max_lat_OU)
norm_1<-1/(-hessian(ml))
normal_approx_1<-dnorm(reg_aux,sigma_max_lat_OU,sqrt(norm_1))

valc1_1<-qnorm(0.005,sigma_max_lat_OU,sqrt(norm_1))
valc2_1<-qnorm(0.995,sigma_max_lat_OU,sqrt(norm_1))
val1<-dnorm(valc1_1,sigma_max_lat_OU,sqrt(norm_1))


plot(reg_aux,normal_approx_1,type="l",main=expression(paste("Profile Likelihood Ratio ",sigma)),xlab=expression(sigma),ylab=expression(paste("Profile Likelihood Ratio ",sigma)),col="blue",lwd=2,cex.main=2,cex.lab=1.2)
lines(c(sigma_max_lat_OU,sigma_max_lat_OU),c(0,3000),col="red",lwd=2)
lines(c(valc1_1,valc1_1),c(0,3000),col="green",lwd=2)
lines(c(valc2_1,valc2_1),c(0,3000),col="green",lwd=2)
lines(reg_aux,exp(val_perfil_aux-max_OU_lat)*(1/(sqrt(norm_1*2*pi))),col="black",lwd=2)


#----End longitude-----------


#------------Simulaciones-------------


par(mfrow=c(2,1))

plot(c(0,index/2),lat,type="l")



T<-100
n<-1000
D<-T/n
t<-seq(D,T,length=n)

sigma<-1
beta<--0.4
alpha<--0.93


#-------------Parameters---------------#
K_1<-cov_wfbm_2(t,beta)
mu_i_1<-0
mu_1<-rep(0,n)
for(i in 1:n){
 mu_1[i]<-mu_i_1
}
sim_lon<-mvrnorm(  1, mu_1, K_1 )
plot(c(0,t),c(mu_i_1,sim_lon), type="l",xlab="t",ylab=expression(zeta["t,f"]),lwd=2,cex.main=2,cex.lab=1.2)



