#!/usr/bin/env Rscript

# =============================================================================
# FAMILIA 2: mu_H(t), posicion integrada de un fOU fraccionario (ifOU)
# =============================================================================
# Se estiman beta, sigma y H. Si el perfil de beta crece hasta BETA_MAX,
# se aplica UNICAMENTE la regla solicitada:
#   beta* = primer beta tal que ell_p(BETA_MAX)-ell_p(beta*) = 0.01.
#
# PARAMETROS QUE PUEDES CAMBIAR:
BETA_BOUNDS <- c(0.01, 400)
H_BOUNDS <- c(0.05, 0.98)
BETA_PLATEAU_TOLERANCE <- 1e-2
DELTA <- 0.5
N_QUAD <- 128L
# =============================================================================

if (!requireNamespace("readxl", quietly = TRUE))
  stop("Falta 'readxl'. Instala una vez con install.packages('readxl').")
args <- commandArgs(trailingOnly = TRUE)
DATA_FILE <- if (length(args)) args[[1L]] else file.choose()

read_ten_series <- function(path) {
  raw <- as.data.frame(readxl::read_excel(path,sheet="in",col_names=FALSE))
  width <- if(ncol(raw)>=20L)4L else 3L;out<-list();k<-0L
  for(bat in 1:5){j<-width*(bat-1L)+1L;d<-raw[,j:(j+2L),drop=FALSE]
    names(d)<-c("time","Longitude","Latitude")
    d[]<-lapply(d,function(x)suppressWarnings(as.numeric(x)));d<-d[complete.cases(d),,drop=FALSE]
    for(coordinate in c("Latitude","Longitude")){k<-k+1L;x<-d[[coordinate]]
      out[[k]]<-list(bat=bat,coordinate=coordinate,t=d$time[-1L],y=x[-1L]-x[1L],
        delta=DELTA,n_original=length(x),n_used=length(x)-1L)}}
  out
}

profile_scale_loglik <- function(y,K){R<-tryCatch(chol((K+t(K))/2),error=function(e)NULL)
  if(is.null(R))return(list(ok=FALSE,logLik=-Inf,scale=NA_real_))
  z<-backsolve(R,y,transpose=TRUE);sigma2<-sum(z^2)/length(y)
  if(!is.finite(sigma2)||sigma2<=0)return(list(ok=FALSE,logLik=-Inf,scale=NA_real_))
  ll<--.5*(length(y)*(log(2*pi)+1+log(sigma2))+2*sum(log(diag(R))))
  list(ok=is.finite(ll),logLik=ll,scale=sigma2)}

gauss_legendre <- local({
  cache <- new.env(parent=emptyenv())
  function(n){key<-as.character(n);if(exists(key,cache,inherits=FALSE))return(get(key,cache))
    i<-1:(n-1L);off<-i/sqrt(4*i^2-1);J<-matrix(0,n,n)
    J[cbind(i,i+1L)]<-off;J[cbind(i+1L,i)]<-off;e<-eigen(J,symmetric=TRUE)
    o<-order(e$values);ans<-list(nodes=e$values[o],weights=2*e$vectors[1L,o]^2)
    assign(key,ans,cache);ans}
})

# Misma covarianza de tu funcion cov_OUH_zeta_2 y la misma recursion S_B.
# La integral doble se reduce algebraicamente a una integral unidimensional;
# esto permite ajustar las diez trayectorias sin cambiar el modelo.
cov_ifou_grid <- function(N,delta,beta,H,n_quad=N_QUAD){
  if(beta<=0||H<=0||H>=1)stop("Parametros fuera del dominio.")
  q<-gauss_legendre(n_quad);r<-delta*(q$nodes+1)/2;wr<-delta*q$weights/2;a<-2*H
  lag<-delta*(0:(N-1L))
  weight_r<-exp(-beta*r)*(-expm1(-2*beta*(delta-r)))/(2*beta)
  int2<-vapply(lag,function(L)sum(wr*weight_r*((L+r)^a+abs(L-r)^a)),numeric(1L))
  kd<-delta*(1:N)
  int1<-vapply(kd,function(x)sum(wr*exp(-beta*r)*(x-r)^a),numeric(1L))
  A<--expm1(-beta*delta)/beta
  id<-abs(outer(0:(N-1L),0:(N-1L),"-"))+1L
  innovation<-.5*A*(outer(int1,rep(1,N))+outer(rep(1,N),int1))-.5*matrix(int2[id],N,N)
  ar<-exp(-beta*delta);K<-innovation
  if(N>=2L){for(i in 2:N)K[i,]<-K[i,]+ar*K[i-1L,]
    for(j in 2:N)K[,j]<-K[,j]+ar*K[,j-1L]}
  (K+t(K))/2
}

cov_observed <- function(s,beta,H){
  index<-as.integer(round(s$t/s$delta))
  if(max(abs(s$t-index*s$delta))>1e-10)stop("Los tiempos no estan en la malla DELTA.")
  K<-cov_ifou_grid(max(index),s$delta,beta,H)
  K[index,index,drop=FALSE]
}

fit_one <- function(s){
  lb<-log(BETA_BOUNDS);hb<-H_BOUNDS
  objective<-function(par){beta<-exp(par[1L]);H<-par[2L]
    ans<-tryCatch(profile_scale_loglik(s$y,cov_observed(s,beta,H)),
      error=function(e)list(ok=FALSE,logLik=-Inf))
    if(ans$ok)-ans$logLik else 1e100}
  beta_starts<-pmin(pmax(c(.3,2,10,60,200,390),BETA_BOUNDS[1]),BETA_BOUNDS[2])
  H_starts<-pmin(pmax(c(.20,.55,.85,.94),hb[1]),hb[2])
  starts<-unique(as.matrix(expand.grid(log_beta=log(beta_starts),H=H_starts)))
  fits<-lapply(seq_len(nrow(starts)),function(i)optim(starts[i,],objective,
    method="L-BFGS-B",lower=c(lb[1],hb[1]),upper=c(lb[2],hb[2]),
    control=list(maxit=350L,factr=1e7,pgtol=1e-7)))
  best<-fits[[which.min(vapply(fits,`[[`,numeric(1L),"value"))]]
  beta_raw<-exp(best$par[1L]);H_raw<-best$par[2L]

  cache<-new.env(parent=emptyenv())
  profile_beta<-function(beta){key<-sprintf("%.10f",beta)
    if(exists(key,cache,inherits=FALSE))return(get(key,cache))
    objH<-function(H){ans<-tryCatch(profile_scale_loglik(s$y,cov_observed(s,beta,H)),
      error=function(e)list(ok=FALSE,logLik=-Inf));if(ans$ok)-ans$logLik else 1e100}
    op<-optimize(objH,hb,tol=2e-8);H<-op$minimum
    K<-cov_observed(s,beta,H);p<-profile_scale_loglik(s$y,K)
    ans<-list(beta=beta,H=H,sigma=sqrt(p$scale),logLik=p$logLik,K=K)
    assign(key,ans,cache);ans}

  endpoint<-profile_beta(BETA_BOUNDS[2])
  raw<-profile_beta(beta_raw)
  plateau_case<-beta_raw>=.90*BETA_BOUNDS[2] || endpoint$logLik>=raw$logLik-1e-5
  if(!plateau_case){selected<-raw;type<-"interior";gap<-NA_real_}else{
    gapfun<-function(beta)endpoint$logLik-profile_beta(beta)$logLik-BETA_PLATEAU_TOLERANCE
    scan<-sort(unique(c(exp(seq(log(BETA_BOUNDS[1]),log(BETA_BOUNDS[2]),length.out=100L)),
      seq(BETA_BOUNDS[1],BETA_BOUNDS[2],length.out=140L))))
    val<-vapply(scan,gapfun,numeric(1L));cross<-which(val<=0 & c(TRUE,head(val,-1)>0));cross<-cross[cross>1L]
    if(!length(cross))stop("No se encontro beta* para Bat ",s$bat," ",s$coordinate)
    i<-cross[1L];bstar<-uniroot(gapfun,c(scan[i-1L],scan[i]),tol=5e-7)$root
    selected<-profile_beta(bstar);type<-"plateau_beta_star"
    gap<-endpoint$logLik-selected$logLik
  }
  k<-3L
  data.frame(bat=s$bat,coordinate=s$coordinate,n_original=s$n_original,n_used=s$n_used,
    model="mu_H_ifOU",k=k,logLik_selected=selected$logLik,AIC=-2*selected$logLik+2*k,
    beta=selected$beta,sigma=selected$sigma,H=selected$H,estimate_type=type,
    logLik_at_beta_max=endpoint$logLik,plateau_gap=gap,
    maximum_at_parameter_limit=type!="interior"||selected$H<=hb[1]+1e-5||selected$H>=hb[2]-1e-5,
    convergence=best$convergence)
}

series<-read_ten_series(DATA_FILE)
results<-do.call(rbind,lapply(series,function(s){
  cat(sprintf("Ajustando mu_H: Bat %d, %s ...\n",s$bat,s$coordinate));fit_one(s)}))
results<-results[order(results$bat,match(results$coordinate,c("Latitude","Longitude"))),]
OUTPUT_FILE<-file.path(dirname(normalizePath(DATA_FILE)),"resultados_mu_H_ifOU.csv")
write.csv(results,OUTPUT_FILE,row.names=FALSE);print(results,row.names=FALSE,digits=10)
cat("\nArchivo escrito en:",OUTPUT_FILE,"\n")

