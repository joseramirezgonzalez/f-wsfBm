#!/usr/bin/env Rscript

# =============================================================================
# FAMILIA 3: mu_{1/2}(t), caso iOU con H fijado en 1/2
# PARAMETROS QUE PUEDES CAMBIAR:
BETA_BOUNDS <- c(0.01,400)
DELTA <- 0.5
N_QUAD <- 128L
# =============================================================================

if(!requireNamespace("readxl",quietly=TRUE))stop("Falta 'readxl'.")
args<-commandArgs(trailingOnly=TRUE);DATA_FILE<-if(length(args))args[[1L]] else file.choose()

read_ten_series<-function(path){raw<-as.data.frame(readxl::read_excel(path,sheet="in",col_names=FALSE))
  width<-if(ncol(raw)>=20L)4L else 3L;out<-list();k<-0L
  for(bat in 1:5){j<-width*(bat-1L)+1L;d<-raw[,j:(j+2L),drop=FALSE]
    names(d)<-c("time","Longitude","Latitude");d[]<-lapply(d,function(x)suppressWarnings(as.numeric(x)))
    d<-d[complete.cases(d),,drop=FALSE]
    for(coordinate in c("Latitude","Longitude")){k<-k+1L;x<-d[[coordinate]]
      out[[k]]<-list(bat=bat,coordinate=coordinate,t=d$time[-1L],y=x[-1L]-x[1L],delta=DELTA,
        n_original=length(x),n_used=length(x)-1L)}};out}

profile_scale_loglik<-function(y,K){R<-tryCatch(chol((K+t(K))/2),error=function(e)NULL)
  if(is.null(R))return(list(ok=FALSE,logLik=-Inf,scale=NA_real_));z<-backsolve(R,y,transpose=TRUE)
  s2<-sum(z^2)/length(y);ll<--.5*(length(y)*(log(2*pi)+1+log(s2))+2*sum(log(diag(R))))
  list(ok=is.finite(ll),logLik=ll,scale=s2)}

gauss_legendre<-local({cache<-new.env(parent=emptyenv());function(n){key<-as.character(n)
  if(exists(key,cache,inherits=FALSE))return(get(key,cache));i<-1:(n-1L);off<-i/sqrt(4*i^2-1)
  J<-matrix(0,n,n);J[cbind(i,i+1L)]<-off;J[cbind(i+1L,i)]<-off;e<-eigen(J,symmetric=TRUE);o<-order(e$values)
  ans<-list(nodes=e$values[o],weights=2*e$vectors[1L,o]^2);assign(key,ans,cache);ans}})

cov_iou_grid<-function(N,delta,beta){H<-.5;a<-1;q<-gauss_legendre(N_QUAD)
  r<-delta*(q$nodes+1)/2;wr<-delta*q$weights/2;lag<-delta*(0:(N-1L))
  weight_r<-exp(-beta*r)*(-expm1(-2*beta*(delta-r)))/(2*beta)
  int2<-vapply(lag,function(L)sum(wr*weight_r*((L+r)^a+abs(L-r)^a)),numeric(1L))
  kd<-delta*(1:N);int1<-vapply(kd,function(x)sum(wr*exp(-beta*r)*(x-r)^a),numeric(1L))
  A<--expm1(-beta*delta)/beta;id<-abs(outer(0:(N-1L),0:(N-1L),"-"))+1L
  K<-.5*A*(outer(int1,rep(1,N))+outer(rep(1,N),int1))-.5*matrix(int2[id],N,N)
  ar<-exp(-beta*delta);if(N>=2L){for(i in 2:N)K[i,]<-K[i,]+ar*K[i-1L,]
    for(j in 2:N)K[,j]<-K[,j]+ar*K[,j-1L]};(K+t(K))/2}

cov_observed<-function(s,beta){id<-as.integer(round(s$t/s$delta));K<-cov_iou_grid(max(id),s$delta,beta)
  K[id,id,drop=FALSE]}

fit_one<-function(s){objective<-function(log_beta){ans<-tryCatch(
    profile_scale_loglik(s$y,cov_observed(s,exp(log_beta))),error=function(e)list(ok=FALSE,logLik=-Inf))
    if(ans$ok)-ans$logLik else 1e100}
  grid<-seq(log(BETA_BOUNDS[1]),log(BETA_BOUNDS[2]),length.out=70L);v<-vapply(grid,objective,numeric(1L))
  ids<-order(v)[1:min(5L,length(v))]
  fits<-lapply(ids,function(i)optimize(objective,c(grid[max(1L,i-1L)],grid[min(length(grid),i+1L)]),tol=1e-9))
  best<-fits[[which.min(vapply(fits,`[[`,numeric(1L),"objective"))]];beta<-exp(best$minimum)
  prof<-profile_scale_loglik(s$y,cov_observed(s,beta));k<-2L
  data.frame(bat=s$bat,coordinate=s$coordinate,n_original=s$n_original,n_used=s$n_used,
    model="mu_half_iOU",k=k,logLik_max=prof$logLik,AIC=-2*prof$logLik+2*k,
    beta=beta,sigma=sqrt(prof$scale),H_fixed=.5,
    maximum_at_parameter_limit=beta<=1.001*BETA_BOUNDS[1]||beta>=.999*BETA_BOUNDS[2],convergence=0L)}

series<-read_ten_series(DATA_FILE);results<-do.call(rbind,lapply(series,function(s){
  cat(sprintf("Ajustando mu_1/2: Bat %d, %s ...\n",s$bat,s$coordinate));fit_one(s)}))
results<-results[order(results$bat,match(results$coordinate,c("Latitude","Longitude"))),]
OUTPUT_FILE<-file.path(dirname(normalizePath(DATA_FILE)),"resultados_mu_half_iOU.csv")
write.csv(results,OUTPUT_FILE,row.names=FALSE);print(results,row.names=FALSE,digits=10)
cat("\nArchivo escrito en:",OUTPUT_FILE,"\n")

