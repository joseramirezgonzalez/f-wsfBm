#!/usr/bin/env Rscript

# =============================================================================
# FAMILIA 7: Stein S^(3), Cov(h)=sigma^2 (1+lambda |h|)^(-3/2)
# PARAMETROS QUE PUEDES CAMBIAR:
LAMBDA_BOUNDS <- c(1e-8, 1)
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
        n_original=length(x),n_used=length(x)-1L)}}
  out
}

profile_scale_loglik <- function(y,K){R<-tryCatch(chol((K+t(K))/2),error=function(e)NULL)
  if(is.null(R))return(list(ok=FALSE,logLik=-Inf,scale=NA_real_))
  z<-backsolve(R,y,transpose=TRUE);sigma2<-sum(z^2)/length(y)
  ll<--.5*(length(y)*(log(2*pi)+1+log(sigma2))+2*sum(log(diag(R))))
  list(ok=is.finite(ll),logLik=ll,scale=sigma2)}

cov_stein_unit <- function(t,lambda)(1+lambda*abs(outer(t,t,"-")))^(-3/2)

fit_one <- function(s){objective<-function(log_lambda){
    ans<-profile_scale_loglik(s$y,cov_stein_unit(s$t,exp(log_lambda)))
    if(ans$ok)-ans$logLik else 1e100}
  opt<-optimize(objective,log(LAMBDA_BOUNDS),tol=1e-10);lambda<-exp(opt$minimum)
  prof<-profile_scale_loglik(s$y,cov_stein_unit(s$t,lambda));k<-2L
  data.frame(bat=s$bat,coordinate=s$coordinate,n_original=s$n_original,n_used=s$n_used,
    model="Stein_S3",k=k,logLik_max=prof$logLik,AIC=-2*prof$logLik+2*k,
    sigma2=prof$scale,lambda=lambda,
    maximum_at_parameter_limit=lambda<=1.001*LAMBDA_BOUNDS[1]||lambda>=.999*LAMBDA_BOUNDS[2],
    convergence=0L)}

series<-read_ten_series(DATA_FILE)
results<-do.call(rbind,lapply(series,function(s){
  cat(sprintf("Ajustando Stein S3: Bat %d, %s ...\n",s$bat,s$coordinate));fit_one(s)}))
results<-results[order(results$bat,match(results$coordinate,c("Latitude","Longitude"))),]
OUTPUT_FILE<-file.path(dirname(normalizePath(DATA_FILE)),"resultados_Stein_S3.csv")
write.csv(results,OUTPUT_FILE,row.names=FALSE);print(results,row.names=FALSE,digits=10)
cat("\nArchivo escrito en:",OUTPUT_FILE,"\n")

