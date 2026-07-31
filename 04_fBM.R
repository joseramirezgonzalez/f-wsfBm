#!/usr/bin/env Rscript

# =============================================================================
# FAMILIA 4: movimiento browniano fraccionario sigma W_H(t)
# =============================================================================
# Este es el punto que estaba mal en el codigo conjunto: aqui la verosimilitud
# y el AIC se calculan EXCLUSIVAMENTE con la covarianza de fBM.
#
# PARAMETROS QUE PUEDES CAMBIAR:
H_BOUNDS <- c(0.01, 0.99)
# =============================================================================

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Falta 'readxl'. Instala una vez con install.packages('readxl').")
}
args <- commandArgs(trailingOnly = TRUE)
DATA_FILE <- if (length(args)) args[[1L]] else file.choose()
if (!file.exists(DATA_FILE)) stop("No existe el archivo seleccionado.")

read_ten_series <- function(path) {
  raw <- as.data.frame(readxl::read_excel(path, sheet = "in", col_names = FALSE),
                       check.names = FALSE)
  width <- if (ncol(raw) >= 20L) 4L else 3L
  if (ncol(raw) < 5L * width) stop("El Excel no contiene cinco bloques de datos.")
  out <- list(); k <- 0L
  for (bat in 1:5) {
    j <- width * (bat - 1L) + 1L
    d <- raw[, j:(j + 2L), drop = FALSE]
    names(d) <- c("time", "Longitude", "Latitude")
    d[] <- lapply(d, function(x) suppressWarnings(as.numeric(x)))
    d <- d[complete.cases(d), , drop = FALSE]
    for (coordinate in c("Latitude", "Longitude")) {
      k <- k + 1L; x <- d[[coordinate]]
      out[[k]] <- list(bat = bat, coordinate = coordinate,
        t = d$time[-1L], y = x[-1L] - x[1L],
        n_original = length(x), n_used = length(x) - 1L)
    }
  }
  out
}

profile_scale_loglik <- function(y, K) {
  K <- (K + t(K)) / 2
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) return(list(ok = FALSE, logLik = -Inf, scale = NA_real_))
  z <- backsolve(R, y, transpose = TRUE)
  sigma2 <- sum(z^2) / length(y)
  if (!is.finite(sigma2) || sigma2 <= 0)
    return(list(ok = FALSE, logLik = -Inf, scale = NA_real_))
  logdet <- 2 * sum(log(diag(R)))
  ll <- -0.5 * (length(y) * (log(2 * pi) + 1 + log(sigma2)) + logdet)
  list(ok = TRUE, logLik = ll, scale = sigma2)
}

cov_fbm_unit <- function(t, H) {
  if (!is.finite(H) || H <= 0 || H >= 1) stop("H debe estar en (0,1).")
  p <- t^(2 * H)
  0.5 * (outer(p, rep(1, length(t))) +
           outer(rep(1, length(t)), p) - abs(outer(t, t, "-"))^(2 * H))
}

fit_one <- function(s) {
  objective <- function(H) {
    ans <- profile_scale_loglik(s$y, cov_fbm_unit(s$t, H))
    if (ans$ok) -ans$logLik else 1e100
  }
  opt <- optimize(objective, H_BOUNDS, tol = 1e-9)
  H <- opt$minimum
  prof <- profile_scale_loglik(s$y, cov_fbm_unit(s$t, H))
  k <- 2L  # H y sigma
  data.frame(
    bat = s$bat, coordinate = s$coordinate,
    n_original = s$n_original, n_used = s$n_used,
    model = "fBM", k = k, logLik_max = prof$logLik,
    AIC = -2 * prof$logLik + 2 * k,
    H = H, sigma = sqrt(prof$scale), sigma2 = prof$scale,
    maximum_at_parameter_limit =
      H <= H_BOUNDS[1] + 1e-5 || H >= H_BOUNDS[2] - 1e-5,
    convergence = 0L
  )
}

series <- read_ten_series(DATA_FILE)
results <- do.call(rbind, lapply(series, function(s) {
  cat(sprintf("Ajustando fBM: Bat %d, %s ...\n", s$bat, s$coordinate))
  fit_one(s)
}))
results <- results[order(results$bat, match(results$coordinate,
                                            c("Latitude", "Longitude"))), ]
OUTPUT_FILE <- file.path(dirname(normalizePath(DATA_FILE)), "resultados_fBM.csv")
write.csv(results, OUTPUT_FILE, row.names = FALSE)
print(results, row.names = FALSE, digits = 10)
cat("\nArchivo escrito en:", OUTPUT_FILE, "\n")

