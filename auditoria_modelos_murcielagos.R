#!/usr/bin/env Rscript

# =============================================================================
# AUDITORIA REPRODUCIBLE DE LOS SIETE MODELOS DE COVARIANZA (TABLA 1)
# =============================================================================
#
# Objetivo
# --------
# Este programa vuelve a ajustar, para Latitud y Longitud de cinco
# murcielagos, las siete familias gaussianas comparadas en la tabla:
#
#   1. zeta_{t,f_{sigma^2,beta}}
#   2. ifOU: mu_H(t)
#   3. iOU:  mu_{1/2}(t)
#   4. fBM:  sigma W_H(t)
#   5. Confluent: A_{sigma^2,eta,alpha,beta}(t)
#   6. Stein S^{(2)}_{sigma^2,lambda}(t)
#   7. Stein S^{(3)}_{sigma^2,lambda}(t)
#
# El archivo esta comentado deliberadamente con mucho detalle. Cada funcion
# de covarianza indica que proceso representa, que parametros usa y que
# cantidad devuelve. No se anade "jitter" silencioso a las matrices: si una
# matriz no es definida positiva, el ajuste se rechaza y queda documentado.
#
# Correcciones metodologicas incorporadas
# ----------------------------------------
# * Todos los modelos usan EXACTAMENTE el mismo vector de n-1 observaciones:
#       y_i = X(t_i) - X(0),  i = 1,...,n-1.
#   Se excluye el cero determinista y_0 = X(0)-X(0) = 0.
# * fBM se evalua con SU PROPIA covarianza/verosimilitud. El codigo publico
#   original calculaba su AIC llamando por error a la verosimilitud iOU.
# * El factor 2 de la covarianza zeta se incluye de acuerdo con la formula del
#   manuscrito. El factor cambia sigma^2, aunque no cambia la log-verosimilitud
#   perfilada ni el AIC.
# * Para ifOU se conserva la regla de meseta solicitada en el articulo:
#       ell_p(400) - ell_p(beta*) = 0.01,
#   y beta* es el primer cruce de la meseta. No se sustituye beta* por 400.
# * Confluent y Stein se reajustan sobre las mismas n-1 observaciones que los
#   modelos no estacionarios.
# * La frontera beta=0 de zeta se compara usando k=1, pues bajo la convencion
#   declarada solo se estima sigma^2 en ese caso.
# * En Confluent se impone alpha <= 15. El programa marca expresamente cuando
#   la solucion cae en esa frontera; no la llama MLE interior.
#
# Modelo probabilistico y verosimilitud
# -------------------------------------
# Para cada familia se supone
#       y ~ N(0, c K(theta)),
# donde K(theta) es una matriz de covarianza con escala unitaria y c es la
# escala. La escala se perfila analiticamente:
#       c_hat = y' K(theta)^(-1) y / m,
# con m=n-1. Asi la optimizacion numerica solo busca los parametros de forma.
# La log-verosimilitud gaussiana maximizada se calcula mediante Cholesky:
#       ell = -1/2 [m log(2 pi) + log|c K| + y'(c K)^(-1)y].
# Finalmente,
#       AIC = -2 ell + 2 k,
# donde k cuenta TODOS los parametros estimados, incluida la escala.
#
# Dependencias
# ------------
# install.packages(c("readxl", "gsl", "DEoptim"))
#
# Ejecucion (desde una terminal)
# -----------------------------
# Rscript auditoria_modelos_murcielagos.R \
#   --data=bat_set.xlsx \
#   --out=resultados_auditoria
#
# Opciones utiles:
#   --models=fBM,zeta,ifOU   ajusta solo esas familias.
#   --quick                 reduce cuadraturas/iteraciones para una prueba.
#                           NO usar --quick para cifras de publicacion.
#   --seed=6501             fija la semilla de DEoptim.
#
# Archivos producidos
# -------------------
# resultados_todos_los_modelos.csv : parametros, logLik, AIC y diagnosticos.
# ganadores_por_serie.csv           : modelo con AIC minimo por serie.
# comparacion_con_referencia.csv    : diferencia contra 70 AIC auditados.
# resumen_verificacion.txt          : lectura humana de las comprobaciones.
# sessionInfo.txt                   : version exacta de R y paquetes.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# ----------------------------------------------------------------------------
# 0. Argumentos, configuracion y dependencias
# ----------------------------------------------------------------------------

script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (!length(file_arg)) return(getwd())
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
}

get_argument <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

ARGS <- commandArgs(trailingOnly = TRUE)
SCRIPT_DIR <- script_directory()
DATA_PATH <- get_argument(ARGS, "data", file.path(SCRIPT_DIR, "bat_set.xlsx"))
OUT_DIR <- get_argument(ARGS, "out", file.path(SCRIPT_DIR, "resultados_auditoria"))
SEED <- as.integer(get_argument(ARGS, "seed", "6501"))
QUICK <- "--quick" %in% ARGS

ALL_MODELS <- c("zeta", "ifOU", "iOU", "fBM", "Confluent", "Stein_S2", "Stein_S3")
MODEL_ARGUMENT <- get_argument(ARGS, "models", paste(ALL_MODELS, collapse = ","))
MODELS_TO_RUN <- trimws(strsplit(MODEL_ARGUMENT, ",", fixed = TRUE)[[1L]])
unknown_models <- setdiff(MODELS_TO_RUN, ALL_MODELS)
if (length(unknown_models)) {
  stop("Modelo(s) desconocido(s): ", paste(unknown_models, collapse = ", "))
}

required_packages <- c("readxl", "gsl", "DEoptim")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Faltan paquetes de R: ", paste(missing_packages, collapse = ", "), "\n",
    "Instalelos con:\n",
    "install.packages(c(\"readxl\", \"gsl\", \"DEoptim\"))"
  )
}

CFG <- list(
  delta = 0.5,               # paso de la malla temporal, en minutos
  beta_max_ifou = 400,
  plateau_tolerance = 1e-2,  # diferencia de logLik que define beta*
  h_bounds = c(0.05, 0.98),
  beta_bounds_ifou = c(0.01, 400),
  beta_bounds_iou = c(0.01, 400),
  beta_bounds_zeta = c(1e-7, 2),
  lambda_bounds = c(1e-8, 1),
  ch_eta_bounds = c(0.005, 5),
  ch_alpha_bounds = c(0.01, 15),
  ch_beta_bounds = c(0.001, 50000),
  n_gl_ifou = if (QUICK) 32L else 72L,
  n_gl_zeta_coarse = if (QUICK) 40L else 64L,
  n_gl_zeta_final = if (QUICK) 80L else 160L,
  de_iter = if (QUICK) 25L else 90L,
  de_population = if (QUICK) 24L else 42L
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

# ----------------------------------------------------------------------------
# 1. Lectura y preparacion de las diez series
# ----------------------------------------------------------------------------

read_bat_series <- function(path, delta = 0.5) {
  if (!file.exists(path)) stop("No se encontro el archivo de datos: ", path)

  raw <- as.data.frame(
    readxl::read_excel(path, sheet = "in", col_names = FALSE),
    check.names = FALSE
  )

  if (ncol(raw) < 15L) {
    stop("La hoja 'in' debe contener cinco bloques de tres columnas (15 columnas).")
  }

  result <- list()
  cursor <- 1L
  for (bat in 1:5) {
    columns <- (3L * (bat - 1L) + 1L):(3L * bat)
    block <- raw[, columns, drop = FALSE]
    keep <- complete.cases(block)
    block <- block[keep, , drop = FALSE]
    block[] <- lapply(block, as.numeric)

    time_all <- block[[1L]]
    longitude_all <- block[[2L]]
    latitude_all <- block[[3L]]

    if (time_all[[1L]] != 0) stop("Bat ", bat, ": el primer tiempo no es cero.")
    if (is.unsorted(time_all, strictly = TRUE)) {
      stop("Bat ", bat, ": los tiempos no son estrictamente crecientes.")
    }
    grid_index <- round(time_all / delta)
    if (max(abs(time_all - delta * grid_index)) > 1e-10) {
      stop("Bat ", bat, ": existe un tiempo fuera de la malla de ", delta, " minutos.")
    }

    # IMPORTANTE: se elimina la primera observacion, que seria exactamente 0.
    # Tanto las familias estacionarias como las no estacionarias reciben este
    # mismo vector y, por tanto, sus AIC se basan en la misma muestra.
    for (coordinate in c("Latitude", "Longitude")) {
      x_all <- if (coordinate == "Latitude") latitude_all else longitude_all
      result[[cursor]] <- list(
        bat = bat,
        coordinate = coordinate,
        time_all = time_all,
        x_all = x_all,
        t = time_all[-1L],
        y = x_all[-1L] - x_all[[1L]],
        n_original = length(time_all),
        n_used = length(time_all) - 1L,
        delta = delta
      )
      cursor <- cursor + 1L
    }
  }

  # Esta comprobacion detecta de inmediato si se cargo otra version del libro.
  expected_n <- c(`1` = 47L, `2` = 45L, `3` = 62L, `4` = 241L, `5` = 165L)
  observed_n <- vapply(result[c(TRUE, FALSE)], `[[`, integer(1L), "n_original")
  names(observed_n) <- as.character(1:5)
  if (!identical(unname(observed_n), unname(expected_n))) {
    warning(
      "Los tamanos de muestra no coinciden con bat_set.xlsx auditado. Observados: ",
      paste(observed_n, collapse = ", ")
    )
  }
  result
}

# ----------------------------------------------------------------------------
# 2. Log-verosimilitud gaussiana con escala perfilada
# ----------------------------------------------------------------------------

profile_scale_loglik <- function(y, K) {
  # K es la covarianza con escala unitaria. La simetrizacion solo elimina
  # error de redondeo de maquina; no altera la formula matematica.
  K <- (K + t(K)) / 2

  if (any(!is.finite(K))) {
    return(list(ok = FALSE, logLik = -Inf, scale = NA_real_, chol = NULL))
  }

  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) {
    return(list(ok = FALSE, logLik = -Inf, scale = NA_real_, chol = NULL))
  }

  # R es triangular superior y satisface K = R'R.
  # z = R'^(-1)y, por lo que z'z = y'K^(-1)y.
  z <- backsolve(R, y, transpose = TRUE)
  m <- length(y)
  scale_hat <- sum(z * z) / m

  if (!is.finite(scale_hat) || scale_hat <= 0) {
    return(list(ok = FALSE, logLik = -Inf, scale = NA_real_, chol = NULL))
  }

  log_det_K <- 2 * sum(log(diag(R)))
  loglik <- -0.5 * (
    m * (log(2 * pi) + 1 + log(scale_hat)) + log_det_K
  )

  list(ok = TRUE, logLik = as.numeric(loglik), scale = scale_hat, chol = R)
}

covariance_diagnostics <- function(K) {
  K <- (K + t(K)) / 2
  R <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(R)) {
    return(list(pd = FALSE, min_chol_diagonal = NA_real_, reciprocal_condition = 0))
  }
  list(
    pd = TRUE,
    min_chol_diagonal = min(diag(R)),
    # rcond cercano a cero indica una matriz numericamente mal condicionada.
    reciprocal_condition = rcond(K)
  )
}

# ----------------------------------------------------------------------------
# 3. Cuadratura de Gauss-Legendre (implementada en R base)
# ----------------------------------------------------------------------------

gauss_legendre <- local({
  cache <- new.env(parent = emptyenv())
  function(n) {
    key <- as.character(as.integer(n))
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    n <- as.integer(n)
    if (n < 2L) stop("La cuadratura necesita al menos dos nodos.")
    i <- 1:(n - 1L)
    off_diagonal <- i / sqrt(4 * i^2 - 1)
    J <- matrix(0, n, n)
    J[cbind(i, i + 1L)] <- off_diagonal
    J[cbind(i + 1L, i)] <- off_diagonal
    eig <- eigen(J, symmetric = TRUE)
    order_nodes <- order(eig$values)
    answer <- list(
      nodes = eig$values[order_nodes],
      weights = 2 * eig$vectors[1L, order_nodes]^2
    )
    assign(key, answer, envir = cache)
    answer
  }
})

x_log_x <- function(x) {
  # Extension continua de x log(x) en x=0: lim_{x->0+} x log(x)=0.
  answer <- numeric(length(x))
  dim(answer) <- dim(x)
  take <- x > 0
  answer[take] <- x[take] * log(x[take])
  answer
}

phi_x_log_x <- function(x) {
  # Primitiva de x log(x), tambien extendida continuamente en cero.
  answer <- numeric(length(x))
  dim(answer) <- dim(x)
  take <- x > 0
  answer[take] <- 0.5 * x[take]^2 * log(x[take]) - 0.25 * x[take]^2
  answer
}

# ----------------------------------------------------------------------------
# 4. Las siete matrices de covarianza
# ----------------------------------------------------------------------------

cov_zeta_beta_zero <- function(t) {
  # CASO DE FRONTERA beta=0.
  #
  # Si f_{sigma^2,0}(u)=sigma^2, la integral de la familia zeta posee esta
  # expresion cerrada. La funcion devuelve K para sigma^2=1; luego la
  # verosimilitud estima sigma^2 por perfil.
  ti <- outer(t, rep(1, length(t)))
  tj <- t(ti)
  K <- phi_x_log_x(ti + tj) + phi_x_log_x(abs(ti - tj)) -
    2 * phi_x_log_x(ti) - 2 * phi_x_log_x(tj)
  (K + t(K)) / 2
}

cov_zeta <- function(t, beta, n_quad = CFG$n_gl_zeta_final) {
  # COVARIANZA zeta_{t,f_{sigma^2,beta}}.
  #
  # Para sigma^2=1 calcula
  #
  #  K(s,t) = 2 integral_0^{min(s,t)} exp(-beta u) [
  #      (s+t-2u)log(s+t-2u)
  #      -(s-u)log(s-u) -(t-u)log(t-u)] du.
  #
  # beta controla el decaimiento de f(u)=sigma^2 exp(-beta u). El factor 2
  # es parte de la formula del manuscrito. La integral se evalua de forma
  # vectorizada con Gauss-Legendre en [0,min(s,t)].
  if (beta == 0) return(cov_zeta_beta_zero(t))
  if (!is.finite(beta) || beta < 0) stop("beta de zeta debe ser no negativo.")

  n <- length(t)
  ti <- outer(t, rep(1, n))
  tj <- t(ti)
  upper <- pmin(ti, tj)
  total <- ti + tj
  quad <- gauss_legendre(n_quad)
  K <- matrix(0, n, n)

  for (q in seq_along(quad$nodes)) {
    # Transformacion [-1,1] -> [0,upper].
    u <- upper * (quad$nodes[[q]] + 1) / 2
    weight <- upper * quad$weights[[q]] / 2
    integrand <- x_log_x(total - 2 * u) -
      x_log_x(ti - u) - x_log_x(tj - u)
    K <- K + 2 * weight * exp(-beta * u) * integrand
  }
  (K + t(K)) / 2
}

cov_fbm <- function(t, H) {
  # COVARIANZA DE MOVIMIENTO BROWNIANO FRACCIONARIO (fBM).
  #
  # Para sigma=1 calcula
  #   K(s,t) = 1/2 [s^(2H) + t^(2H) - |t-s|^(2H)].
  #
  # H in (0,1) es el parametro de Hurst:
  #   H=1/2  : Browniano ordinario;
  #   H>1/2  : incrementos persistentemente correlacionados;
  #   H<1/2  : incrementos antipersistentes.
  # La escala perfilada es sigma^2; en la salida se reporta sigma=sqrt(scale).
  if (!is.finite(H) || H <= 0 || H >= 1) stop("H debe pertenecer a (0,1).")
  p <- t^(2 * H)
  0.5 * (outer(p, rep(1, length(t))) + outer(rep(1, length(t)), p) -
    abs(outer(t, t, "-"))^(2 * H))
}

cov_stein <- function(t, lambda, dimension) {
  # COVARIANZA DE STEIN.
  #
  # Para sigma^2=1 calcula
  #   K(s,t) = (1 + lambda |t-s|)^(-d/2),
  # con d=2 para S y d=3 para S^(3).
  #
  # lambda>0 controla la rapidez del decaimiento temporal. La escala perfilada
  # se reporta como sigma^2. Es una covarianza estacionaria: depende solo del
  # rezago |t-s|, no de los tiempos absolutos.
  if (!is.finite(lambda) || lambda <= 0) stop("lambda debe ser positivo.")
  if (!dimension %in% c(2L, 3L)) stop("dimension debe ser 2 o 3.")
  (1 + lambda * abs(outer(t, t, "-")))^(-dimension / 2)
}

cov_confluent <- function(t, eta, alpha, beta) {
  # COVARIANZA HIPERGEOMETRICA CONFLUENTE (CH).
  #
  # Para sigma^2=1 calcula
  #   K(s,t) = Gamma(eta+alpha)/Gamma(eta)
  #            * U(alpha, 1-eta, eta (|t-s|/beta)^2),
  # donde U es la funcion hipergeometrica confluente de Tricomi.
  #
  # eta, alpha y beta son positivos. beta es una escala temporal; eta y alpha
  # controlan forma, regularidad y cola. En rezago cero se usa exactamente
  # K(t,t)=1 para evitar evaluar U en z=0.
  if (any(!is.finite(c(eta, alpha, beta))) || any(c(eta, alpha, beta) <= 0)) {
    stop("eta, alpha y beta deben ser positivos.")
  }

  lag <- abs(outer(t, t, "-"))
  unique_lag <- sort(unique(as.vector(lag)))
  z <- eta * (unique_lag / beta)^2
  value <- rep(NA_real_, length(z))
  value[z == 0] <- 1
  positive <- z > 0

  if (any(positive)) {
    prefactor <- exp(lgamma(eta + alpha) - lgamma(eta))
    # strict=FALSE permite recuperar el valor numerico junto a situaciones en
    # que GSL solo emite un estado de precision; la finitud y la positividad se
    # comprueban inmediatamente despues y Cholesky valida la matriz completa.
    u_value <- gsl::hyperg_U(alpha, 1 - eta, z[positive], strict = FALSE)
    value[positive] <- prefactor * u_value
  }

  if (any(!is.finite(value)) || any(value <= 0)) {
    stop("La covarianza Confluent produjo valores no finitos o no positivos.")
  }
  K <- matrix(value[match(as.vector(lag), unique_lag)], nrow(lag), ncol(lag))
  (K + t(K)) / 2
}

cov_ifou_grid <- function(max_grid_index, delta, beta, H, n_quad = CFG$n_gl_ifou) {
  # COVARIANZA DEL PROCESO ifOU EN UNA MALLA REGULAR.
  #
  # El proceso de posicion es
  #   mu_H(t) = integral_0^t exp[-beta(t-u)] W_H(u) du,
  # y, para sigma=1,
  #   Cov(mu_H(t),mu_H(s))
  #     = integral_0^t integral_0^s exp(-beta v)
  #       c_H(t-v,s-u) exp(-beta u) du dv,
  #   c_H(x,y)=1/2[x^(2H)+y^(2H)-|x-y|^(2H)].
  #
  # Evaluar esa doble integral para cada par (s,t) seria muy costoso. Esta
  # funcion hace una reduccion algebraica exacta a integrales unidimensionales
  # dentro de cada intervalo de longitud delta y despues aplica la recursion
  # OU en ambas direcciones de la matriz. La aproximacion numerica restante es
  # solo la cuadratura Gauss-Legendre unidimensional.
  #
  # iOU es exactamente el caso H=1/2 y usa esta misma funcion; de ese modo no
  # puede evaluarse accidentalmente con la covarianza de otra familia.
  if (beta <= 0 || H <= 0 || H >= 1) stop("Parametros ifOU fuera de dominio.")
  N <- as.integer(max_grid_index)
  a <- 2 * H
  quad <- gauss_legendre(n_quad)
  r <- delta * (quad$nodes + 1) / 2
  wr <- delta * quad$weights / 2
  lag <- delta * (0:(N - 1L))

  # integral doble reducida en r=v-u. La expresion entre parentesis combina
  # los rezagos positivos y negativos.
  weight_r <- exp(-beta * r) * (-expm1(-2 * beta * (delta - r))) / (2 * beta)
  int2 <- vapply(lag, function(L) {
    sum(wr * weight_r * ((L + r)^a + abs(L - r)^a))
  }, numeric(1L))

  # Primera parte de c_H: integral exp(-beta u)(k*delta-u)^(2H) du.
  k_delta <- delta * (1:N)
  int1 <- vapply(k_delta, function(kd) {
    sum(wr * exp(-beta * r) * (kd - r)^a)
  }, numeric(1L))

  A <- -expm1(-beta * delta) / beta
  lag_index <- abs(outer(0:(N - 1L), 0:(N - 1L), "-")) + 1L
  innovation <- 0.5 * A * (
    outer(int1, rep(1, N)) + outer(rep(1, N), int1)
  ) - 0.5 * matrix(int2[lag_index], N, N)

  # Si ar=exp(-beta*delta), filtrar por filas y columnas equivale a acumular
  # las contribuciones de todos los intervalos pasados de la integral OU.
  ar <- exp(-beta * delta)
  K <- innovation
  if (N >= 2L) {
    for (i in 2:N) K[i, ] <- K[i, ] + ar * K[i - 1L, ]
    for (j in 2:N) K[, j] <- K[, j] + ar * K[, j - 1L]
  }
  (K + t(K)) / 2
}

cov_ifou_observed <- function(series, beta, H, n_quad = CFG$n_gl_ifou) {
  # Construye primero la malla completa de 0.5 min y luego extrae solo los
  # tiempos realmente observados. Esto respeta los huecos irregulares de cada
  # trayectoria sin fingir observaciones ausentes.
  index <- as.integer(round(series$t / series$delta))
  if (max(abs(series$t - series$delta * index)) > 1e-10) {
    stop("Los tiempos no pertenecen a la malla declarada.")
  }
  full <- cov_ifou_grid(max(index), series$delta, beta, H, n_quad)
  full[index, index, drop = FALSE]
}

# ----------------------------------------------------------------------------
# 5. Utilidades comunes para devolver y auditar cada ajuste
# ----------------------------------------------------------------------------

empty_parameter_list <- function() {
  list(
    beta = NA_real_, sigma = NA_real_, sigma2 = NA_real_, H = NA_real_,
    eta = NA_real_, alpha = NA_real_, lambda = NA_real_
  )
}

new_result <- function(series, model, k, likelihood, parameters, K,
                       estimate_type = "interior", convergence = 0L,
                       note = "") {
  pars <- modifyList(empty_parameter_list(), parameters)
  diagnostics <- covariance_diagnostics(K)
  data.frame(
    bat = series$bat,
    coordinate = series$coordinate,
    n_original = series$n_original,
    n_used = series$n_used,
    model = model,
    k = as.integer(k),
    logLik = likelihood,
    AIC = -2 * likelihood + 2 * k,
    beta = pars$beta,
    sigma = pars$sigma,
    sigma2 = pars$sigma2,
    H = pars$H,
    eta = pars$eta,
    alpha = pars$alpha,
    lambda = pars$lambda,
    estimate_type = estimate_type,
    convergence = as.integer(convergence),
    covariance_pd = diagnostics$pd,
    min_chol_diagonal = diagnostics$min_chol_diagonal,
    reciprocal_condition = diagnostics$reciprocal_condition,
    note = note,
    check.names = FALSE
  )
}

safe_objective <- function(y, K_function) {
  tryCatch({
    likelihood <- profile_scale_loglik(y, K_function())
    if (!likelihood$ok || !is.finite(likelihood$logLik)) 1e50 else -likelihood$logLik
  }, error = function(e) 1e50)
}

# ----------------------------------------------------------------------------
# 6. Ajuste de fBM
# ----------------------------------------------------------------------------

fit_fbm <- function(series) {
  objective <- function(H) {
    safe_objective(series$y, function() cov_fbm(series$t, H))
  }
  optimum <- optimize(objective, interval = CFG$h_bounds, tol = 1e-9)
  H <- optimum$minimum
  K <- cov_fbm(series$t, H)
  likelihood <- profile_scale_loglik(series$y, K)
  sigma <- sqrt(likelihood$scale)
  new_result(
    series, "fBM", 2L, likelihood$logLik,
    list(H = H, sigma = sigma), K,
    note = "Cov=sigma^2*K_fBM(H); AIC calculado con la verosimilitud fBM."
  )
}

# ----------------------------------------------------------------------------
# 7. Ajuste de Stein d=2 y d=3
# ----------------------------------------------------------------------------

fit_stein <- function(series, dimension) {
  bounds <- log(CFG$lambda_bounds)
  objective <- function(log_lambda) {
    lambda <- exp(log_lambda)
    safe_objective(series$y, function() cov_stein(series$t, lambda, dimension))
  }
  optimum <- optimize(objective, interval = bounds, tol = 1e-10)
  lambda <- exp(optimum$minimum)
  K <- cov_stein(series$t, lambda, dimension)
  likelihood <- profile_scale_loglik(series$y, K)
  model <- paste0("Stein_S", dimension)
  new_result(
    series, model, 2L, likelihood$logLik,
    list(sigma2 = likelihood$scale, lambda = lambda), K,
    note = paste0("Cov=sigma^2*(1+lambda*|h|)^(-", dimension, "/2).")
  )
}

# ----------------------------------------------------------------------------
# 8. Ajuste de iOU (caso H=1/2 de ifOU)
# ----------------------------------------------------------------------------

fit_iou <- function(series) {
  lower <- log(CFG$beta_bounds_iou[[1L]])
  upper <- log(CFG$beta_bounds_iou[[2L]])
  grid <- seq(lower, upper, length.out = if (QUICK) 35L else 70L)

  objective <- function(log_beta) {
    beta <- exp(log_beta)
    safe_objective(series$y, function() {
      cov_ifou_observed(series, beta, H = 0.5, n_quad = CFG$n_gl_ifou)
    })
  }

  values <- vapply(grid, objective, numeric(1L))
  candidates <- order(values)[seq_len(min(5L, length(values)))]
  fits <- lapply(candidates, function(i) {
    interval <- c(grid[max(1L, i - 1L)], grid[min(length(grid), i + 1L)])
    optimize(objective, interval = interval, tol = 1e-8)
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1L), "objective"))]]
  beta <- exp(best$minimum)
  K <- cov_ifou_observed(series, beta, H = 0.5, n_quad = CFG$n_gl_ifou)
  likelihood <- profile_scale_loglik(series$y, K)

  new_result(
    series, "iOU", 2L, likelihood$logLik,
    list(beta = beta, sigma = sqrt(likelihood$scale)), K,
    note = "iOU es ifOU con H fijado en 1/2; se estiman beta y sigma."
  )
}

# ----------------------------------------------------------------------------
# 9. Ajuste ifOU y seleccion de beta* en perfiles con meseta
# ----------------------------------------------------------------------------

fit_ifou <- function(series) {
  log_beta_bounds <- log(CFG$beta_bounds_ifou)
  h_bounds <- CFG$h_bounds

  objective <- function(par) {
    beta <- exp(par[[1L]])
    H <- par[[2L]]
    if (par[[1L]] < log_beta_bounds[[1L]] || par[[1L]] > log_beta_bounds[[2L]] ||
        H < h_bounds[[1L]] || H > h_bounds[[2L]]) return(1e50)
    safe_objective(series$y, function() {
      cov_ifou_observed(series, beta, H, n_quad = CFG$n_gl_ifou)
    })
  }

  # La malla inicial evita depender de un unico punto inicial en una superficie
  # que puede ser muy plana para beta grande.
  starting_grid <- expand.grid(
    beta = c(0.3, 2, 10, 60, 200, 390),
    H = c(0.20, 0.55, 0.85, 0.94)
  )
  starting_grid$value <- mapply(function(beta, H) objective(c(log(beta), H)),
    starting_grid$beta, starting_grid$H)
  starting_grid <- starting_grid[order(starting_grid$value), , drop = FALSE]
  number_starts <- if (QUICK) 3L else 7L

  fits <- lapply(seq_len(min(number_starts, nrow(starting_grid))), function(i) {
    optim(
      par = c(log(starting_grid$beta[[i]]), starting_grid$H[[i]]),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(log_beta_bounds[[1L]], h_bounds[[1L]]),
      upper = c(log_beta_bounds[[2L]], h_bounds[[2L]]),
      control = list(maxit = if (QUICK) 120L else 350L, factr = 1e7, pgtol = 1e-7)
    )
  })
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1L), "value"))]]
  beta_raw <- exp(best$par[[1L]])
  H_raw <- best$par[[2L]]

  # Si el maximo restringido esta junto a beta=400, el parametro no queda
  # identificado de manera util. Se aplica la regla beta* de la meseta.
  plateau_case <- beta_raw >= 0.90 * CFG$beta_max_ifou

  if (!plateau_case) {
    K <- cov_ifou_observed(series, beta_raw, H_raw, n_quad = CFG$n_gl_ifou)
    likelihood <- profile_scale_loglik(series$y, K)
    return(new_result(
      series, "ifOU", 3L, likelihood$logLik,
      list(beta = beta_raw, sigma = sqrt(likelihood$scale), H = H_raw), K,
      convergence = best$convergence,
      note = "Maximo interior de (beta,H); sigma se perfila analiticamente."
    ))
  }

  # Perfil ell_p(beta): para cada beta se vuelve a maximizar en H y sigma.
  cache <- new.env(parent = emptyenv())
  profile_beta <- function(beta) {
    key <- sprintf("%.10f", beta)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    objective_H <- function(H) {
      safe_objective(series$y, function() {
        cov_ifou_observed(series, beta, H, n_quad = CFG$n_gl_ifou)
      })
    }
    H_fit <- optimize(objective_H, interval = h_bounds, tol = 2e-8)$minimum
    K <- cov_ifou_observed(series, beta, H_fit, n_quad = CFG$n_gl_ifou)
    likelihood <- profile_scale_loglik(series$y, K)
    answer <- list(
      beta = beta, H = H_fit, sigma = sqrt(likelihood$scale),
      logLik = likelihood$logLik, K = K
    )
    assign(key, answer, envir = cache)
    answer
  }

  endpoint <- profile_beta(CFG$beta_max_ifou)
  root_function <- function(beta) {
    endpoint$logLik - profile_beta(beta)$logLik - CFG$plateau_tolerance
  }

  # Se busca el PRIMER cruce, no simplemente un valor cercano a 400.
  scan_beta <- sort(unique(c(
    exp(seq(log(CFG$beta_bounds_ifou[[1L]]), log(CFG$beta_max_ifou), length.out = 90L)),
    seq(1, CFG$beta_max_ifou, length.out = 120L)
  )))
  scan_value <- vapply(scan_beta, root_function, numeric(1L))
  crossing <- which(scan_value <= 0 & c(TRUE, head(scan_value, -1L) > 0))
  crossing <- crossing[crossing > 1L]
  if (!length(crossing)) {
    stop("No se encontro el primer cruce de la meseta ifOU para Bat ",
      series$bat, " ", series$coordinate)
  }
  i <- crossing[[1L]]
  beta_star <- uniroot(
    root_function,
    interval = c(scan_beta[[i - 1L]], scan_beta[[i]]),
    tol = 5e-7
  )$root
  selected <- profile_beta(beta_star)
  plateau_gap <- endpoint$logLik - selected$logLik

  new_result(
    series, "ifOU", 3L, selected$logLik,
    list(beta = beta_star, sigma = selected$sigma, H = selected$H), selected$K,
    estimate_type = "plateau_beta_star",
    convergence = best$convergence,
    note = sprintf(
      "beta*=primer cruce; ell_p(400)-ell_p(beta*)=%.8f (objetivo=0.01).",
      plateau_gap
    )
  )
}

# ----------------------------------------------------------------------------
# 10. Ajuste zeta, incluido el caso identificable beta=0
# ----------------------------------------------------------------------------

fit_zeta <- function(series) {
  log_bounds <- log(CFG$beta_bounds_zeta)

  objective_coarse <- function(log_beta) {
    beta <- exp(log_beta)
    safe_objective(series$y, function() {
      cov_zeta(series$t, beta, n_quad = CFG$n_gl_zeta_coarse)
    })
  }
  grid <- seq(log_bounds[[1L]], log_bounds[[2L]], length.out = if (QUICK) 32L else 55L)
  values <- vapply(grid, objective_coarse, numeric(1L))
  best_index <- order(values)[seq_len(min(5L, length(values)))]
  coarse_fits <- lapply(best_index, function(i) {
    interval <- c(grid[max(1L, i - 1L)], grid[min(length(grid), i + 1L)])
    optimize(objective_coarse, interval = interval, tol = 2e-7)
  })
  coarse <- coarse_fits[[which.min(vapply(coarse_fits, `[[`, numeric(1L), "objective"))]]

  # La cuadratura fina se usa tanto para refinar beta como para el AIC final.
  objective_final <- function(log_beta) {
    beta <- exp(log_beta)
    safe_objective(series$y, function() {
      cov_zeta(series$t, beta, n_quad = CFG$n_gl_zeta_final)
    })
  }
  final_interval <- pmax(log_bounds[[1L]], coarse$minimum - 0.8)
  final_upper <- pmin(log_bounds[[2L]], coarse$minimum + 0.8)
  refined <- optimize(objective_final, interval = c(final_interval, final_upper), tol = 2e-8)
  beta <- exp(refined$minimum)
  K_interior <- cov_zeta(series$t, beta, n_quad = CFG$n_gl_zeta_final)
  likelihood_interior <- profile_scale_loglik(series$y, K_interior)
  aic_interior <- -2 * likelihood_interior$logLik + 2 * 2L

  # La frontera tiene una formula cerrada y, segun la convencion de la tabla,
  # un unico parametro estimado: sigma^2.
  K_boundary <- cov_zeta_beta_zero(series$t)
  likelihood_boundary <- profile_scale_loglik(series$y, K_boundary)
  aic_boundary <- -2 * likelihood_boundary$logLik + 2 * 1L

  if (aic_boundary < aic_interior) {
    return(new_result(
      series, "zeta", 1L, likelihood_boundary$logLik,
      list(beta = 0, sigma2 = likelihood_boundary$scale), K_boundary,
      estimate_type = "boundary_beta_zero",
      note = "Frontera beta=0: f_{sigma^2,0}(u)=sigma^2; solo sigma^2 se estima."
    ))
  }

  new_result(
    series, "zeta", 2L, likelihood_interior$logLik,
    list(beta = beta, sigma2 = likelihood_interior$scale), K_interior,
    note = "Covarianza zeta con el factor 2 del manuscrito."
  )
}

# ----------------------------------------------------------------------------
# 11. Ajuste Confluent con busqueda global y control de frontera alpha=15
# ----------------------------------------------------------------------------

fit_confluent <- function(series) {
  lower <- log(c(
    CFG$ch_eta_bounds[[1L]], CFG$ch_alpha_bounds[[1L]], CFG$ch_beta_bounds[[1L]]
  ))
  upper <- log(c(
    CFG$ch_eta_bounds[[2L]], CFG$ch_alpha_bounds[[2L]], CFG$ch_beta_bounds[[2L]]
  ))

  objective <- function(log_parameters) {
    if (any(log_parameters < lower) || any(log_parameters > upper)) return(1e50)
    p <- exp(log_parameters)
    safe_objective(series$y, function() {
      cov_confluent(series$t, eta = p[[1L]], alpha = p[[2L]], beta = p[[3L]])
    })
  }

  # DEoptim reduce el riesgo de declarar como maximo un optimo local. La
  # semilla depende de serie/coordenada y hace reproducible la busqueda.
  set.seed(SEED + 10L * series$bat + if (series$coordinate == "Longitude") 1L else 0L)
  global <- DEoptim::DEoptim(
    fn = objective,
    lower = lower,
    upper = upper,
    control = DEoptim::DEoptim.control(
      NP = CFG$de_population,
      itermax = CFG$de_iter,
      trace = FALSE,
      reltol = 1e-9,
      steptol = if (QUICK) 10L else 35L
    )
  )
  global_start <- global$optim$bestmem

  generic_starts <- rbind(
    global_start,
    log(c(0.05, 0.5, 10)),
    log(c(0.1, 5, 100)),
    log(c(0.5, 1, 1000)),
    log(c(1, 8, 100)),
    log(c(2, 3, 1000)),
    log(c(0.1, 10, 10000))
  )

  local_fits <- lapply(seq_len(nrow(generic_starts)), function(i) {
    fit_lbfgs <- optim(
      generic_starts[i, ], objective, method = "L-BFGS-B",
      lower = lower, upper = upper,
      control = list(maxit = if (QUICK) 250L else 1200L, factr = 1e5, pgtol = 1e-8)
    )
    # Nelder-Mead es util cerca de una frontera cuando las diferencias finitas
    # de L-BFGS-B resultan demasiado pequenas para mover los parametros.
    fit_nm <- optim(
      fit_lbfgs$par, objective, method = "Nelder-Mead",
      control = list(maxit = if (QUICK) 300L else 1600L, reltol = 1e-11)
    )
    if (fit_nm$value < fit_lbfgs$value) fit_nm else fit_lbfgs
  })
  best <- local_fits[[which.min(vapply(local_fits, `[[`, numeric(1L), "value"))]]

  # Auditoria explicita de alpha=15: no basta con esperar que un optimizador
  # continuo llegue numericamente exactamente a la cota.
  alpha_max <- CFG$ch_alpha_bounds[[2L]]
  lower_boundary <- lower[c(1L, 3L)]
  upper_boundary <- upper[c(1L, 3L)]
  objective_boundary <- function(log_eta_beta) {
    if (any(log_eta_beta < lower_boundary) || any(log_eta_beta > upper_boundary)) return(1e50)
    eta <- exp(log_eta_beta[[1L]])
    beta <- exp(log_eta_beta[[2L]])
    safe_objective(series$y, function() {
      cov_confluent(series$t, eta = eta, alpha = alpha_max, beta = beta)
    })
  }
  boundary_start <- best$par[c(1L, 3L)]
  boundary_lbfgs <- optim(
    boundary_start, objective_boundary, method = "L-BFGS-B",
    lower = lower_boundary, upper = upper_boundary,
    control = list(maxit = if (QUICK) 250L else 1000L, factr = 1e5, pgtol = 1e-9)
  )
  boundary_nm <- optim(
    boundary_lbfgs$par, objective_boundary, method = "Nelder-Mead",
    control = list(maxit = if (QUICK) 300L else 1600L, reltol = 1e-12)
  )
  boundary_best <- if (boundary_nm$value < boundary_lbfgs$value) boundary_nm else boundary_lbfgs

  if (boundary_best$value <= best$value + 1e-7) {
    p <- c(exp(boundary_best$par[[1L]]), alpha_max, exp(boundary_best$par[[2L]]))
    estimate_type <- "boundary_alpha_max"
    convergence <- boundary_best$convergence
    note <- paste0(
      "alpha=15 es la cota impuesta; solucion de frontera, no MLE interior. ",
      "El AIC depende de alpha_max=15."
    )
  } else {
    p <- exp(best$par)
    estimate_type <- "interior"
    convergence <- best$convergence
    note <- "Maximo interior bajo los limites numericos declarados."
  }

  K <- cov_confluent(series$t, eta = p[[1L]], alpha = p[[2L]], beta = p[[3L]])
  likelihood <- profile_scale_loglik(series$y, K)
  new_result(
    series, "Confluent", 4L, likelihood$logLik,
    list(sigma2 = likelihood$scale, eta = p[[1L]], alpha = p[[2L]], beta = p[[3L]]),
    K, estimate_type = estimate_type, convergence = convergence, note = note
  )
}

# ----------------------------------------------------------------------------
# 12. Pruebas unitarias de formulas y convenciones
# ----------------------------------------------------------------------------

run_unit_tests <- function(series_list) {
  test_t <- c(0.5, 1.0, 2.0)

  K_fbm <- cov_fbm(test_t, 0.7)
  stopifnot(max(abs(diag(K_fbm) - test_t^1.4)) < 1e-12)
  stopifnot(max(abs(K_fbm - t(K_fbm))) < 1e-12)

  K_s2 <- cov_stein(test_t, 0.1, 2L)
  K_s3 <- cov_stein(test_t, 0.1, 3L)
  stopifnot(max(abs(diag(K_s2) - 1)) < 1e-14)
  stopifnot(max(abs(diag(K_s3) - 1)) < 1e-14)

  K_ch <- cov_confluent(test_t, eta = 0.8, alpha = 1.3, beta = 20)
  stopifnot(max(abs(diag(K_ch) - 1)) < 1e-12)

  K_zeta0 <- cov_zeta_beta_zero(test_t)
  stopifnot(max(abs(K_zeta0 - t(K_zeta0))) < 1e-12)
  stopifnot(!is.null(tryCatch(chol(K_zeta0), error = function(e) NULL)))

  K_ifou <- cov_ifou_grid(4L, 0.5, beta = 2, H = 0.7, n_quad = max(32L, CFG$n_gl_ifou))
  stopifnot(max(abs(K_ifou - t(K_ifou))) < 1e-10)
  stopifnot(!is.null(tryCatch(chol(K_ifou), error = function(e) NULL)))

  # Todos los modelos deben recibir n-1 valores y ninguno puede recibir el
  # cero artificial inicial.
  stopifnot(all(vapply(series_list, function(s) length(s$y) == s$n_original - 1L, logical(1L))))
  stopifnot(all(vapply(series_list, function(s) length(s$t) == length(s$y), logical(1L))))

  TRUE
}

# ----------------------------------------------------------------------------
# 13. Los 70 AIC de referencia de la auditoria independiente
# ----------------------------------------------------------------------------

reference_aic <- function() {
  # Estas cifras se obtuvieron con una segunda implementacion numerica
  # independiente, Cholesky exacto y las mismas formulas/muestra. No se usan
  # para ajustar: sirven unicamente como prueba posterior de reproduccion.
  read.csv(text = "bat,coordinate,model,AIC_reference
1,Latitude,zeta,-263.701076537
1,Latitude,ifOU,-282.661438446
1,Latitude,iOU,-253.478751996
1,Latitude,fBM,-284.704242512
1,Latitude,Confluent,-281.393186160
1,Latitude,Stein_S2,-246.901251560
1,Latitude,Stein_S3,-247.010015560
1,Longitude,zeta,-231.464821427
1,Longitude,ifOU,-245.337382764
1,Longitude,iOU,-230.073170623
1,Longitude,fBM,-247.362608553
1,Longitude,Confluent,-238.114206470
1,Longitude,Stein_S2,-220.468538070
1,Longitude,Stein_S3,-220.586619390
2,Latitude,zeta,-277.520167903
2,Latitude,ifOU,-284.920943657
2,Latitude,iOU,-272.593669633
2,Latitude,fBM,-286.942742246
2,Latitude,Confluent,-272.878065520
2,Latitude,Stein_S2,-261.650478710
2,Latitude,Stein_S3,-261.766212940
2,Longitude,zeta,-252.475202362
2,Longitude,ifOU,-266.697225571
2,Longitude,iOU,-257.457582029
2,Longitude,fBM,-268.720143183
2,Longitude,Confluent,-254.920928840
2,Longitude,Stein_S2,-247.873063750
2,Longitude,Stein_S3,-247.985421880
3,Latitude,zeta,-583.553861991
3,Latitude,ifOU,-586.563527330
3,Latitude,iOU,-585.752890945
3,Latitude,fBM,-580.228267988
3,Latitude,Confluent,-568.230233040
3,Latitude,Stein_S2,-536.697890500
3,Latitude,Stein_S3,-536.812727030
3,Longitude,zeta,-564.410087985
3,Longitude,ifOU,-559.912475124
3,Longitude,iOU,-533.218420657
3,Longitude,fBM,-561.926479655
3,Longitude,Confluent,-545.173007740
3,Longitude,Stein_S2,-467.699620820
3,Longitude,Stein_S3,-467.815707420
4,Latitude,zeta,-2580.945958430
4,Latitude,ifOU,-2534.840871530
4,Latitude,iOU,-2473.366688040
4,Latitude,fBM,-2536.626687130
4,Latitude,Confluent,-2519.239253260
4,Latitude,Stein_S2,-2298.239587900
4,Latitude,Stein_S3,-2298.356694130
4,Longitude,zeta,-2244.660254580
4,Longitude,ifOU,-2259.060635570
4,Longitude,iOU,-2146.934895560
4,Longitude,fBM,-2261.089635750
4,Longitude,Confluent,-2240.079860790
4,Longitude,Stein_S2,-1913.007319920
4,Longitude,Stein_S3,-1913.124974590
5,Latitude,zeta,-1648.925360670
5,Latitude,ifOU,-1649.603825190
5,Latitude,iOU,-1590.779785110
5,Latitude,fBM,-1651.633442950
5,Latitude,Confluent,-1640.839433270
5,Latitude,Stein_S2,-1488.965544750
5,Latitude,Stein_S3,-1489.182217480
5,Longitude,zeta,-1585.845254060
5,Longitude,ifOU,-1581.289654290
5,Longitude,iOU,-1535.582065500
5,Longitude,fBM,-1583.310030550
5,Longitude,Confluent,-1573.800992070
5,Longitude,Stein_S2,-1392.852906620
5,Longitude,Stein_S3,-1393.075461410")
}

reference_tolerance <- function(model) {
  # Confluent usa optimizacion global y una funcion especial; se permite una
  # tolerancia algo mayor. Para zeta, la tolerancia cubre la cuadratura finita.
  unname(c(
    zeta = 0.03, ifOU = 0.03, iOU = 0.02, fBM = 0.005,
    Confluent = 0.15, Stein_S2 = 0.005, Stein_S3 = 0.005
  )[model])
}

# ----------------------------------------------------------------------------
# 14. Ejecucion completa y escritura de resultados
# ----------------------------------------------------------------------------

fit_one_model <- function(series, model) {
  switch(model,
    zeta = fit_zeta(series),
    ifOU = fit_ifou(series),
    iOU = fit_iou(series),
    fBM = fit_fbm(series),
    Confluent = fit_confluent(series),
    Stein_S2 = fit_stein(series, 2L),
    Stein_S3 = fit_stein(series, 3L),
    stop("Modelo no implementado: ", model)
  )
}

error_result <- function(series, model, message) {
  data.frame(
    bat = series$bat, coordinate = series$coordinate,
    n_original = series$n_original, n_used = series$n_used,
    model = model, k = NA_integer_, logLik = NA_real_, AIC = NA_real_,
    beta = NA_real_, sigma = NA_real_, sigma2 = NA_real_, H = NA_real_,
    eta = NA_real_, alpha = NA_real_, lambda = NA_real_,
    estimate_type = "ERROR", convergence = 999L, covariance_pd = FALSE,
    min_chol_diagonal = NA_real_, reciprocal_condition = NA_real_,
    note = message, check.names = FALSE
  )
}

main <- function() {
  cat("\n=== Auditoria de siete modelos de covarianza ===\n")
  cat("Datos :", normalizePath(DATA_PATH), "\n")
  cat("Salida:", normalizePath(OUT_DIR, mustWork = FALSE), "\n")
  cat("Modo  :", if (QUICK) "QUICK (no publicacion)" else "AUDITORIA", "\n")
  cat("Modelos:", paste(MODELS_TO_RUN, collapse = ", "), "\n\n")

  series_list <- read_bat_series(DATA_PATH, CFG$delta)
  run_unit_tests(series_list)
  cat("Pruebas unitarias de formulas y muestra: OK\n\n")

  rows <- list()
  cursor <- 1L
  for (series in series_list) {
    for (model in MODELS_TO_RUN) {
      label <- sprintf("Bat %d %-9s %-10s", series$bat, series$coordinate, model)
      cat("Ajustando", label, "... ")
      start_time <- proc.time()[[3L]]
      fit <- tryCatch(
        fit_one_model(series, model),
        error = function(e) error_result(series, model, conditionMessage(e))
      )
      elapsed <- proc.time()[[3L]] - start_time
      if (is.finite(fit$AIC[[1L]])) {
        cat(sprintf("AIC = %.6f  (%.1f s)\n", fit$AIC[[1L]], elapsed))
      } else {
        cat(sprintf("ERROR: %s  (%.1f s)\n", fit$note[[1L]], elapsed))
      }
      rows[[cursor]] <- fit
      cursor <- cursor + 1L
    }
  }
  results <- do.call(rbind, rows)
  results <- results[order(results$bat, results$coordinate, match(results$model, ALL_MODELS)), ]

  # Ganador por AIC dentro de los modelos efectivamente ejecutados.
  valid <- results[is.finite(results$AIC), , drop = FALSE]
  split_series <- split(valid, interaction(valid$bat, valid$coordinate, drop = TRUE))
  winners <- do.call(rbind, lapply(split_series, function(d) d[which.min(d$AIC), , drop = FALSE]))
  winners <- winners[order(winners$bat, winners$coordinate), ]

  # Comparacion independiente de los 70 valores. Si se ejecuta un subconjunto,
  # la tabla contiene solamente las filas disponibles.
  reference <- reference_aic()
  comparison <- merge(
    results[, c("bat", "coordinate", "model", "AIC", "estimate_type", "note")],
    reference,
    by = c("bat", "coordinate", "model"),
    all.x = TRUE,
    sort = FALSE
  )
  comparison$absolute_difference <- abs(comparison$AIC - comparison$AIC_reference)
  comparison$tolerance <- vapply(comparison$model, reference_tolerance, numeric(1L))
  comparison$passes_reference <- comparison$absolute_difference <= comparison$tolerance
  comparison <- comparison[order(
    comparison$bat, comparison$coordinate, match(comparison$model, ALL_MODELS)
  ), ]

  write.csv(results, file.path(OUT_DIR, "resultados_todos_los_modelos.csv"), row.names = FALSE)
  write.csv(winners, file.path(OUT_DIR, "ganadores_por_serie.csv"), row.names = FALSE)
  write.csv(comparison, file.path(OUT_DIR, "comparacion_con_referencia.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

  failed_reference <- comparison[
    is.na(comparison$passes_reference) | !comparison$passes_reference, , drop = FALSE
  ]
  boundary_ch <- results[
    results$model == "Confluent" & results$estimate_type == "boundary_alpha_max", , drop = FALSE
  ]

  summary_lines <- c(
    "AUDITORIA DE LA TABLA DE MODELOS",
    "================================",
    "",
    paste0("Archivo de datos: ", normalizePath(DATA_PATH)),
    paste0("Muestra comun: y_i=X(t_i)-X(0), i=1,...,n-1."),
    paste0("Cero inicial excluido para las siete familias."),
    paste0("Modo quick: ", QUICK),
    "",
    "Ganadores por AIC:",
    paste0(
      "  Bat ", winners$bat, " ", winners$coordinate, ": ", winners$model,
      " (AIC=", sprintf("%.6f", winners$AIC), ")"
    ),
    "",
    paste0("Filas comparadas con referencia: ", nrow(comparison)),
    paste0("Filas dentro de tolerancia: ", sum(comparison$passes_reference, na.rm = TRUE)),
    paste0("Filas fuera de tolerancia: ", nrow(failed_reference)),
    "",
    paste0("Ajustes Confluent en alpha_max=15: ", nrow(boundary_ch)),
    if (nrow(boundary_ch)) paste0(
      "  Bat ", boundary_ch$bat, " ", boundary_ch$coordinate,
      ": alpha=", boundary_ch$alpha
    ) else "  Ninguno.",
    "",
    "Nota: una solucion Confluent en alpha=15 es una solucion de frontera",
    "dependiente de la cota numerica; no debe llamarse MLE interior."
  )
  writeLines(summary_lines, file.path(OUT_DIR, "resumen_verificacion.txt"), useBytes = TRUE)

  cat("\n=== Ganadores ===\n")
  print(winners[, c("bat", "coordinate", "model", "AIC")], row.names = FALSE)
  cat("\nComparaciones fuera de tolerancia:", nrow(failed_reference), "\n")
  if (nrow(failed_reference)) {
    print(failed_reference[, c(
      "bat", "coordinate", "model", "AIC", "AIC_reference",
      "absolute_difference", "tolerance"
    )], row.names = FALSE)
  }
  cat("\nResultados escritos en:", normalizePath(OUT_DIR), "\n")
  invisible(results)
}

if (sys.nframe() == 0L) main()
