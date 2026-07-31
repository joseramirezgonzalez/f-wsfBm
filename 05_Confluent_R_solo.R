#!/usr/bin/env Rscript

# =============================================================================
# MODELO A: covarianza hipergeometrica confluyente (CH)
# Ajuste de los 10 casos usando EXCLUSIVAMENTE R.
#
# Cov(A(t), A(s)) = sigma2 * Gamma(eta + alpha) / Gamma(eta) *
#   U(alpha, 1 - eta, eta * (abs(t-s) / beta)^2)
#
# GPBayes usa
#   Gamma(nu + tail) / Gamma(nu) *
#   U(tail, 1 - nu, (h / range)^2).
# Por tanto, la correspondencia EXACTA es:
#   nu   = eta
#   tail = alpha
#   range = beta / sqrt(eta)
#
# PARAMETROS EDITABLES (los mismos limites se usan en los 10 casos)
# =============================================================================

ETA_BOUNDS   <- c(0.005, 5)
ALPHA_BOUNDS <- c(0.01, 15)
BETA_BOUNDS  <- c(0.001, 50000)

# Controles de la busqueda global. Para una comprobacion mas intensa se puede
# aumentar DE_ITERMAX o GLOBAL_RUNS. Los valores siguientes ya realizan dos
# busquedas globales independientes y refinaciones locales/acotadas.
GLOBAL_RUNS <- 2L
DE_NP       <- 60L
DE_ITERMAX  <- 100L
LOCAL_STARTS <- 6L
RANDOM_SEED <- 20260731L

# Estabilizacion NUMERICA fija; no es un nugget estimado y no aumenta k.
JITTER_SEQUENCE <- c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5)

# Si faltan paquetes, el script intenta instalarlos desde CRAN.
AUTO_INSTALL_PACKAGES <- TRUE

# =============================================================================

required_packages <- c("readxl", "GPBayes", "DEoptim")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]

if (length(missing_packages)) {
  if (!AUTO_INSTALL_PACKAGES) {
    stop(
      "Faltan paquetes de R: ", paste(missing_packages, collapse = ", "),
      ". Instala con install.packages(c(",
      paste(sprintf("'%s'", missing_packages), collapse = ", "), "))."
    )
  }
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(still_missing)) {
  stop(
    "No se pudieron cargar los paquetes: ",
    paste(still_missing, collapse = ", "),
    ". Revisa los mensajes de instalacion anteriores."
  )
}

options(digits = 12)

args <- commandArgs(trailingOnly = TRUE)
DATA_FILE <- if (length(args)) args[[1L]] else file.choose()
DATA_FILE <- normalizePath(DATA_FILE, mustWork = TRUE)
OUTPUT_FILE <- file.path(dirname(DATA_FILE), "resultados_Confluent_GPBayes.csv")


# Los valores actuales de la tabla se incluyen UNICAMENTE como puntos iniciales.
# DEoptim tambien genera una poblacion aleatoria en todo el dominio permitido;
# ningun AIC de referencia interviene en la funcion que se maximiza.
table_starts <- data.frame(
  bat = rep(1:5, each = 2L),
  coordinate = rep(c("Latitude", "Longitude"), 5L),
  eta = c(0.791, 0.737, 0.690, 0.670, 0.955,
          0.893, 0.868, 0.850, 0.899, 0.971),
  alpha = c(8.770, 9.476, 15, 15, 1.319,
            9.107, 9.132, 9.148, 8.954, 9.783),
  beta = c(11300.714, 8278.885, 7207.818, 5923.863, 135.603,
           1974.214, 3336.100, 16400.040, 420.182, 311.460),
  stringsAsFactors = FALSE
)


read_ten_series <- function(path) {
  raw <- as.data.frame(
    readxl::read_excel(path, sheet = "in", col_names = FALSE)
  )

  # Archivo 2D: tres columnas por murcielago.
  # Archivo con altitud: cuatro columnas por murcielago.
  width <- if (ncol(raw) >= 20L) 4L else 3L
  out <- vector("list", 10L)
  k <- 0L

  for (bat in 1:5) {
    first_col <- width * (bat - 1L) + 1L
    if (first_col + 2L > ncol(raw)) {
      stop("El Excel no contiene las tres columnas esperadas para Bat #", bat, ".")
    }

    d <- raw[, first_col:(first_col + 2L), drop = FALSE]
    names(d) <- c("time", "Longitude", "Latitude")
    d[] <- lapply(d, function(x) suppressWarnings(as.numeric(x)))
    d <- d[complete.cases(d), , drop = FALSE]

    if (nrow(d) < 3L) {
      stop("Hay menos de tres observaciones completas para Bat #", bat, ".")
    }

    for (coordinate in c("Latitude", "Longitude")) {
      k <- k + 1L
      x <- d[[coordinate]]
      out[[k]] <- list(
        bat = bat,
        coordinate = coordinate,
        t = d$time[-1L],
        y = x[-1L] - x[1L]
      )
    }
  }

  out
}


make_context <- function(series) {
  lag_matrix <- abs(outer(series$t, series$t, "-"))
  lag_values <- sort(unique(as.numeric(lag_matrix)))

  list(
    bat = series$bat,
    coordinate = series$coordinate,
    y = as.numeric(series$y),
    n = length(series$y),
    lag_values = lag_values,
    lag_index = match(as.numeric(lag_matrix), lag_values),
    cache = new.env(hash = TRUE, parent = emptyenv())
  )
}


profile_scale_loglik <- function(y, K) {
  K <- (K + t(K)) / 2
  n <- length(y)
  base <- max(mean(diag(K)), 1e-12)

  for (jitter in JITTER_SEQUENCE) {
    R <- tryCatch(
      chol(K + diag(base * jitter, n)),
      error = function(e) NULL
    )

    if (is.null(R)) next

    z <- tryCatch(
      backsolve(R, y, transpose = TRUE),
      error = function(e) NULL
    )
    if (is.null(z)) next

    sigma2 <- sum(z^2) / n
    if (!is.finite(sigma2) || sigma2 <= 0) next

    logLik <- -0.5 * (
      n * (log(2 * pi) + 1 + log(sigma2)) +
        2 * sum(log(diag(R)))
    )

    if (is.finite(logLik)) {
      return(list(ok = TRUE, logLik = logLik, sigma2 = sigma2, jitter = jitter))
    }
  }

  list(ok = FALSE, logLik = -Inf, sigma2 = NA_real_, jitter = NA_real_)
}


covariance_from_log_parameters <- function(log_parameters, context) {
  natural <- exp(log_parameters)
  eta <- natural[1L]
  alpha <- natural[2L]
  beta <- natural[3L]

  # Transformacion exacta entre la parametrizacion del articulo y GPBayes.
  gp_range <- beta / sqrt(eta)

  rho <- suppressWarnings(
    GPBayes::CH(
      matrix(context$lag_values, ncol = 1L),
      range = gp_range,
      tail = alpha,
      nu = eta
    )
  )
  rho <- as.numeric(rho)

  if (length(rho) != length(context$lag_values) ||
      any(!is.finite(rho)) || any(rho < 0)) {
    stop("GPBayes::CH devolvio una covarianza no evaluable.")
  }

  # En h=0 la correlacion teorica es exactamente uno.
  rho[context$lag_values == 0] <- 1
  K <- matrix(rho[context$lag_index], nrow = context$n, ncol = context$n)
  (K + t(K)) / 2
}


make_objective <- function(context) {
  function(log_parameters) {
    key <- paste(sprintf("%.13g", log_parameters), collapse = "|")
    if (exists(key, envir = context$cache, inherits = FALSE)) {
      return(get(key, envir = context$cache, inherits = FALSE))
    }

    value <- tryCatch({
      K <- covariance_from_log_parameters(log_parameters, context)
      fit <- profile_scale_loglik(context$y, K)
      if (fit$ok) -fit$logLik else 1e100
    }, error = function(e) 1e100)

    if (!is.finite(value)) value <- 1e100
    assign(key, value, envir = context$cache)
    value
  }
}


clip_natural <- function(x) {
  pmin(
    pmax(x, c(ETA_BOUNDS[1L], ALPHA_BOUNDS[1L], BETA_BOUNDS[1L])),
    c(ETA_BOUNDS[2L], ALPHA_BOUNDS[2L], BETA_BOUNDS[2L])
  )
}


starting_values <- function(bat, coordinate) {
  current <- table_starts[
    table_starts$bat == bat & table_starts$coordinate == coordinate,
    c("eta", "alpha", "beta"),
    drop = FALSE
  ]

  generic <- rbind(
    c(0.05, 0.5, 10),
    c(0.1, 5, 100),
    c(0.5, 1, 1000),
    c(1, 8, 100),
    c(2, 3, 1000),
    c(0.1, 10, 10000),
    c(0.8, ALPHA_BOUNDS[2L], 500),
    c(0.8, ALPHA_BOUNDS[2L], 10000)
  )

  if (nrow(current)) {
    current <- as.numeric(current[1L, ])
    around_current <- rbind(
      current,
      current * c(0.8, 1, 0.5),
      current * c(1.2, 1, 2),
      c(current[1L], ALPHA_BOUNDS[2L], current[3L])
    )
    generic <- rbind(around_current, generic)
  }

  generic <- t(apply(generic, 1L, clip_natural))
  unique(log(generic))
}


make_initial_population <- function(np, lower, upper, seeds) {
  population <- sapply(
    seq_along(lower),
    function(j) runif(np, min = lower[j], max = upper[j])
  )
  population <- matrix(population, nrow = np, ncol = length(lower))

  n_seed <- min(nrow(seeds), np)
  if (n_seed > 0L) population[seq_len(n_seed), ] <- seeds[seq_len(n_seed), ]
  population
}


unique_rows <- function(x, digits = 10L) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  keys <- apply(round(x, digits), 1L, paste, collapse = "|")
  x[!duplicated(keys), , drop = FALSE]
}


safe_nlminb <- function(start, objective, lower, upper) {
  start <- pmin(pmax(start, lower), upper)
  tryCatch(
    nlminb(
      start = start,
      objective = objective,
      lower = lower,
      upper = upper,
      control = list(
        eval.max = 2500L,
        iter.max = 1200L,
        rel.tol = 1e-11,
        x.tol = 1e-9
      )
    ),
    error = function(e) list(
      par = start,
      objective = objective(start),
      convergence = 999L,
      message = conditionMessage(e)
    )
  )
}


refine_boundary_face <- function(best, fixed_index, fixed_value,
                                 objective, lower, upper) {
  free_index <- setdiff(1:3, fixed_index)
  start_free <- best[free_index]

  face_objective <- function(free_parameters) {
    full <- best
    full[fixed_index] <- fixed_value
    full[free_index] <- free_parameters
    objective(full)
  }

  fit <- safe_nlminb(
    start = start_free,
    objective = face_objective,
    lower = lower[free_index],
    upper = upper[free_index]
  )

  full <- best
  full[fixed_index] <- fixed_value
  full[free_index] <- fit$par
  list(par = full, objective = objective(full))
}


fit_one_series <- function(series) {
  context <- make_context(series)
  objective <- make_objective(context)

  lower <- log(c(ETA_BOUNDS[1L], ALPHA_BOUNDS[1L], BETA_BOUNDS[1L]))
  upper <- log(c(ETA_BOUNDS[2L], ALPHA_BOUNDS[2L], BETA_BOUNDS[2L]))
  seeds <- starting_values(series$bat, series$coordinate)

  de_fits <- vector("list", GLOBAL_RUNS)
  candidates <- seeds

  for (run in seq_len(GLOBAL_RUNS)) {
    set.seed(
      RANDOM_SEED + 1000L * series$bat +
        10L * match(series$coordinate, c("Latitude", "Longitude")) + run
    )

    initial_population <- make_initial_population(DE_NP, lower, upper, seeds)
    de_fits[[run]] <- DEoptim::DEoptim(
      fn = objective,
      lower = lower,
      upper = upper,
      control = DEoptim::DEoptim.control(
        NP = DE_NP,
        itermax = DE_ITERMAX,
        strategy = 6L,
        p = 0.2,
        CR = 0.9,
        F = 0.8,
        trace = FALSE,
        initialpop = initial_population,
        reltol = 1e-9,
        steptol = DE_ITERMAX
      )
    )

    candidates <- rbind(
      candidates,
      de_fits[[run]]$optim$bestmem,
      de_fits[[run]]$member$pop
    )
  }

  candidates <- unique_rows(candidates)
  candidate_values <- apply(candidates, 1L, objective)
  candidates <- candidates[order(candidate_values), , drop = FALSE]
  local_starts <- candidates[seq_len(min(LOCAL_STARTS, nrow(candidates))), , drop = FALSE]

  local_fits <- lapply(
    seq_len(nrow(local_starts)),
    function(i) safe_nlminb(local_starts[i, ], objective, lower, upper)
  )

  local_parameters <- do.call(rbind, lapply(local_fits, `[[`, "par"))
  local_values <- vapply(local_fits, `[[`, numeric(1L), "objective")
  best_index <- which.min(local_values)
  best <- as.numeric(local_parameters[best_index, ])
  best_value <- objective(best)

  # Comprobacion explicita de las seis caras del dominio. Esto es importante
  # porque alpha=15 es una solucion de frontera en varios de los diez casos.
  face_fits <- list()
  face_counter <- 0L
  for (fixed_index in 1:3) {
    for (fixed_value in c(lower[fixed_index], upper[fixed_index])) {
      face_counter <- face_counter + 1L
      face_fits[[face_counter]] <- refine_boundary_face(
        best, fixed_index, fixed_value, objective, lower, upper
      )
    }
  }

  for (face in face_fits) {
    if (is.finite(face$objective) && face$objective < best_value) {
      best <- face$par
      best_value <- face$objective
    }
  }

  natural <- exp(best)
  natural <- clip_natural(natural)
  best <- log(natural)
  K <- covariance_from_log_parameters(best, context)
  profile <- profile_scale_loglik(context$y, K)

  if (!profile$ok) {
    stop(
      "No se pudo evaluar la solucion final de Bat #", series$bat,
      " - ", series$coordinate, "."
    )
  }

  data.frame(
    bat = series$bat,
    coordinate = series$coordinate,
    sigma2 = profile$sigma2,
    eta = natural[1L],
    alpha = natural[2L],
    beta = natural[3L],
    AIC = -2 * profile$logLik + 2 * 4L,
    stringsAsFactors = FALSE
  )
}


series <- read_ten_series(DATA_FILE)
results_list <- vector("list", length(series))

for (i in seq_along(series)) {
  s <- series[[i]]
  cat(sprintf(
    "\n[%d/10] Ajustando A: Bat #%d - %s ...\n",
    i, s$bat, s$coordinate
  ))

  results_list[[i]] <- fit_one_series(s)
  partial <- do.call(rbind, results_list[seq_len(i)])
  partial <- partial[
    order(partial$bat, match(partial$coordinate, c("Latitude", "Longitude"))),
    , drop = FALSE
  ]
  write.csv(partial, OUTPUT_FILE, row.names = FALSE)
  print(results_list[[i]], row.names = FALSE, digits = 10)
}

results <- do.call(rbind, results_list)
results <- results[
  order(results$bat, match(results$coordinate, c("Latitude", "Longitude"))),
  , drop = FALSE
]

cat("\nRESULTADOS FINALES DEL MODELO A (k = 4)\n")
print(results, row.names = FALSE, digits = 10)
cat("\nArchivo escrito en:\n", OUTPUT_FILE, "\n", sep = "")

