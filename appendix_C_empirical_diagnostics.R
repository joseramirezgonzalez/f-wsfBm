#!/usr/bin/env Rscript

# =============================================================================
# APENDICE C: DIAGNOSTICOS EMPIRICOS DE LAS TRAYECTORIAS DE MURCIELAGOS
# =============================================================================
#
# Este archivo reproduce, desde el libro bat_set.xlsx, todos los calculos y
# todas las figuras del Apendice C del material suplementario:
#
#   1. medias moviles (ventana centrada de 10 observaciones);
#   2. varianzas moviles (ventana centrada de 10 observaciones);
#   3. ACF de las series originales;
#   4. ACF de las primeras diferencias;
#   5. ACF de las segundas diferencias;
#   6. pruebas ADF y KPSS para niveles y primeras diferencias (Tabla C.1);
#   7. detrended fluctuation analysis, DFA (Figura C.9 y Tabla C.2);
#   8. resumen grafico de los exponentes DFA (Figura C.10).
#
# El codigo no ajusta los modelos de la Tabla 1 del manuscrito principal. Su
# objetivo exclusivo es reproducir el Apendice C. El ajuste de fBM y de los
# demas modelos se encuentra en el archivo de auditoria de modelos entregado
# por separado.
#
# EJECUCION
# ---------
#
#   Rscript appendix_C_empirical_diagnostics.R \
#     --data=bat_set.xlsx \
#     --out=resultados_apendice_C
#
# DEPENDENCIA
# -----------
#
# Solo se necesita el paquete "readxl" para leer el libro de Excel. Tambien se
# acepta una copia CSV de las mismas 15 columnas, en cuyo caso no se requieren
# paquetes. Todos los calculos estadisticos y todas las figuras se construyen
# con R base. Esto evita que cambios entre versiones alteren las cifras.
#
# CONVENCIONES IMPORTANTES
# ------------------------
#
# * Las coordenadas se analizan en el orden en que aparecen en el libro:
#   tiempo, longitud y latitud, en cinco bloques consecutivos de tres columnas.
# * Las pruebas de estacionariedad se aplican a las observaciones ordenadas.
# * La ADF incluye constante y tendencia lineal, como tseries::adf.test.
# * La KPSS contrasta estacionariedad alrededor de un nivel y utiliza la regla
#   corta de Newey--West de tseries::kpss.test (lshort = TRUE).
# * Los valores p de ADF se acotan a [0.01, 0.99] y los de KPSS a [0.01, 0.10],
#   tal como hace tseries. Por eso 0.0100 y 0.1000 son limites de reporte.
# * DFA usa el perfil acumulado, segmentos no superpuestos que comienzan en la
#   primera observacion y una tendencia lineal dentro de cada segmento.
# * Los seis tamanos de ventana mas grandes del grid de cada murcielago forman
#   el rango de regresion mostrado en rojo en la Figura C.9.
#
# INTERPRETACION
# --------------
#
# La estacionariedad de una primera diferencia no significa independencia:
# una sucesion estacionaria puede conservar autocorrelacion. Por ello, que la
# ACF de las primeras diferencias tenga barras significativas no contradice,
# por si solo, una descripcion con incrementos estacionarios como fBM.
# Asimismo, las pruebas ADF/KPSS son diagnosticos finito-muestrales y no
# sustituyen la comparacion de verosimilitudes de los modelos completos.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)


# -----------------------------------------------------------------------------
# 1. Argumentos de linea de comandos y dependencias
# -----------------------------------------------------------------------------

parse_arguments <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    data = if (nzchar(Sys.getenv("APPENDIX_C_DATA"))) {
      Sys.getenv("APPENDIX_C_DATA")
    } else {
      "bat_set.xlsx"
    },
    out = if (nzchar(Sys.getenv("APPENDIX_C_OUT"))) {
      Sys.getenv("APPENDIX_C_OUT")
    } else {
      "resultados_apendice_C"
    },
    skip_figures = identical(Sys.getenv("APPENDIX_C_SKIP_FIGURES"), "1")
  )
  for (argument in arguments) {
    if (grepl("^--data=", argument)) {
      defaults$data <- sub("^--data=", "", argument)
    } else if (grepl("^--out=", argument)) {
      defaults$out <- sub("^--out=", "", argument)
    } else if (argument == "--skip-figures") {
      # Opcion util para una comprobacion numerica rapida en entornos sin
      # dispositivo grafico. La ejecucion normal NO debe usar esta opcion.
      defaults$skip_figures <- TRUE
    } else if (argument %in% c("--help", "-h")) {
      cat(
        "Uso:\n",
        "  Rscript appendix_C_empirical_diagnostics.R ",
        "--data=bat_set.xlsx --out=resultados_apendice_C ",
        "[--skip-figures]\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Argumento no reconocido: ", argument, call. = FALSE)
    }
  }
  defaults
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Falta el paquete '", package, "'. Instalelo con ",
      "install.packages(\"", package, "\").",
      call. = FALSE
    )
  }
}

ARGS <- parse_arguments()

DATA_PATH <- normalizePath(ARGS$data, mustWork = TRUE)
if (tolower(tools::file_ext(DATA_PATH)) %in% c("xlsx", "xls")) {
  require_package("readxl")
}
OUT_DIR <- ARGS$out
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
DATA_DIR <- file.path(OUT_DIR, "figure_data")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. Lectura y auditoria estructural del libro de datos
# -----------------------------------------------------------------------------

read_bat_workbook <- function(path) {
  # El libro contiene una hoja llamada "in". Cada murcielago ocupa tres
  # columnas consecutivas: tiempo, longitud y latitud. Las filas vacias al
  # final de cada bloque se eliminan solo para ese murcielago.
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("xlsx", "xls")) {
    raw <- as.data.frame(
      readxl::read_excel(path, sheet = "in", col_names = FALSE),
      check.names = FALSE
    )
  } else if (extension == "csv") {
    raw <- utils::read.csv(
      path, header = FALSE, check.names = FALSE,
      na.strings = c("", "NA"), strip.white = TRUE
    )
  } else {
    stop("El archivo de datos debe ser .xlsx, .xls o .csv.", call. = FALSE)
  }

  if (ncol(raw) < 15L) {
    stop("La hoja 'in' debe contener al menos 15 columnas.", call. = FALSE)
  }

  bats <- vector("list", 5L)
  for (bat in seq_len(5L)) {
    columns <- (3L * (bat - 1L) + 1L):(3L * bat)
    block <- raw[, columns, drop = FALSE]
    keep <- complete.cases(block)
    block <- block[keep, , drop = FALSE]

    bats[[bat]] <- data.frame(
      time = as.numeric(block[[1L]]),
      Longitude = as.numeric(block[[2L]]),
      Latitude = as.numeric(block[[3L]])
    )

    if (any(!is.finite(as.matrix(bats[[bat]])))) {
      stop("Hay valores no finitos en el bloque del murcielago ", bat, ".")
    }
    if (is.unsorted(bats[[bat]]$time, strictly = TRUE)) {
      stop("Los tiempos del murcielago ", bat, " no son crecientes.")
    }
  }

  names(bats) <- paste0("Bat", seq_len(5L))
  bats
}

BATS <- read_bat_workbook(DATA_PATH)

# Estos tamanos identifican de forma inequívoca el libro usado para las tablas
# del articulo. El control evita ejecutar inadvertidamente el script sobre una
# version re-muestreada o truncada de los datos.
EXPECTED_N <- c(47L, 45L, 62L, 241L, 165L)
OBSERVED_N <- unname(vapply(BATS, nrow, integer(1L)))
if (!identical(OBSERVED_N, EXPECTED_N)) {
  stop(
    "Tamanos de muestra inesperados. Esperados: ",
    paste(EXPECTED_N, collapse = ", "), "; observados: ",
    paste(OBSERVED_N, collapse = ", "), ".",
    call. = FALSE
  )
}

data_audit <- do.call(rbind, lapply(seq_along(BATS), function(bat) {
  data.frame(
    Bat = bat,
    n = nrow(BATS[[bat]]),
    time_min = min(BATS[[bat]]$time),
    time_max = max(BATS[[bat]]$time),
    median_time_gap = median(diff(BATS[[bat]]$time)),
    min_time_gap = min(diff(BATS[[bat]]$time)),
    max_time_gap = max(diff(BATS[[bat]]$time))
  )
}))
write.csv(data_audit, file.path(TAB_DIR, "data_audit.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 3. Medias y varianzas moviles
# -----------------------------------------------------------------------------

centered_rolling_apply <- function(x, width = 10L, FUN) {
  # Para width=10, la primera ventana x[1:10] se asigna a la sexta
  # observacion. Esta es la convencion centrada de stats::filter para una
  # ventana par y reproduce las Figuras C.4 y C.5.
  x <- as.numeric(x)
  n <- length(x)
  if (width < 2L || width > n) stop("Ancho de ventana no valido.")

  result <- rep(NA_real_, n)
  left <- floor(width / 2L)
  right <- width - left - 1L
  centers <- (left + 1L):(n - right)

  for (center in centers) {
    index <- (center - left):(center + right)
    result[center] <- FUN(x[index])
  }
  result
}

rolling_rows <- list()
rolling_index <- 0L
for (bat in seq_along(BATS)) {
  for (coordinate in c("Longitude", "Latitude")) {
    x <- BATS[[bat]][[coordinate]]
    rolling_index <- rolling_index + 1L
    rolling_rows[[rolling_index]] <- data.frame(
      Bat = bat,
      Coordinate = coordinate,
      Observation = seq_along(x) - 1L,
      Time = BATS[[bat]]$time,
      Value = x,
      Rolling_mean = centered_rolling_apply(x, 10L, mean),
      # stats::var usa el denominador width-1, es decir, varianza muestral.
      Rolling_variance = centered_rolling_apply(x, 10L, stats::var)
    )
  }
}
rolling_data <- do.call(rbind, rolling_rows)
write.csv(
  rolling_data,
  file.path(DATA_DIR, "rolling_mean_variance_values.csv"),
  row.names = FALSE
)

plot_rolling_panels <- function(variable, y_text, filename) {
  grDevices::pdf(filename, width = 11, height = 8.5, useDingbats = FALSE)
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); grDevices::dev.off()}, add = TRUE)
  par(mfrow = c(5, 2), mar = c(2.2, 3.4, 2.2, 1.0), oma = c(1.5, 1.5, 0, 0))

  for (bat in seq_along(BATS)) {
    for (coordinate in c("Longitude", "Latitude")) {
      take <- rolling_data$Bat == bat & rolling_data$Coordinate == coordinate
      panel <- rolling_data[take, ]
      good <- is.finite(panel[[variable]])
      plot(
        panel$Observation[good], panel[[variable]][good],
        type = "l", col = "black", lwd = 1,
        xlab = if (bat == 5L) "Observation index" else "",
        ylab = y_text,
        main = paste0("Bat ", bat, " - ", coordinate, ": ",
                      if (variable == "Rolling_mean") "Rolling Mean" else
                        "Rolling Variance"),
        cex.main = 0.78, cex.axis = 0.72, cex.lab = 0.72
      )
    }
  }
}

if (!ARGS$skip_figures) {
  plot_rolling_panels(
    "Rolling_mean", "Rolling mean",
    file.path(FIG_DIR, "Rmean.pdf")
  )
  plot_rolling_panels(
    "Rolling_variance", "Rolling variance",
    file.path(FIG_DIR, "Rvariance.pdf")
  )
}


# -----------------------------------------------------------------------------
# 4. Funciones de autocorrelacion: nivel, primera y segunda diferencia
# -----------------------------------------------------------------------------

default_acf_lag <- function(n) {
  # Esta es la regla por defecto usada por stats::acf para una serie univariada.
  min(n - 1L, floor(10 * log10(n)))
}

acf_numeric <- function(x, lag.max = default_acf_lag(length(x))) {
  # stats::acf calcula la autocorrelacion muestral con el mismo denominador n
  # para todas las autocovarianzas. Se conserva el rezago cero en la salida.
  as.numeric(
    stats::acf(
      as.numeric(x), lag.max = lag.max, type = "correlation",
      plot = FALSE, demean = TRUE, na.action = na.fail
    )$acf
  )
}

acf_rows <- list()
acf_index <- 0L
for (bat in seq_along(BATS)) {
  for (coordinate in c("Longitude", "Latitude")) {
    x0 <- BATS[[bat]][[coordinate]]
    for (difference_order in 0:2) {
      x <- if (difference_order == 0L) x0 else
        diff(x0, differences = difference_order)
      values <- acf_numeric(x)
      bound <- stats::qnorm(0.975) / sqrt(length(x))
      acf_index <- acf_index + 1L
      acf_rows[[acf_index]] <- data.frame(
        Bat = bat,
        Coordinate = coordinate,
        Difference_order = difference_order,
        n = length(x),
        Lag = seq_along(values) - 1L,
        ACF = values,
        Lower_95 = -bound,
        Upper_95 = bound,
        Outside_95 = seq_along(values) > 1L & abs(values) > bound
      )
    }
  }
}
acf_data <- do.call(rbind, acf_rows)
write.csv(acf_data, file.path(DATA_DIR, "acf_values.csv"), row.names = FALSE)

acf_diagnostic_summary <- do.call(rbind, lapply(
  split(
    acf_data,
    interaction(
      acf_data$Bat, acf_data$Coordinate, acf_data$Difference_order,
      drop = TRUE
    )
  ),
  function(panel) {
    outside <- panel$Lag[panel$Outside_95]
    data.frame(
      Bat = panel$Bat[1L],
      Coordinate = panel$Coordinate[1L],
      Difference_order = panel$Difference_order[1L],
      n = panel$n[1L],
      Lag1_ACF = panel$ACF[panel$Lag == 1L],
      Bound_95 = panel$Upper_95[1L],
      Lags_outside_95 = if (length(outside)) paste(outside, collapse = ",") else ""
    )
  }
))
acf_diagnostic_summary <- acf_diagnostic_summary[
  order(
    acf_diagnostic_summary$Difference_order,
    acf_diagnostic_summary$Bat,
    match(acf_diagnostic_summary$Coordinate, c("Longitude", "Latitude"))
  ),
]
write.csv(
  acf_diagnostic_summary,
  file.path(TAB_DIR, "acf_diagnostic_summary.csv"),
  row.names = FALSE
)

plot_acf_grid <- function(difference_order, filename, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(filename, width = 11, height = 9.5, useDingbats = FALSE)
  } else {
    grDevices::png(
      filename, width = 1800, height = 2400, res = 220,
      type = if (capabilities("cairo")) "cairo" else getOption("bitmapType")
    )
  }

  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); grDevices::dev.off()}, add = TRUE)
  par(mfrow = c(5, 2), mar = c(2.3, 3.2, 2.8, 0.8), oma = c(1.3, 1.3, 0, 0))

  prefix <- c("", "first-differenced ", "second-differenced ")[
    difference_order + 1L
  ]

  for (bat in seq_along(BATS)) {
    for (coordinate in c("Longitude", "Latitude")) {
      x0 <- BATS[[bat]][[coordinate]]
      x <- if (difference_order == 0L) x0 else
        diff(x0, differences = difference_order)
      values <- acf_numeric(x)
      lags <- seq_along(values) - 1L
      bound <- stats::qnorm(.975) / sqrt(length(x))
      y_limits <- range(c(values, -bound, bound), finite = TRUE)
      padding <- .06 * diff(y_limits)

      plot(
        lags, values,
        type = "h", lwd = 1, col = "black",
        ylim = y_limits + c(-padding, padding),
        xlab = if (bat == 5L) "Lag" else "",
        ylab = "ACF", main = "",
        cex.axis = .72, cex.lab = .72
      )
      abline(h = 0, col = "black", lwd = .8)
      abline(
        h = c(-bound, bound), col = "blue",
        lty = if (difference_order == 2L) 1 else 2, lwd = .9
      )
      title(
        main = paste0("ACF Bat ", bat, " - ", prefix, coordinate),
        line = .8, cex.main = .74
      )
    }
  }
}

if (!ARGS$skip_figures) {
  plot_acf_grid(0L, file.path(FIG_DIR, "acf.pdf"), "pdf")
  plot_acf_grid(1L, file.path(FIG_DIR, "ACFD.pdf"), "pdf")
  plot_acf_grid(
    2L,
    file.path(FIG_DIR, "Figure_C8_ACF_second_differences.png"),
    "png"
  )
}

# Control especifico de la afirmacion de la Figura C.8: en las diez series la
# ACF de rezago uno debe caer por debajo del limite inferior aproximado.
second_difference_summary <- subset(acf_diagnostic_summary, Difference_order == 2L)
if (!all(second_difference_summary$Lag1_ACF < -second_difference_summary$Bound_95)) {
  stop("No se reproduce la afirmacion de rezago uno de la Figura C.8.")
}


# -----------------------------------------------------------------------------
# 5. Prueba ADF: nula de raiz unitaria / no estacionariedad
# -----------------------------------------------------------------------------

adf_test <- function(x, k = NULL) {
  # Regresion ADF usada por tseries::adf.test:
  #
  #   Delta x_t = a + b*x_(t-1) + c*t
  #               + sum_j gamma_j Delta x_(t-j) + error_t.
  #
  # La estadistica es el cociente t asociado con b. La nula es b=0 (raiz
  # unitaria); un valor p menor que 0.05 se etiqueta "Stat." en la Tabla C.1.
  x <- as.numeric(stats::na.omit(x))
  N <- length(x)
  if (N < 8L) return(c(statistic = NA, p.value = NA, lag = NA))
  if (is.null(k)) k <- floor((N - 1L)^(1 / 3))

  dx <- diff(x - x[1L])
  n <- length(dx)
  tt <- (k + 1L):n
  response <- dx[tt]
  design <- cbind(intercept = 1, lagged_level = x[tt], trend = tt)
  if (k > 0L) {
    lagged_differences <- sapply(seq_len(k), function(j) dx[tt - j])
    design <- cbind(design, lagged_differences)
  }

  fit <- lm.fit(design, response)
  residuals <- fit$residuals
  degrees_freedom <- length(response) - ncol(design)
  sigma2 <- sum(residuals^2) / degrees_freedom
  covariance <- sigma2 * solve(crossprod(design))
  statistic <- unname(
    fit$coefficients["lagged_level"] /
      sqrt(covariance["lagged_level", "lagged_level"])
  )

  # Tabla de superficies de respuesta utilizada por tseries. Las columnas
  # corresponden a tamanos de muestra y las filas a probabilidades.
  probabilities <- c(.01, .025, .05, .10, .90, .95, .975, .99)
  sample_sizes <- c(25, 50, 100, 250, 500, 100000)
  response_table <- rbind(
    c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
    c(-3.95, -3.80, -3.73, -3.69, -3.68, -3.66),
    c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41),
    c(-3.24, -3.18, -3.15, -3.13, -3.13, -3.12),
    c(-1.14, -1.19, -1.22, -1.23, -1.24, -1.25),
    c(-0.80, -0.87, -0.90, -0.92, -0.93, -0.94),
    c(-0.50, -0.58, -0.62, -0.64, -0.65, -0.66),
    c( 0.00, -0.01, -0.01, -0.02, -0.02, -0.02)
  )

  # tseries interpola con n = length(diff(x)) = N-1.
  critical_values <- apply(
    response_table, 1L,
    function(row) approx(
      sample_sizes, row, xout = N - 1L, rule = 2
    )$y
  )
  p_value <- approx(
    critical_values, probabilities, xout = statistic,
    rule = 2, ties = "ordered"
  )$y
  p_value <- max(.01, min(.99, p_value))

  c(statistic = statistic, p.value = p_value, lag = k)
}


# -----------------------------------------------------------------------------
# 6. Prueba KPSS: nula de estacionariedad alrededor de un nivel
# -----------------------------------------------------------------------------

kpss_level_test <- function(x, lags = NULL) {
  # Bajo la nula, x_t = nivel constante + proceso estacionario. Se restan la
  # media y se forman sumas parciales S_t. La estadistica es
  #
  #   eta / LRV,  eta = n^(-2) sum_t S_t^2,
  #
  # donde LRV es la varianza de largo plazo estimada con pesos de Bartlett.
  # Un valor p menor que 0.05 se etiqueta "Non-stat." en la Tabla C.1.
  x <- as.numeric(stats::na.omit(x))
  n <- length(x)
  if (n < 6L) return(c(statistic = NA, p.value = NA, lag = NA))

  residuals <- x - mean(x)
  if (is.null(lags)) {
    # Regla corta de tseries::kpss.test(..., lshort=TRUE).
    lags <- min(n - 2L, floor(4 * (n / 100)^.25))
  }

  eta <- sum(cumsum(residuals)^2) / n^2
  long_run_variance <- sum(residuals^2) / n
  if (lags > 0L) {
    for (lag in seq_len(lags)) {
      autocovariance <- sum(
        residuals[(lag + 1L):n] * residuals[seq_len(n - lag)]
      ) / n
      bartlett_weight <- 1 - lag / (lags + 1L)
      long_run_variance <- long_run_variance +
        2 * bartlett_weight * autocovariance
    }
  }

  statistic <- eta / long_run_variance
  critical_values <- c(.347, .463, .574, .739)
  probabilities <- c(.10, .05, .025, .01)
  p_value <- approx(
    critical_values, probabilities, xout = statistic, rule = 2
  )$y
  p_value <- max(.01, min(.10, p_value))

  c(statistic = statistic, p.value = p_value, lag = lags)
}


# -----------------------------------------------------------------------------
# 7. Tabla C.1: ADF y KPSS para niveles y primeras diferencias
# -----------------------------------------------------------------------------

stationarity_rows <- list()
stationarity_index <- 0L
for (difference_order in 0:1) {
  for (bat in seq_along(BATS)) {
    for (coordinate in c("Longitude", "Latitude")) {
      x0 <- BATS[[bat]][[coordinate]]
      x <- if (difference_order == 0L) x0 else diff(x0)
      adf <- adf_test(x)
      kpss <- kpss_level_test(x)

      stationarity_index <- stationarity_index + 1L
      stationarity_rows[[stationarity_index]] <- data.frame(
        Bat = bat,
        Series = if (difference_order == 0L) coordinate else
          paste0("Delta ", coordinate),
        Difference_order = difference_order,
        n = length(x),
        ADF_statistic = unname(adf["statistic"]),
        ADF_p = unname(adf["p.value"]),
        ADF_lag = as.integer(adf["lag"]),
        KPSS_statistic = unname(kpss["statistic"]),
        KPSS_p = unname(kpss["p.value"]),
        KPSS_lag = as.integer(kpss["lag"]),
        ADF_decision = if (adf["p.value"] < .05) "Stat." else "Non-stat.",
        KPSS_decision = if (kpss["p.value"] < .05) "Non-stat." else "Stat."
      )
    }
  }
}
stationarity_table <- do.call(rbind, stationarity_rows)
write.csv(
  stationarity_table,
  file.path(TAB_DIR, "table_C1_adf_kpss.csv"),
  row.names = FALSE
)

# Valores p publicados, en el mismo orden del bucle anterior. Esta comparacion
# constituye una prueba automatica de reproduccion de la Tabla C.1.
EXPECTED_ADF_P <- c(
  .9724, .4993, .8021, .8010, .0100, .3952, .0415, .4264, .5748, .4934,
  .3341, .0100, .3980, .1770, .0100, .0100, .0100, .0262, .3486, .2016
)
EXPECTED_KPSS_P <- c(
  .0100, .0100, .0100, .0100, .0100, .0100, .0100, .0100, .0100, .0100,
  .1000, .1000, .1000, .0178, .0991, .1000, .0100, .1000, .0100, .0100
)

if (!identical(round(stationarity_table$ADF_p, 4), EXPECTED_ADF_P)) {
  stop("Los valores p de ADF no reproducen la Tabla C.1 al cuarto decimal.")
}
if (!identical(round(stationarity_table$KPSS_p, 4), EXPECTED_KPSS_P)) {
  stop("Los valores p de KPSS no reproducen la Tabla C.1 al cuarto decimal.")
}

latex_escape_delta <- function(series) {
  if (grepl("^Delta ", series)) {
    paste0("$\\Delta$ ", sub("^Delta ", "", series))
  } else {
    series
  }
}

write_table_C1_latex <- function(table, path) {
  connection <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  writeLines(c(
    "% Generado automaticamente por appendix_C_empirical_diagnostics.R",
    "\\begin{tabular}{@{}cccccc@{}}",
    "\\toprule",
    "\\textbf{Bat} & \\textbf{Series} & \\textbf{ADF $p$} & \\textbf{KPSS $p$} & \\textbf{ADF} & \\textbf{KPSS} \\\\",
    "\\midrule",
    "\\multicolumn{6}{l}{\\textit{Panel A: Raw series}} \\\\",
    "\\midrule"
  ), connection)

  for (panel in 0:1) {
    if (panel == 1L) {
      writeLines(c(
        "\\midrule",
        "\\multicolumn{6}{l}{\\textit{Panel B: First differences}} \\\\",
        "\\midrule"
      ), connection)
    }
    rows <- table[table$Difference_order == panel, ]
    for (i in seq_len(nrow(rows))) {
      writeLines(sprintf(
        "%d & %s & %.4f & %.4f & %s & %s \\\\",
        rows$Bat[i], latex_escape_delta(rows$Series[i]),
        rows$ADF_p[i], rows$KPSS_p[i],
        rows$ADF_decision[i], rows$KPSS_decision[i]
      ), connection)
    }
  }
  writeLines(c("\\bottomrule", "\\end{tabular}"), connection)
}

write_table_C1_latex(
  stationarity_table,
  file.path(TAB_DIR, "table_C1_adf_kpss.tex")
)


# -----------------------------------------------------------------------------
# 8. Detrended fluctuation analysis (DFA de orden 1)
# -----------------------------------------------------------------------------

dfa_scale_grid <- function(n) {
  # Para las tres series cortas se usan los 12 enteros 4,...,15. Para las dos
  # series largas se usan 12 escalas logaritmicamente espaciadas desde 8 hasta
  # floor(n/4), redondeadas al entero mas cercano. Esta regla reproduce los
  # puntos negros de la Figura C.9:
  #
  # Bat 4: 8,10,12,14,17,20,24,29,35,42,50,60
  # Bat 5: 8, 9,11,12,14,17,20,23,26,30,35,41.
  if (n <= 62L) return(4:15)
  unique(as.integer(round(exp(seq(log(8), log(floor(n / 4)), length.out = 12L)))))
}

dfa_fluctuation <- function(x, scale) {
  # Paso 1: perfil acumulado de la serie centrada.
  profile <- cumsum(as.numeric(x) - mean(x))

  # Paso 2: segmentos consecutivos no superpuestos, comenzando en el origen.
  number_of_segments <- floor(length(profile) / scale)
  if (number_of_segments < 2L) return(NA_real_)

  squared_residuals <- numeric(0)
  local_time <- seq_len(scale)
  design <- cbind(intercept = 1, local_time = local_time)

  # Paso 3: en cada segmento se elimina una tendencia lineal (DFA1).
  for (segment in 0:(number_of_segments - 1L)) {
    index <- segment * scale + seq_len(scale)
    local_profile <- profile[index]
    coefficients <- qr.solve(design, local_profile)
    residuals <- local_profile - as.numeric(design %*% coefficients)
    squared_residuals <- c(squared_residuals, residuals^2)
  }

  # Paso 4: raiz de la media de todos los residuos cuadrados.
  sqrt(mean(squared_residuals))
}

dfa_fit <- function(x, bat) {
  scales <- dfa_scale_grid(length(x))
  fluctuation <- vapply(scales, function(s) dfa_fluctuation(x, s), numeric(1L))
  if (any(!is.finite(fluctuation) | fluctuation <= 0)) {
    stop("DFA produjo una fluctuacion no positiva para Bat ", bat, ".")
  }

  # Los seis puntos de mayor escala son los puntos rojos de la Figura C.9.
  selected <- tail(seq_along(scales), 6L)
  regression <- lm(log(fluctuation[selected]) ~ log(scales[selected]))
  interval <- confint(regression, level = .95)[2L, ]

  list(
    scales = scales,
    fluctuation = fluctuation,
    selected = selected,
    regression = regression,
    alpha = unname(coef(regression)[2L]),
    ci_low = unname(interval[1L]),
    ci_high = unname(interval[2L]),
    r_squared = summary(regression)$r.squared
  )
}

dfa_results <- vector("list", 10L)
dfa_rows <- list()
dfa_curve_rows <- list()
dfa_index <- 0L
for (bat in seq_along(BATS)) {
  for (coordinate in c("Longitude", "Latitude")) {
    dfa_index <- dfa_index + 1L
    fit <- dfa_fit(BATS[[bat]][[coordinate]], bat)
    dfa_results[[dfa_index]] <- c(list(Bat = bat, Coordinate = coordinate), fit)

    fitted_fluctuation <- rep(NA_real_, length(fit$scales))
    fitted_fluctuation[fit$selected] <- exp(predict(fit$regression))

    dfa_rows[[dfa_index]] <- data.frame(
      Bat = bat,
      Coordinate = coordinate,
      alpha = fit$alpha,
      CI_low = fit$ci_low,
      CI_high = fit$ci_high,
      R_squared = fit$r_squared,
      regression_scales = paste(fit$scales[fit$selected], collapse = ",")
    )
    dfa_curve_rows[[dfa_index]] <- data.frame(
      Bat = bat,
      Coordinate = coordinate,
      Scale = fit$scales,
      Fluctuation = fit$fluctuation,
      Selected_for_regression = seq_along(fit$scales) %in% fit$selected,
      Fitted_fluctuation = fitted_fluctuation
    )
  }
}
dfa_table_long <- do.call(rbind, dfa_rows)
dfa_curve_data <- do.call(rbind, dfa_curve_rows)
write.csv(
  dfa_table_long,
  file.path(TAB_DIR, "table_C2_dfa_long.csv"),
  row.names = FALSE
)
write.csv(
  dfa_curve_data,
  file.path(DATA_DIR, "dfa_fluctuation_values.csv"),
  row.names = FALSE
)

# Tabla ancha con una fila por murcielago, equivalente a la Tabla C.2.
dfa_table_wide <- do.call(rbind, lapply(seq_len(5L), function(bat) {
  longitude <- subset(dfa_table_long, Bat == bat & Coordinate == "Longitude")
  latitude <- subset(dfa_table_long, Bat == bat & Coordinate == "Latitude")
  data.frame(
    Bat = bat,
    Longitude_alpha = longitude$alpha,
    Longitude_CI_low = longitude$CI_low,
    Longitude_CI_high = longitude$CI_high,
    Longitude_R_squared = longitude$R_squared,
    Latitude_alpha = latitude$alpha,
    Latitude_CI_low = latitude$CI_low,
    Latitude_CI_high = latitude$CI_high,
    Latitude_R_squared = latitude$R_squared
  )
}))
write.csv(
  dfa_table_wide,
  file.path(TAB_DIR, "table_C2_dfa.csv"),
  row.names = FALSE
)

# Prueba automatica de reproduccion de la Tabla C.2.
EXPECTED_DFA <- data.frame(
  Bat = 1:5,
  Longitude_alpha = c(2.349, 1.727, 1.103, 1.907, 2.084),
  Longitude_CI_low = c(1.839, .963, .830, 1.799, 1.305),
  Longitude_CI_high = c(2.858, 2.490, 1.376, 2.016, 2.863),
  Longitude_R_squared = c(.976, .908, .969, .998, .932),
  Latitude_alpha = c(1.933, 1.836, 1.577, 1.962, 1.982),
  Latitude_CI_low = c(1.688, 1.603, .539, 1.827, 1.320),
  Latitude_CI_high = c(2.179, 2.069, 2.615, 2.098, 2.644),
  Latitude_R_squared = c(.992, .992, .816, .998, .945)
)

dfa_columns <- setdiff(names(EXPECTED_DFA), "Bat")
if (!isTRUE(all.equal(
  round(dfa_table_wide[dfa_columns], 3),
  EXPECTED_DFA[dfa_columns],
  check.attributes = FALSE,
  tolerance = 0
))) {
  stop("Los resultados DFA no reproducen la Tabla C.2 al tercer decimal.")
}

write_table_C2_latex <- function(table, path) {
  connection <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  writeLines(c(
    "% Generado automaticamente por appendix_C_empirical_diagnostics.R",
    "\\begin{tabular}{|c|c|c|c|c|}",
    "\\hline",
    paste0(
      "\\textbf{Bat ID} & \\textbf{Longitude: $\\widehat\\alpha$ [95\\% CI]} & ",
      "\\textbf{$R^2$} & \\textbf{Latitude: $\\widehat\\alpha$ [95\\% CI]} & ",
      "\\textbf{$R^2$} \\\\"
    ),
    "\\hline"
  ), connection)

  for (i in seq_len(nrow(table))) {
    writeLines(sprintf(
      paste0(
        "%d & %.3f [%.3f, %.3f] & %.3f & ",
        "%.3f [%.3f, %.3f] & %.3f \\\\"
      ),
      table$Bat[i],
      table$Longitude_alpha[i], table$Longitude_CI_low[i],
      table$Longitude_CI_high[i], table$Longitude_R_squared[i],
      table$Latitude_alpha[i], table$Latitude_CI_low[i],
      table$Latitude_CI_high[i], table$Latitude_R_squared[i]
    ), connection)
  }
  writeLines(c("\\hline", "\\end{tabular}"), connection)
}

write_table_C2_latex(
  dfa_table_wide,
  file.path(TAB_DIR, "table_C2_dfa.tex")
)


# -----------------------------------------------------------------------------
# 9. Figuras DFA
# -----------------------------------------------------------------------------

plot_dfa_curves <- function(path) {
  grDevices::pdf(path, width = 12, height = 5.8, useDingbats = FALSE)
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); grDevices::dev.off()}, add = TRUE)
  par(mfrow = c(2, 5), mar = c(3.0, 3.2, 2.2, .8), oma = c(1.0, 1.2, 2.0, .5))

  # La figura se organiza por coordenada: Longitud arriba, Latitud abajo.
  for (coordinate in c("Longitude", "Latitude")) {
    for (bat in seq_len(5L)) {
      result <- dfa_results[[which(vapply(
        dfa_results,
        function(z) z$Bat == bat && z$Coordinate == coordinate,
        logical(1L)
      ))]]
      selected <- result$selected

      plot(
        result$scales, result$fluctuation,
        log = "xy", type = "o", pch = 16, cex = .55,
        col = "black", lwd = 1,
        xlab = if (coordinate == "Latitude") "Window size (scale s)" else "",
        ylab = if (bat == 1L) "Fluctuation F(s)" else "",
        main = paste("Bat", bat), cex.main = .85,
        cex.axis = .72, cex.lab = .72
      )
      points(
        result$scales[selected], result$fluctuation[selected],
        pch = 16, col = "#C43C2C", cex = .7
      )
      lines(
        result$scales[selected], exp(predict(result$regression)),
        col = "#C43C2C", lwd = 1.7
      )
      if (bat == 5L) {
        mtext(coordinate, side = 4, line = .15, cex = .75)
      }
    }
  }
  mtext(
    "Detrended Fluctuation Analysis (DFA)",
    side = 3, outer = TRUE, line = .3, font = 2
  )
}

plot_dfa_summary <- function(path) {
  grDevices::pdf(path, width = 8.5, height = 8, useDingbats = FALSE)
  old_par <- par(no.readonly = TRUE)
  on.exit({par(old_par); grDevices::dev.off()}, add = TRUE)
  par(mfrow = c(2, 1), mar = c(3.0, 4.0, 2.2, 1.0), oma = c(1.4, .6, 2.0, .5))

  y_limits <- range(dfa_table_long[, c("CI_low", "CI_high")], finite = TRUE)
  y_limits <- c(min(.40, y_limits[1L]), max(2.9, y_limits[2L]))

  for (coordinate in c("Longitude", "Latitude")) {
    panel <- subset(dfa_table_long, Coordinate == coordinate)
    symbol <- if (coordinate == "Longitude") 16 else 17
    plot(
      panel$Bat, panel$alpha,
      ylim = y_limits, xlim = c(.7, 5.3), xaxt = "n",
      pch = symbol, cex = .9, xlab = "", ylab = "DFA exponent",
      main = coordinate, cex.main = .9
    )
    axis(1, at = 1:5, labels = paste("Bat ID", 1:5), cex.axis = .75)
    segments(
      panel$Bat, panel$CI_low, panel$Bat, panel$CI_high,
      lwd = 1.2
    )
    segments(
      panel$Bat - .08, panel$CI_low, panel$Bat + .08, panel$CI_low,
      lwd = 1.2
    )
    segments(
      panel$Bat - .08, panel$CI_high, panel$Bat + .08, panel$CI_high,
      lwd = 1.2
    )
    abline(h = .5, lty = 2, col = "grey30")
    abline(h = 1, lty = 4, col = "grey30")
  }
  mtext(
    "Estimated DFA exponents across bats and coordinates",
    side = 3, outer = TRUE, line = .3, font = 2
  )
}

if (!ARGS$skip_figures) {
  plot_dfa_curves(file.path(FIG_DIR, "dfa_curves_english.pdf"))
  plot_dfa_summary(file.path(FIG_DIR, "dfa_alpha_summary_english.pdf"))
}


# -----------------------------------------------------------------------------
# 10. Resumen reproducible y controles finales
# -----------------------------------------------------------------------------

first_difference_tests <- subset(stationarity_table, Difference_order == 1L)
joint_stationary <- with(
  first_difference_tests,
  ADF_decision == "Stat." & KPSS_decision == "Stat."
)
kpss_stationary <- first_difference_tests$KPSS_decision == "Stat."

if (sum(joint_stationary) != 4L) {
  stop("No se reproducen los cuatro casos conjuntamente estacionarios.")
}
if (sum(kpss_stationary) != 6L) {
  stop("No se reproducen los seis casos sin rechazo KPSS en diferencias.")
}

validation_lines <- c(
  "VALIDACION DEL APENDICE C",
  "==========================",
  paste("Archivo de datos:", DATA_PATH),
  paste("MD5 del archivo:", unname(tools::md5sum(DATA_PATH))),
  paste("Tamanos por murcielago:", paste(OBSERVED_N, collapse = ", ")),
  "",
  "Tabla C.1:",
  "  - Los 20 valores p de ADF coinciden al cuarto decimal: SI.",
  "  - Los 20 valores p de KPSS coinciden al cuarto decimal: SI.",
  paste0(
    "  - Primeras diferencias clasificadas como estacionarias por ambas ",
    "pruebas: ", sum(joint_stationary), " de 10."
  ),
  paste0(
    "  - Primeras diferencias para las que KPSS no rechaza estacionariedad: ",
    sum(kpss_stationary), " de 10."
  ),
  "",
  "Figura C.8:",
  "  - ACF(1) por debajo del limite inferior en 10 de 10 segundas diferencias: SI.",
  "",
  "Tabla C.2:",
  "  - Los 10 exponentes, 20 limites y 10 R^2 coinciden al tercer decimal: SI.",
  "  - Los 10 exponentes puntuales son mayores que 1: SI.",
  "",
  "Interpretacion:",
  paste0(
    "  La evidencia contra estacionariedad de los niveles es general, pero ",
    "la evidencia sobre estacionariedad de las primeras diferencias es mixta."
  ),
  paste0(
    "  Los cuatro casos apoyados simultaneamente por ADF y KPSS, y los dos ",
    "casos adicionales sin rechazo KPSS, hacen plausible considerar modelos ",
    "con incrementos estacionarios (como fBM) junto con alternativas cuyos ",
    "incrementos no son estacionarios."
  ),
  paste0(
    "  Esta lectura es compatible con una comparacion por AIC, pero no ",
    "sustituye el ajuste de la covarianza completa ni implica que las diez ",
    "trayectorias tengan incrementos estacionarios."
  )
)
writeLines(validation_lines, file.path(OUT_DIR, "APPENDIX_C_validation.txt"))

capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

cat(
  paste0(
    "Apendice C reproducido correctamente.\n",
    "Resultados: ", normalizePath(OUT_DIR, mustWork = TRUE), "\n",
    "ADF/KPSS: coincidencia exacta al redondeo publicado.\n",
    "DFA: coincidencia exacta al redondeo publicado.\n"
  )
)
