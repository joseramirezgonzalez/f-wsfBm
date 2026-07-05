library(shiny)
library(MASS)
library(leaflet)
library(readxl)
library(mvtnorm)

# Function to convert R-style expressions to LaTeX
toLatex <- function(expr) {
  expr <- gsub("sqrt\\(([^\\)]+)\\)", "\\\\sqrt{\\1}", expr)
  expr <- gsub("sin\\(([^\\)]+)\\)", "\\\\sin{\\1}", expr)
  expr <- gsub("cos\\(([^\\)]+)\\)", "\\\\cos{\\1}", expr)
  expr <- gsub("tan\\(([^\\)]+)\\)", "\\\\tan{\\1}", expr)
  expr <- gsub("log\\(([^\\)]+)\\)", "\\\\log{\\1}", expr)
  expr <- gsub("exp\\(([^\\)]+)\\)", "\\\\exp{\\1}", expr)
  expr <- gsub("\\*", " \\\\cdot ", expr)
  expr <- gsub("\\^([a-zA-Z0-9]+)", "^{\\1}", expr)
  return(expr)
}

# Wrap functions
wrap_longitude <- function(lon) {
  ((lon + 180) %% 360) - 180
}

wrap_latitude <- function(lat) {
  ((lat + 90) %% 180) - 90
}

# Exponential class used in the paper:
# f_{sigma^2,beta}(u) = sigma^2 exp(-beta u).
# If beta = 0, this reduces to the one-parameter boundary model
# f_{sigma^2,0}(u) = sigma^2.
make_exp_weight <- function(sigma2, beta) {
  sigma2 <- as.numeric(sigma2)
  beta <- as.numeric(beta)
  if (!is.finite(sigma2) || sigma2 <= 0) {
    stop("sigma^2 must be positive.")
  }
  if (!is.finite(beta) || beta < 0) {
    stop("beta must be nonnegative.")
  }
  if (abs(beta) < .Machine$double.eps) {
    return(function(u) rep(sigma2, length(u)))
  }
  function(u) sigma2 * exp(-beta * u)
}

# Covariance kernel construction (simulation tabs)
cov_f_log <- function(f, u, s, t) {
  return(2 * f(u) * (-xlogx_safe(t - u) - xlogx_safe(s - u) + xlogx_safe(t + s - 2 * u)))
}

cov_f_log_w <- function(f, s, t) {
  g <- function(z) cov_f_log(f, z, s, t)
  result <- tryCatch(integrate(g, lower = 0, upper = min(s, t))$value, error = function(e) NA)
  return(result)
}

cov_f_weighted <- function(f, t) {
  n <- length(t)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:i) {
      K[i, j] <- cov_f_log_w(f, t[j], t[i])
      K[j, i] <- K[i, j]
    }
  }
  return(K)
}

# ------------------------------------------------------------
# Inference (practical likelihood) for the exponential class
# ------------------------------------------------------------

xlogx_safe <- function(x) {
  # Convention: x * log(x) = 0 at x = 0
  ifelse(x <= 0, 0, x * log(x))
}

Phi_xlogx <- function(x) {
  # Antiderivative of x log(x), with the continuous value 0 at x=0.
  y <- numeric(length(x))
  ind <- x > 0
  y[ind] <- 0.5 * x[ind]^2 * log(x[ind]) - 0.25 * x[ind]^2
  y
}

cov_2_beta0_matrix <- function(t) {
  # Closed-form unit-scale covariance for beta=0, i.e. f_{sigma^2,0}(u)=sigma^2.
  # This implementation avoids building several auxiliary n x n matrices.
  # K(s,t)=Phi(s+t)+Phi(|s-t|)-2Phi(s)-2Phi(t), where Phi'(x)=x log(x).
  n <- length(t)
  K <- matrix(0, n, n)
  Pt <- Phi_xlogx(t)
  for (j in seq_len(n)) {
    s <- t[j]
    tv <- t[j:n]
    vals <- Phi_xlogx(s + tv) + Phi_xlogx(abs(s - tv)) - 2 * Pt[j] - 2 * Pt[j:n]
    K[j:n, j] <- vals
    K[j, j:n] <- vals
  }
  K
}

D_2 <- function(u, s, t, b) {
  exp(-b * u) * ( -xlogx_safe(t - u) - xlogx_safe(s - u) + xlogx_safe(t + s - 2 * u) )
}

cov_2 <- function(s, t, b) {
  y <- min(s, t)
  if (y <= 0) return(0)
  g <- function(z) D_2(z, s = s, t = t, b = b)
  out <- tryCatch(integrate(g, lower = 0, upper = y, rel.tol = 1e-6)$value,
                  error = function(e) NA_real_)
  return(2 * out)
}

cov_wfbm_2 <- function(t, b) {
  b <- as.numeric(b)
  if (!is.finite(b) || b < 0) stop("beta must be nonnegative.")
  if (abs(b) < .Machine$double.eps) {
    return(cov_2_beta0_matrix(t))
  }

  n <- length(t)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:i) {
      K[i, j] <- cov_2(t[j], t[i], b = b)
      K[j, i] <- K[i, j]
    }
  }
  return(K)
}

log_verosimilitude_practical <- function(T, b, mu, index) {
  # mu: vector of observations (lon or lat), length = k+1 (including time 0)
  # index: integer vector of observed time indices (in steps), length = k
  # b: beta. The value b=0 corresponds to the boundary model f_{sigma^2,0}(u)=sigma^2.
  k <- length(mu) - 1
  n <- index[k]
  D <- T / n
  t_obs <- index * D

  M <- cov_wfbm_2(t_obs, b = b)
  delta_mu <- as.vector(mu[2:(k + 1)] - mu[1])

  diag_mean <- mean(diag(M), na.rm = TRUE)
  if (!is.finite(diag_mean) || diag_mean <= 0) diag_mean <- 1
  diag(M) <- diag(M) + 1e-10 * diag_mean

  R <- chol(M)
  z <- forwardsolve(t(R), delta_mu)
  Q <- sum(z^2)
  logdetM <- 2 * sum(log(diag(R)))

  # Under f_{sigma^2,beta}(u)=sigma^2 exp(-beta u), M is the covariance
  # matrix obtained with sigma^2=1. Hence the profiled MLE is Q_beta/k.
  sigma2_hat <- Q / k

  ll <- -0.5 * (k * log(2 * pi) + k * log(sigma2_hat) + logdetM + k)
  return(c(ll, sigma2_hat))
}

infer_beta_sigma <- function(T, mu, index, beta_lower, beta_upper,
                             grid_n = 60, compute_profile = TRUE,
                             use_beta0 = FALSE) {
  log_vero_pract <- function(x) {
    -log_verosimilitude_practical(T, x, mu, index)[1]
  }
  log_vero_pract_2 <- function(x) {
    log_verosimilitude_practical(T, x, mu, index)
  }

  res <- list()

  if (isTRUE(use_beta0)) {
    aux <- log_vero_pract_2(0)
    res$beta_hat <- 0
    res$beta_fixed <- TRUE
    res$sigma2_hat <- aux[2]
    res$logLik_hat <- aux[1]
    res$npar <- 1
    res$aic_hat <- 2 - 2 * aux[1]
    res$model_label <- "beta = 0"
    res$grid_beta <- NULL
    res$profile_ll <- NULL
    return(res)
  }

  # Profile likelihood in beta.
  if (compute_profile) {
    reg_b <- seq(beta_lower, beta_upper, length.out = grid_n)
    profile_ll <- rep(NA_real_, grid_n)
    for (i in seq_along(reg_b)) {
      aux <- log_vero_pract_2(reg_b[i])
      profile_ll[i] <- aux[1]
    }
    res$grid_beta <- reg_b
    res$profile_ll <- profile_ll
  }

  beta_hat <- optim((beta_lower + beta_upper) / 2,
                    fn = log_vero_pract,
                    lower = beta_lower, upper = beta_upper,
                    method = "L-BFGS-B")$par
  aux_hat <- log_vero_pract_2(beta_hat)
  sigma2_hat <- aux_hat[2]
  ll_hat <- aux_hat[1]

  res$beta_hat <- beta_hat
  res$beta_fixed <- FALSE
  res$sigma2_hat <- sigma2_hat
  res$logLik_hat <- ll_hat
  res$npar <- 2
  res$aic_hat <- 4 - 2 * ll_hat
  res$model_label <- "beta > 0"
  return(res)
}

# ------------------------------------------------------------
# Prediction (Gaussian conditional mean / simulations)
# ------------------------------------------------------------
# Given observed positions mu at times index*D (with mu[1] = initial position),
# and MLEs (beta_hat, sigma2_hat), we simulate future points using Gaussian
# conditional laws and return the mean predicted trajectory.
predict_future_conditional <- function(mu, index, D, beta_hat, sigma2_hat, sim_step) {
  sim_step <- as.integer(sim_step)
  if (sim_step <= 0) return(NULL)

  n_obs <- length(index)
  if (length(mu) != (n_obs + 1)) {
    stop("mu must have length = length(index) + 1 (including time 0).")
  }

  mu0 <- as.numeric(mu[1])
  y_obs <- as.numeric(mu[2:(n_obs + 1)])  # observed positions at t_obs

  t_obs <- as.numeric(index) * D
  t_future <- (max(as.numeric(index)) + seq_len(sim_step)) * D

  # Build the unit-scale covariance matrix on observed+future time points.
  t_all <- c(t_obs, t_future)
  M_all <- cov_wfbm_2(t_all, b = beta_hat)

  M_oo <- M_all[1:n_obs, 1:n_obs, drop = FALSE]
  M_of <- M_all[1:n_obs, (n_obs + 1):(n_obs + sim_step), drop = FALSE]
  M_fo <- t(M_of)
  M_ff <- M_all[(n_obs + 1):(n_obs + sim_step), (n_obs + 1):(n_obs + sim_step), drop = FALSE]

  # Jitter for numerical stability
  diag_mean <- mean(diag(M_oo))
  if (!is.finite(diag_mean) || diag_mean <= 0) diag_mean <- 1
  eps <- 1e-10 * diag_mean
  M_oo_j <- M_oo
  diag(M_oo_j) <- diag(M_oo_j) + eps

  diff_obs <- y_obs - rep(mu0, n_obs)

  # Conditional mean: scale cancels, so use M instead of Sigma.
  A <- tryCatch(
    solve(M_oo_j, diff_obs),
    error = function(e) {
      M_tmp <- M_oo_j
      diag(M_tmp) <- diag(M_tmp) + 1e-6 * diag_mean
      solve(M_tmp, diff_obs)
    }
  )
  mean_future <- as.numeric(rep(mu0, sim_step) + M_fo %*% A)

  # Conditional covariance: Sigma = sigma2_hat * M.
  B <- tryCatch(
    solve(M_oo_j, M_of),
    error = function(e) {
      M_tmp <- M_oo_j
      diag(M_tmp) <- diag(M_tmp) + 1e-6 * diag_mean
      solve(M_tmp, M_of)
    }
  )
  cov_base <- M_ff - M_fo %*% B
  cov_base <- (cov_base + t(cov_base)) / 2

  Sigma_future <- as.numeric(sigma2_hat) * cov_base
  Sigma_future <- (Sigma_future + t(Sigma_future)) / 2

  # Clamp numerical negatives
  v <- diag(Sigma_future)
  v[!is.finite(v)] <- NA_real_
  v[v < 0] <- 0
  sd_future <- sqrt(v)

  return(list(
    t_future = t_future,
    mean = mean_future,
    sd = sd_future,
    Sigma = Sigma_future
  ))
}


make_psd <- function(Sigma, eps = 1e-10) {
  # Symmetrize and ensure positive semidefinite by eigenvalue clipping.
  Sigma <- (Sigma + t(Sigma)) / 2
  diag_mean <- mean(diag(Sigma), na.rm = TRUE)
  if (!is.finite(diag_mean) || diag_mean <= 0) diag_mean <- 1
  diag(Sigma) <- diag(Sigma) + eps * diag_mean
  ee <- eigen(Sigma, symmetric = TRUE)
  vals <- ee$values
  vals[vals < 0] <- 0
  Sigma_psd <- ee$vectors %*% diag(vals, length(vals)) %*% t(ee$vectors)
  Sigma_psd <- (Sigma_psd + t(Sigma_psd)) / 2
  return(Sigma_psd)
}

ui <- fluidPage(
  withMathJax(titlePanel(HTML("Process $$\\zeta_{t,f}$$"))),
  div(style = "background-color: #e8f5e9; padding: 15px; border-radius: 8px; margin-bottom: 20px; font-size: 16px; font-family: 'Courier New', monospace;",
      withMathJax(HTML("Covariance function:<br>$$K_f(s,t) := 2 \\int_0^{s \\wedge t} f(r)\\left[(s+t-2r)\\log(s+t-2r) - (s-r)\\log(s-r) - (t-r)\\log(t-r)\\right]dr$$"))
  ),
  tags$script(HTML("
    function toggleBetaControls() {
      var lon0 = $('#beta0_lon').is(':checked');
      var lat0 = $('#beta0_lat').is(':checked');
      $('#beta_low_lon').prop('disabled', lon0);
      $('#beta_high_lon').prop('disabled', lon0);
      $('#beta_low_lat').prop('disabled', lat0);
      $('#beta_high_lat').prop('disabled', lat0);
      $('#grid_block').prop('disabled', lon0 && lat0);
      $('#compute_profile').prop('disabled', lon0 && lat0);
    }
    $(document).on('shiny:connected', function() {
      toggleBetaControls();
      $(document).on('change', '#beta0_lon, #beta0_lat', toggleBetaControls);
    });
  ")),
  tabsetPanel(
    tabPanel("Telemetry Simulation",
      sidebarLayout(
        sidebarPanel(
          numericInput("T", "Total Time (T):", value = 10),
          numericInput("n", "Number of Steps (n):", value = 100),
          numericInput("zeta1_0", HTML("Initial value &zeta;<sub>1</sub>(0):"), value = 0),
          numericInput("zeta2_0", HTML("Initial value &zeta;<sub>2</sub>(0):"), value = 0),
          numericInput("beta1", HTML("&beta;<sub>1</sub>"), value = 0.5, min = 0),
          numericInput("sigma1", HTML("&sigma;<sub>1</sub><sup>2</sup>"), value = 1.2, min = 1e-8),
          numericInput("beta2", HTML("&beta;<sub>2</sub>"), value = 0.45, min = 0),
          numericInput("sigma2", HTML("&sigma;<sub>2</sub><sup>2</sup>"), value = 0.95, min = 1e-8),
          actionButton("simular", "Simulate Trajectory")
          ,
        ),
        mainPanel(
          leafletOutput("mapPlot"),
          plotOutput("flon_flat_plot")
        )
      )
    ),
    tabPanel(withMathJax(HTML("Simulated trajectory $$\\zeta_{t,f}$$")),
      sidebarLayout(
        sidebarPanel(
          textInput("user_expr", "Enter the function f(x):", value = "x*sqrt(x)*(1+sin(x))"),
          numericInput("T_sim", "Total Time T:", value = 10),
          numericInput("n_sim", "Number of Steps n:", value = 100),
          numericInput("mu_i_1", HTML("Initial value $$\\zeta_{0,f}$$:"), value = 0),
          actionButton("evaluar", "Save f(x)"),
          actionButton("simulate_traj", "Simulate Trajectory", class = "btn btn-primary")
        ),
        mainPanel(
          withMathJax(HTML("<h3>Simulated trajectory $$\\zeta_{t,f}$$</h3>")),
          div(style = "background-color: #fdf6e3; padding: 20px; border-radius: 10px; margin-bottom: 20px; font-size: 16px; overflow-x: auto; white-space: normal; word-break: break-word;",
            strong("Condition:"),
            withMathJax(HTML("$$f:\\mathbb{R}_{+}\\to \\mathbb{R}_{+} \\text{ is a measurable function such that, for any } \\delta > 0,$$ $$\\int_0^\\delta f(u)\\, du < \\infty$$"))
          ),
          div(style = "background-color: #e7f0fd; padding: 20px; border-radius: 10px; margin-bottom: 20px; font-size: 16px; overflow-x: auto; white-space: normal; word-break: break-word;",
            uiOutput("output_funcion")
          ),
          plotOutput("traj_plot"),
          plotOutput("output_fx_plot")
        )
      )
    ),
    tabPanel("Inference (Practical Likelihood)",
      sidebarLayout(
        sidebarPanel(
          helpText("Upload the telemetry Excel file (same 3-column structure per animal as in the example table): time index, longitude, latitude."),
          fileInput("data_file", "Telemetry data (Excel .xlsx):", accept = c(".xlsx", ".xls")),
          actionButton("instructions_btn", "Instructions", class = "btn-info"),
          tags$hr(),
          numericInput("animal_number", "Animal number:", value = 4, min = 1, step = 1),
          numericInput("delta", HTML("&Delta; (time step D):"), value = 0.5, min = 1e-6, step = 0.1),

          tags$hr(),

          h4("Beta search intervals / boundary model"),
          checkboxInput("beta0_lon", HTML("Longitude: use boundary model \\(\\beta=0\\)"), value = FALSE),
          checkboxInput("beta0_lat", HTML("Latitude: use boundary model \\(\\beta=0\\)"), value = FALSE),
          fluidRow(
            column(6, numericInput("beta_low_lon", "Longitude: beta lower", value = 0.02, min = 1e-6, step = 0.01)),
            column(6, numericInput("beta_high_lon", "Longitude: beta upper", value = 0.08, min = 1e-6, step = 0.01))
          ),
          fluidRow(
            column(6, numericInput("beta_low_lat", "Latitude: beta lower", value = 0.02, min = 1e-6, step = 0.01)),
            column(6, numericInput("beta_high_lat", "Latitude: beta upper", value = 0.08, min = 1e-6, step = 0.01))
          ),

          tags$hr(),

          numericInput("grid_block", "Profile grid size:", value = 60, min = 10, step = 10),
          checkboxInput("compute_profile", "Compute profile likelihood plots", value = TRUE),

          tags$hr(),
          h4("Prediction"),
          numericInput("Sim_m", "Sim_m (number of simulations):", value = 200, min = 0, step = 50),
          numericInput("sim_step", "sim_step (steps ahead):", value = 20, min = 0, step = 1),
          numericInput("alpha", "alpha (for 1-alpha confidence rectangles):", value = 0.05, min = 1e-6, max = 0.99, step = 0.01),
          helpText("We simulate sim_step future points with step size Delta and display only the mean predicted trajectory."),

          actionButton("run_inference", "Run inference")
        ),
        mainPanel(
          h3("Estimated parameters"),
          uiOutput("inference_table_ui"),
          tags$hr(),
          fluidRow(
            column(6,
              h4("Longitude"),
              plotOutput("profile_lon_plot", height = "300px")
            ),
            column(6,
              h4("Latitude"),
              plotOutput("profile_lat_plot", height = "300px")
            )
          ),
          tags$hr(),
          h4("Observed trajectory + predicted means + (1-alpha) confidence rectangles"),
          leafletOutput("mapPlot_animal", height = "460px")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # ---------------------------
  # Instructions (Inference tab)
  # ---------------------------
  observeEvent(input$instructions_btn, {
    showModal(modalDialog(
      title = "Inference instructions (practical likelihood)",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close")
      ),
      tags$div(
        style = "line-height: 1.35; font-size: 14px;",
        tags$p("This window explains how to run the inference and how to interpret each output shown in the Inference tab."),
        tags$h4("1) What you need"),
        tags$ul(
          tags$li(tags$b("Telemetry Excel file (.xlsx)"), " with the expected 3-column block structure per animal:"),
          tags$ul(
            tags$li(tags$b("Column 1 of the block:"), " time (or time index / day information)"),
            tags$li(tags$b("Column 2 of the block:"), " longitude observations"),
            tags$li(tags$b("Column 3 of the block:"), " latitude observations")
          ),
          tags$li("If the file contains multiple animals, they must be stored as consecutive 3-column blocks: (time, lon, lat) for animal 1, then (time, lon, lat) for animal 2, etc."),
          tags$li("The first row is assumed to correspond to the initial time (time 0). Missing entries should be standard NA/blank cells.")
        ),
        tags$h4("2) Step-by-step procedure"),
        tags$ol(
          tags$li(tags$b("Upload the Excel file"), " using the file selector (no file path typing is required)."),
          tags$li(tags$b("Select the Animal number"), " you want to analyze. The app will read the corresponding (time, lon, lat) block."),
          tags$li(tags$b("Set Δ (Delta)"), " the time step size used to convert the time column into integer indices. The app computes the observation indices as round(time / Δ)."),
          tags$li(tags$b("Choose beta search intervals"), " separately for Longitude and Latitude when fitting the two-parameter model with β>0."),
          tags$li(tags$b("Boundary model β=0"), " may be selected separately for Longitude and Latitude. In that case the app fits the one-parameter model ", tags$span(HTML("f<sub>σ²,0</sub>(u)=σ²")), "."),
          tags$li(tags$b("Profile grid settings (optional)"), ":"),
          tags$ul(
            tags$li(tags$b("Profile grid size"), " controls the resolution of the profile likelihood curve in β."),
            tags$li(tags$b("Compute profile likelihood plots"), " must be enabled if you want the β profile plots.")
          ),
          tags$li(tags$b("Prediction settings"), ":"),
          tags$ul(
            tags$li(tags$b("sim_step"), " = how many steps ahead you want to predict. Each step has size Δ."),
            tags$li(tags$b("alpha"), " sets the confidence level 1-α for the theoretical conditional rectangles."),
            tags$li(tags$b("Sim_m"), " = how many Monte Carlo simulations are used to compute the Monte Carlo conditional mean (yellow). If Sim_m = 0, the yellow curve is not computed.")
          ),
          tags$li(tags$b("Click 'Run inference'"), " and wait until the progress bar completes.")
        ),
        tags$h4("3) What the app computes"),
        tags$ul(
          tags$li(tags$b("β-hat and σ²-hat (MLE)"), " are computed separately for Longitude and Latitude under ", tags$span(HTML("f<sub>σ²,β</sub>(u)=σ² exp(-βu)")), "."),
          tags$li(tags$b("AIC"), " is reported as AIC = 4 - 2 logLik for the two-parameter model and AIC = 2 - 2 logLik for the boundary model β=0."),
          tags$li(tags$b("Theoretical conditional prediction"), " is obtained using Gaussian conditional distributions for future observations given the observed trajectory and the fitted parameters.")
        ),
        tags$h4("4) How to read the outputs"),
        tags$h5("A) Estimated parameters table (top)"),
        tags$ul(
          tags$li(tags$b("β-hat"), ": estimated β for each coordinate."),
          tags$li(tags$b("σ²-hat"), ": estimated scale parameter for each coordinate."),
          tags$li(tags$b("AIC"), ": model selection criterion; smaller is better when comparing models on the same data.")
        ),
        tags$h5("B) Profile likelihood plots (middle)"),
        tags$ul(
          tags$li("Each curve shows the log-likelihood as a function of β over your chosen interval."),
          tags$li("The vertical ", tags$b(tags$span(style="color:red;", "red dashed line")), " marks the MLE β-hat.")
        ),
        tags$h5("C) Map (bottom)"),
        tags$p("The map shows the observed path and the fitted/predicted components. The legend at the bottom-right summarizes the colors:"),
        tags$ul(
          tags$li(tags$b("Black"), ": observed trajectory from the uploaded Excel file (selected animal)."),
          tags$li(tags$b(tags$span(style="color:red;", "Red")), ": theoretical conditional mean path (future) based on Gaussian conditional expectation."),
          tags$li(tags$b(tags$span(style="color:blue;", "Blue")), ": theoretical conditional (1-α) confidence rectangles at each future step (marginal intervals for lon and lat)."),
          tags$li(tags$b(tags$span(style="color:gold;", "Yellow")), ": Monte Carlo conditional mean (average of Sim_m future simulations)."),
          tags$li("The map also marks the end of the theoretical mean and the end of the Monte Carlo mean with labels.")
        ),
        tags$h4("5) Common issues and fixes"),
        tags$ul(
          tags$li(tags$b("No plot / no map"), ": make sure you uploaded a valid Excel file and clicked 'Run inference'."),
          tags$li(tags$b("Animal number exceeds columns"), ": your file may not have that many 3-column blocks. Reduce Animal number."),
          tags$li(tags$b("Time indices not strictly increasing"), ": check the time column or adjust Δ so that round(time/Δ) increases strictly."),
          tags$li(tags$b("Numerical errors in covariance inversion"), ": try a smaller sim_step or adjust the β interval (sometimes extremely small/large β can cause instability).")
        )
      )
    ))
  })

  funcion_guardada <- reactiveVal(NULL)
  latex_guardado <- reactiveVal("")
  inference_res <- reactiveVal(NULL)

  observeEvent(input$run_inference, {
    req(input$data_file)
    validate(
      need(isTRUE(input$beta0_lon) || input$beta_low_lon < input$beta_high_lon, "Longitude interval: beta lower must be < beta upper."),
      need(isTRUE(input$beta0_lat) || input$beta_low_lat < input$beta_high_lat, "Latitude interval: beta lower must be < beta upper."),
      need(isTRUE(input$beta0_lon) || input$beta_low_lon > 0, "Longitude beta lower must be positive unless the boundary model beta=0 is selected."),
      need(isTRUE(input$beta0_lat) || input$beta_low_lat > 0, "Latitude beta lower must be positive unless the boundary model beta=0 is selected.")
    )

    withProgress(message = "Running inference...", value = 0, {
      incProgress(0.05, detail = "Reading data")

      datos <- data.frame(read_excel(input$data_file$datapath, col_names = FALSE))
      number <- as.integer(input$animal_number)

      col_day <- 3 * (number - 1) + 1
      col_lon <- 3 * (number - 1) + 2
      col_lat <- 3 * number

      validate(
        need(ncol(datos) >= col_lat, "The selected animal number exceeds the number of columns in the uploaded file.")
      )

      lon_col <- datos[[col_lon]]
      ok <- !(is.na(lon_col) | lon_col == "NA" | lon_col == "")
      validate(need(any(ok), "No valid longitude data found for the selected animal number."))

      k <- max(which(ok))

      lon <- as.numeric(lon_col[1:k])
      lat <- as.numeric(datos[[col_lat]][1:k])

      D <- as.numeric(input$delta)
      validate(need(D > 0, "Delta (D) must be positive."))

      days_raw <- as.numeric(datos[[col_day]][2:k])  # exclude time 0
      validate(need(all(!is.na(days_raw)), "Time index column contains missing values for the selected animal."))

      index <- as.integer(round(days_raw / D))
      validate(need(length(index) == (k - 1), "Index length mismatch (check the input file format)."))
      validate(need(all(index >= 1), "All inferred indices must be >= 1 (check Delta and time column)."))
      validate(need(all(diff(index) > 0), "Time indices must be strictly increasing (check the time column and Delta)."))

      T <- index[length(index)] * D

      incProgress(0.15, detail = "Estimating longitude parameters")
      res_lon <- infer_beta_sigma(
        T = T, mu = lon, index = index,
        beta_lower = input$beta_low_lon, beta_upper = input$beta_high_lon,
        grid_n = input$grid_block, compute_profile = input$compute_profile,
        use_beta0 = input$beta0_lon
      )

      incProgress(0.55, detail = "Estimating latitude parameters")
      res_lat <- infer_beta_sigma(
        T = T, mu = lat, index = index,
        beta_lower = input$beta_low_lat, beta_upper = input$beta_high_lat,
        grid_n = input$grid_block, compute_profile = input$compute_profile,
        use_beta0 = input$beta0_lat
      )

      # Prediction (future conditional mean + confidence bounds)
      sim_step <- as.integer(input$sim_step)
      Sim_m <- as.integer(input$Sim_m)  # kept for UI compatibility (not required for the mean/CIs)
      alpha <- as.numeric(input$alpha)
      validate(need(is.finite(alpha) && alpha > 0 && alpha < 1, "alpha must be in (0,1)."))

      pred_lon <- NULL
      pred_lat <- NULL
      if (sim_step > 0) {
        incProgress(0.80, detail = "Computing conditional mean + confidence rectangles")
        pred_lon <- tryCatch(
          predict_future_conditional(mu = lon, index = index, D = D,
                                     beta_hat = res_lon$beta_hat, sigma2_hat = res_lon$sigma2_hat,
                                     sim_step = sim_step),
          error = function(e) NULL
        )
        pred_lat <- tryCatch(
          predict_future_conditional(mu = lat, index = index, D = D,
                                     beta_hat = res_lat$beta_hat, sigma2_hat = res_lat$sigma2_hat,
                                     sim_step = sim_step),
          error = function(e) NULL
        )

        # Add marginal (1-alpha) intervals per future step: mean +/- z * sd
        z <- qnorm(1 - alpha / 2)
        if (!is.null(pred_lon) && !is.null(pred_lon$mean) && !is.null(pred_lon$sd)) {
          pred_lon$lower <- pred_lon$mean - z * pred_lon$sd
          pred_lon$upper <- pred_lon$mean + z * pred_lon$sd
        }
        if (!is.null(pred_lat) && !is.null(pred_lat$mean) && !is.null(pred_lat$sd)) {
          pred_lat$lower <- pred_lat$mean - z * pred_lat$sd
          pred_lat$upper <- pred_lat$mean + z * pred_lat$sd

        # Monte Carlo predicted mean (Sim_m simulations) from the conditional law
        mc_lon_mean <- NULL
        mc_lat_mean <- NULL
        if (Sim_m > 0 && !is.null(pred_lon) && !is.null(pred_lat) &&
            !is.null(pred_lon$Sigma) && !is.null(pred_lat$Sigma)) {

          Sigma_lon <- tryCatch(make_psd(pred_lon$Sigma), error = function(e) NULL)
          Sigma_lat <- tryCatch(make_psd(pred_lat$Sigma), error = function(e) NULL)

          if (!is.null(Sigma_lon)) {
            sims_lon <- tryCatch(
              rmvnorm(Sim_m, mean = pred_lon$mean, sigma = Sigma_lon),
              error = function(e) NULL
            )
            if (!is.null(sims_lon)) mc_lon_mean <- as.numeric(colMeans(sims_lon))
          }

          if (!is.null(Sigma_lat)) {
            sims_lat <- tryCatch(
              rmvnorm(Sim_m, mean = pred_lat$mean, sigma = Sigma_lat),
              error = function(e) NULL
            )
            if (!is.null(sims_lat)) mc_lat_mean <- as.numeric(colMeans(sims_lat))
          }

          if (!is.null(mc_lon_mean)) pred_lon$mc_mean <- mc_lon_mean
          if (!is.null(mc_lat_mean)) pred_lat$mc_mean <- mc_lat_mean
        }
        }
      }

      inference_res(list(
        T = T, D = D, k = k,
        lon = list(mu = lon, res = res_lon),
        lat = list(mu = lat, res = res_lat),
        pred = list(sim_step = sim_step, Sim_m = Sim_m, alpha = alpha, lon = pred_lon, lat = pred_lat)
      ))

      incProgress(1, detail = "Done")
    })
  })

  output$inference_table_ui <- renderUI({
    r <- inference_res()
    req(r)

    fmt <- function(x) formatC(x, digits = 6, format = "f")
    fmt_beta <- function(x, fixed) {
      if (isTRUE(fixed)) return("0 (fixed)")
      formatC(x, digits = 6, format = "f")
    }

    rows <- list(
      list("Longitude", r$lon$res$model_label, fmt_beta(r$lon$res$beta_hat, r$lon$res$beta_fixed),
           r$lon$res$sigma2_hat, r$lon$res$logLik_hat, r$lon$res$aic_hat),
      list("Latitude", r$lat$res$model_label, fmt_beta(r$lat$res$beta_hat, r$lat$res$beta_fixed),
           r$lat$res$sigma2_hat, r$lat$res$logLik_hat, r$lat$res$aic_hat)
    )

    tbl <- tags$table(
      class = "table table-striped",
      tags$thead(
        tags$tr(
          tags$th("Coordinate"),
          tags$th("Model"),
          tags$th(HTML("\\(\\hat{\\beta}\\)")),
          tags$th(HTML("\\(\\widehat{\\sigma}^{\\,2}\\)")),
          tags$th("logLik"),
          tags$th("AIC")
        )
      ),
      tags$tbody(
        lapply(rows, function(rr) {
          tags$tr(
            tags$td(rr[[1]]),
            tags$td(rr[[2]]),
            tags$td(rr[[3]]),
            tags$td(fmt(rr[[4]])),
            tags$td(fmt(rr[[5]])),
            tags$td(fmt(rr[[6]]))
          )
        })
      )
    )

    withMathJax(tagList(
      tbl,
      tags$script(HTML("if (window.MathJax) { if (MathJax.typesetPromise) { MathJax.typesetPromise(); } else if (MathJax.Hub && MathJax.Hub.Queue) { MathJax.Hub.Queue(['Typeset', MathJax.Hub]); } }"))
    ))
  })

output$profile_lon_plot <- renderPlot({
    r <- inference_res()
    req(r)

    if (is.null(r$lon$res$grid_beta)) {
      plot.new()
      msg <- if (isTRUE(r$lon$res$beta_fixed)) "Boundary model beta = 0: no beta profile" else "Profile likelihood plot was not computed"
      text(0.5, 0.5, msg, cex = 1.1)
      return(invisible(NULL))
    }

    plot(r$lon$res$grid_beta, r$lon$res$profile_ll, type = "l",
         xlab = expression(beta), ylab = "profile log-likelihood",
         main = "Profile Likelihood (Longitude)")
    abline(v = r$lon$res$beta_hat, lty = 2, col = "red", lwd = 2)
    grid(lty = 3, col = "gray80")
    box()
  })

  output$profile_lat_plot <- renderPlot({
    r <- inference_res()
    req(r)

    if (is.null(r$lat$res$grid_beta)) {
      plot.new()
      msg <- if (isTRUE(r$lat$res$beta_fixed)) "Boundary model beta = 0: no beta profile" else "Profile likelihood plot was not computed"
      text(0.5, 0.5, msg, cex = 1.1)
      return(invisible(NULL))
    }

    plot(r$lat$res$grid_beta, r$lat$res$profile_ll, type = "l",
         xlab = expression(beta), ylab = "profile log-likelihood",
         main = "Profile Likelihood (Latitude)")
    abline(v = r$lat$res$beta_hat, lty = 2, col = "red", lwd = 2)
    grid(lty = 3, col = "gray80")
    box()
  })

  observeEvent(input$evaluar, {
    req(input$user_expr)
    f_string <- paste0("function(x) { ", input$user_expr, " }")
    f <- eval(parse(text = f_string))
    funcion_guardada(f)
    latex_guardado(toLatex(input$user_expr))
  })

  output$output_funcion <- renderUI({
    req(latex_guardado() != "")
    withMathJax(HTML(paste0("$$\\text{Saved function:}\\quad f(x) = ", latex_guardado(), "$$")))
  })

  observeEvent(input$simular, {
    T <- input$T
    n <- input$n
    D <- T / n
    t <- seq(D, T, length.out = n)

    f_lon <- make_exp_weight(input$sigma1, input$beta1)
    f_lat <- make_exp_weight(input$sigma2, input$beta2)

    x_vals <- seq(0, 50, by = 0.1)
    y_flon <- f_lon(x_vals)
    y_flat <- f_lat(x_vals)

    K1 <- cov_f_weighted(f_lon, t)
    K2 <- cov_f_weighted(f_lat, t)

    mu1 <- rep(input$zeta1_0, n)
    mu2 <- rep(input$zeta2_0, n)

    lon <- mvrnorm(1, mu1, K1)
    lat <- mvrnorm(1, mu2, K2)

    df <- data.frame(
      lon = wrap_longitude(lon),
      lat = wrap_latitude(lat)
    )

    output$mapPlot <- renderLeaflet({
      leaflet(df) %>%
        addTiles() %>%
        addPolylines(~lon, ~lat, color = "#2c3e50", weight = 3) %>%
        addCircleMarkers(~lon, ~lat, radius = 3, color = "#2c3e50", fillOpacity = 0.9) %>%
        addMarkers(lng = df$lon[1], lat = df$lat[1], label = "Start", labelOptions = labelOptions(noHide = TRUE)) %>%
        addMarkers(lng = df$lon[n], lat = df$lat[n], label = "End", labelOptions = labelOptions(noHide = TRUE))
    })

    output$flon_flat_plot <- renderPlot({
      plot(x_vals, y_flon, type = "l", col = "#1f77b4", lwd = 3, ylim = range(c(y_flon, y_flat)),
           xlab = "u", ylab = "f(u)", main = "Covariance Kernel Functions",
           cex.lab = 1.4, cex.main = 1.6)
      lines(x_vals, y_flat, col = "#d62728", lwd = 3)
      legend("topright",
             legend = c(expression(f[lon](u)), expression(f[lat](u))),
             col = c("#1f77b4", "#d62728"), lty = 1, lwd = 3, bty = "n")
      grid(lty = 3, col = "gray80")
      box()
    })
  })

  observeEvent(input$simulate_traj, {
    req(funcion_guardada())
    f <- funcion_guardada()
    T <- input$T_sim
    n <- input$n_sim
    D <- T / n
    t <- seq(D, T, length.out = n)
    K <- cov_f_weighted(f, t)
    mu_i_1 <- input$mu_i_1
    mu <- rep(mu_i_1, n)
    sim <- mvrnorm(1, mu, K)

    output$traj_plot <- renderPlot({
      plot(c(0, t), c(mu_i_1, sim), type = "l", lwd = 3, col = "#2c3e50",
           xlab = "t", ylab = expression(zeta["t,f"]),
           main = expression("Simulated trajectory " * zeta["t,f"]),
           cex.lab = 1.6, cex.main = 1.8, las = 1)
      points(0, mu_i_1, pch = 16, col = "green", cex = 1.5)
      points(t[n], sim[n], pch = 16, col = "red", cex = 1.5)
      legend("topleft", legend = c("Start", "End"), col = c("green", "red"), pch = 16, bty = "n")
      grid(lty = 3, col = "gray70")
      box()
    })

    output$output_fx_plot <- renderPlot({
      x_vals <- seq(0, 50, by = 0.1)
      y_vals <- f(x_vals)
      plot(x_vals, y_vals, type = "l", col = "#1b7837", lwd = 2,
           xlab = "x", ylab = "f(x)", cex.lab = 1.4)
      grid(lty = 3, col = "gray80")
      box()
    })
  })
  # ---------------------------
  # Map of selected animal trajectory (Tab 3 - Inference)
  # ---------------------------
  output$mapPlot_animal <- renderLeaflet({
    r <- inference_res()
    req(r)

    lon_obs <- r$lon$mu
    lat_obs <- r$lat$mu

    df_obs <- data.frame(
      lon = wrap_longitude(as.numeric(lon_obs)),
      lat = wrap_latitude(as.numeric(lat_obs))
    )
    df_obs <- df_obs[complete.cases(df_obs), , drop = FALSE]
    validate(need(nrow(df_obs) >= 2, "Not enough valid (lon, lat) points to draw the trajectory."))

    m <- leaflet(df_obs) %>%
      addTiles() %>%
      addPolylines(~lon, ~lat, weight = 3, color = "#2c3e50", opacity = 0.9) %>%
      addCircleMarkers(~lon, ~lat, radius = 3, color = "#2c3e50", fillOpacity = 0.9) %>%
      addMarkers(lng = df_obs$lon[1], lat = df_obs$lat[1], label = "Start",
                 labelOptions = labelOptions(noHide = TRUE)) %>%
      addMarkers(lng = df_obs$lon[nrow(df_obs)], lat = df_obs$lat[nrow(df_obs)], label = "End (observed)",
                 labelOptions = labelOptions(noHide = TRUE))

    # Add mean predicted continuation (conditional mean) + confidence rectangles
    if (!is.null(r$pred) && !is.null(r$pred$lon) && !is.null(r$pred$lat)) {
      lon_mean <- r$pred$lon$mean
      lat_mean <- r$pred$lat$mean

      if (is.numeric(lon_mean) && is.numeric(lat_mean) &&
          length(lon_mean) == length(lat_mean) &&
          length(lon_mean) >= 1 &&
          all(is.finite(lon_mean)) && all(is.finite(lat_mean))) {

        # Mean predicted points (red), starting from the last observed point
        df_pred_pts <- data.frame(
          lon = wrap_longitude(as.numeric(lon_mean)),
          lat = wrap_latitude(as.numeric(lat_mean))
        )
        df_pred_pts <- df_pred_pts[complete.cases(df_pred_pts), , drop = FALSE]

        # Polyline from last observed to predicted mean path
        df_pred_line <- data.frame(
          lon = wrap_longitude(c(df_obs$lon[nrow(df_obs)], as.numeric(lon_mean))),
          lat = wrap_latitude(c(df_obs$lat[nrow(df_obs)], as.numeric(lat_mean)))
        )
        df_pred_line <- df_pred_line[complete.cases(df_pred_line), , drop = FALSE]

        if (nrow(df_pred_line) >= 2) {
          m <- m %>%
            addPolylines(data = df_pred_line, ~lon, ~lat, weight = 4, color = "red", opacity = 0.9)
        }

        if (nrow(df_pred_pts) >= 1) {
          m <- m %>%
            addCircleMarkers(data = df_pred_pts, ~lon, ~lat,
                             radius = 5, color = "red", fillOpacity = 0.9)

        # Monte Carlo predicted mean (yellow), if available
        if (!is.null(r$pred$lon$mc_mean) && !is.null(r$pred$lat$mc_mean)) {
          lon_mc <- r$pred$lon$mc_mean
          lat_mc <- r$pred$lat$mc_mean

          if (is.numeric(lon_mc) && is.numeric(lat_mc) &&
              length(lon_mc) == length(lat_mc) &&
              length(lon_mc) >= 1 &&
              all(is.finite(lon_mc)) && all(is.finite(lat_mc))) {

            df_mc_pts <- data.frame(
              lon = wrap_longitude(as.numeric(lon_mc)),
              lat = wrap_latitude(as.numeric(lat_mc))
            )
            df_mc_pts <- df_mc_pts[complete.cases(df_mc_pts), , drop = FALSE]

            df_mc_line <- data.frame(
              lon = wrap_longitude(c(df_obs$lon[nrow(df_obs)], as.numeric(lon_mc))),
              lat = wrap_latitude(c(df_obs$lat[nrow(df_obs)], as.numeric(lat_mc)))
            )
            df_mc_line <- df_mc_line[complete.cases(df_mc_line), , drop = FALSE]

            if (nrow(df_mc_line) >= 2) {
              m <- m %>%
                addPolylines(data = df_mc_line, ~lon, ~lat, weight = 4,
                            color = "yellow", opacity = 0.9)
            }
            if (nrow(df_mc_pts) >= 1) {
              m <- m %>%
                addCircleMarkers(data = df_mc_pts, ~lon, ~lat,
                                 radius = 5, color = "yellow", fillOpacity = 0.9)
            # Mark end of Monte Carlo mean trajectory
            m <- m %>%
              addMarkers(lng = df_mc_pts$lon[nrow(df_mc_pts)], lat = df_mc_pts$lat[nrow(df_mc_pts)],
                         label = "End (Monte Carlo mean)",
                         labelOptions = labelOptions(noHide = TRUE))
            }
          }
        }

        }

        # Confidence rectangles per future step (blue), based on marginal intervals
        if (!is.null(r$pred$lon$lower) && !is.null(r$pred$lon$upper) &&
            !is.null(r$pred$lat$lower) && !is.null(r$pred$lat$upper)) {

          lon_lo <- wrap_longitude(as.numeric(r$pred$lon$lower))
          lon_hi <- wrap_longitude(as.numeric(r$pred$lon$upper))
          lat_lo <- wrap_latitude(as.numeric(r$pred$lat$lower))
          lat_hi <- wrap_latitude(as.numeric(r$pred$lat$upper))

          # Build rectangles; ensure bounds order
          lng1 <- pmin(lon_lo, lon_hi)
          lng2 <- pmax(lon_lo, lon_hi)
          lt1  <- pmin(lat_lo, lat_hi)
          lt2  <- pmax(lat_lo, lat_hi)

          okr <- is.finite(lng1) & is.finite(lng2) & is.finite(lt1) & is.finite(lt2)
          if (any(okr)) {
            steps <- which(okr)
            rect_df <- data.frame(
              lng1 = lng1[okr],
              lng2 = lng2[okr],
              lat1 = lt1[okr],
              lat2 = lt2[okr],
              step = steps
            )
            lab <- paste0("Step +", rect_df$step, " (1-α CI, α=", round(r$pred$alpha, 4), ")")
            m <- m %>%
              addRectangles(data = rect_df,
                            lng1 = ~lng1, lat1 = ~lat1, lng2 = ~lng2, lat2 = ~lat2,
                            fill = FALSE, color = "blue", weight = 2, opacity = 0.8,
                            label = lab)
          }
        }

        # Mark last predicted mean point
        if (nrow(df_pred_pts) >= 1) {
          m <- m %>%
            addMarkers(lng = df_pred_pts$lon[nrow(df_pred_pts)], lat = df_pred_pts$lat[nrow(df_pred_pts)],
                       label = "End (predicted mean)",
                       labelOptions = labelOptions(noHide = TRUE))
        }
      }
    }
    # Legend / labels
    alpha_txt <- if (!is.null(r$pred) && !is.null(r$pred$alpha)) round(r$pred$alpha, 4) else NA
    Sim_m_txt <- if (!is.null(r$pred) && !is.null(r$pred$Sim_m)) r$pred$Sim_m else NA
    m <- m %>% addControl(
      html = paste0(
        "<div style='background: rgba(255,255,255,0.95); padding:10px; border-radius:8px; font-size:13px;'>",
        "<b>Legend</b><br>",
        "<span style='color:black; font-weight:600;'>Black</span>: observed trajectory<br>",
        "<span style='color:red; font-weight:600;'>Red</span>: theoretical conditional mean<br>",
        "<span style='color:blue; font-weight:600;'>Blue</span>: theoretical (1-α) confidence rectangles (α=", alpha_txt, ")<br>",
        "<span style='color:gold; font-weight:600;'>Yellow</span>: Monte Carlo predicted mean (Sim_m=", Sim_m_txt, ")",
        "</div>"
      ),
      position = "bottomright"
    )


    m
  })

}

shinyApp(ui, server)
