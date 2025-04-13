library(shiny)
library(MASS)
library(leaflet)

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

# Covariance kernel construction
cov_f_log <- function(f, u, s, t) {
  return(2 * f(u) * (-(t - u) * log(t - u) - (s - u) * log(s - u) + (t + s - 2 * u) * log(t + s - 2 * u)))
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

ui <- fluidPage(
  withMathJax(titlePanel(HTML("Process $$\\zeta_{t,f}$$"))),
  div(style = "background-color: #e8f5e9; padding: 15px; border-radius: 8px; margin-bottom: 20px; font-size: 16px; font-family: 'Courier New', monospace;",
      withMathJax(HTML("Covariance function:<br>$$K_f(s,t) := 2 \\int_0^{s \\wedge t} f(r)\\left[(s+t-2r)\\log(s+t-2r) - (s-r)\\log(s-r) - (t-r)\\log(t-r)\\right]dr$$"))
  ),
  tabsetPanel(
    tabPanel("Telemetry Simulation",
      sidebarLayout(
        sidebarPanel(
          numericInput("T", "Total Time (T):", value = 10),
          numericInput("n", "Number of Steps (n):", value = 100),
          numericInput("zeta1_0", HTML("Initial value &zeta;<sub>1</sub>(0):"), value = 0),
          numericInput("zeta2_0", HTML("Initial value &zeta;<sub>2</sub>(0):"), value = 0),
          numericInput("beta1", HTML("&beta;<sub>1</sub>"), value = 0.5),
          numericInput("sigma1", HTML("&sigma;<sub>1</sub>"), value = 1.2),
          numericInput("beta2", HTML("&beta;<sub>2</sub>"), value = 0.45),
          numericInput("sigma2", HTML("&sigma;<sub>2</sub>"), value = 0.95),
          actionButton("simular", "Simulate Trajectory")
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
    )
  )
)

server <- function(input, output, session) {
  funcion_guardada <- reactiveVal(NULL)
  latex_guardado <- reactiveVal("")

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

    f_lon <- function(u) (input$sigma1^2 / (2 * input$beta1)) * exp(-input$beta1 * u)
    f_lat <- function(u) (input$sigma2^2 / (2 * input$beta2)) * exp(-input$beta2 * u)

    x_vals <- seq(0, 50, by = 0.1)
    y_flon <- f_lon(x_vals)
    y_flat <- f_lat(x_vals)

    K1 <- cov_f_weighted(f_lon, t)
    K2 <- cov_f_weighted(f_lat, t)

    mu1 <- rep(input$zeta1_0, n)
    mu2 <- rep(input$zeta2_0, n)

    lon <- mvrnorm(1, mu1, K1)
    lat <- mvrnorm(1, mu2, K2)

    df <- data.frame(lon = lon, lat = lat)

    output$mapPlot <- renderLeaflet({
      leaflet(df) %>%
        addTiles() %>%
        addPolylines(~lon, ~lat, color = "#2c3e50", weight = 3) %>%
        addCircleMarkers(~lon, ~lat, radius = 3, color = "#2c3e50", fillOpacity = 0.9) %>%
        addMarkers(lng = lon[1], lat = lat[1], label = "Start", labelOptions = labelOptions(noHide = TRUE)) %>%
        addMarkers(lng = lon[n], lat = lat[n], label = "End", labelOptions = labelOptions(noHide = TRUE))
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
}

shinyApp(ui, server)
