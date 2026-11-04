# Prototype created by Steven Van Aken
# All demand, supplier volume, capacity and cost figures must use the
# same planning period. Shipment sizes and FTL use the same quantity unit.
#
# The Shiny process receives uploaded data. Unlike a browser-only tool,
# deploying this app remotely sends the data to the Shiny server.
setwd
library(shiny)
library(leaflet)
library(DT)
library(geosphere)
library(lpSolve)

# ---------------------------------------------------------------------
# Demonstration inputs: replace with your own data.
# Warehouse_IDs: blank = any warehouse; otherwise IDs separated by "|".
# Move_range_km: 0 = fixed; Inf = unrestricted.
# Supplier volumes are allocated in proportion to warehouse throughput.
# ---------------------------------------------------------------------

demo_customers <- data.frame(
  Customer_ID = c("Brussels", "Antwerp", "Ghent", "Liege", "Lille", "Eindhoven"),
  Latitude = c(50.8503, 51.2194, 51.0543, 50.6326, 50.6292, 51.4416),
  Longitude = c(4.3517, 4.4025, 3.7174, 5.5797, 3.0573, 5.4697),
  Demand = c(1000, 1400, 800, 700, 900, 1100),
  Warehouse_IDs = "",
  Country = c("BE", "BE", "BE", "BE", "FR", "NL"),
  Group = "",
  Shipment_size = c(10, 15, 8, 10, 12, 14),
  stringsAsFactors = FALSE
)

demo_warehouses <- data.frame(
  Warehouse_ID = "Antwerp_DC",
  Latitude = 51.2194,
  Longitude = 4.4025,
  Move_range_km = 0,
  Capacity = Inf,
  stringsAsFactors = FALSE
)

demo_suppliers <- data.frame(
  Supplier_ID = "Liege",
  Latitude = 50.6326,
  Longitude = 5.5797,
  Volume = 5900,
  Shipment_size = 24,
  stringsAsFactors = FALSE
)

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

distance_matrix <- function(a, b, circuity = 1) {
  geosphere::distm(
    as.matrix(a[, c("Longitude", "Latitude")]),
    as.matrix(b[, c("Longitude", "Latitude")]),
    fun = geosphere::distHaversine
  ) / 1000 * circuity
}

# Cost fraction for a shipment occupying x of a full truck:
# cubic Bezier with fixed endpoints (0,0), (1,1).
# Control points must be monotonic for stable inversion.
curve_fraction <- function(x, p) {
  t <- seq(0, 1, length.out = 2001)
  bx <- 3 * (1-t)^2 * t * p$x1 +
    3 * (1-t) * t^2 * p$x2 + t^3
  by <- 3 * (1-t)^2 * t * p$y1 +
    3 * (1-t) * t^2 * p$y2 + t^3
  approx(bx, by, xout = pmin(1, pmax(0, x)), ties = "ordered")$y
}

# Cost per quantity per km.
# Shipment sizes above FTL are rejected during validation.
unit_rate <- function(size, multiplier, p) {
  p$ftl_rate * multiplier * curve_fraction(size / p$ftl, p) / size
}

prepare_data <- function(c, w, s, p) {
  for (nm in c("Latitude", "Longitude", "Demand", "Shipment_size"))
    c[[nm]] <- num(c[[nm]])
  for (nm in c("Latitude", "Longitude", "Move_range_km", "Capacity"))
    w[[nm]] <- num(w[[nm]])
  for (nm in c("Latitude", "Longitude", "Volume", "Shipment_size"))
    s[[nm]] <- num(s[[nm]])
  
  assert(nrow(c) > 0, "At least one customer is required.")
  assert(!anyDuplicated(c$Customer_ID), "Customer IDs must be unique.")
  assert(!anyDuplicated(w$Warehouse_ID), "Warehouse IDs must be unique.")
  assert(all(is.finite(c$Demand) & c$Demand > 0),
         "Customer demand must be finite and positive.")
  
  check_coords <- function(d) {
    all(is.finite(d$Latitude) & abs(d$Latitude) <= 85 &
          is.finite(d$Longitude) & abs(d$Longitude) <= 180)
  }
  assert(check_coords(c) && check_coords(w) && check_coords(s),
         "Invalid coordinates. This prototype supports latitudes -85 to 85.")
  assert(all(!is.na(w$Capacity) & w$Capacity > 0),
         "Warehouse capacities must be positive or Inf.")
  assert(all(!is.na(w$Move_range_km) & w$Move_range_km >= 0),
         "Move ranges must be nonnegative or Inf.")
  assert(p$k >= nrow(w), "Warehouse count is below predefined count.")
  assert(p$x1 <= p$x2 && p$y1 <= p$y2,
         "Bezier controls must satisfy X1 <= X2 and Y1 <= Y2.")
  
  if (p$auto_ship) c$Shipment_size <- p$customer_size
  assert(all(is.finite(c$Shipment_size) &
               c$Shipment_size > 0 & c$Shipment_size <= p$ftl),
         "Customer shipment sizes must be > 0 and <= FTL.")
  
  if (p$suppliers) {
    assert(nrow(s) > 0, "Supply constraint is enabled but suppliers are empty.")
    assert(all(is.finite(s$Volume) & s$Volume > 0),
           "Supplier volumes must be finite and positive.")
    assert(all(is.finite(s$Shipment_size) &
                 s$Shipment_size > 0 & s$Shipment_size <= p$ftl),
           "Supplier shipment sizes must be > 0 and <= FTL.")
  }
  
  set.seed(p$seed)
  if (p$random_pct > 0) {
    a <- p$random_pct / 100
    factor <- if (p$distribution == "Uniform") {
      runif(nrow(c), 1-a, 1+a)
    } else {
      # Truncated normal: +/- percentage is the allowed range.
      pmin(1+a, pmax(1-a, rnorm(nrow(c), 1, a/3)))
    }
    c$Demand <- c$Demand * factor
  }
  
  # Only observed points define the free-location bounding box.
  # User-entered predefined warehouses can remain outside that box.
  list(customers = c, warehouses = w, suppliers = s)
}

# ---------------------------------------------------------------------
# Solver: binary customer/group assignment + continuous relocation.
# This is a multi-start alternating heuristic, not a global certificate.
# Capacity-constrained assignments are solved as binary programs.
# ---------------------------------------------------------------------

solve_network <- function(c, w, s, p) {
  set.seed(p$seed)
  
  groups <- if (p$groups) {
    ifelse(is.na(c$Group) | trimws(c$Group) == "",
           paste0("customer:", c$Customer_ID),
           paste0("group:", c$Group))
  } else paste0("customer:", c$Customer_ID)
  
  members <- split(seq_len(nrow(c)), groups)
  gdemand <- vapply(members, function(ix) sum(c$Demand[ix]), numeric(1))
  
  all_ids <- c(w$Warehouse_ID,
               if (p$k > nrow(w))
                 paste0("CoG_", seq_len(p$k - nrow(w))) else character())
  assert(!anyDuplicated(all_ids), "Generated and predefined IDs overlap.")
  
  allowed <- matrix(TRUE, nrow(c), p$k)
  if (p$assignments) {
    for (i in seq_len(nrow(c))) {
      text <- c$Warehouse_IDs[i]
      if (!is.na(text) && nzchar(trimws(text))) {
        ids <- trimws(strsplit(text, "|", fixed = TRUE)[[1]])
        assert(all(ids %in% all_ids),
               paste("Unknown Warehouse_IDs for", c$Customer_ID[i]))
        allowed[i, ] <- all_ids %in% ids
      }
    }
  }
  
  group_allowed <- do.call(rbind, lapply(members, function(ix) {
    apply(allowed[ix, , drop = FALSE], 2, all)
  }))
  assert(all(rowSums(group_allowed) > 0),
         "A customer/group has no commonly allowed warehouse.")
  
  cap <- c(w$Capacity, rep(Inf, p$k - nrow(w)))
  customer_rate <- unit_rate(c$Shipment_size, p$customer_ratio, p)
  
  supplier_rate <- if (p$suppliers)
    unit_rate(s$Shipment_size, p$supplier_ratio, p) else numeric()
  
  arc_cost <- function(loc) {
    d <- distance_matrix(c, loc, p$circuity)
    costs <- sweep(d, 1, c$Demand * customer_rate, "*")
    if (p$improve_service)
      costs <- costs + sweep(pmax(d - p$lead_km, 0), 1,
                             c$Demand * p$service_penalty, "*")
    
    if (p$suppliers) {
      inbound <- colSums(sweep(
        distance_matrix(s, loc, p$circuity), 1,
        s$Volume * supplier_rate, "*"
      )) / sum(c$Demand)
      costs <- costs + outer(c$Demand, inbound)
    }
    costs
  }
  
  assign_customers <- function(loc) {
    costs <- arc_cost(loc)
    gcost <- do.call(rbind, lapply(members, function(ix)
      colSums(costs[ix, , drop = FALSE])))
    ng <- length(members)
    nv <- ng * p$k
    
    # Variable order: group changes fastest, warehouse slowest.
    eq <- matrix(0, ng, nv)
    for (g in seq_len(ng))
      eq[g, g + (seq_len(p$k)-1L)*ng] <- 1
    
    mat <- eq
    dirs <- rep("=", ng)
    rhs <- rep(1, ng)
    
    if (p$capacities) {
      for (j in which(is.finite(cap))) {
        row <- numeric(nv)
        row[seq_len(ng) + (j-1L)*ng] <- gdemand
        mat <- rbind(mat, row)
        dirs <- c(dirs, "<=")
        rhs <- c(rhs, cap[j])
      }
    }
    
    forbidden <- which(!as.vector(group_allowed))
    if (length(forbidden)) {
      row <- numeric(nv)
      row[forbidden] <- 1
      mat <- rbind(mat, row)
      dirs <- c(dirs, "=")
      rhs <- c(rhs, 0)
    }
    
    lp <- lpSolve::lp("min", as.vector(gcost), mat, dirs, rhs,
                      all.bin = TRUE)
    assert(lp$status == 0,
           "No feasible assignment, or assignment solver failed.")
    z <- matrix(lp$solution, ng, p$k)
    ga <- max.col(z, ties.method = "first")
    a <- integer(nrow(c))
    for (g in seq_along(members)) a[members[[g]]] <- ga[g]
    a
  }
  
  lower <- c(min(c$Longitude), min(c$Latitude))
  upper <- c(max(c$Longitude), max(c$Latitude))
  # Avoid zero-width optim() bounds.
  upper <- pmax(upper, lower + 1e-5)
  
  enforce_radius <- function(x, j) {
    if (j > nrow(w) || !is.finite(w$Move_range_km[j])) return(x)
    origin <- c(w$Longitude[j], w$Latitude[j])
    radius <- w$Move_range_km[j]
    if (radius == 0) return(origin)
    if (distHaversine(origin, x)/1000 <= radius) return(x)
    destPoint(origin, bearing(origin, x), radius * 1000)[1, ]
  }
  
  score <- function(loc, a) {
    sum(arc_cost(loc)[cbind(seq_len(nrow(c)), a)])
  }
  
  best <- NULL
  for (run in seq_len(p$starts)) {
    # Independent random starts from observed customers.
    ix <- sample(seq_len(nrow(c)), p$k, replace = p$k > nrow(c))
    loc <- c[, c("Longitude", "Latitude")][ix, , drop = FALSE]
    rownames(loc) <- NULL
    if (nrow(w))
      loc[seq_len(nrow(w)), ] <- w[, c("Longitude", "Latitude")]
    
    previous <- Inf
    a <- assign_customers(loc)
    
    for (iteration in seq_len(p$iterations)) {
      for (j in seq_len(p$k)) {
        if (j <= nrow(w) && w$Move_range_km[j] == 0) next
        if (!any(a == j)) next
        
        objective <- function(x) {
          candidate <- loc
          candidate[j, ] <- enforce_radius(x, j)
          score(candidate, a)
        }
        
        bounds_low <- lower
        bounds_high <- upper
        if (j <= nrow(w)) {
          bounds_low <- pmin(lower, c(w$Longitude[j], w$Latitude[j]))
          bounds_high <- pmax(upper, c(w$Longitude[j], w$Latitude[j]))
        }
        
        fit <- optim(
          pmin(bounds_high, pmax(bounds_low, as.numeric(loc[j, ]))),
          objective, method = "L-BFGS-B",
          lower = bounds_low, upper = bounds_high,
          control = list(maxit = 100)
        )
        proposed <- loc
        proposed[j, ] <- enforce_radius(fit$par, j)
        if (score(proposed, a) <= score(loc, a))
          loc <- proposed
      }
      
      a <- assign_customers(loc)
      value <- score(loc, a)
      if (abs(previous - value) <= p$tolerance * max(1, abs(value))) break
      previous <- value
    }
    
    if (is.null(best) || value < best$objective)
      best <- list(loc = loc, assignment = a, objective = value,
                   iteration = iteration)
  }
  
  loc <- best$loc
  a <- best$assignment
  d <- distance_matrix(c, loc, p$circuity)[cbind(seq_len(nrow(c)), a)]
  throughput <- vapply(seq_len(p$k), function(j)
    sum(c$Demand[a == j]), numeric(1))
  
  out_cost <- sum(c$Demand * customer_rate * d)
  in_cost <- if (p$suppliers) sum(sweep(
    distance_matrix(s, loc, p$circuity), 1,
    s$Volume * supplier_rate, "*"
  ) * rep(throughput / sum(c$Demand), each = nrow(s))) else 0
  
  wh <- data.frame(
    Warehouse_ID = all_ids,
    Longitude = loc$Longitude,
    Latitude = loc$Latitude,
    Throughput = throughput,
    Capacity = cap
  )
  
  assignments <- c
  assignments$Warehouse_ID <- all_ids[a]
  assignments$Distance_km <- d
  assignments$Within_lead_time <- d <= p$lead_km
  assignments$Outbound_cost <- c$Demand * customer_rate * d
  
  list(warehouses = wh, customers = assignments, suppliers = s,
       p = p, objective = best$objective,
       outbound = out_cost, inbound = in_cost,
       iteration = best$iteration)
}

# ---------------------------------------------------------------------
# User interface
# ---------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Warehouse Centre of Gravity — Prototype by Steven Van Aken"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      numericInput("k", "Total number of warehouses", 2, min = 1),
      actionButton("run", "Calculate", class = "btn-primary"),
      verbatimTextOutput("calculation_status"),
      hr(),
      tabsetPanel(
        tabPanel("Parameters",
                 selectInput("units", "Display distance", "Kilometers"),
                 numericInput("circuity", "Distance circuity factor", 1.25, min = 1),
                 numericInput("lead", "Lead-time distance (display units)", 200, min = 0),
                 checkboxInput("improve", "Improve service level", FALSE),
                 numericInput("penalty", "Penalty / qty / km above limit", 0.1, min = 0),
                 checkboxInput("supply", "Include supplier costs: option 3 approximation",
                               FALSE),
                 numericInput("random", "Randomize demand +/- (%)", 0, 0, 99),
                 selectInput("distribution", "Demand distribution",
                             c("Uniform", "Truncated normal")),
                 numericInput("seed", "Random seed", 123, min = 1),
                 numericInput("starts", "Independent starts", 3, 1, 30),
                 numericInput("iterations", "Maximum alternating iterations", 20, 1),
                 numericInput("tolerance", "Relative convergence tolerance", 1e-5,
                              min = 1e-9)
        ),
        tabPanel("Rates",
                 numericInput("ftl", "Full truck load quantity", 24, min = 0.01),
                 numericInput("rate", "FTL cost / km", 1.5, min = 0),
                 checkboxInput("auto", "Override individual customer shipment sizes",
                               FALSE),
                 numericInput("csize", "Default customer shipment size", 10, min = .01),
                 numericInput("cratio", "Customer rate multiplier", 1, min = 0),
                 numericInput("sratio", "Supplier rate multiplier", 1, min = 0),
                 helpText("Bezier endpoints: (0,0) and (100%,100%).",
                          "Prototype defaults, not verified Stelling defaults."),
                 sliderInput("x1", "Start control X (%)", 0, 100, 10),
                 sliderInput("y1", "Start control Y (%)", 0, 100, 40),
                 sliderInput("x2", "End control X (%)", 0, 100, 60),
                 sliderInput("y2", "End control Y (%)", 0, 100, 90)
        ),
        tabPanel("Constraints",
                 checkboxInput("predefined", "Use predefined warehouses", TRUE),
                 checkboxInput("assignments", "Fixed / subset customer assignments", FALSE),
                 checkboxInput("capacity", "Hard capacity limits", FALSE),
                 checkboxInput("groups", "Single-source customer groups", FALSE),
                 helpText("Capacity assignment is binary, not splittable LP.",
                          "Supplier volume is not a product-level balance model.")
        ),
        tabPanel("Map",
                 checkboxInput("flows", "Display customer flows", TRUE),
                 checkboxInput("sflows", "Display supplier flows", TRUE),
                 checkboxInput("ids", "Display customer IDs", FALSE),
                 checkboxInput("black", "Black customers beyond lead-time limit", TRUE),
                 checkboxInput("rings", "Display lead-time distance circles", FALSE),
                 checkboxInput("throughput", "Display warehouse throughput", TRUE),
                 checkboxInput("fixedsize", "Fixed-size warehouse circles", FALSE),
                 sliderInput("csize_map", "Customer size factor", 1, 12, 4),
                 sliderInput("wsize_map", "CoG size factor", 2, 20, 10),
                 sliderInput("opacity", "Customer opacity", 0, 1, .75),
                 sliderInput("wopacity", "CoG opacity", 0, 1, .9),
                 sliderInput("flow_size", "Flow size factor", .2, 5, 1),
                 numericInput("map_lat", "Map centre latitude", 50.85),
                 numericInput("map_lon", "Map centre longitude", 4.35),
                 numericInput("map_zoom", "Map zoom", 6, 1, 18),
                 actionButton("recenter", "Apply map centre")
        ),
        tabPanel("Future release",
                 helpText("Visible backlog only. These controls are disabled."),
                 tags$fieldset(
                   disabled = "disabled",
                   selectInput("supply_future", "Additional supply options",
                               c("1: closest warehouse", "2: LP auto-sized suppliers",
                                 "4: via closest warehouse")),
                   checkboxInput("shortage_future", "Forbid supply shortage", FALSE),
                   checkboxInput("border_future", "Border crossing points", FALSE),
                   checkboxInput("satellite_future", "Main / satellite warehouses", FALSE),
                   checkboxInput("city_future", "Display closest city", FALSE),
                   numericInput("population_future", "Minimum city population", 15000),
                   textInput("country_future", "Closest-city country filter", ""),
                   sliderInput("sensitivity_future", "Cost-increase area (%)", 1, 10, 5),
                   checkboxInput("enumerate_future", "Enumerate open/close networks", FALSE),
                   checkboxInput("fleet_future", "Own-transport estimator", FALSE)
                 )
        )
      )
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Map",
                 fluidRow(
                   column(4, selectInput("country", "Customer country filter", "All")),
                   column(4, selectInput("warehouse_filter", "Warehouse filter", "All")),
                   column(4, selectInput("distance_filter", "Distance filter",
                                         c("All", "Within limit", "Beyond limit")))
                 ),
                 leafletOutput("map", height = 650),
                 verbatimTextOutput("status"),
                 helpText("Filters affect display and service metrics only;",
                          "they never silently remove customers from optimisation.")
        ),
        tabPanel("Customers",
                 fileInput("customers_csv", "CSV upload: header required"),
                 helpText("Exact column names must match the table. UTF-8, comma-separated.",
                          "Warehouse_IDs uses | between allowed warehouse IDs."),
                 DTOutput("customers_table")
        ),
        tabPanel("Warehouses",
                 fileInput("warehouses_csv", "CSV upload: header required"),
                 helpText("Use Inf for unlimited capacity or movement.",
                          "Move_range_km is always measured in kilometres."),
                 DTOutput("warehouses_table")
        ),
        tabPanel("Suppliers",
                 fileInput("suppliers_csv", "CSV upload: header required"),
                 DTOutput("suppliers_table")
        ),
        tabPanel("Results - Disabled",
                 h4("Network totals"),
                 DTOutput("totals"),
                 h4("Warehouse report"),
                 DTOutput("warehouse_report"),
                 h4("Filtered service metrics"),
                 DTOutput("metrics"),
                 h4("Customer assignments"),
                 DTOutput("assignment_report")
        ),
        tabPanel("Transport curve", plotOutput("curve"))
      )
    )
  )
)

# ---------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------

server <- function(input, output, session) {
  calculation_status <- reactiveVal(
    "Ready. Click Calculate to start."
  )
  
  output$calculation_status <- renderText({
    calculation_status()
  })
  
  tables <- reactiveValues(
    customers = demo_customers,
    warehouses = demo_warehouses,
    suppliers = demo_suppliers
  )
  result <- reactiveVal(NULL)
  
  for (kind in c("customers", "warehouses", "suppliers")) local({
    key <- kind
    schema <- names(isolate(tables[[key]]))
    
    output[[paste0(key, "_table")]] <- renderDT({
      datatable(tables[[key]], editable = TRUE, rownames = FALSE,
                options = list(scrollX = TRUE, pageLength = 10))
    })
    
    observeEvent(input[[paste0(key, "_table_cell_edit")]], {
      e <- input[[paste0(key, "_table_cell_edit")]]
      d <- tables[[key]]
      # DT uses zero-based column indices when rownames = FALSE.
      j <- e$col + 1L
      d[e$row, j] <- DT::coerceValue(e$value, d[e$row, j])
      tables[[key]] <- d
      result(NULL)
    })
    
    observeEvent(input[[paste0(key, "_csv")]], {
      tryCatch({
        file <- input[[paste0(key, "_csv")]]
        d <- read.csv(file$datapath, stringsAsFactors = FALSE,
                      check.names = FALSE, na.strings = "NA")
        assert(all(schema %in% names(d)),
               paste("Required columns:", paste(schema, collapse = ", ")))
        tables[[key]] <- d[, schema, drop = FALSE]
        result(NULL)
      }, error = function(e) showNotification(conditionMessage(e),
                                              type = "error", duration = NULL))
    })
  })
  
  params <- reactive({
    list(
      k = as.integer(input$k),
      circuity = input$circuity,
      lead_km = input$lead * if (input$units == "Miles") 1.609344 else 1,
      improve_service = input$improve,
      service_penalty = input$penalty,
      suppliers = input$supply,
      random_pct = input$random,
      distribution = input$distribution,
      seed = as.integer(input$seed),
      starts = as.integer(input$starts),
      iterations = as.integer(input$iterations),
      tolerance = input$tolerance,
      ftl = input$ftl,
      ftl_rate = input$rate,
      auto_ship = input$auto,
      customer_size = input$csize,
      customer_ratio = input$cratio,
      supplier_ratio = input$sratio,
      x1 = input$x1/100, y1 = input$y1/100,
      x2 = input$x2/100, y2 = input$y2/100,
      assignments = input$assignments,
      capacities = input$capacity,
      groups = input$groups,
      predefined = input$predefined
    )
  })
  
  # Prevent showing stale results after a model parameter changes.
  observeEvent(params(), result(NULL), ignoreInit = TRUE)
  
  observeEvent(input$run, {
    result(NULL)
    tryCatch({
      p <- params()
      w <- tables$warehouses
      if (!p$predefined) w <- w[FALSE, , drop = FALSE]
      d <- prepare_data(tables$customers, w, tables$suppliers, p)
      withProgress(message = "Optimising network", value = 0, {
        r <- solve_network(d$customers, d$warehouses, d$suppliers, p)
        result(r)
      })
    }, error = function(e) showNotification(
      conditionMessage(e), type = "error", duration = NULL
    ))
  })
  
  observeEvent(result(), {
    r <- result()
    if (is.null(r)) return()
    updateSelectInput(session, "country",
                      choices = c("All", sort(unique(r$customers$Country))),
                      selected = "All")
    updateSelectInput(session, "warehouse_filter",
                      choices = c("All", r$warehouses$Warehouse_ID),
                      selected = "All")
  })
  
  filtered <- reactive({
    r <- req(result())
    d <- r$customers
    if (input$country != "All")
      d <- d[d$Country == input$country, , drop = FALSE]
    if (input$warehouse_filter != "All")
      d <- d[d$Warehouse_ID == input$warehouse_filter, , drop = FALSE]
    if (input$distance_filter == "Within limit")
      d <- d[d$Within_lead_time, , drop = FALSE]
    if (input$distance_filter == "Beyond limit")
      d <- d[!d$Within_lead_time, , drop = FALSE]
    d
  })
  
  output$map <- renderLeaflet({
    leaflet() |>
      addTiles(
        urlTemplate = "https://tile.openstreetmap.org/{z}/{x}/{y}.png",
        attribution = '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
      ) |>
      setView(4.35, 50.85, 6)
  })
  
  observeEvent(input$recenter, {observeEvent(input$run, {
    result(NULL)
    
    stage <- "Reading parameters"
    calculation_status("Calculate received. Preparing the model...")
    message("[CoG] Calculate received")
    
    tryCatch(
      {
        p <- params()
        
        stage <- "Preparing input tables"
        message("[CoG] ", stage)
        
        w <- tables$warehouses
        
        if (!p$predefined) {
          w <- w[FALSE, , drop = FALSE]
        }
        
        stage <- "Validating input data"
        message("[CoG] ", stage)
        
        d <- prepare_data(
          tables$customers,
          w,
          tables$suppliers,
          p
        )
        
        stage <- "Running the optimisation solver"
        message("[CoG] ", stage)
        
        calculation_status(
          paste(
            "Running optimisation:",
            nrow(d$customers), "customers;",
            p$k, "warehouses;",
            p$starts, "start(s)."
          )
        )
        
        r <- withProgress(
          message = "Optimising warehouse network",
          value = 0,
          {
            incProgress(0.1, detail = "Starting solver")
            
            solution <- solve_network(
              d$customers,
              d$warehouses,
              d$suppliers,
              p
            )
            
            incProgress(0.9, detail = "Calculation finished")
            solution
          }
        )
        
        stage <- "Publishing the results"
        result(r)
        
        calculation_status(
          sprintf(
            paste0(
              "Calculation completed.\n",
              "Objective: %.2f\n",
              "Warehouses: %d\n",
              "Check the Map and Results tabs."
            ),
            r$objective,
            nrow(r$warehouses)
          )
        )
        
        message("[CoG] Calculation completed")
      },
      
      error = function(e) {
        error_text <- paste0(
          "Calculation failed.\n",
          "Stage: ", stage, "\n",
          "Error: ", conditionMessage(e)
        )
        
        calculation_status(error_text)
        message("[CoG] ", error_text)
        
        showNotification(
          error_text,
          type = "error",
          duration = NULL
        )
      }
    )
  }, ignoreInit = TRUE)})
  
  observe({
    r <- result()
    proxy <- leafletProxy("map") |>
      clearMarkers() |> clearShapes() |> clearControls()
    if (is.null(r)) return()
    
    d <- filtered()
    w <- r$warehouses
    pal <- colorFactor("Dark2", domain = w$Warehouse_ID)
    colours <- pal(d$Warehouse_ID)
    if (input$black)
      colours[!d$Within_lead_time] <- "#000000"
    
    if (input$flows && nrow(d)) {
      for (i in seq_len(nrow(d))) {
        j <- match(d$Warehouse_ID[i], w$Warehouse_ID)
        proxy <- proxy |> addPolylines(
          lng = c(w$Longitude[j], d$Longitude[i]),
          lat = c(w$Latitude[j], d$Latitude[i]),
          color = pal(w$Warehouse_ID[j]),
          weight = input$flow_size, opacity = .35
        )
      }
    }
    
    if (r$p$suppliers && input$sflows) {
      s <- r$suppliers
      for (i in seq_len(nrow(s))) for (j in seq_len(nrow(w))) {
        if (w$Throughput[j] <= 0) next
        volume <- s$Volume[i] * w$Throughput[j] / sum(w$Throughput)
        proxy <- proxy |> addPolylines(
          lng = c(s$Longitude[i], w$Longitude[j]),
          lat = c(s$Latitude[i], w$Latitude[j]),
          color = "#7b3294", dashArray = "5,5",
          weight = input$flow_size,
          label = paste("Inbound quantity:", round(volume, 1))
        )
      }
      proxy <- proxy |> addCircleMarkers(
        data = s, lng = ~Longitude, lat = ~Latitude,
        radius = 7, color = "#7b3294", fillOpacity = .9,
        label = ~as.character(Supplier_ID)
      )
    }
    
    if (nrow(d)) {
      labels <- if (input$ids) as.character(d$Customer_ID) else NULL
      proxy <- proxy |> addCircleMarkers(
        data = d, lng = ~Longitude, lat = ~Latitude,
        radius = input$csize_map * sqrt(d$Demand / mean(d$Demand)),
        color = colours, fillOpacity = input$opacity,
        weight = 1, label = labels
      )
    }
    
    labels <- if (input$throughput)
      paste(w$Warehouse_ID, "| throughput:", round(w$Throughput, 1))
    else w$Warehouse_ID
    
    radii <- if (input$fixedsize) rep(input$wsize_map, nrow(w)) else
      input$wsize_map * sqrt(pmax(w$Throughput, 1) /
                               max(1, mean(w$Throughput)))
    
    proxy <- proxy |> addCircleMarkers(
      data = w, lng = ~Longitude, lat = ~Latitude,
      radius = radii, color = pal(w$Warehouse_ID),
      fillColor = "white", weight = 3,
      fillOpacity = input$wopacity, label = labels
    )
    
    if (input$rings)
      proxy <- proxy |> addCircles(
        data = w, lng = ~Longitude, lat = ~Latitude,
        # Lead time is circuity-adjusted: convert to geographic radius.
        radius = r$p$lead_km / r$p$circuity * 1000,
        color = pal(w$Warehouse_ID), fill = FALSE, weight = 1
      )
    
    proxy |> addLegend(position = "bottomright", pal = pal,
                       values = w$Warehouse_ID, title = "Warehouse")
  })
  
  observeEvent(result(), {
    r <- result()
    if (is.null(r)) return()
    coords <- rbind(
      r$customers[, c("Longitude", "Latitude")],
      r$warehouses[, c("Longitude", "Latitude")]
    )
    leafletProxy("map") |> fitBounds(
      min(coords$Longitude), min(coords$Latitude),
      max(coords$Longitude), max(coords$Latitude)
    )
  })
  
  output$status <- renderText({
    r <- result()
    if (is.null(r)) return("No current solution. Click Calculate.")
    paste(
      "Feasible heuristic solution; no global optimality guarantee.",
      "\nObjective includes transport costs and any service penalty.",
      "\nAlternating iterations in best run:", r$iteration
    )
  })
  
  output$totals <- renderDT({
    r <- req(result())
    datatable(data.frame(
      Metric = c("Outbound transportation cost",
                 "Inbound transportation cost",
                 "Service penalty",
                 "Total objective"),
      Value = c(r$outbound, r$inbound,
                r$objective-r$outbound-r$inbound, r$objective)
    ), rownames = FALSE, options = list(dom = "t"))
  })
  
  output$warehouse_report <- renderDT({
    datatable(req(result())$warehouses, rownames = FALSE,
              options = list(scrollX = TRUE))
  })
  
  output$assignment_report <- renderDT({
    d <- filtered()
    if (input$units == "Miles") {
      d$Distance_miles <- d$Distance_km / 1.609344
      d$Distance_km <- NULL
    }
    datatable(d, rownames = FALSE, options = list(scrollX = TRUE))
  })
  
  output$metrics <- renderDT({
    d <- filtered()
    if (!nrow(d)) return(datatable(data.frame(Message = "No matching customers")))
    divisor <- if (input$units == "Miles") 1.609344 else 1
    x <- d$Distance_km / divisor
    datatable(data.frame(
      Metric = c("Customers included", "Demand included",
                 "Demand-weighted average distance", "Minimum distance",
                 "Average distance", "Maximum distance",
                 "Customers within lead-time limit (%)",
                 "Demand within lead-time limit (%)"),
      Value = c(nrow(d), sum(d$Demand),
                weighted.mean(x, d$Demand), min(x), mean(x), max(x),
                100 * mean(d$Within_lead_time),
                100 * sum(d$Demand[d$Within_lead_time])/sum(d$Demand))
    ), rownames = FALSE, options = list(dom = "t"))
  })
  
  output$curve <- renderPlot({
    p <- params()
    req(p$x1 <= p$x2, p$y1 <= p$y2)
    x <- seq(0, 1, length.out = 200)
    plot(x * 100, curve_fraction(x, p) * 100,
         type = "l", lwd = 2, col = "#2166ac",
         xlab = "Shipment size (% FTL)",
         ylab = "Transport cost (% FTL cost)",
         main = "Cubic Bezier transport-rate curve (how shipment cost changes with shipment size)")
    abline(0, 1, lty = 2, col = "grey60")
    grid()
  })
}

shinyApp(ui, server)