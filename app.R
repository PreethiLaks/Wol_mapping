library(shiny)
library(terra)
library(sf)

options(terra.plot.legend = TRUE)
terraOptions(memfrac = 0.8, progress = 0)

# ----------------------------
# File paths
# ----------------------------
paths <- list(
  ara  = file.path("data", "An_arabiensis.tif"),
  col  = file.path("data", "An_coluzzii.geotiff"),
  mou  = file.path("data", "An_moucheti.tif"),
  ste  = file.path("data", "An_stephensi.tif"),
  pfpr = file.path("data", "PfPR_mean.tif"),
  itn  = file.path("data", "ITN_use_rate.tif"),
  inc  = file.path("data", "Pf_Incidence_mean_2000.tif"),
  pop  = file.path("data", "POP_MEAN_2015_2025.tif"),
  tt   = file.path("data", "Global_travel_time.tif"),
  gadm = file.path("data", "gadm_410.gpkg")
)

for (nm in names(paths)) {
  if (!file.exists(paths[[nm]])) stop("File not found: ", paths[[nm]])
}

# ----------------------------
# Helper functions
# ----------------------------
to01 <- function(r) {
  mx <- suppressWarnings(as.numeric(global(r, "max", na.rm = TRUE)))
  r1 <- if (is.finite(mx) && mx > 1) r / 100 else r
  clamp(r1, 0, 1, values = TRUE)
}

safe_div <- function(num, den, eps = 1e-9) num / (den + eps)

fast_align <- function(r, template, method = "bilinear") {
  if (!identical(crs(r), crs(template))) {
    project(r, template, method = method)
  } else {
    resample(r, template, method = method)
  }
}

tier_mask_value <- function(tier_r, k) {
  m <- tier_r
  m[m != k] <- NA
  m[m == k] <- 1
  m
}

tier_any_mask <- function(tier_r) {
  m <- tier_r
  m[!(m %in% 1:3)] <- NA
  m[m %in% 1:3] <- 1
  m
}

make_tier <- function(x, valid_mask, breaks, labels) {
  out <- rast(x)
  values(out) <- NA
  
  out[!is.na(valid_mask) & x >= breaks[1] & x < breaks[2]] <- labels[1]
  out[!is.na(valid_mask) & x >= breaks[2] & x < breaks[3]] <- labels[2]
  out[!is.na(valid_mask) & x >= breaks[3]] <- labels[3]
  
  out
}

coverage_table_countries <- function(tier_r, v0) {
  count_countries <- function(mask_r, v0) {
    ex <- terra::extract(mask_r, v0, fun = sum, na.rm = TRUE, weights = TRUE)
    if (is.null(ex) || nrow(ex) == 0) return(0)
    vals <- ex[, ncol(ex)]
    vals[is.na(vals)] <- 0
    sum(vals > 0)
  }
  
  m1   <- mask(tier_mask_value(tier_r, 1), v0)
  m2   <- mask(tier_mask_value(tier_r, 2), v0)
  m3   <- mask(tier_mask_value(tier_r, 3), v0)
  m123 <- mask(tier_any_mask(tier_r), v0)
  
  data.frame(
    Tier = c("Tier 1", "Tier 2", "Tier 3", "Tier 1-3 (total)"),
    Countries_covered = c(
      count_countries(m1, v0),
      count_countries(m2, v0),
      count_countries(m3, v0),
      count_countries(m123, v0)
    ),
    stringsAsFactors = FALSE
  )
}

impact_table <- function(tier_r, layers, burden_r, pop_r) {
  bur_valid <- mask(burden_r, layers$burden_mask)
  pop_valid <- mask(pop_r, layers$burden_mask)
  
  m1   <- tier_mask_value(tier_r, 1)
  m2   <- tier_mask_value(tier_r, 2)
  m3   <- tier_mask_value(tier_r, 3)
  m123 <- tier_any_mask(tier_r)
  
  # Population by tier
  p1 <- as.numeric(global(mask(pop_valid, m1), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(p1)) p1 <- 0
  p2 <- as.numeric(global(mask(pop_valid, m2), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(p2)) p2 <- 0
  p3 <- as.numeric(global(mask(pop_valid, m3), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(p3)) p3 <- 0
  p123 <- as.numeric(global(mask(pop_valid, m123), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(p123)) p123 <- 0
  
  # Mean annual cases by tier
  c1 <- as.numeric(global(mask(bur_valid, m1), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(c1)) c1 <- 0
  c2 <- as.numeric(global(mask(bur_valid, m2), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(c2)) c2 <- 0
  c3 <- as.numeric(global(mask(bur_valid, m3), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(c3)) c3 <- 0
  c123 <- as.numeric(global(mask(bur_valid, m123), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(c123)) c123 <- 0
  
  # Total WHO AFRO mean annual cases
  afro_cases_all <- as.numeric(global(mask(bur_valid, layers$v0), "sum", na.rm = TRUE)[1, 1])
  if (!is.finite(afro_cases_all) || afro_cases_all <= 0) afro_cases_all <- NA_real_
  
  fmt_cases_pct <- function(x, total) {
    if (is.na(total)) return(paste0(format(round(x), big.mark = ","), " (NA)"))
    pct <- 100 * x / total
    paste0(
      format(round(x), big.mark = ","),
      " (",
      formatC(pct, format = "f", digits = 2),
      "%)"
    )
  }
  
  data.frame(
    Metric = c(
      "Population benefiting - Tier 1",
      "Population benefiting - Tier 2",
      "Population benefiting - Tier 3",
      "Population benefiting - Tier 1-3 combined",
      "Estimated mean annual cases - Tier 1 (% WHO AFRO)",
      "Estimated mean annual cases - Tier 2 (% WHO AFRO)",
      "Estimated mean annual cases - Tier 3 (% WHO AFRO)",
      "Estimated mean annual cases - Tier 1-3 combined (% WHO AFRO)"
    ),
    Value = c(
      format(round(p1), big.mark = ","),
      format(round(p2), big.mark = ","),
      format(round(p3), big.mark = ","),
      format(round(p123), big.mark = ","),
      fmt_cases_pct(c1, afro_cases_all),
      fmt_cases_pct(c2, afro_cases_all),
      fmt_cases_pct(c3, afro_cases_all),
      fmt_cases_pct(c123, afro_cases_all)
    ),
    stringsAsFactors = FALSE
  )
}

compute_zoom <- function(r, v0b, auto_zoom, zoom_pad) {
  xlim <- ylim <- NULL
  
  if (isTRUE(auto_zoom)) {
    idx <- which(!is.na(values(r)))
    if (length(idx) >= 1) {
      xy <- xyFromCell(r, idx)
      xr <- range(xy[, 1], na.rm = TRUE)
      yr <- range(xy[, 2], na.rm = TRUE)
      
      pad <- zoom_pad / 100
      dx <- diff(xr)
      dy <- diff(yr)
      rx <- res(r)[1]
      ry <- res(r)[2]
      
      if (!is.finite(dx) || dx == 0) dx <- 2 * rx
      if (!is.finite(dy) || dy == 0) dy <- 2 * ry
      
      xlim <- c(xr[1] - dx * pad, xr[2] + dx * pad)
      ylim <- c(yr[1] - dy * pad, yr[2] + dy * pad)
      
      xlim <- c(max(xlim[1], xmin(v0b)), min(xlim[2], xmax(v0b)))
      ylim <- c(max(ylim[1], ymin(v0b)), min(ylim[2], ymax(v0b)))
    }
  }
  
  list(xlim = xlim, ylim = ylim)
}

plot_tier_map <- function(tier_r, layers, title_text, auto_zoom, zoom_pad) {
  tier_cols <- c("1" = "#1a9850", "2" = "#f16913", "3" = "#d73027")
  tier_labels <- c("1" = "Overall Tier 1", "2" = "Overall Tier 2", "3" = "Overall Tier 3")
  
  vals <- sort(unique(na.omit(as.vector(values(tier_r)))))
  vals <- vals[vals %in% 1:3]
  
  if (length(vals) == 0) {
    plot.new()
    title(paste0(title_text, "\n(No eligible cells)"))
    return(invisible(NULL))
  }
  
  z <- compute_zoom(tier_r, layers$v0, auto_zoom, zoom_pad)
  af_ext <- ext(layers$v0)
  
  xlim <- if (is.null(z$xlim)) c(xmin(af_ext), xmax(af_ext)) else z$xlim
  ylim <- if (is.null(z$ylim)) c(ymin(af_ext), ymax(af_ext)) else z$ylim
  
  e_view <- ext(xlim[1], xlim[2], ylim[1], ylim[2])
  v0_view <- crop(layers$v0, e_view)
  v1_view <- crop(layers$v1, e_view)
  
  plot(v0_view, xlim = xlim, ylim = ylim, col = "white", border = NA, axes = TRUE, main = title_text)
  plot(tier_r, add = TRUE, col = tier_cols, breaks = c(0.5, 1.5, 2.5, 3.5), legend = FALSE, colNA = NA)
  plot(v1_view, add = TRUE, lwd = 0.6, col = NA, border = "#BFD7FF")
  plot(v0_view, add = TRUE, lwd = 1.2, col = NA, border = "#6996E8")
  
  legend(
    "bottomleft",
    legend = tier_labels[as.character(vals)],
    fill   = tier_cols[as.character(vals)],
    border = NA,
    bty    = "n",
    cex    = 0.95,
    title  = "Overall Tier"
  )
}

# ----------------------------
# Data loading
# ----------------------------
PFPR_MEAN <- rast(paths$pfpr)
ITN_MEAN  <- rast(paths$itn)
INC_MEAN  <- rast(paths$inc)
TT_MEAN <- rast(paths$tt)
POP_MEAN  <- rast(paths$pop)

species_rasters <- list(
  arabiensis = rast(paths$ara),
  coluzzii   = rast(paths$col),
  moucheti   = rast(paths$mou),
  stephensi  = rast(paths$ste)
)

afro_iso <- c(
  "DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM","COG","CIV","COD","GNQ","ERI",
  "SWZ","ETH","GAB","GMB","GHA","GIN","GNB","KEN","LSO","LBR","MDG","MWI","MLI","MRT","MUS","MOZ",
  "NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","ZAF","SSD","TGO","UGA","TZA","ZMB","ZWE"
)

a0 <- st_read(paths$gadm, layer = "ADM_0", quiet = TRUE)
a1 <- st_read(paths$gadm, layer = "ADM_1", quiet = TRUE)

ADMIN0 <- vect(a0[a0$GID_0 %in% afro_iso, ])
ADMIN1 <- vect(a1[a1$GID_0 %in% afro_iso, ])

# ----------------------------
# Precompute common template and aligned layers
# ----------------------------
ADMIN0_PFPR <- if (!identical(crs(ADMIN0), crs(PFPR_MEAN))) project(ADMIN0, crs(PFPR_MEAN)) else ADMIN0
TEMPLATE <- crop(PFPR_MEAN, ADMIN0_PFPR, snap = "out")

ADMIN0_ALIGNED <- if (!identical(crs(ADMIN0), crs(TEMPLATE))) project(ADMIN0, crs(TEMPLATE)) else ADMIN0
ADMIN1_ALIGNED <- if (!identical(crs(ADMIN1), crs(TEMPLATE))) project(ADMIN1, crs(TEMPLATE)) else ADMIN1

PFPR_ALIGNED <- to01(fast_align(PFPR_MEAN, TEMPLATE))
ITN_ALIGNED  <- to01(fast_align(ITN_MEAN,  TEMPLATE))
INC_ALIGNED  <- fast_align(INC_MEAN, TEMPLATE, method = "bilinear")
TT_ALIGNED <- fast_align(TT_MEAN, TEMPLATE, method = "bilinear")
POP_ALIGNED  <- fast_align(POP_MEAN, TEMPLATE, method = "near")

SPECIES_ALIGNED <- lapply(species_rasters, function(r) fast_align(to01(r), TEMPLATE))

CELL_AREA_KM2 <- cellSize(TEMPLATE, unit = "km")
names(CELL_AREA_KM2) <- "cell_km2"
# ----------------------------
# Use-case rules
# ----------------------------
use_cases <- list(
  "1a" = list(
    itn_breaks = c(-Inf, 0.60, 0.80, Inf),
    itn_labels = c(1, 2, 3),
    pfpr_breaks = c(-Inf, 0.15, 0.40, Inf),
    pfpr_labels = c(1, 2, 3)
  ),
  "1b" = list(
    itn_breaks = c(-Inf, 0.70, 0.80, Inf),
    itn_labels = c(1, 2, 3),
    pfpr_breaks = c(-Inf, 0.15, 0.40, Inf),
    pfpr_labels = c(3, 2, 1)
  ),
  "2" = list(
    itn_breaks = c(-Inf, 0.60, 0.80, Inf),
    itn_labels = c(3, 2, 1),
    pfpr_breaks = c(-Inf, 0.15, 0.40, Inf),
    pfpr_labels = c(3, 2, 1)
  ),
  "3" = list(
    itn_breaks = c(-Inf, 0.60, 0.80, Inf),
    itn_labels = c(3, 2, 1),
    pfpr_breaks = c(-Inf, 0.05, 0.15, 0.40),
    pfpr_labels = c(1, 2, 3)
  )
)

# ----------------------------
# UI
# ----------------------------
ui <- fluidPage(
  titlePanel("Wolbachia mapping"),
  sidebarLayout(
    sidebarPanel(
      selectInput("species", "Species mask",
                  c("arabiensis", "coluzzii", "moucheti", "stephensi"),
                  selected = "moucheti"
      ),
      hr(),
      checkboxInput("auto_zoom", "Auto-zoom to eligible area", TRUE),
      sliderInput("zoom_pad", "Zoom padding (%)", min = 0, max = 25, value = 10, step = 1)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Use case 1a", plotOutput("map_1a", height = 650), h4("Impact Summary"), tableOutput("impact_1a"), h4("Coverage"), tableOutput("cov_1a")),
        tabPanel("Use case 1b", plotOutput("map_1b", height = 650), h4("Impact Summary"), tableOutput("impact_1b"), h4("Coverage"), tableOutput("cov_1b")),
        tabPanel("Use case 2",  plotOutput("map_2",  height = 650), h4("Impact Summary"), tableOutput("impact_2"),  h4("Coverage"), tableOutput("cov_2")),
        tabPanel("Use case 3",  plotOutput("map_3",  height = 650), h4("Impact Summary"), tableOutput("impact_3"),  h4("Coverage"), tableOutput("cov_3")),
        tabPanel("Use case 4", plotOutput("map_4", height = 650), h4("Impact Summary"), tableOutput("impact_4"), h4("Coverage"), tableOutput("cov_4"))
      )
    )
  )
)

# ----------------------------
# Server
# ----------------------------
server <- function(input, output, session) {
  
  afro_template <- reactive(TEMPLATE)
  
  base_layers <- reactive({
    tmpl <- TEMPLATE
    v0 <- ADMIN0_ALIGNED
    v1 <- ADMIN1_ALIGNED
    
    species_aligned <- SPECIES_ALIGNED
    base_p <- species_aligned[[input$species]]
    pres_mask <- ifel(!is.na(base_p) & base_p > 0, 1, NA)
    
    base_p0 <- base_p
    base_p0[is.na(base_p0)] <- 0
    
    other_names <- setdiff(names(species_aligned), input$species)
    others_stack <- rast(lapply(other_names, function(nm) {
      rr <- species_aligned[[nm]]
      rr[is.na(rr)] <- 0
      rr
    }))
    
    others_sum <- app(others_stack, sum)
    dom_share <- safe_div(base_p0, base_p0 + others_sum)
    
    pfpr01 <- PFPR_ALIGNED
    itn01  <- ITN_ALIGNED
    inc_al <- INC_ALIGNED
    pop_al <- POP_ALIGNED
    tt_al  <- TT_ALIGNED
    
    data_mask <- ifel(!is.na(pfpr01) & !is.na(itn01) & !is.na(base_p) & (base_p > 0), 1, NA)
    burden_mask <- ifel(!is.na(inc_al) & !is.na(pop_al) & (pop_al > 0) & (inc_al > 1e-12), 1, NA)
    
    list(
      tmpl = tmpl,
      v0 = v0,
      v1 = v1,
      pfpr01 = pfpr01,
      itn01 = itn01,
      pfpr_afro = mask(pfpr01, v0),
      pres_mask = pres_mask,
      dom_share = dom_share,
      data_mask = data_mask,
      burden_mask = burden_mask,
      inc_al = inc_al,
      pop_al = pop_al,
      tt_al = tt_al
    )
  })
  
  burden_aligned <- reactive({
    b <- base_layers()$pop_al * (base_layers()$inc_al / 1000)
    names(b) <- "cases_year"
    b
  })
  
  tier_vector <- function(layers) {
    ds <- layers$dom_share
    out <- rast(layers$tmpl)
    values(out) <- NA
    
    out[!is.na(layers$data_mask) & !is.na(layers$pres_mask) & ds >= 0.70] <- 1
    out[!is.na(layers$data_mask) & !is.na(layers$pres_mask) & ds >= 0.40 & ds < 0.70] <- 2
    out[!is.na(layers$data_mask) & !is.na(layers$pres_mask) & ds < 0.40] <- 3
    out
  }
  
  tier_itn <- function(layers, case_name) {
    cfg <- use_cases[[case_name]]
    make_tier(layers$itn01, layers$data_mask, cfg$itn_breaks, cfg$itn_labels)
  }
  
  tier_pfpr_overall <- function(layers, case_name) {
    cfg <- use_cases[[case_name]]
    make_tier(layers$pfpr01, layers$data_mask, cfg$pfpr_breaks, cfg$pfpr_labels)
  }
  
  tier_pfpr_denom <- function(layers, case_name) {
    cfg <- use_cases[[case_name]]
    out <- make_tier(layers$pfpr_afro, layers$pfpr_afro, cfg$pfpr_breaks, cfg$pfpr_labels)
    mask(out, layers$v0)
  }
  
  overall_tier <- function(layers, case_name) {
    tv <- tier_vector(layers)
    ti <- tier_itn(layers, case_name)
    tp <- tier_pfpr_overall(layers, case_name)
    
    out <- rast(layers$tmpl)
    values(out) <- NA
    valid <- !is.na(layers$data_mask) & !is.na(layers$pres_mask)
    
    out[(tv == 1) & (ti == 1) & (tp == 1) & valid] <- 1
    out[(tv == 2) & (ti == 2) & (tp == 2) & valid] <- 2
    out[(tv == 3) & (ti == 3) & (tp == 3) & valid] <- 3
    
    mask(out, layers$v0)
  }
  
  overall_tier_uc4 <- function(layers) {
    ste <- SPECIES_ALIGNED$stephensi
    ste <- mask(ste, layers$v0)
    tt  <- mask(layers$tt_al, layers$v0)
    pf  <- mask(layers$pfpr01, layers$v0)
    
    ste_tier <- ifel(ste >= 0.75, 1,
                     ifel(ste >= 0.15, 2,
                          ifel(ste >= 0.05, 3, NA)
                     )
    )
    
    tt_tier <- ifel(tt <= 60, 1,
                    ifel(tt <= 180, 2,
                         ifel(tt <= 360, 3, NA)
                    )
    )
    
    pfpr_tier <- ifel(pf >= 0.40, 1,
                      ifel(pf >= 0.15, 2,
                           ifel(pf >= 0, 3, NA)
                      )
    )
    
    out <- rast(layers$tmpl)
    values(out) <- NA
    
    out[ste_tier == 1 & tt_tier == 1 & pfpr_tier == 1] <- 1
    
    out[
      ste_tier %in% c(1, 2) &
        tt_tier %in% c(1, 2) &
        pfpr_tier %in% c(1, 2) &
        is.na(out)
    ] <- 2
    
    out[
      !is.na(ste_tier) &
        !is.na(tt_tier) &
        !is.na(pfpr_tier) &
        is.na(out)
    ] <- 3
    
    mask(out, layers$v0)
  }
  
  make_outputs <- function(case_name) {
    output[[paste0("map_", case_name)]] <- renderPlot({
      layers <- base_layers()
      plot_tier_map(
        overall_tier(layers, case_name),
        layers,
        paste0("Use case ", case_name, " - ", input$species),
        input$auto_zoom,
        input$zoom_pad
      )
    })
    
    output[[paste0("impact_", case_name)]] <- renderTable({
      layers <- base_layers()
      impact_table(overall_tier(layers, case_name), layers, burden_aligned(), layers$pop_al)
    })
    
    output[[paste0("cov_", case_name)]] <- renderTable({
      layers <- base_layers()
      coverage_table_countries(
        overall_tier(layers, case_name),
        layers$v0
      )
    })
  }
  
  for (nm in c("1a", "1b", "2", "3")) {
    local({
      case_name <- nm
      make_outputs(case_name)
    })
  }
  
  output$map_4 <- renderPlot({
    layers <- base_layers()
    plot_tier_map(
      overall_tier_uc4(layers),
      layers,
      "Use case 4 - An. stephensi urban/peri-urban vulnerability",
      input$auto_zoom,
      input$zoom_pad
    )
  })
  
  output$impact_4 <- renderTable({
    layers <- base_layers()
    impact_table(overall_tier_uc4(layers), layers, burden_aligned(), layers$pop_al)
  })
  
  output$cov_4 <- renderTable({
    layers <- base_layers()
    coverage_table_countries(
      overall_tier_uc4(layers),
      layers$v0
    )
  })
}

shinyApp(ui, server)