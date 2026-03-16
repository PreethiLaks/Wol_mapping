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

coverage_table_pfpr_denom_afro <- function(tier_r, v0, pfpr_band_afro_r) {
  A <- cellSize(tier_r, unit = "km")
  names(A) <- "cell_km2"
  
  count_countries <- function(mask_r, v0) {
    ex <- terra::extract(mask_r, v0, fun = sum, na.rm = TRUE, weights = TRUE)
    if (is.null(ex) || nrow(ex) == 0) return(0)
    vals <- ex[, ncol(ex)]
    vals[is.na(vals)] <- 0
    sum(vals > 0)
  }
  
  get_denom <- function(k) {
    mk <- tier_mask_value(pfpr_band_afro_r, k)
    mk <- mask(mk, v0)
    val <- as.numeric(global(mask(A, mk), "sum", na.rm = TRUE)[1, 1])
    if (!is.finite(val)) val <- 0
    val
  }
  
  get_stats <- function(k) {
    mk <- if (identical(k, "any")) tier_any_mask(tier_r) else tier_mask_value(tier_r, k)
    mk <- mask(mk, v0)
    
    area_km2 <- as.numeric(global(mask(A, mk), "sum", na.rm = TRUE)[1, 1])
    pixels   <- as.numeric(global(mk, "sum", na.rm = TRUE)[1, 1])
    
    if (!is.finite(area_km2)) area_km2 <- 0
    if (!is.finite(pixels)) pixels <- 0
    
    countries <- count_countries(mk, v0)
    
    c(area_km2 = area_km2, pixels = pixels, countries = countries)
  }
  
  s1 <- get_stats(1)
  s2 <- get_stats(2)
  s3 <- get_stats(3)
  sT <- get_stats("any")
  
  d1 <- get_denom(1)
  d2 <- get_denom(2)
  d3 <- get_denom(3)
  dT <- d1 + d2 + d3
  
  data.frame(
    Tier = c("Tier 1", "Tier 2", "Tier 3", "Tier 1-3 (total)"),
    Area_km2 = c(s1["area_km2"], s2["area_km2"], s3["area_km2"], sT["area_km2"]),
    Pct_of_AFRO_PfPR_band_km2 = c(
      ifelse(d1 > 0, 100 * s1["area_km2"] / d1, NA_real_),
      ifelse(d2 > 0, 100 * s2["area_km2"] / d2, NA_real_),
      ifelse(d3 > 0, 100 * s3["area_km2"] / d3, NA_real_),
      ifelse(dT > 0, 100 * sT["area_km2"] / dT, NA_real_)
    ),
    Pixels = c(s1["pixels"], s2["pixels"], s3["pixels"], sT["pixels"]),
    Countries_covered = c(s1["countries"], s2["countries"], s3["countries"], sT["countries"]),
    stringsAsFactors = FALSE
  )
}

impact_table <- function(tier_r, layers, burden_r, pop_r) {
  bur_valid <- mask(burden_r, layers$burden_mask)
  pop_valid <- mask(pop_r, layers$burden_mask)
  
  m1 <- tier_mask_value(tier_r, 1)
  m2 <- tier_mask_value(tier_r, 2)
  m3 <- tier_mask_value(tier_r, 3)
  many <- tier_any_mask(tier_r)
  
  t1_cases <- as.numeric(global(mask(bur_valid, m1), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(t1_cases)) t1_cases <- 0
  t2_cases <- as.numeric(global(mask(bur_valid, m2), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(t2_cases)) t2_cases <- 0
  t3_cases <- as.numeric(global(mask(bur_valid, m3), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(t3_cases)) t3_cases <- 0
  t123_cases <- as.numeric(global(mask(bur_valid, many), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(t123_cases)) t123_cases <- 0
  
  t1_pop <- as.numeric(global(mask(pop_valid, m1), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(t1_pop)) t1_pop <- 0
  t123_pop <- as.numeric(global(mask(pop_valid, many), "sum", na.rm = TRUE)[1, 1]); if (!is.finite(t123_pop)) t123_pop <- 0
  
  cases_per_1000_t1 <- ifelse(t1_pop > 0, 1000 * t1_cases / t1_pop, NA_real_)
  
  afro_cases_all <- as.numeric(global(mask(bur_valid, layers$v0), "sum", na.rm = TRUE)[1, 1])
  if (!is.finite(afro_cases_all) || afro_cases_all <= 0) afro_cases_all <- NA_real_
  
  pct_all_afro_t1   <- ifelse(is.na(afro_cases_all), NA_real_, 100 * t1_cases / afro_cases_all)
  pct_all_afro_t123 <- ifelse(is.na(afro_cases_all), NA_real_, 100 * t123_cases / afro_cases_all)
  
  data.frame(
    Metric = c(
      "Tier-1 population benefiting",
      "Tier-1 burden (cases/year)",
      "Tier-2 burden (cases/year)",
      "Tier-3 burden (cases/year)",
      "Burden (Tier 1-3 combined) (cases/year)",
      "Tier-1 burden per 1000 pop",
      "Population benefiting (Tier 1-3 combined)",
      "% of ALL-WHO-AFRO burden in Tier 1 (%)",
      "% of ALL-WHO-AFRO burden in Tier 1-3 combined (%)"
    ),
    Value = c(
      format(round(t1_pop), big.mark = ","),
      format(round(t1_cases), big.mark = ","),
      format(round(t2_cases), big.mark = ","),
      format(round(t3_cases), big.mark = ","),
      format(round(t123_cases), big.mark = ","),
      ifelse(is.na(cases_per_1000_t1), "NA", round(cases_per_1000_t1, 4)),
      format(round(t123_pop), big.mark = ","),
      ifelse(is.na(pct_all_afro_t1), "NA", formatC(pct_all_afro_t1, format = "f", digits = 6)),
      ifelse(is.na(pct_all_afro_t123), "NA", formatC(pct_all_afro_t123, format = "f", digits = 6))
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
        tabPanel("Use case 3",  plotOutput("map_3",  height = 650), h4("Impact Summary"), tableOutput("impact_3"),  h4("Coverage"), tableOutput("cov_3"))
      )
    )
  )
)

# ----------------------------
# Server
# ----------------------------
server <- function(input, output, session) {
  
  afro_template <- reactive({
    v0_pf <- if (!identical(crs(ADMIN0), crs(PFPR_MEAN))) project(ADMIN0, crs(PFPR_MEAN)) else ADMIN0
    crop(PFPR_MEAN, v0_pf, snap = "out")
  })
  
  base_layers <- reactive({
    tmpl <- afro_template()
    
    v0 <- if (!identical(crs(ADMIN0), crs(tmpl))) project(ADMIN0, crs(tmpl)) else ADMIN0
    v1 <- if (!identical(crs(ADMIN1), crs(tmpl))) project(ADMIN1, crs(tmpl)) else ADMIN1
    
    species_aligned <- lapply(species_rasters, function(r) fast_align(to01(r), tmpl))
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
    
    pfpr01 <- to01(fast_align(PFPR_MEAN, tmpl))
    itn01  <- to01(fast_align(ITN_MEAN, tmpl))
    inc_al <- fast_align(INC_MEAN, tmpl, method = "bilinear")
    pop_al <- fast_align(POP_MEAN, tmpl, method = "near")
    
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
      pop_al = pop_al
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
      coverage_table_pfpr_denom_afro(
        overall_tier(layers, case_name),
        layers$v0,
        tier_pfpr_denom(layers, case_name)
      )
    })
  }
  
  for (nm in c("1a", "1b", "2", "3")) {
    local({
      case_name <- nm
      make_outputs(case_name)
    })
  }
}

shinyApp(ui, server)