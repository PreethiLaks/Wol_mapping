library(shiny)
library(terra)
library(sf)

options(terra.plot.legend = TRUE)
terraOptions(memfrac = 0.8, progress = 0, tempdir = tempdir())

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
  gadm = file.path("data", "gadm_410.gpkg"),
  # ── NEW: gambiae included for dominance calculation only, not mapped ────────
  gam  = file.path("data", "2017_Anopheles_gambiae.Mean_Decompressed.geotiff")
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
  out[!is.na(valid_mask) & x >= breaks[3]]                  <- labels[3]
  out
}

# FIX: count countries using terra::relate instead of repeated extract()
coverage_table_countries <- function(tier_r, v0) {
  count_countries_fast <- function(mask_r, v0) {
    # rasterize v0 to get a country-ID raster, then check overlap with mask
    cid <- rasterize(v0, mask_r, field = seq_len(nrow(v0)))
    hit <- !is.na(cid) & !is.na(mask_r)
    ids <- unique(na.omit(values(ifel(hit, cid, NA))))
    length(ids)
  }
  
  m1   <- tier_mask_value(tier_r, 1)
  m2   <- tier_mask_value(tier_r, 2)
  m3   <- tier_mask_value(tier_r, 3)
  m123 <- tier_any_mask(tier_r)
  
  data.frame(
    Tier = c("Tier 1", "Tier 2", "Tier 3", "Tier 1-3 (total)"),
    Countries_covered = c(
      count_countries_fast(m1,   v0),
      count_countries_fast(m2,   v0),
      count_countries_fast(m3,   v0),
      count_countries_fast(m123, v0)
    ),
    stringsAsFactors = FALSE
  )
}

impact_table <- function(tier_r, layers, burden_r, pop_r) {
  bur_valid <- mask(burden_r, layers$burden_mask)
  pop_valid <- mask(pop_r,    layers$burden_mask)
  
  # Stack all tier masks and do one global() call per variable
  m1   <- tier_mask_value(tier_r, 1)
  m2   <- tier_mask_value(tier_r, 2)
  m3   <- tier_mask_value(tier_r, 3)
  m123 <- tier_any_mask(tier_r)
  
  pop_stack <- c(
    mask(pop_valid, m1),
    mask(pop_valid, m2),
    mask(pop_valid, m3),
    mask(pop_valid, m123)
  )
  bur_stack <- c(
    mask(bur_valid, m1),
    mask(bur_valid, m2),
    mask(bur_valid, m3),
    mask(bur_valid, m123)
  )
  
  # One global() call each for pop and burden stacks
  pop_sums <- as.numeric(global(pop_stack, "sum", na.rm = TRUE)[, 1])
  bur_sums <- as.numeric(global(bur_stack, "sum", na.rm = TRUE)[, 1])
  pop_sums[!is.finite(pop_sums)] <- 0
  bur_sums[!is.finite(bur_sums)] <- 0
  
  # ── DENOMINATOR 1: WHO AFRO (original) ──────────────────────────────────────
  afro_cases_all <- as.numeric(global(mask(bur_valid, layers$v0), "sum", na.rm = TRUE)[1, 1])
  if (!is.finite(afro_cases_all) || afro_cases_all <= 0) afro_cases_all <- NA_real_
  
  # ── DENOMINATOR 2 (NEW): cases within species-present range (pres_mask) ────
  # Rationale: pres_mask (base_p > 0) defines the outer boundary of where
  # Wolbachia deployment is biologically relevant for the selected species.
  # Dominance is now computed over all 5 vectors (including An. gambiae),
  # so dom_share already reflects the correct competitive context.
  # This answers "of all cases in the species range, what % are captured?"
  species_range_denom <- as.numeric(
    global(mask(bur_valid, layers$pres_mask), "sum", na.rm = TRUE)[1, 1]
  )
  if (!is.finite(species_range_denom) || species_range_denom <= 0) species_range_denom <- NA_real_
  # ────────────────────────────────────────────────────────────────────────────
  
  fmt_pct <- function(x, total) {
    if (is.na(total)) return("NA")
    formatC(100 * x / total, format = "f", digits = 2)
  }
  
  fmt_cases_row <- function(x) {
    paste0(
      format(round(x), big.mark = ","),
      "  |  ",
      fmt_pct(x, afro_cases_all), "% (AFRO)",
      "  /  ",
      fmt_pct(x, species_range_denom), "% (species range)"
    )
  }
  
  data.frame(
    Metric = c(
      "Population benefiting - Tier 1",
      "Population benefiting - Tier 2",
      "Population benefiting - Tier 3",
      "Population benefiting - Tier 1-3 combined",
      "Estimated mean annual cases - Tier 1",
      "Estimated mean annual cases - Tier 2",
      "Estimated mean annual cases - Tier 3",
      "Estimated mean annual cases - Tier 1-3 combined"
    ),
    Value = c(
      format(round(pop_sums[1]), big.mark = ","),
      format(round(pop_sums[2]), big.mark = ","),
      format(round(pop_sums[3]), big.mark = ","),
      format(round(pop_sums[4]), big.mark = ","),
      fmt_cases_row(bur_sums[1]),   # <- shows both % side by side
      fmt_cases_row(bur_sums[2]),
      fmt_cases_row(bur_sums[3]),
      fmt_cases_row(bur_sums[4])
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
      dx <- diff(xr); dy <- diff(yr)
      rx <- res(r)[1]; ry <- res(r)[2]
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
  tier_cols   <- c("1" = "#1a9850", "2" = "#f16913", "3" = "#d73027")
  tier_labels <- c("1" = "Overall Tier 1", "2" = "Overall Tier 2", "3" = "Overall Tier 3")
  
  vals <- sort(unique(na.omit(as.vector(values(tier_r)))))
  vals <- vals[vals %in% 1:3]
  
  if (length(vals) == 0) {
    plot.new(); title(paste0(title_text, "\n(No eligible cells)")); return(invisible(NULL))
  }
  
  z      <- compute_zoom(tier_r, layers$v0, auto_zoom, zoom_pad)
  af_ext <- ext(layers$v0)
  xlim   <- if (is.null(z$xlim)) c(xmin(af_ext), xmax(af_ext)) else z$xlim
  ylim   <- if (is.null(z$ylim)) c(ymin(af_ext), ymax(af_ext)) else z$ylim
  e_view <- ext(xlim[1], xlim[2], ylim[1], ylim[2])
  
  v0_view <- crop(layers$v0, e_view)
  v1_view <- crop(layers$v1, e_view)
  
  plot(v0_view, xlim = xlim, ylim = ylim, col = "white", border = NA, axes = TRUE, main = title_text)
  plot(tier_r,  add = TRUE, col = tier_cols, breaks = c(0.5, 1.5, 2.5, 3.5), legend = FALSE, colNA = NA)
  plot(v1_view, add = TRUE, lwd = 0.6, col = NA, border = "#BFD7FF")
  plot(v0_view, add = TRUE, lwd = 1.2, col = NA, border = "#6996E8")
  
  legend("bottomleft",
         legend = tier_labels[as.character(vals)],
         fill   = tier_cols[as.character(vals)],
         border = NA, bty = "n", cex = 0.95, title = "Overall Tier"
  )
}

# ----------------------------
# Data loading  (runs ONCE at startup)
# ----------------------------
PFPR_MEAN       <- rast(paths$pfpr)
ITN_MEAN        <- rast(paths$itn)
INC_MEAN        <- rast(paths$inc)
TT_MEAN         <- rast(paths$tt)
POP_MEAN        <- rast(paths$pop)
species_rasters <- list(
  arabiensis = rast(paths$ara),
  coluzzii   = rast(paths$col),
  moucheti   = rast(paths$mou),
  stephensi  = rast(paths$ste),
  gambiae    = rast(paths$gam)   # ── NEW: added to enable dropdown mapping
)

afro_iso <- c(
  "DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM","COG","CIV","COD","GNQ","ERI",
  "SWZ","ETH","GAB","GMB","GHA","GIN","GNB","KEN","LSO","LBR","MDG","MWI","MLI","MRT","MUS","MOZ",
  "NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","ZAF","SSD","TGO","UGA","TZA","ZMB","ZWE"
)

a0     <- st_read(paths$gadm, layer = "ADM_0", quiet = TRUE)
a1     <- st_read(paths$gadm, layer = "ADM_1", quiet = TRUE)
ADMIN0 <- vect(a0[a0$GID_0 %in% afro_iso, ])
ADMIN1 <- vect(a1[a1$GID_0 %in% afro_iso, ])

# ----------------------------
# Precompute ALL aligned layers at startup  (runs ONCE)
# ----------------------------
ADMIN0_PFPR    <- if (!identical(crs(ADMIN0), crs(PFPR_MEAN))) project(ADMIN0, crs(PFPR_MEAN)) else ADMIN0
TEMPLATE       <- crop(PFPR_MEAN, ADMIN0_PFPR, snap = "out")
ADMIN0_ALIGNED <- if (!identical(crs(ADMIN0), crs(TEMPLATE))) project(ADMIN0, crs(TEMPLATE)) else ADMIN0
ADMIN1_ALIGNED <- if (!identical(crs(ADMIN1), crs(TEMPLATE))) project(ADMIN1, crs(TEMPLATE)) else ADMIN1

PFPR_ALIGNED   <- to01(fast_align(PFPR_MEAN, TEMPLATE))
ITN_ALIGNED    <- to01(fast_align(ITN_MEAN,  TEMPLATE))
INC_ALIGNED    <- fast_align(INC_MEAN, TEMPLATE, method = "bilinear")
TT_ALIGNED     <- fast_align(TT_MEAN,  TEMPLATE, method = "bilinear")
POP_ALIGNED    <- fast_align(POP_MEAN, TEMPLATE, method = "near")
# ── All 5 species aligned at startup including gambiae ──────────────────────
# gambiae is now in SPECIES_ALIGNED so it appears in the dropdown and
# its presence probability is always included in dominance calculations
# via setdiff() — when gambiae is selected it competes against the other 4,
# and when another species is selected gambiae is in others_sum automatically.
SPECIES_ALIGNED <- lapply(species_rasters, function(r) fast_align(to01(r), TEMPLATE))
# ────────────────────────────────────────────────────────────────────────────

# Precompute burden raster globally (pop * incidence / 1000) — species-independent
BURDEN_GLOBAL  <- POP_ALIGNED * (INC_ALIGNED / 1000)
names(BURDEN_GLOBAL) <- "cases_year"

# ----------------------------
# Use-case rules
# ----------------------------
use_cases <- list(
  "1a" = list(
    itn_breaks  = c(-Inf, 0.60, 0.80, Inf), itn_labels  = c(1, 2, 3),
    pfpr_breaks = c(-Inf, 0.15, 0.40, Inf), pfpr_labels = c(1, 2, 3)
  ),
  "1b" = list(
    itn_breaks  = c(-Inf, 0.70, 0.80, Inf), itn_labels  = c(1, 2, 3),
    pfpr_breaks = c(-Inf, 0.15, 0.40, Inf), pfpr_labels = c(3, 2, 1)
  ),
  "2"  = list(
    itn_breaks  = c(-Inf, 0.60, 0.80, Inf), itn_labels  = c(3, 2, 1),
    pfpr_breaks = c(-Inf, 0.15, 0.40, Inf), pfpr_labels = c(3, 2, 1)
  ),
  "3"  = list(
    itn_breaks  = c(-Inf, 0.60, 0.80, Inf), itn_labels  = c(3, 2, 1),
    pfpr_breaks = c(-Inf, 0.05, 0.15, 0.40), pfpr_labels = c(1, 2, 3)
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
                  c("arabiensis", "coluzzii", "moucheti", "stephensi", "gambiae"),
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
        tabPanel("Use case 4",  plotOutput("map_4",  height = 650), h4("Impact Summary"), tableOutput("impact_4"),  h4("Coverage"), tableOutput("cov_4"))
      )
    )
  )
)

# ----------------------------
# Server
# ----------------------------
server <- function(input, output, session) {
  
  # ── OPTIMIZATION 1: species-level layers ──────────────────────────────────
  # Recomputes ONLY when species changes, not on zoom/padding tweaks.
  species_layers <- reactive({
    sp       <- input$species
    base_p   <- SPECIES_ALIGNED[[sp]]
    base_p0  <- base_p; base_p0[is.na(base_p0)] <- 0
    
    others_stack <- rast(lapply(
      setdiff(names(SPECIES_ALIGNED), sp),   # all 5 minus selected — gambiae included automatically
      function(nm) { r <- SPECIES_ALIGNED[[nm]]; r[is.na(r)] <- 0; r }
    ))
    others_sum <- app(others_stack, sum)
    dom_share  <- safe_div(base_p0, base_p0 + others_sum)
    
    pres_mask   <- ifel(!is.na(base_p) & base_p > 0, 1, NA)
    data_mask   <- ifel(!is.na(PFPR_ALIGNED) & !is.na(ITN_ALIGNED) & !is.na(base_p) & (base_p > 0), 1, NA)
    burden_mask <- ifel(!is.na(INC_ALIGNED) & !is.na(POP_ALIGNED) & (POP_ALIGNED > 0) & (INC_ALIGNED > 1e-12), 1, NA)
    
    list(
      v0          = ADMIN0_ALIGNED,
      v1          = ADMIN1_ALIGNED,
      pfpr01      = PFPR_ALIGNED,
      itn01       = ITN_ALIGNED,
      pfpr_afro   = mask(PFPR_ALIGNED, ADMIN0_ALIGNED),
      pres_mask   = pres_mask,
      dom_share   = dom_share,
      data_mask   = data_mask,
      burden_mask = burden_mask,
      inc_al      = INC_ALIGNED,
      pop_al      = POP_ALIGNED,
      tt_al       = TT_ALIGNED
    )
  })
  
  # ── OPTIMIZATION 2: burden is species-independent, computed once ──────────
  # (BURDEN_GLOBAL already computed at startup above)
  
  # ── OPTIMIZATION 3: memoised tier rasters — one per use-case per species ──
  # bindCache() caches the result keyed on (species, case_name).
  # The tier raster is computed once and reused for map + impact + coverage.
  
  tier_for <- function(case_name) {
    reactive({
      layers <- species_layers()
      
      tv <- {
        # Proportional dominance share with minimum presence filter ───────────
        # dom_share = selected species / all 5 species combined (per pixel).
        # Tiers reflect how much of the local vector pool the species represents.
        # base_p >= 0.10 floor excludes pixels where MAP model has <10%
        # confidence the species is present — below this, deployment is not
        # meaningful regardless of dominance share.
        ds <- layers$dom_share
        bp <- SPECIES_ALIGNED[[input$species]]; bp[is.na(bp)] <- 0
        valid_dom <- !is.na(layers$data_mask) & !is.na(layers$pres_mask) & bp >= 0.10
        out <- rast(TEMPLATE); values(out) <- NA
        out[valid_dom & ds >= 0.70] <- 1
        out[valid_dom & ds >= 0.40 & ds < 0.70] <- 2
        out[valid_dom & ds <  0.40] <- 3
        # ────────────────────────────────────────────────────────────────────
        out
      }
      
      cfg <- use_cases[[case_name]]
      ti  <- make_tier(layers$itn01,  layers$data_mask, cfg$itn_breaks,  cfg$itn_labels)
      tp  <- make_tier(layers$pfpr01, layers$data_mask, cfg$pfpr_breaks, cfg$pfpr_labels)
      
      out <- rast(TEMPLATE); values(out) <- NA
      valid <- !is.na(layers$data_mask) & !is.na(layers$pres_mask)
      out[(tv == 1) & (ti == 1) & (tp == 1) & valid] <- 1
      out[(tv == 2) & (ti == 2) & (tp == 2) & valid] <- 2
      out[(tv == 3) & (ti == 3) & (tp == 3) & valid] <- 3
      mask(out, layers$v0)
    }) |> bindCache(input$species, case_name)
  }
  
  tier_1a <- tier_for("1a")
  tier_1b <- tier_for("1b")
  tier_2  <- tier_for("2")
  tier_3  <- tier_for("3")
  
  # UC4 tier — always stephensi-based, no species dropdown dependency
  tier_4 <- reactive({
    layers <- species_layers()   # still needed for v0, pfpr01, tt_al
    ste <- mask(SPECIES_ALIGNED$stephensi, layers$v0)
    tt  <- mask(layers$tt_al,  layers$v0)
    pf  <- mask(layers$pfpr01, layers$v0)
    
    ste_tier  <- ifel(ste >= 0.75, 1, ifel(ste >= 0.15, 2, ifel(ste >= 0.05, 3, NA)))
    tt_tier   <- ifel(tt  <= 60,  1, ifel(tt  <= 180,  2, ifel(tt  <= 360,  3, NA)))
    pfpr_tier <- ifel(pf  >= 0.40, 1, ifel(pf  >= 0.15, 2, ifel(pf  >= 0,   3, NA)))
    
    out <- rast(TEMPLATE); values(out) <- NA
    out[ste_tier == 1 & tt_tier == 1 & pfpr_tier == 1] <- 1
    out[ste_tier %in% c(1,2) & tt_tier %in% c(1,2) & pfpr_tier %in% c(1,2) & is.na(out)] <- 2
    out[!is.na(ste_tier) & !is.na(tt_tier) & !is.na(pfpr_tier) & is.na(out)] <- 3
    mask(out, layers$v0)
  }) |> bindCache(input$species)   # UC4 result also invalidates if species changes (v0 mask path)
  
  # ── OPTIMIZATION 4: zoom is a DISPLAY-only concern ─────────────────────────
  # Separate zoom-reactive so it doesn't re-trigger heavy raster computations.
  zoom_params <- reactive(list(auto = input$auto_zoom, pad = input$zoom_pad))
  
  # ── Wire outputs ────────────────────────────────────────────────────────────
  make_outputs <- function(case_name, tier_r_reactive) {
    output[[paste0("map_", case_name)]] <- renderPlot({
      layers <- species_layers()
      zp     <- zoom_params()
      plot_tier_map(
        tier_r_reactive(),
        layers,
        paste0("Use case ", case_name, " - ", input$species),
        zp$auto, zp$pad
      )
    })
    
    output[[paste0("impact_", case_name)]] <- renderTable({
      layers <- species_layers()
      impact_table(tier_r_reactive(), layers, BURDEN_GLOBAL, layers$pop_al)
    })
    
    output[[paste0("cov_", case_name)]] <- renderTable({
      layers <- species_layers()
      coverage_table_countries(tier_r_reactive(), layers$v0)
    })
  }
  
  make_outputs("1a", tier_1a)
  make_outputs("1b", tier_1b)
  make_outputs("2",  tier_2)
  make_outputs("3",  tier_3)
  
  output$map_4 <- renderPlot({
    layers <- species_layers()
    zp     <- zoom_params()
    plot_tier_map(
      tier_4(),
      layers,
      "Use case 4 - An. stephensi urban/peri-urban vulnerability",
      zp$auto, zp$pad
    )
  })
  
  output$impact_4 <- renderTable({
    layers <- species_layers()
    impact_table(tier_4(), layers, BURDEN_GLOBAL, layers$pop_al)
  })
  
  output$cov_4 <- renderTable({
    layers <- species_layers()
    coverage_table_countries(tier_4(), layers$v0)
  })
}

shinyApp(ui, server)