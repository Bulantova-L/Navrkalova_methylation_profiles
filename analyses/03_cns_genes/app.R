# Interaktivní prohlížeč metylačních profilů (CNS geny)
#
# Ocekava strukturu:
#   results/<sample>/<sample>_cns_mapped.bed
# se sloupci (bez hlavicky): chr, start, end, gene, methylation, coverage
#
# Spusteni: shiny::runApp("cesta/k/teto/slozce")

library(shiny)
library(data.table)
library(ggplot2)
library(reshape2)

# ---- Paleta ----
CAT_COLORS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
DIVERGING_LOW  <- "#184f95"
DIVERGING_MID  <- "#f0efec"
DIVERGING_HIGH <- "#b23b3b"
GRID_COLOR   <- "#e1e0d9"
AXIS_COLOR   <- "#c3c2b7"
MUTED_COLOR  <- "#898781"
INK_COLOR    <- "#0b0b0b"

# ---- Nacteni dat ----

read_cns_file <- function(file_path) {
  sample_id <- sub("_cns_mapped\\.bed$", "", basename(file_path))

  dt <- fread(file_path, header = FALSE)
  if (ncol(dt) != 6) {
    stop("Neocekavany pocet sloupcu v souboru: ", file_path, " (", ncol(dt), ")")
  }
  setnames(dt, c("chr", "start", "end", "gene", "methylation", "coverage"))

  dt[, sample := sample_id]
  dt[, methylation := as.numeric(methylation)]
  dt[, coverage := as.numeric(coverage)]
  dt[, probe_id := paste0(chr, ":", start, "-", end)]
  dt[, pos := (start + end) / 2]

  dt
}

# runApp() nastavi pracovni adresar na slozku aplikace (analyses/03_cns_genes) -> results/ je o 2 urovne vys
results_dir <- if (dir.exists("results")) "results" else "../../results"
bed_files <- list.files(results_dir, pattern = "_cns_mapped\\.bed$", recursive = TRUE, full.names = TRUE)

if (length(bed_files) == 0) {
  stop("Ve slozce '", results_dir, "' nebyl nalezen zadny soubor '*_cns_mapped.bed'.")
}

all_data <- rbindlist(lapply(bed_files, read_cns_file), fill = TRUE)

genes   <- sort(unique(all_data$gene))
samples <- sort(unique(all_data$sample))

sample_colors <- setNames(CAT_COLORS[seq_along(samples)], samples)

# ---- UI ----

ui <- fluidPage(
  titlePanel("Metylacni profily vzorku - CNS geny"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", "Gen:", choices = genes, selected = genes[1]),
      hr(),
      helpText("Korelacni heatmapa porovnava metylacni hodnoty jednotlivych probu vybraneho genu mezi vzorky."),
      helpText("Metylacni profil zobrazuje metylaci podle genomicke pozice pro vsechny vzorky."),
      hr(),
      tags$b("Vzorky v datech:"),
      tags$ul(lapply(samples, tags$li))
    ),
    mainPanel(
      h4("Korelace mezi vzorky"),
      plotOutput("corr_heatmap", height = "420px"),
      hr(),
      h4("Metylacni profil genu"),
      plotOutput("gene_profile", height = "700px")
    )
  )
)

# ---- Server ----

server <- function(input, output, session) {

  gene_data <- reactive({
    all_data[gene == input$gene]
  })

  gene_wide <- reactive({
    d <- gene_data()
    dcast(d, probe_id + pos ~ sample, value.var = "methylation")
  })

  output$corr_heatmap <- renderPlot({
    w <- gene_wide()
    present_samples <- intersect(samples, names(w))

    validate(need(length(present_samples) >= 2,
                   "Pro vypocet korelace jsou potreba alespon 2 vzorky s daty pro tento gen."))

    mat <- as.matrix(w[, present_samples, drop = FALSE])
    cor_mat <- suppressWarnings(cor(mat, use = "pairwise.complete.obs"))

    cor_df <- melt(cor_mat, varnames = c("sample_x", "sample_y"), value.name = "correlation")
    cor_df$sample_x <- factor(cor_df$sample_x, levels = present_samples)
    cor_df$sample_y <- factor(cor_df$sample_y, levels = rev(present_samples))

    ggplot(cor_df, aes(x = sample_x, y = sample_y, fill = correlation)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.2f", correlation)), color = INK_COLOR, size = 4.2) +
      scale_fill_gradient2(
        low = DIVERGING_LOW, mid = DIVERGING_MID, high = DIVERGING_HIGH,
        midpoint = 0, limits = c(-1, 1), name = "Korelace"
      ) +
      labs(x = NULL, y = NULL, title = paste0("Korelace metylace mezi vzorky - ", input$gene)) +
      coord_fixed() +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = INK_COLOR),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(color = INK_COLOR)
      )
  })

    output$gene_profile <- renderPlot({
    d <- gene_data()
    validate(need(nrow(d) > 0, "Pro tento gen nejsou k dispozici zadna data."))

    d$sample <- factor(d$sample, levels = samples)

    ggplot(d, aes(x = pos, xend = pos, y = 0, yend = methylation, color = sample)) +
      geom_segment(linewidth = 1) +
      facet_wrap(~ sample, ncol = 1, strip.position = "left") +
      scale_color_manual(values = sample_colors, guide = "none") +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
      labs(
        x = "Genomicka pozice (chr, stred probu)",
        y = "Metylace",
        title = paste0("Metylacni profil - ", input$gene)
      ) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = GRID_COLOR),
        axis.line = element_line(color = AXIS_COLOR),
        axis.text = element_text(color = MUTED_COLOR),
        plot.title = element_text(color = INK_COLOR),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(color = INK_COLOR, face = "bold", hjust = 0)
      )
  })
}
shinyApp(ui, server)
