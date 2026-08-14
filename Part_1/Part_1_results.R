## LIFE748 Assessment 2 - Part 1: Results visualisation

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)



dir.create("Figures", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)


## Read QUAST assembly comparison
quast <- read.delim(
  "Results/QUAST_Comparison/report.tsv",
  check.names = FALSE
)

head(quast)


## Select key assembly metrics for the report
assembly_metrics <- quast %>%
  filter(
    Assembly %in% c(
      "# contigs",
      "Largest contig",
      "Total length",
      "GC (%)",
      "N50",
      "L50",
      "# N's per 100 kbp"
    )
  )

assembly_metrics



## Extract runtime and memory information

extract_runtime <- function(file) {
  
  lines <- readLines(file)
  
  elapsed_line <- grep(
    "Elapsed \\(wall clock\\)",
    lines,
    value = TRUE
  )
  
  memory_line <- grep(
    "Maximum resident set size",
    lines,
    value = TRUE
  )
  
  elapsed <- sub(
    ".*: ",
    "",
    elapsed_line
  )
  
  memory_kb <- as.numeric(
    sub(
      ".*: ",
      "",
      memory_line
    )
  )
  
  list(
    elapsed = elapsed,
    memory_gb = memory_kb / 1000000
  )
}


flye_runtime <- extract_runtime(
  "Results/flye_runtime.txt"
)

spades_runtime <- extract_runtime(
  "Results/spades_runtime.txt"
)

prokka_runtime <- extract_runtime(
  "Results/prokka_runtime.txt"
)

bakta_runtime <- extract_runtime(
  "Results/bakta_runtime.txt"
)



##Convert elapsed time from mm:ss to decimal minutes

time_to_minutes <- function(x) {
  
  parts <- strsplit(x, ":")[[1]]
  
  if (length(parts) == 2) {
    
    as.numeric(parts[1]) +
      as.numeric(parts[2]) / 60
    
  } else {
    
    as.numeric(parts[1]) * 60 +
      as.numeric(parts[2]) +
      as.numeric(parts[3]) / 60
  }
}


flye_minutes <- time_to_minutes(
  flye_runtime$elapsed
)

spades_minutes <- time_to_minutes(
  spades_runtime$elapsed
)

prokka_minutes <- time_to_minutes(
  prokka_runtime$elapsed
)

bakta_minutes <- time_to_minutes(
  bakta_runtime$elapsed
)


## Assembly benchmarking table

assembly_table <- data.frame(
  Metric = c(
    "Contigs >=500 bp",
    "Total assembly length (bp)",
    "Largest contig (bp)",
    "N50 (bp)",
    "L50",
    "GC content (%)",
    "Ns per 100 kbp",
    "Runtime (min)",
    "Peak memory (GB)"
  ),
  
  Flye = c(
    quast$Flye[quast$Assembly == "# contigs"],
    quast$Flye[quast$Assembly == "Total length"],
    quast$Flye[quast$Assembly == "Largest contig"],
    quast$Flye[quast$Assembly == "N50"],
    quast$Flye[quast$Assembly == "L50"],
    quast$Flye[quast$Assembly == "GC (%)"],
    quast$Flye[quast$Assembly == "# N's per 100 kbp"],
    flye_minutes,
    flye_runtime$memory_gb
  ),
  
  SPAdes = c(
    quast$SPAdes[quast$Assembly == "# contigs"],
    quast$SPAdes[quast$Assembly == "Total length"],
    quast$SPAdes[quast$Assembly == "Largest contig"],
    quast$SPAdes[quast$Assembly == "N50"],
    quast$SPAdes[quast$Assembly == "L50"],
    quast$SPAdes[quast$Assembly == "GC (%)"],
    quast$SPAdes[quast$Assembly == "# N's per 100 kbp"],
    spades_minutes,
    spades_runtime$memory_gb
  )
)

assembly_table

write.csv(
  assembly_table,
  "Tables/assembly_comparison.csv",
  row.names = FALSE
)



##Parse annotation summaries

prokka_lines <- readLines(
  "Annotation/Prokka/GN3_Prokka.txt"
)

get_prokka_value <- function(name) {
  
  line <- grep(
    paste0("^", name, ":"),
    prokka_lines,
    value = TRUE
  )
  
  if (length(line) == 0) {
    return(NA)
  }
  
  as.numeric(
    trimws(
      sub(
        paste0("^", name, ":"),
        "",
        line
      )
    )
  )
}


bakta_lines <- readLines(
  "Annotation/Bakta/GN3_Bakta.txt"
)

get_bakta_value <- function(name) {
  
  line <- grep(
    paste0("^", name, ":"),
    bakta_lines,
    value = TRUE
  )
  
  if (length(line) == 0) {
    return(NA)
  }
  
  as.numeric(
    trimws(
      sub(
        paste0("^", name, ":"),
        "",
        line
      )
    )
  )
}


## Annotation comparison table

annotation_table <- data.frame(
  Feature = c(
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "ncRNA regions",
    "CRISPR arrays",
    "sORFs",
    "oriC",
    "oriV",
    "oriT",
    "Runtime (min)",
    "Peak memory (GB)"
  ),
  
  Prokka = c(
    get_prokka_value("CDS"),
    get_prokka_value("rRNA"),
    get_prokka_value("tRNA"),
    get_prokka_value("tmRNA"),
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    prokka_minutes,
    prokka_runtime$memory_gb
  ),
  
  Bakta = c(
    get_bakta_value("CDSs"),
    get_bakta_value("rRNAs"),
    get_bakta_value("tRNAs"),
    get_bakta_value("tmRNAs"),
    get_bakta_value("ncRNAs"),
    get_bakta_value("ncRNA regions"),
    get_bakta_value("CRISPR arrays"),
    get_bakta_value("sORFs"),
    get_bakta_value("oriCs"),
    get_bakta_value("oriVs"),
    get_bakta_value("oriTs"),
    bakta_minutes,
    bakta_runtime$memory_gb
  )
)

annotation_table

write.csv(
  annotation_table,
  "Tables/annotation_comparison.csv",
  row.names = FALSE
)



## Assembly performance figure


assembly_plot_data <- data.frame(
  Assembler = factor(
    c("Flye", "SPAdes"),
    levels = c("Flye", "SPAdes")
  ),
  
  Contigs = c(
    quast$Flye[quast$Assembly == "# contigs"],
    quast$SPAdes[quast$Assembly == "# contigs"]
  ),
  
  Runtime = c(
    flye_minutes,
    spades_minutes
  ),
  
  Memory = c(
    flye_runtime$memory_gb,
    spades_runtime$memory_gb
  )
)

## Consistent colours for both assemblers
assembler_cols <- c(
  "Flye" = "#00BFC4",
  "SPAdes" = "#F8766D"
)


## Panel A - contig count
p1 <- ggplot(
  assembly_plot_data,
  aes(
    x = Assembler,
    y = Contigs,
    fill = Assembler
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = Contigs),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = assembler_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "A",
    x = NULL,
    y = "Contigs ≥500 bp"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )


## Panel B - runtime
p2 <- ggplot(
  assembly_plot_data,
  aes(
    x = Assembler,
    y = Runtime,
    fill = Assembler
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f", Runtime)),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = assembler_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "B",
    x = NULL,
    y = "Runtime (min)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )


## Panel C -peak memory
p3 <- ggplot(
  assembly_plot_data,
  aes(
    x = Assembler,
    y = Memory,
    fill = Assembler
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f", Memory)),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = assembler_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "C",
    x = NULL,
    y = "Peak memory (GB)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )


assembly_figure <- p1 + p2 + p3 +
  plot_annotation(
    title = "Genome assembly performance",
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold"
      )
    )
  )

assembly_figure


## SAVE
ggsave(
  "Figures/assembly_performance.png",
  assembly_figure,
  width = 10,
  height = 4,
  dpi = 300
)