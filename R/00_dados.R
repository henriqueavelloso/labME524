# Leitura comum: CSV, perfis e log-verossimilhança do ZIP.
# Os outros scripts começam com source_irmao("00_dados.R").

achar_lab2 <- function() {
  dir <- normalizePath(getwd(), mustWork = TRUE)
  for (i in seq_len(6)) {
    if (file.exists(file.path(dir, "visitas_medicas.csv"))) {
      return(dir)
    }
    candidato <- file.path(dir, "Lab2", "visitas_medicas.csv")
    if (file.exists(candidato)) {
      return(normalizePath(file.path(dir, "Lab2")))
    }
    pai <- dirname(dir)
    if (pai == dir) break
    dir <- pai
  }
  stop(
    "Nao achei visitas_medicas.csv. ",
    "No Posit, deixe o working directory na pasta Lab2 ou na pasta do projeto."
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Instale o pacote ggplot2 antes de rodar os scripts.")
}
library(ggplot2)

lab2 <- achar_lab2()
dir_fig <- file.path(lab2, "relatorio", "figures")
dir_out <- file.path(lab2, "relatorio", "output")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

dados <- read.csv(
  file.path(lab2, "visitas_medicas.csv"),
  stringsAsFactors = FALSE
)
dados$perfil <- factor(dados$perfil, levels = c("baixo", "medio", "alto"))
perfis <- split(dados$visitas, dados$perfil)

# Log-verossimilhança observada, sem o fatorial (não depende de p nem de lambda).
loglik_zip <- function(y, p, lambda) {
  n0 <- sum(y == 0)
  ypos <- y[y > 0]
  n0 * log(p + (1 - p) * exp(-lambda)) +
    sum(log(1 - p) - lambda + ypos * log(lambda))
}

# Histograma pedido na seção de dados.
g_hist <- ggplot(dados, aes(x = visitas)) +
  geom_histogram(binwidth = 1, boundary = -0.5, fill = "grey35", color = "white") +
  facet_wrap(~ perfil, nrow = 1) +
  labs(
    title = "Consultas anuais por perfil de risco",
    x = "Número de consultas",
    y = "Clientes"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(dir_fig, "hist_visitas.png"),
  g_hist,
  width = 8, height = 3.2, dpi = 130
)
