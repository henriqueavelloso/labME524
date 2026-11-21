# Roda as quatro partes na ordem e grava figuras e CSV.
#
# No Posit Cloud, deixe o working directory na pasta Lab2
# (a que contém visitas_medicas.csv) ou na pasta do projeto, e então:
#   source("R/run_all.R")
# ou, se o projeto abrir na pasta acima de Lab2:
#   source("Lab2/R/run_all.R")
#
# Saída: relatorio/figures/ e relatorio/output/.
# Baixe essas duas pastas de volta para preencher o relatório.

source_irmao <- function(nome) {
  candidatos <- c(
    nome,
    file.path("R", nome),
    file.path("Lab2", "R", nome)
  )
  existe <- candidatos[file.exists(candidatos)]
  if (!length(existe)) {
    stop("Nao achei ", nome, ". Working directory: ", getwd())
  }
  source(existe[1], local = FALSE)
}

source_irmao("01_optim_nr.R")
source_irmao("02_em.R")
source_irmao("03_gibbs.R")
source_irmao("04_tabelas.R")
