# Parte 4 e a comparação da Questão 8.
# Não reestima nada: lê os CSV gravados pelos scripts anteriores.
# Prêmio anual: C = 250 * lambda * (1 - p) + 100.

source_irmao <- function(nome) {
  candidatos <- c(nome, file.path("R", nome), file.path("Lab2", "R", nome))
  existe <- candidatos[file.exists(candidatos)]
  if (!length(existe)) {
    stop("Nao achei ", nome, ". Working directory: ", getwd())
  }
  source(existe[1], local = FALSE)
}

source_irmao("00_dados.R")

ler <- function(arquivo) {
  caminho <- file.path(dir_out, arquivo)
  if (!file.exists(caminho)) {
    stop("Falta ", caminho, ". Rode antes o script que grava esse CSV.")
  }
  read.csv(caminho, stringsAsFactors = FALSE)
}

premio <- function(p, lambda) 250 * lambda * (1 - p) + 100

optim_df <- ler("optim.csv")
em_df <- ler("em.csv")
gibbs_df <- ler("gibbs.csv")

# Os três chutes de optim() e do EM concordam; fica o Ponto 1.
optim_p1 <- optim_df[optim_df$inicializacao == "Ponto 1", ]
em_p1 <- em_df[em_df$inicializacao == "Ponto 1", ]

comparacao <- data.frame(
  perfil = optim_p1$perfil,
  p_optim = optim_p1$p_estimado,
  lambda_optim = optim_p1$lambda_estimado,
  p_em = em_p1$p[match(optim_p1$perfil, em_p1$perfil)],
  lambda_em = em_p1$lambda[match(optim_p1$perfil, em_p1$perfil)],
  p_gibbs = gibbs_df$mediana_p[match(optim_p1$perfil, gibbs_df$perfil)],
  lambda_gibbs = gibbs_df$mediana_lambda[match(optim_p1$perfil, gibbs_df$perfil)],
  ic_p_inf = gibbs_df$ic_p_inf[match(optim_p1$perfil, gibbs_df$perfil)],
  ic_p_sup = gibbs_df$ic_p_sup[match(optim_p1$perfil, gibbs_df$perfil)],
  ic_lambda_inf = gibbs_df$ic_lambda_inf[match(optim_p1$perfil, gibbs_df$perfil)],
  ic_lambda_sup = gibbs_df$ic_lambda_sup[match(optim_p1$perfil, gibbs_df$perfil)],
  stringsAsFactors = FALSE
)
comparacao$premio_optim <- premio(comparacao$p_optim, comparacao$lambda_optim)
comparacao$premio_em <- premio(comparacao$p_em, comparacao$lambda_em)
comparacao$premio_gibbs <- premio(comparacao$p_gibbs, comparacao$lambda_gibbs)

# EM parado na primeira atualização, para a Questão 11.
comparacao$premio_em_iter1 <- premio(
  em_p1$p_iter1[match(comparacao$perfil, em_p1$perfil)],
  em_p1$lambda_iter1[match(comparacao$perfil, em_p1$perfil)]
)

print(comparacao)
write.csv(comparacao, file.path(dir_out, "comparacao.csv"), row.names = FALSE)
