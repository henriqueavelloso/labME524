# Parte 1. optim() e 10 iterações de Newton-Raphson.
# Mesmos três chutes do relatório: alpha = -2, -0.5, 1 e beta = 0, 1.5, 2.

source_irmao <- function(nome) {
  candidatos <- c(nome, file.path("R", nome), file.path("Lab2", "R", nome))
  existe <- candidatos[file.exists(candidatos)]
  if (!length(existe)) {
    stop("Nao achei ", nome, ". Working directory: ", getwd())
  }
  source(existe[1], local = FALSE)
}

source_irmao("00_dados.R")
set.seed(248214)

p_func <- function(alpha) 1 / (1 + exp(-alpha))
lambda_func <- function(beta) exp(beta)

negative_logv <- function(params, y) {
  -loglik_zip(y, p_func(params[1]), lambda_func(params[2]))
}

pontos_iniciais <- data.frame(
  inicializacao = c("Ponto 1", "Ponto 2", "Ponto 3"),
  alpha = c(-2, -0.5, 1),
  beta = c(0, 1.5, 2)
)

resultados <- do.call(rbind, lapply(names(perfis), function(perfil) {
  y <- perfis[[perfil]]
  do.call(rbind, lapply(seq_len(nrow(pontos_iniciais)), function(i) {
    alpha_0 <- pontos_iniciais$alpha[i]
    beta_0 <- pontos_iniciais$beta[i]
    ajuste <- optim(par = c(alpha_0, beta_0), fn = negative_logv, y = y)
    data.frame(
      perfil = perfil,
      inicializacao = pontos_iniciais$inicializacao[i],
      p_inicial = p_func(alpha_0),
      lambda_inicial = lambda_func(beta_0),
      p_estimado = p_func(ajuste$par[1]),
      lambda_estimado = lambda_func(ajuste$par[2]),
      logv = -ajuste$value,
      convergencia = ajuste$convergence,
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(resultados) <- NULL
print(resultados)
write.csv(resultados, file.path(dir_out, "optim.csv"), row.names = FALSE)

g_optim <- ggplot(resultados) +
  geom_segment(
    aes(
      x = p_inicial, y = lambda_inicial,
      xend = p_estimado, yend = lambda_estimado,
      color = inicializacao
    ),
    arrow = arrow(length = grid::unit(0.15, "cm")),
    alpha = 0.7
  ) +
  geom_point(aes(x = p_inicial, y = lambda_inicial, color = inicializacao),
             size = 3, shape = 1, stroke = 1.2) +
  geom_point(aes(x = p_estimado, y = lambda_estimado, color = inicializacao),
             size = 3) +
  facet_wrap(~ perfil, nrow = 1) +
  labs(
    title = "Estimação por máxima verossimilhança",
    subtitle = "Três pontos iniciais por perfil",
    x = "p — probabilidade de zero estrutural",
    y = expression(lambda ~ "— taxa da Poisson"),
    color = "Inicialização"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(dir_fig, "grafico_optim.png"), g_optim,
       width = 10, height = 4, dpi = 130)

score_zip <- function(p, lambda, y) {
  n_0 <- sum(y == 0)
  n_pos <- sum(y > 0)
  S_pos <- sum(y[y > 0])
  A <- p + (1 - p) * exp(-lambda)
  c(
    n_0 * (1 - exp(-lambda)) / A - n_pos / (1 - p),
    n_0 * (-exp(-lambda) + p * exp(-lambda)) / A - n_pos + S_pos / lambda
  )
}

hessiana_zip <- function(p, lambda, y) {
  n_0 <- sum(y == 0)
  n_pos <- sum(y > 0)
  S_pos <- sum(y[y > 0])
  A <- p + (1 - p) * exp(-lambda)
  matrix(
    c(
      -n_0 * (1 - exp(-lambda))^2 / A^2 - n_pos / (1 - p)^2,
      n_0 * exp(-lambda) / A^2,
      n_0 * exp(-lambda) / A^2,
      n_0 * exp(-lambda) * p * (1 - p) / A^2 - S_pos / lambda^2
    ),
    nrow = 2, byrow = TRUE
  )
}

newton_raphson <- function(theta, y, B) {
  historico <- matrix(NA_real_, nrow = B + 1, ncol = 2)
  colnames(historico) <- c("p", "lambda")
  historico[1, ] <- theta
  for (i in seq_len(B)) {
    # Dez passos da fórmula, mesmo se p sair de (0, 1).
    # É esse registro, e não um break, que mostra a falha do método.
    U <- score_zip(theta[1], theta[2], y)
    H <- hessiana_zip(theta[1], theta[2], y)
    theta <- theta - solve(H, U)
    historico[i + 1, ] <- theta
  }
  historico
}

B <- 10
chutes_nr <- data.frame(
  inicializacao = paste("Ponto", 1:3),
  p_inicial = p_func(c(-2, -0.5, 1)),
  lambda_inicial = lambda_func(c(0, 1.5, 2))
)

resultados_nr <- data.frame()
historico_nr <- data.frame()
for (perfil in names(perfis)) {
  for (j in seq_len(nrow(chutes_nr))) {
    ponto <- chutes_nr[j, ]
    ajuste <- newton_raphson(
      c(ponto$p_inicial, ponto$lambda_inicial),
      perfis[[perfil]],
      B
    )
    resultados_nr <- rbind(resultados_nr, data.frame(
      perfil = perfil,
      inicializacao = ponto$inicializacao,
      p_inicial = ponto$p_inicial,
      lambda_inicial = ponto$lambda_inicial,
      p_estimado = ajuste[B + 1, 1],
      lambda_estimado = ajuste[B + 1, 2],
      stringsAsFactors = FALSE
    ))
    historico_nr <- rbind(historico_nr, data.frame(
      perfil = perfil,
      inicializacao = ponto$inicializacao,
      iteracao = 0:B,
      p = ajuste[, 1],
      lambda = ajuste[, 2],
      stringsAsFactors = FALSE
    ))
  }
}
print(resultados_nr)
write.csv(resultados_nr, file.path(dir_out, "newton.csv"), row.names = FALSE)

historico_long <- rbind(
  data.frame(
    historico_nr[c("perfil", "inicializacao", "iteracao")],
    parametro = "p", valor = historico_nr$p
  ),
  data.frame(
    historico_nr[c("perfil", "inicializacao", "iteracao")],
    parametro = "lambda", valor = historico_nr$lambda
  )
)
g_nr <- ggplot(
  historico_long,
  aes(x = iteracao, y = valor, color = inicializacao, group = inicializacao)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  facet_grid(parametro ~ perfil, scales = "free_y") +
  scale_x_continuous(breaks = 0:B) +
  labs(
    title = "Trajetória do Newton-Raphson",
    x = "Iteração", y = "Valor do parâmetro", color = "Ponto inicial"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(dir_fig, "newton_raphson_iteracoes.png"), g_nr,
       width = 10, height = 5.5, dpi = 130)
