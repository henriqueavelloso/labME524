# Parte 2. EM do ZIP, no formato da mistura da aula:
# passo E (peso w_i), passo M (duas médias), guarda a log-verossimilhança.

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

em_zip <- function(y, p0, lambda0, tol = 1e-10, max_iter = 500) {
  p <- p0
  lambda <- lambda0
  ll_hist <- numeric(max_iter)
  p_hist <- numeric(max_iter)
  lambda_hist <- numeric(max_iter)

  for (t in seq_len(max_iter)) {
    # Passo E: probabilidade de zero estrutural dado y_i.
    # Quem consultou não pode ser estrutural.
    w <- numeric(length(y))
    zeros <- y == 0
    w[zeros] <- p / (p + (1 - p) * exp(-lambda))

    # Passo M: p é a fração esperada de estruturais;
    # lambda é a média de consultas só na parte Poisson (peso 1 - w).
    p <- mean(w)
    lambda <- sum((1 - w) * y) / sum(1 - w)

    ll_hist[t] <- loglik_zip(y, p, lambda)
    p_hist[t] <- p
    lambda_hist[t] <- lambda

    if (t > 1 && abs(ll_hist[t] - ll_hist[t - 1]) < tol) {
      return(list(
        p = p, lambda = lambda,
        loglik = ll_hist[t], iter = t,
        ll_hist = ll_hist[seq_len(t)],
        p_hist = p_hist[seq_len(t)],
        lambda_hist = lambda_hist[seq_len(t)]
      ))
    }
  }

  list(
    p = p, lambda = lambda,
    loglik = ll_hist[max_iter], iter = max_iter,
    ll_hist = ll_hist, p_hist = p_hist, lambda_hist = lambda_hist
  )
}

chutes <- data.frame(
  inicializacao = c("Ponto 1", "Ponto 2", "Ponto 3"),
  p0 = plogis(c(-2, -0.5, 1)),
  lambda0 = exp(c(0, 1.5, 2))
)

resultados_em <- data.frame()
historico_em <- data.frame()

for (perfil in names(perfis)) {
  y <- perfis[[perfil]]
  for (j in seq_len(nrow(chutes))) {
    chute <- chutes[j, ]
    fit <- em_zip(y, chute$p0, chute$lambda0)
    resultados_em <- rbind(resultados_em, data.frame(
      perfil = perfil,
      inicializacao = chute$inicializacao,
      p = fit$p,
      lambda = fit$lambda,
      loglik = fit$loglik,
      iter = fit$iter,
      p_iter1 = fit$p_hist[1],
      lambda_iter1 = fit$lambda_hist[1],
      stringsAsFactors = FALSE
    ))
    historico_em <- rbind(historico_em, data.frame(
      perfil = perfil,
      inicializacao = chute$inicializacao,
      iteracao = seq_along(fit$ll_hist),
      loglik = fit$ll_hist,
      stringsAsFactors = FALSE
    ))
  }
}

print(resultados_em)
write.csv(resultados_em, file.path(dir_out, "em.csv"), row.names = FALSE)
write.csv(historico_em, file.path(dir_out, "em_historico.csv"), row.names = FALSE)

# Sobe a cada iteração? (propriedade do EM)
monotonia <- tapply(
  historico_em$loglik,
  list(historico_em$perfil, historico_em$inicializacao),
  function(v) all(diff(v) >= -1e-8)
)
print(monotonia)

g_em <- ggplot(
  historico_em,
  aes(x = iteracao, y = loglik, color = inicializacao)
) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ perfil, scales = "free_y", nrow = 1) +
  labs(
    title = "Log-verossimilhança observada ao longo do EM",
    x = "Iteração",
    y = expression(ell(p, lambda)),
    color = "Ponto inicial"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(dir_fig, "em_loglik.png"), g_em,
       width = 9, height = 3.4, dpi = 130)
