# Parte 3. Gibbs do ZIP com aumento de dados.
# Mesma ordem da aula da mistura: sorteia a latente, depois os parâmetros
# nas condicionais conjugadas. B e burn-in são os do enunciado.

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

gibbs_zip <- function(y, B = 10000, burnin = 2000) {
  n <- length(y)
  zeros <- which(y == 0)
  cadeia_p <- numeric(B)
  cadeia_lambda <- numeric(B)

  # Chute inicial: proporção de zeros acima da Poisson, e lambda que recupera a média.
  media <- mean(y)
  p <- max(mean(y == 0) - exp(-media), 0.05)
  lambda <- media / (1 - p)
  cadeia_p[1] <- p
  cadeia_lambda[1] <- lambda
  z <- rep(0, n)

  for (t in 2:B) {
    # Z_i | resto. Se y_i > 0, Z_i = 0 com probabilidade 1.
    prob_estrutural <- p / (p + (1 - p) * exp(-lambda))
    z[] <- 0
    z[zeros] <- rbinom(length(zeros), size = 1, prob = prob_estrutural)

    # p | z ~ Beta(1 + soma z, 1 + n - soma z)
    s <- sum(z)
    p <- rbeta(1, shape1 = 1 + s, shape2 = 1 + n - s)

    # lambda | y, z ~ Gamma(forma, taxa), só nos clientes do componente Poisson
    poisson <- z == 0
    lambda <- rgamma(
      1,
      shape = 1 + sum(y[poisson]),
      rate = 0.01 + sum(poisson)
    )

    cadeia_p[t] <- p
    cadeia_lambda[t] <- lambda
  }

  idx <- (burnin + 1):B
  p_post <- cadeia_p[idx]
  l_post <- cadeia_lambda[idx]
  list(
    p = cadeia_p,
    lambda = cadeia_lambda,
    resumo = data.frame(
      mediana_p = median(p_post),
      ic_p_inf = unname(quantile(p_post, 0.025)),
      ic_p_sup = unname(quantile(p_post, 0.975)),
      mediana_lambda = median(l_post),
      ic_lambda_inf = unname(quantile(l_post, 0.025)),
      ic_lambda_sup = unname(quantile(l_post, 0.975)),
      mediana_1a_metade_p = median(p_post[seq_len(length(p_post) / 2)]),
      mediana_2a_metade_p = median(p_post[(length(p_post) / 2 + 1):length(p_post)])
    )
  )
}

B <- 10000
burnin <- 2000
resumo_gibbs <- data.frame()
cadeia_plot <- data.frame()

for (perfil in names(perfis)) {
  fit <- gibbs_zip(perfis[[perfil]], B = B, burnin = burnin)
  resumo_gibbs <- rbind(resumo_gibbs, data.frame(
    perfil = perfil,
    fit$resumo,
    stringsAsFactors = FALSE
  ))
  # Uma a cada dez iterações, só para o gráfico não ficar pesado.
  idx_plot <- seq(1, B, by = 10)
  cadeia_plot <- rbind(cadeia_plot, data.frame(
    perfil = perfil,
    iteracao = idx_plot,
    p = fit$p[idx_plot],
    lambda = fit$lambda[idx_plot],
    stringsAsFactors = FALSE
  ))
}

print(resumo_gibbs)
write.csv(resumo_gibbs, file.path(dir_out, "gibbs.csv"), row.names = FALSE)

cadeia_long <- rbind(
  data.frame(
    cadeia_plot[c("perfil", "iteracao")],
    parametro = "p", valor = cadeia_plot$p
  ),
  data.frame(
    cadeia_plot[c("perfil", "iteracao")],
    parametro = "lambda", valor = cadeia_plot$lambda
  )
)

g_gibbs <- ggplot(cadeia_long, aes(x = iteracao, y = valor)) +
  geom_line(linewidth = 0.3, color = "grey25") +
  geom_vline(xintercept = burnin, linetype = 2, color = "firebrick") +
  facet_grid(parametro ~ perfil, scales = "free_y") +
  labs(
    title = "Trajetória do Gibbs sampler",
    subtitle = "A linha tracejada marca o fim do burn-in (2.000 iterações)",
    x = "Iteração",
    y = "Valor amostrado"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(dir_fig, "gibbs_trajetorias.png"), g_gibbs,
       width = 9, height = 4.6, dpi = 130)
