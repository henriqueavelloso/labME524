############### Parte 1 - Otimização numérica ############### 

# seed

set.seed(248214)

###### pacotes #####

library(tidyverse)
library(ggplot2)

##### leitura dos dados #####

dados <- read_csv("dados/visitas_medicas.csv")

# separando cada perfil

for (perfil in unique(dados$perfil)) {
  assign(perfil, dados[dados$perfil == perfil, ])
}

# análise descritiva das visitas por perfil

dados |>
  group_by(perfil) |>
  summarise(
    n = n(),
    media = mean(visitas),
    mediana = median(visitas),
    maximo = max(visitas),
    proporcao_zeros = mean(visitas == 0),
    .groups = "drop"
  )

##### modelo ZIP #####

# parâmetros

# p - precisamos que 0 < p < 1 ==> vamos usar a curva logistica

p_func <- function(alpha) {
  1 / (1 + exp(-alpha))
}

# lambda - precisamos que lambda > 0 ==> vamos usar exponencial

lambda_func <- function(beta) {
  exp(beta)
}

### log verossimilhança ###

# caso y_i = 0

log_v_zero <- function(p, lambda, y_length){
  log_v <- 0
  for (i in seq_len(y_length)){
    value <- sum(log(p + (1 - p)*exp(-lambda)))
    log_v <- log_v + value
  }
  return(log_v)
}

# caso y_i > 0

log_v_nonzero <- function(y, p, lambda){
  log_v <- 0
  for (y_i in y){
    value <- sum(log(1 - p) - lambda + y_i*log(lambda))
    log_v <- log_v + value
  }
  return(log_v)
}

# calculando log verossimilhança

calculate_logv <- function(y, p, lambda){
  y_zero <- y[y == 0]
  y_length <- length(y_zero)
  
  y_nonzero <- y[y > 0]
  
  value_log_v_zero <- log_v_zero(p, lambda, y_length)
  value_log_v_nonzero <- log_v_nonzero(y_nonzero, p, lambda)
  
  value_log_v <- value_log_v_zero + value_log_v_nonzero
  
  return(value_log_v)
}

# calculando -log verossimilhança

negative_logv <- function(params, y) {
  alpha <- params[1]
  beta <- params[2]
  
  p <- p_func(alpha)
  lambda <- lambda_func(beta)
  
  resultado <- -calculate_logv(y, p, lambda)
  
  return(resultado)
}

# optim()

# pontos iniciais

pontos_iniciais <- data.frame(
  inicializacao = c("Ponto 1", "Ponto 2", "Ponto 3"),
  alpha = c(-2, -0.5, 1),
  beta = c(0, 1.5, 2)
)

# otimizacao para cada perfil

# perfis
perfis <- list(
  baixo = baixo,
  medio = medio,
  alto = alto
)

resultados <- do.call(rbind, lapply(names(perfis), function(perfil) {
  
  dados <- perfis[[perfil]]
  
  do.call(rbind, lapply(seq_len(nrow(pontos_iniciais)), function(i) {
    
    # pontos iniciais
    alpha_0 <- pontos_iniciais$alpha[i]
    beta_0 <- pontos_iniciais$beta[i]
    
    # otimizando
    resultado <- optim(
      par = c(alpha_0, beta_0),
      fn = negative_logv,
      y = dados$visitas
    )
    
    data.frame(
      perfil = perfil,
      inicializacao = pontos_iniciais$inicializacao[i],
      p_inicial = p_func(alpha_0),
      lambda_inicial = lambda_func(beta_0),
      p_estimado = p_func(resultado$par[1]),
      lambda_estimado = lambda_func(resultado$par[2]),
      logv = -resultado$value,
      convergencia = resultado$convergence
    )
  }))
}))

rownames(resultados) <- NULL

### resultado dos estimadores ###

print(resultados)

### visualização -- comparação dos diferentes pontos inicias para cada perfil de risco

ggplot(resultados) +
  
  geom_segment(
    aes(
      x = p_inicial,
      y = lambda_inicial,
      xend = p_estimado,
      yend = lambda_estimado,
      color = inicializacao
    ),
    arrow = arrow(length = grid::unit(0.15, "cm")),
    alpha = 0.7
  ) +
  
  geom_point(
    aes(
      x = p_inicial,
      y = lambda_inicial,
      color = inicializacao
    ),
    size = 3,
    shape = 1,
    stroke = 1.2
  ) +
  
  geom_point(
    aes(
      x = p_estimado,
      y = lambda_estimado,
      color = inicializacao
    ),
    size = 3
  ) +
  
  facet_wrap(~perfil, nrow = 1) +
  
  labs(
    title = "Estimação por Máxima Verossimilhança",
    subtitle = "Comparação entre diferentes pontos iniciais e perfis",
    x = "p — Probabilidade de zero estrutural",
    y = expression(lambda ~ "— Taxa da Poisson"),
    color = "Inicialização"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(shape = 1, fill = NA, size = 3)
    )
  ) +
  
  theme_minimal(base_size = 12)


##### Newton-Raphson #####

# primeiras derivadas em relação a p e lambda

score_zip <- function(p, lambda, y) {
  
  n_0 <- sum(y == 0)
  
  n_pos <- sum(y > 0)
  
  S_pos <- sum(y[y > 0])
  
  score_p <- (n_0*(1 - exp(-lambda)))/(p + (1 - p)*exp(-lambda)) - n_pos/(1 - p)
  
  score_lambda <- (n_0*(-exp(-lambda) + p*exp(-lambda)))/(p + (1 - p)*exp(-lambda)) - n_pos + S_pos/lambda
  
  return(c(score_p, score_lambda))
}

# matriz hessiana

hessiana_zip <- function(p, lambda, y) {
  
  n_0 <- sum(y == 0)
  n_pos <- sum(y > 0)
  S_pos <- sum(y[y > 0])
  
  H_pp <- (-n_0*(1 - exp(-lambda))^2)/(p + (1 - p)*exp(-lambda))^2 - n_pos/(1 - p)^2
  
  H_p_lambda <- (n_0*exp(-lambda))/(p + (1 - p)*exp(-lambda))^2
  
  H_lambda_lambda <- (n_0*exp(-lambda)*p*(1 - p))/(p + (1 - p)*exp(-lambda))^2 - S_pos/lambda^2
  
  H <- matrix(
    c(H_pp, H_p_lambda,
      H_p_lambda, H_lambda_lambda),
    nrow = 2,
    byrow = TRUE
  )
  
  return(H)
}

##### aplicando o newton-raphson #####

pontos_iniciais <- data.frame(
  inicializacao = paste("Ponto", 1:3),
  p_inicial = plogis(c(-2, -0.5, 1)),
  lambda_inicial = exp(c(0, 1.5, 2))
)

# numero de iteracoes
B <- 10

perfis <- list(
  baixo = baixo$visitas,
  medio = medio$visitas,
  alto = alto$visitas
)

newton_raphson <- function(theta, y, B) {
  
  historico <- matrix(
    NA_real_,
    nrow = B + 1,
    ncol = 2
  )
  
  colnames(historico) <- c("p", "lambda")
  historico[1, ] <- theta
  
  for (i in seq_len(B)) {
    # score
    U <- score_zip(theta[1], theta[2], y)
    # hessiana
    H <- hessiana_zip(theta[1], theta[2], y)
    # theta_i+1
    theta <- theta - solve(H, U)
    # atualizand parametros
    historico[i + 1, ] <- theta
  }
  
  return(historico)
}

resultados_nr <- data.frame()
historico_nr <- data.frame()

for (perfil in names(perfis)) {
  for (j in 1:nrow(pontos_iniciais)) {
    
    ponto <- pontos_iniciais[j, ]
    
    theta_inicial <- c(
      ponto$p_inicial,
      ponto$lambda_inicial
    )
    
    # newton-raphson
    ajuste <- newton_raphson(
      theta_inicial, perfis[[perfil]], B
    )
    
    # resultados finais
    resultados_nr <- bind_rows(
      resultados_nr,
      data.frame(
        perfil = perfil,
        inicializacao = ponto$inicializacao,
        p_inicial = ponto$p_inicial,
        lambda_inicial = ponto$lambda_inicial,
        p_estimado = ajuste[B + 1, 1],
        lambda_estimado = ajuste[B + 1, 2]
      )
    )
    
    # histórico das iterações
    historico_nr <- bind_rows(
      historico_nr,
      data.frame(
        perfil = perfil,
        inicializacao = ponto$inicializacao,
        iteracao = 0:B,
        p = ajuste[, 1],
        lambda = ajuste[, 2]
      )
    )
  }
}

print(resultados_nr)

# visualização -- comparação dos diferentes pontos inicias para cada perfil de risco

ggplot(resultados_nr) +
  
  geom_segment(
    aes(
      x = p_inicial,
      y = lambda_inicial,
      xend = p_estimado,
      yend = lambda_estimado,
      color = inicializacao
    ),
    arrow = arrow(length = grid::unit(0.15, "cm")),
    alpha = 0.7
  ) +
  
  geom_point(
    aes(
      x = p_inicial,
      y = lambda_inicial,
      color = inicializacao
    ),
    size = 3,
    shape = 1,
    stroke = 1.2
  ) +
  
  geom_point(
    aes(
      x = p_estimado,
      y = lambda_estimado,
      color = inicializacao
    ),
    size = 3
  ) +
  
  facet_wrap(~perfil, nrow = 1) +
  
  labs(
    title = "Estimação pelo Método de Newton–Raphson",
    subtitle = "Comparação entre diferentes pontos iniciais e perfis",
    x = "p — Probabilidade de zero estrutural",
    y = expression(lambda ~ "— Taxa da Poisson"),
    color = "Inicialização"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(shape = 1, fill = NA, size = 3)
    )
  ) +
  
  theme_minimal(base_size = 12)

### grafico temporal ###

historico_long <- historico_nr %>%
  pivot_longer(
    cols = c(p, lambda),
    names_to = "parametro",
    values_to = "valor"
  )

ggplot(
  historico_long,
  aes(
    x = iteracao,
    y = valor,
    color = inicializacao,
    group = inicializacao
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_grid(parametro ~ perfil, scales = "free_y") +
  scale_x_continuous(breaks = 0:B) +
  labs(
    title = "Trajetória das estimativas pelo método de Newton-Raphson",
    subtitle = "Evolução dos parâmetros ao longo das iterações",
    x = "Iteração",
    y = "Valor do parâmetro",
    color = "Ponto inicial"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15),
    strip.text = element_text(size = 17, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
