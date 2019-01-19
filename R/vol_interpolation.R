## Interpolacao da superficie de volatilidade
# Pacotes para analisar:
# akima
# splines
# MBA
# stats::spline
# chebpol
# NMOF para alguns metodos aplicados a financas

# Verifies if needed packages are installed
pkgs_att <- c("akima",
              "splines",
              "MBA",
              "ks",
              "tidyverse",
              "ggthemes",
              "plotly")
new.pkgs <- pkgs_att[!(pkgs_att %in% installed.packages()[,"Package"])]

# Installs new packages
if (length(new.pkgs)) install.packages(new.pkgs)

# Loads and attaches all needed packages
for (s in pkgs_att) {
  if (!library(s, character.only = TRUE, logical.return = TRUE)) 
    stop("Error loading packages.")
}

rates <- read_csv("input/Interest_rates.csv",
                  col_types = cols(date = col_date(format = "%m/%d/%Y"))) %>% 
  mutate(log_rate = log(1 + rate / 100)) %>% 
  select(date, period, log_rate)

vol_data <- read_csv("input/IVP_Delta_surface.csv",
                     col_types = cols(date = col_date(format = "%m/%d/%Y"))) %>% 
  left_join(rates, by = c("date", "period"))

symbols_list <- unique(vol_data$symbol)
dates_list <- unique(vol_data$date)

symbol_sel <- "IWM"
data_sel <- as.Date("2017-09-21", format = "%Y-%m-%d")

vol_selected <- vol_data %>% 
  filter(date == data_sel, symbol == symbol_sel) %>% 
  select(-c(date, symbol, exchange, moneyness, strike))

delta <- unique(vol_selected$delta)
periodo <- unique(vol_selected$period)
impvol <- matrix(vol_selected$iv, ncol = length(periodo), byrow = FALSE)
row.names(impvol) <- paste("Delta", delta)
colnames(impvol) <- paste("Period", periodo)

# Numero de pontos de interpolacao em ambos eixos
npoints <- max(periodo) - min(periodo) + 1

# Remove variaveis nao necessarias ----------------------------------------

rm(list = c("rates", "new.pkgs", "pkgs_att", "s"))

# Interpolacao Bicubica de Akima -------------------------------------------

dperiod <- 1 # variacao do periodo na superficie interpolada
ddelta <- (max(delta) - min(delta)) / (npoints - 1)  # variacao do delta na superficie interpolada

# superficie retangular
# akima_cub <- bicubic.grid(x, y, z, nx = npoints, ny = npoints)
akima_cub <- bicubic.grid(delta, periodo, impvol, dx = ddelta, dy = dperiod) 

# Interpolacao Bilinear de Akima ------------------------------------------

akima_bil <- bilinear.grid(x, y, z, dx = ddelta, dy = dperiod) 

# Interpolacao Bicubica com dados irregulares ----------------------------
# Não funcionou muito bem, varios NA

akima_int <- interp(x = vol_selected$delta,
                    y = vol_selected$period,
                    z = vol_selected$iv,
                    linear = FALSE,
                    nx = npoints,
                    ny = npoints)

# Avaliacao da VI em pontos especificos Akima -----------------------------

# Quero avaliar varias datas para um mesmo delta = 0.45
xeval <- 45
yeval <- seq(30, 1080, by = 44)
xyeval <- expand.grid(xeval, yeval)
akima_eval <- bicubic(x, y, z, x0 = xyeval$Var1, y0 = xyeval$Var2) %>% 
  as_data_frame()

# B-Splines com pacote MBA ------------------------------------------------

xyz_mba <- vol_selected %>% 
  select(delta, period, iv) %>% 
  rename(x = delta,
         y = period,
         z = iv)
mba_spline <- mba.surf(xyz_mba, 
                       no.X = npoints, 
                       no.Y = npoints, 
                       sp = FALSE)

# Avaliacao da VI em pontos especificos MBA -----------------------------

mba_eval <- mba.points(xyz_mba, xyeval)$xyz.est %>% 
  as_data_frame()

# Plots -------------------------------------------------------------------

# Plot da Superficie gerada - Akima, MBA

img_lst <- akima_cub
# img_lst <- mba_spline$xyz.est

zmin <- min(img_lst$z, na.rm = TRUE)
zmax <- max(img_lst$z, na.rm = TRUE)
breaks <- pretty(c(zmin,zmax),10)
colors <- heat.colors(length(breaks) - 1)
image(img_lst, breaks = breaks, col = colors)
contour(img_lst, levels = breaks,  add = TRUE)


# ggplot - Akima, MBA
# Escolher um modelo para o ggp_tbl

ggp_tbl <- as_tibble(interp2xyz(akima_cub))
# ggp_tbl <- as_tibble(interp2xyz(mba_spline$xyz.est))

ggplot(ggp_tbl, aes(x = x, y = y)) +
  geom_raster(aes(fill = z), interpolate = FALSE) +
  geom_contour(aes(z = z), color = "black", size = 0.1) +
  scale_fill_viridis_c() +
  scale_x_reverse() +
  guides(fill = guide_legend(title = "Vol")) +
  labs(title = "Contorno superfície B-spline",
       x = "Delta",
       y = "Dias") +
  theme_classic()

# ggplot Smile - Akima, MBA

periodos <-  c(44, 51, 58, 65)
ggp_tbl %>% 
  filter(y %in% periodos) %>% 
  mutate(y = as.factor(y)) %>% 
  ggplot(aes(x = x, y = z, group = y)) +
  geom_line(aes(color = y)) +
  scale_color_viridis_d() +
  scale_x_reverse() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = paste("Smiles de volatilidade"),
       x = "Delta",
       y = "Volatilidade Implícita") +
  guides(color = guide_legend(title = "Periodo")) +
  theme_bw()


# Plotly 3D Surface ---------------------------------------------------
# Reversão do eixo em graficos 3D nao eh suportado
# No plotly tem que botar delta no y e periodo no x !!
# Reduz o numero de pontos do grid para ficar mais leve a superficie

idx <- seq(1, npoints, by = 2)

vol_interp <- list(x = akima_cub$x[idx],
                   y = akima_cub$y[idx],
                   z = akima_cub$z[idx, idx])

vol_interp <- list(x = mba_spline$xyz.est$x[idx],
                   y = mba_spline$xyz.est$y[idx],
                   z = mba_spline$xyz.est$z[idx, idx])

# vol_interp <- mba_spline$xyz.est

vol_eval <- akima_eval
# vol_eval <- mba_eval

plot_ly(y = vol_interp$x, 
        x = vol_interp$y, 
        z = vol_interp$z) %>% 
  add_surface() %>% 
  layout(scene = list(xaxis = list(title = "Dias"),
                      yaxis = list(title = "Delta"),
                      camera = list(eye = list(x = -1.55,
                                               y = -1.55,
                                               z = 0.25)))) %>% 
  add_markers(y = vol_selected$delta,
              x = vol_selected$period,
              z = vol_selected$iv,
              size = 0.1) %>% 
  add_markers(y = vol_eval$x,
              x = vol_eval$y,
              z = vol_eval$z)

# Inverter o Delta --------------------------------------------------------

prices <- vol_selected %>% 
  rename(S = stock_price_for_iv) %>% 
  mutate(delta = delta / 100,
         tau = period / 360,
         omega2 = iv^2 * tau,
         log_mF = sqrt(omega2) * qnorm(delta) - (omega2 / 2),
         log_m = log_mF - log_rate * tau,
         K0 = S * exp(-log_mF),
         K = S * exp(-log_m),
         d1 = (log_mF + omega2 / 2) / sqrt(omega2),
         d2 = d1 - sqrt(omega2),
         C0 = S * delta - K0 * pnorm(d2))  

# Criou-se um grid irregular nos strikes que nao eh facilmente
# interpolavel na superficie diretamente.
# Interpolar os smiles, fixando os periodos.

# Strikes usados para avaliar a interpolacao
# Deve ser uma lista de N elementos onde N eh o numero de periodos diferentes
strikes_eval <- list(seq(min(prices$K), 
                         max(prices$K), 
                         length.out = npoints)) %>% 
  rep(length(periodo))

# interp_dist guarda as interpolacoes de preco da call e a distribuicao
# implicita desta

interp_dist <- prices %>% 
  nest(-c(period, log_rate)) %>% 
  mutate(strikes_eval = strikes_eval,
         bspline = map(data, ~interpSpline(.x$K, .x$C0,
                                           bSpline = TRUE)),
         bspline_eval = map2(bspline, strikes_eval, ~predict(.x, .y)$y),
         smooth = map(data, ~smooth.spline(.x$K, .x$C0, spar = 0.5)),
         smooth_eval = map2(smooth, strikes_eval, ~predict(.x, .y)$y),
         smooth_2deriv = map2(smooth, strikes_eval, ~pmax(predict(.x, .y, 
                                                            deriv = 2)$y,
                                                    0))) %>% 
  select(-c(data, bspline, smooth)) %>% 
  unnest() %>% 
  mutate(impdist = exp(log_rate * period / 360) * smooth_2deriv) %>% 
  select(-smooth_2deriv) %>% 
  gather(key = modelo, value = valor, c(bspline_eval, 
                                        smooth_eval,
                                        impdist))


# Plot dos Smiles e das distribuicoes implicitas
# ATENCAO: metodos spline extrapolam de forma linear
# primeira derivada do ultimo ponto 
interp_dist %>% 
  filter(period %in% periodos,
         modelo %in% c("bspline_eval", "smooth_eval")) %>% 
  mutate(period = as.factor(period)) %>% 
  ggplot(aes(x = strikes_eval, y = valor, group = period)) +
  geom_line(aes(color = period)) +
  facet_wrap(~modelo) +
  scale_color_viridis_d() +
  scale_y_continuous(labels = scales::dollar) +
  labs(x = "Strike",
       y = "Preço Call") +
  guides(color = guide_legend(title = "Periodo")) +
  theme_bw()

interp_dist %>% 
  filter(period %in% periodos,
         modelo %in% c("impdist")) %>% 
  mutate(period = as.factor(period)) %>% 
  ggplot(aes(x = strikes_eval, y = valor, group = period)) +
  geom_line(aes(color = period)) +
  scale_color_viridis_d() +
  labs(x = "Strike",
       y = "p(x)") +
  guides(color = guide_legend(title = "Periodo")) +
  theme_bw()

# Checa se a integral da distribuicao soma 1

distribuicao <- interp_dist %>% 
  group_by(period) %>% 
  filter(modelo == "impdist") %>% 
  summarise(CDF = sum(valor * (strikes - lag(strikes)),
                      na.rm = TRUE))
# Nao esta somando 1, longe!!
# Variavel strikes em prices eh dependente do periodo e esta limitada
# a [min, max] K de cada periodo.
# strikes NAO eh necessaria para implementar as interpolacoes, sendo usada
# apenas para a avaliacao destas
# CRIAR uma nova strikes que seja a mesma para todos os periodos e cubra
# valores de strike indo desde muito baixos ate muito altos.
# O maior periodo contem a maior amplitude de strikes

