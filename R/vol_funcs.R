## Interpolacao da superficie de volatilidade
# Pacotes para analisar:
# akima
# splines
# MBA
# base::spline
# chebpol

# Verifies if needed packages are installed
pkgs_att <- c("akima",
              "splines",
              "MBA",
              "ks",
              "rgl",
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

vol_data <- read_csv("input/IVP_Delta_surface.csv",
                     col_types = cols(date = col_date(format = "%m/%d/%Y")))

vol_surface_raw <- vol_data %>% 
  select(-c(exchange, moneyness, strike)) %>% 
  nest(-c(date, symbol))

# Se houver mais de uma data ou ativo, selecionar uma combinacao

# vol_selected <- vol_surface_raw %>% 
#   filter(date == date_selected, symbol == symbol_selected)

vol_selected <- vol_surface_raw %>% 
  unnest()
x <- unique(vol_selected$delta)
y <- unique(vol_selected$period)
z <- matrix(vol_selected$iv, ncol = length(y), byrow = FALSE)
row.names(z) <- paste("Delta", x)
colnames(z) <- paste("Period", y)

# Interpolacao Bicubica de Akima -------------------------------------------

ddelta <- 1 # variacao do delta na superficie interpolada
dperiod <- 7 # variacao do periodo na superficie interpolada
ny <- ((max(y) - min(y)) / dperiod) + 1 # numero de pontos delta e periodo
nx <- ny
#akima <- bicubic.grid(x, y, z, nx = nx, ny = ny)
# superficie retangular
akima_cub <- bicubic.grid(x, y, z, dx = ddelta, dy = dperiod) 

# Interpolacao Bilinear de Akima ------------------------------------------

akima_bil <- bilinear.grid(x, y, z, dx = ddelta, dy = dperiod) 

# Interpolacao Bicubica com dados irregulares ----------------------------
# Não funcionou muito bem, varios NA

npoints <- 100
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
                       no.X = nx, 
                       no.Y = ny, 
                       sp = FALSE)

# Avaliacao da VI em pontos especificos MBA -----------------------------

mba_eval <- mba.points(xyz_mba, xyeval)$xyz.est %>% 
  as_data_frame()

# Plots -------------------------------------------------------------------

# Plot da Superficie gerada - Akima, MBA

# img_lst <- akima_cub
img_lst <- mba_spline$xyz.est

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

periodo <-  c(44, 51, 58, 65)
ggp_tbl %>% 
  filter(y %in% periodo) %>% 
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

vol_interp <- akima_cub
# vol_interp <- mba_spline$xyz.est

vol_eval <- akima_eval
# vol_eval <- mba_eval

plot_ly(y = vol_interp$x, 
        x = vol_interp$y, 
        z = vol_interp$z) %>% 
  add_surface() %>% 
  layout(scene = list(xaxis = list(title = "Dias"),
                      yaxis = list(title = "Delta"))) %>% 
  add_markers(y = vol_selected$delta,
              x = vol_selected$period,
              z = vol_selected$iv,
              size = 0.1) %>% 
  add_markers(y = vol_eval$x,
              x = vol_eval$y,
              z = vol_eval$z)


  


