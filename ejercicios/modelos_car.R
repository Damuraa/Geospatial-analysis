# ============================================================
# MODELOS CAR PARA FRA - CUENCAS
# ============================================================

# ------------------------------------------------------------
# 1. Paquetes
# ------------------------------------------------------------

paquetes <- c(
  "sf",
  "spdep",
  "CARBayes",
  "coda",
  "dplyr"
)

instalar <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]

if (length(instalar) > 0) {
  install.packages(instalar, dependencies = TRUE)
}

library(sf)
library(spdep)
library(CARBayes)
library(coda)
library(dplyr)

# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

ruta <- "carto/cuencas_sar_R.gpkg"

gdf <- st_read(
  ruta,
  layer = "cuencas_sar_R"
)

print(st_crs(gdf))
print(names(gdf))
print(nrow(gdf))

# ------------------------------------------------------------
# 3. Preparar variables
# ------------------------------------------------------------

gdf <- gdf %>%
  mutate(
    FRA = as.numeric(FRA),
    area_km2 = as.numeric(area_km2),
    Beta_media = as.numeric(Beta_media),
    Beta_std = as.numeric(Beta_std)
  ) %>%
  filter(
    !is.na(FRA),
    !is.na(area_km2),
    !is.na(Beta_media),
    !is.na(Beta_std),
    FRA >= 0,
    FRA <= 1
  )

# Escalar predictoras
gdf <- gdf %>%
  mutate(
    area_km2_z = as.numeric(scale(area_km2)),
    Beta_media_z = as.numeric(scale(Beta_media)),
    Beta_std_z = as.numeric(scale(Beta_std))
  )

# ID espacial
gdf$ID <- 1:nrow(gdf)

# DataFrame sin geometría para CARBayes
datos_modelo <- st_drop_geometry(gdf)

formula_base <- FRA ~ area_km2_z + Beta_media_z + Beta_std_z

# ------------------------------------------------------------
# 4. Matriz de vecindad Queen
# ------------------------------------------------------------

queen_nb <- poly2nb(
  gdf,
  queen = TRUE
)

queen_lw <- nb2listw(
  queen_nb,
  style = "W",
  zero.policy = TRUE
)

W_mat <- nb2mat(
  queen_nb,
  style = "B",
  zero.policy = TRUE
)

print(summary(queen_nb))

cat("\nDimensión matriz W:\n")
print(dim(W_mat))

cat("\nNúmero de vecinos por cuenca:\n")
print(table(card(queen_nb)))

# ------------------------------------------------------------
# 5. Parámetros MCMC
# ------------------------------------------------------------

burnin <- 20000
n_sample <- 80000
thin <- 10

# En Windows es más estable usar 1 core
n_chains <- 3
n_cores <- 1

# ------------------------------------------------------------
# 6. Modelo base bayesiano sin CAR
# ------------------------------------------------------------

cat("\n================ MODELO BASE S.glm ================\n")

m_glm_car <- S.glm(
  formula = formula_base,
  data = datos_modelo,
  family = "gaussian",
  burnin = burnin,
  n.sample = n_sample,
  thin = thin,
  n.chains = n_chains,
  n.cores = n_cores
)

print(m_glm_car)
print(m_glm_car$modelfit)

moran_glm_car <- moran.mc(
  residuals(m_glm_car),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos modelo base CARBayes:\n")
print(moran_glm_car)

# ------------------------------------------------------------
# 7. Modelo CAR Leroux
# ------------------------------------------------------------

cat("\n================ CAR LEROUX ================\n")

m_leroux <- S.CARleroux(
  formula = formula_base,
  data = datos_modelo,
  family = "gaussian",
  W = W_mat,
  burnin = burnin,
  n.sample = n_sample,
  thin = thin,
  n.chains = n_chains,
  n.cores = n_cores
)

print(m_leroux)
print(m_leroux$modelfit)
print(summary(m_leroux$samples))

moran_leroux <- moran.mc(
  residuals(m_leroux),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos Leroux:\n")
print(moran_leroux)

# ------------------------------------------------------------
# 8. Modelo ICAR aproximado: Leroux con rho = 1
# ------------------------------------------------------------

cat("\n================ ICAR rho = 1 ================\n")

m_icar <- S.CARleroux(
  formula = formula_base,
  data = datos_modelo,
  family = "gaussian",
  W = W_mat,
  rho = 1,
  burnin = burnin,
  n.sample = n_sample,
  thin = thin,
  n.chains = n_chains,
  n.cores = n_cores
)

print(m_icar)
print(m_icar$modelfit)
print(summary(m_icar$samples))

moran_icar <- moran.mc(
  residuals(m_icar),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos ICAR:\n")
print(moran_icar)

# ============================================================
# 9. COMPARACIÓN FINAL SIN BYM
# ============================================================

# ------------------------------------------------------------
# Función robusta para extraer coeficientes beta
# ------------------------------------------------------------

extraer_beta <- function(modelo, nombre_modelo) {
  
  resumen_beta <- summary(
    modelo$samples$beta,
    quantiles = c(0.025, 0.975)
  )
  
  beta_mean <- resumen_beta$statistics[, "Mean"]
  beta_sd <- resumen_beta$statistics[, "SD"]
  beta_ci <- resumen_beta$quantiles
  
  nombres <- rownames(resumen_beta$statistics)
  
  if (is.null(nombres) || length(nombres) == 0) {
    nombres <- c(
      "(Intercept)",
      "area_km2_z",
      "Beta_media_z",
      "Beta_std_z"
    )[seq_along(beta_mean)]
  }
  
  tabla <- data.frame(
    Modelo = rep(nombre_modelo, length(beta_mean)),
    Variable = nombres,
    Media = as.numeric(beta_mean),
    SD = as.numeric(beta_sd),
    IC_2.5 = as.numeric(beta_ci[, 1]),
    IC_97.5 = as.numeric(beta_ci[, 2])
  )
  
  tabla$Significativo_95 <- !(
    tabla$IC_2.5 <= 0 &
      tabla$IC_97.5 >= 0
  )
  
  return(tabla)
}

# ------------------------------------------------------------
# Extraer coeficientes beta
# ------------------------------------------------------------

beta_glm_car <- extraer_beta(m_glm_car, "Base S.glm")
beta_leroux <- extraer_beta(m_leroux, "Leroux")
beta_icar <- extraer_beta(m_icar, "ICAR")

tabla_betas_car <- bind_rows(
  beta_glm_car,
  beta_leroux,
  beta_icar
)

cat("\n================ COEFICIENTES BETA ================\n")
print(tabla_betas_car)

# ------------------------------------------------------------
# Función robusta para extraer rho
# ------------------------------------------------------------

extraer_rho <- function(modelo) {
  
  if (is.null(modelo$samples$rho)) {
    return(
      list(
        media = NA,
        ic_inf = NA,
        ic_sup = NA,
        resumen = NA
      )
    )
  }
  
  rho_summary <- summary(
    modelo$samples$rho,
    quantiles = c(0.025, 0.975)
  )
  
  if (is.null(dim(rho_summary$statistics))) {
    rho_media <- as.numeric(rho_summary$statistics["Mean"])
  } else {
    rho_media <- as.numeric(rho_summary$statistics[, "Mean"][1])
  }
  
  if (is.null(dim(rho_summary$quantiles))) {
    rho_inf <- as.numeric(rho_summary$quantiles["2.5%"])
    rho_sup <- as.numeric(rho_summary$quantiles["97.5%"])
  } else {
    rho_inf <- as.numeric(rho_summary$quantiles[, "2.5%"][1])
    rho_sup <- as.numeric(rho_summary$quantiles[, "97.5%"][1])
  }
  
  return(
    list(
      media = rho_media,
      ic_inf = rho_inf,
      ic_sup = rho_sup,
      resumen = rho_summary
    )
  )
}

# ------------------------------------------------------------
# Extraer rho de Leroux
# ------------------------------------------------------------

rho_leroux_obj <- extraer_rho(m_leroux)

rho_leroux <- rho_leroux_obj$media
rho_leroux_inf <- rho_leroux_obj$ic_inf
rho_leroux_sup <- rho_leroux_obj$ic_sup

cat("\n================ RHO LEROUX ================\n")
print(rho_leroux_obj$resumen)

cat("\nRho medio Leroux:", rho_leroux, "\n")
cat("IC 95% rho:", rho_leroux_inf, "-", rho_leroux_sup, "\n")

# ------------------------------------------------------------
# Tabla comparativa CAR
# ------------------------------------------------------------

tabla_car <- data.frame(
  Modelo = c(
    "Base S.glm",
    "Leroux",
    "ICAR"
  ),
  DIC = c(
    as.numeric(m_glm_car$modelfit["DIC"]),
    as.numeric(m_leroux$modelfit["DIC"]),
    as.numeric(m_icar$modelfit["DIC"])
  ),
  WAIC = c(
    as.numeric(m_glm_car$modelfit["WAIC"]),
    as.numeric(m_leroux$modelfit["WAIC"]),
    as.numeric(m_icar$modelfit["WAIC"])
  ),
  p_d = c(
    as.numeric(m_glm_car$modelfit["p.d"]),
    as.numeric(m_leroux$modelfit["p.d"]),
    as.numeric(m_icar$modelfit["p.d"])
  ),
  LMPL = c(
    as.numeric(m_glm_car$modelfit["LMPL"]),
    as.numeric(m_leroux$modelfit["LMPL"]),
    as.numeric(m_icar$modelfit["LMPL"])
  ),
  loglikelihood = c(
    as.numeric(m_glm_car$modelfit["loglikelihood"]),
    as.numeric(m_leroux$modelfit["loglikelihood"]),
    as.numeric(m_icar$modelfit["loglikelihood"])
  ),
  Moran_I_residuos = c(
    as.numeric(moran_glm_car$statistic),
    as.numeric(moran_leroux$statistic),
    as.numeric(moran_icar$statistic)
  ),
  p_Moran_residuos = c(
    as.numeric(moran_glm_car$p.value),
    as.numeric(moran_leroux$p.value),
    as.numeric(moran_icar$p.value)
  ),
  rho_media = c(
    NA,
    rho_leroux,
    1
  ),
  rho_IC_2.5 = c(
    NA,
    rho_leroux_inf,
    1
  ),
  rho_IC_97.5 = c(
    NA,
    rho_leroux_sup,
    1
  )
)

cat("\n================ TABLA COMPARATIVA CAR ================\n")
print(tabla_car)

# ------------------------------------------------------------
# Guardar tablas
# ------------------------------------------------------------

write.csv(
  tabla_car,
  "resultados_comparacion_modelos_car.csv",
  row.names = FALSE
)

write.csv(
  tabla_betas_car,
  "resultados_betas_modelos_car.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Guardar predicciones y residuos
# ------------------------------------------------------------

gdf$pred_CAR_base <- m_glm_car$fitted.values
gdf$res_CAR_base <- residuals(m_glm_car)

gdf$pred_Leroux <- m_leroux$fitted.values
gdf$res_Leroux <- residuals(m_leroux)

gdf$pred_ICAR <- m_icar$fitted.values
gdf$res_ICAR <- residuals(m_icar)

st_write(
  gdf,
  "cuencas_modelos_car_resultados.gpkg",
  layer = "cuencas_modelos_car_resultados",
  delete_layer = TRUE
)

# ------------------------------------------------------------
# Gráficos de convergencia
# ------------------------------------------------------------

png(
  "traza_rho_leroux.png",
  width = 1200,
  height = 700
)

plot(
  m_leroux$samples$rho,
  main = "Traza MCMC - rho Leroux",
  ylab = "rho",
  xlab = "Iteración"
)

dev.off()

png(
  "trazas_beta_leroux.png",
  width = 1200,
  height = 900
)

plot(
  m_leroux$samples$beta,
  main = "Trazas MCMC - coeficientes beta Leroux"
)

dev.off()

# ------------------------------------------------------------
# Mensaje final
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ANÁLISIS CAR FINALIZADO CORRECTAMENTE\n")
cat("Archivos generados:\n")
cat("- resultados_comparacion_modelos_car.csv\n")
cat("- resultados_betas_modelos_car.csv\n")
cat("- cuencas_modelos_car_resultados.gpkg\n")
cat("- traza_rho_leroux.png\n")
cat("- trazas_beta_leroux.png\n")
cat("============================================================\n")
# ============================================================
# PSEUDO R2 Y MÉTRICAS DE DESEMPEÑO PARA MODELOS CAR
# ============================================================

# Variable observada
y_obs <- datos_modelo$FRA

# ------------------------------------------------------------
# Función para calcular métricas
# ------------------------------------------------------------

calcular_metricas_car <- function(modelo, nombre_modelo, y_obs) {
  
  # Predicción posterior media del modelo
  y_pred <- as.numeric(modelo$fitted.values)
  
  # Residuos
  residuos <- y_obs - y_pred
  
  # Suma de cuadrados residual
  sse <- sum(residuos^2, na.rm = TRUE)
  
  # Suma de cuadrados total
  sst <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
  
  # Pseudo R2
  pseudo_r2 <- 1 - (sse / sst)
  
  # R2 por correlación observado-predicho
  r2_cor <- cor(y_obs, y_pred, use = "complete.obs")^2
  
  # Error medio absoluto
  mae <- mean(abs(residuos), na.rm = TRUE)
  
  # Raíz del error cuadrático medio
  rmse <- sqrt(mean(residuos^2, na.rm = TRUE))
  
  data.frame(
    Modelo = nombre_modelo,
    Pseudo_R2 = pseudo_r2,
    R2_cor = r2_cor,
    MAE = mae,
    RMSE = rmse
  )
}

# ------------------------------------------------------------
# Calcular métricas para cada modelo CAR
# ------------------------------------------------------------

metricas_car <- bind_rows(
  calcular_metricas_car(m_glm_car, "Base S.glm", y_obs),
  calcular_metricas_car(m_leroux, "Leroux", y_obs),
  calcular_metricas_car(m_icar, "ICAR", y_obs)
)

cat("\n================ MÉTRICAS CAR ================\n")
print(metricas_car)

# Guardar tabla
write.csv(
  metricas_car,
  "metricas_modelos_car.csv",
  row.names = FALSE
)

# ============================================================
# EXPORTAR CAPA FINAL CON RESULTADOS CAR LEROUX
# ============================================================

# Verificar que existan los objetos necesarios
if (!exists("gdf")) {
  stop("No existe el objeto gdf. Debes correr primero la carga de datos.")
}

if (!exists("m_leroux")) {
  stop("No existe el modelo m_leroux. Debes correr primero el modelo CAR Leroux.")
}

# ------------------------------------------------------------
# 1. Agregar predicción y residuos del modelo Leroux
# ------------------------------------------------------------

gdf$FRA_pred_CAR_Leroux <- as.numeric(m_leroux$fitted.values)

gdf$res_CAR_Leroux <- gdf$FRA - gdf$FRA_pred_CAR_Leroux

gdf$res_abs_CAR_Leroux <- abs(gdf$res_CAR_Leroux)

# ------------------------------------------------------------
# 2. Asegurar variables estandarizadas
# ------------------------------------------------------------

if (!"area_km2_z" %in% names(gdf)) {
  gdf$area_km2_z <- as.numeric(scale(gdf$area_km2))
}

if (!"Beta_media_z" %in% names(gdf)) {
  gdf$Beta_media_z <- as.numeric(scale(gdf$Beta_media))
}

if (!"Beta_std_z" %in% names(gdf)) {
  gdf$Beta_std_z <- as.numeric(scale(gdf$Beta_std))
}

# ------------------------------------------------------------
# 3. Crear capa limpia
# ------------------------------------------------------------

# OJO:
# No se selecciona "geometry" porque en tu objeto la geometría
# se llama "geom" y sf la conserva automáticamente.

gdf_car_leroux <- gdf %>%
  dplyr::select(
    dplyr::any_of(c(
      "HYBAS_ID",
      "FRA",
      "FRA_pred_CAR_Leroux",
      "res_CAR_Leroux",
      "res_abs_CAR_Leroux",
      "area_km2",
      "area_km2_z",
      "Beta_media",
      "Beta_media_z",
      "Beta_std",
      "Beta_std_z"
    ))
  )

# Verificar que sigue siendo sf
print(class(gdf_car_leroux))
print(st_geometry_type(gdf_car_leroux))
print(st_crs(gdf_car_leroux))

# ------------------------------------------------------------
# 4. Guardar GeoPackage
# ------------------------------------------------------------

st_write(
  gdf_car_leroux,
  "carto/cuencas_CAR_Leroux.gpkg",
  layer = "cuencas_CAR_Leroux",
  delete_layer = TRUE
)

cat("\n============================================================\n")
cat("Capa CAR Leroux exportada correctamente:\n")
cat("carto/cuencas_CAR_Leroux.gpkg\n")
cat("Layer: cuencas_CAR_Leroux\n")
cat("============================================================\n")