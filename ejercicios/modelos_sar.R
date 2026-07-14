# ============================================================
# MODELOS SAR PARA FRA - CUENCAS
# ============================================================

# ============================================================
# 1. Configurar librería personal e instalar paquetes
# ============================================================

# Ver librerías disponibles
print(.libPaths())

# Paquetes necesarios
paquetes <- c(
  "jsonlite",
  "rlang",
  "languageserver",
  "sf",
  "spdep",
  "spatialreg",
  "dplyr"
)

# Instalar solo los que falten
instalar <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]

if (length(instalar) > 0) {
  install.packages(
    instalar,
    dependencies = TRUE
  )
}

# Cargar paquetes
library(sf)
library(spdep)
library(spatialreg)
library(dplyr)

# ------------------------------------------------------------
# 2. Cargar datos
# ------------------------------------------------------------

ruta <- "carto/cuencas_sar_R.gpkg"

gdf <- st_read(
  ruta,
  layer = "cuencas_sar_R"
)

# Revisar
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

# Escalar variables predictoras
gdf <- gdf %>%
  mutate(
    area_km2_z = as.numeric(scale(area_km2)),
    Beta_media_z = as.numeric(scale(Beta_media)),
    Beta_std_z = as.numeric(scale(Beta_std))
  )

# Fórmula base
formula_base <- FRA ~ area_km2_z + Beta_media_z + Beta_std_z

# ------------------------------------------------------------
# 4. Matriz Queen
# ------------------------------------------------------------

queen_nb <- poly2nb(gdf, queen = TRUE)

queen_lw <- nb2listw(
  queen_nb,
  style = "W",
  zero.policy = TRUE
)

print(summary(queen_nb))

# ------------------------------------------------------------
# 5. OLS base
# ------------------------------------------------------------

m_ols <- lm(
  formula_base,
  data = gdf
)

cat("\n================ OLS ================\n")
print(summary(m_ols))
print(AIC(m_ols))

moran_ols <- moran.mc(
  residuals(m_ols),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos OLS:\n")
print(moran_ols)

# ------------------------------------------------------------
# 6. SLX - Spatially Lagged X
# ------------------------------------------------------------

m_slx <- lmSLX(
  formula_base,
  data = gdf,
  listw = queen_lw,
  zero.policy = TRUE
)

cat("\n================ SLX ================\n")
print(summary(m_slx))
print(AIC(m_slx))

moran_slx <- moran.mc(
  residuals(m_slx),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos SLX:\n")
print(moran_slx)


# ============================================================
# SLX REDUCIDO: solo rezago espacial del área
# ============================================================

gdf$lag_area_km2_z <- lag.listw(
  queen_lw,
  gdf$area_km2_z,
  zero.policy = TRUE
)

m_slx_red <- lm(
  FRA ~ area_km2_z + Beta_media_z + Beta_std_z + lag_area_km2_z,
  data = gdf
)

cat("\n================ SLX REDUCIDO ================\n")
print(summary(m_slx_red))
print(AIC(m_slx_red))

moran_slx_red <- moran.mc(
  residuals(m_slx_red),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos SLX reducido:\n")
print(moran_slx_red)
# ------------------------------------------------------------
# 7. SEM - Spatial Error Model
# ------------------------------------------------------------

m_sem <- errorsarlm(
  formula_base,
  data = gdf,
  listw = queen_lw,
  zero.policy = TRUE
)

cat("\n================ SEM ================\n")
print(summary(m_sem))
print(AIC(m_sem))

moran_sem <- moran.mc(
  residuals(m_sem),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos SEM:\n")
print(moran_sem)

# ------------------------------------------------------------
# 8. SAR-Lag - Spatial Lag Model
# ------------------------------------------------------------

m_sar <- lagsarlm(
  formula_base,
  data = gdf,
  listw = queen_lw,
  type = "lag",
  zero.policy = TRUE
)

cat("\n================ SAR-LAG ================\n")
print(summary(m_sar))
print(AIC(m_sar))

moran_sar <- moran.mc(
  residuals(m_sar),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos SAR:\n")
print(moran_sar)

# Impactos SAR
cat("\nImpactos SAR:\n")
imp_sar <- impacts(
  m_sar,
  listw = queen_lw,
  R = 200,
  zero.policy = TRUE
)

print(summary(imp_sar, zstats = TRUE, short = TRUE))

# ------------------------------------------------------------
# 9. SDM - Spatial Durbin Model
# ------------------------------------------------------------

m_sdm <- lagsarlm(
  formula_base,
  data = gdf,
  listw = queen_lw,
  type = "mixed",
  zero.policy = TRUE
)

cat("\n================ SDM ================\n")
print(summary(m_sdm))
print(AIC(m_sdm))

moran_sdm <- moran.mc(
  residuals(m_sdm),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos SDM:\n")
print(moran_sdm)

# Impactos SDM
cat("\nImpactos SDM:\n")
imp_sdm <- impacts(
  m_sdm,
  listw = queen_lw,
  R = 200,
  zero.policy = TRUE
)

print(summary(imp_sdm, zstats = TRUE, short = TRUE))

# ------------------------------------------------------------
# 10. SDEM - Spatial Durbin Error Model
# ------------------------------------------------------------

m_sdem <- errorsarlm(
  formula_base,
  data = gdf,
  listw = queen_lw,
  etype = "emixed",
  zero.policy = TRUE
)

cat("\n================ SDEM ================\n")
print(summary(m_sdem))
print(AIC(m_sdem))

moran_sdem <- moran.mc(
  residuals(m_sdem),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos SDEM:\n")
print(moran_sdem)

# ------------------------------------------------------------
# 11. SAC / SARAR - Spatial Lag + Spatial Error
# ------------------------------------------------------------

m_sac <- sacsarlm(
  formula_base,
  data = gdf,
  listw = queen_lw,
  type = "sac",
  zero.policy = TRUE
)

cat("\n================ SAC / SARAR ================\n")
print(summary(m_sac))
print(AIC(m_sac))

moran_sac <- moran.mc(
  residuals(m_sac),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos SAC:\n")
print(moran_sac)

# ------------------------------------------------------------
# 12. Manski / GNS - SAR + WX + Error espacial
# ------------------------------------------------------------

m_gns <- sacsarlm(
  formula_base,
  data = gdf,
  listw = queen_lw,
  type = "sacmixed",
  zero.policy = TRUE
)

cat("\n================ Manski / GNS ================\n")
print(summary(m_gns))
print(AIC(m_gns))

moran_gns <- moran.mc(
  residuals(m_gns),
  listw = queen_lw,
  nsim = 999,
  zero.policy = TRUE
)

cat("\nMoran residuos GNS:\n")
print(moran_gns)

# ------------------------------------------------------------
# 13. Tabla resumen
# ------------------------------------------------------------

tabla_final <- data.frame(
  Modelo = c(
    "OLS",
    "SLX completo",
    "SLX reducido",
    "SEM",
    "SAR",
    "SDM",
    "SDEM",
    "SAC",
    "GNS"
  ),
  AIC = c(
    AIC(m_ols),
    AIC(m_slx),
    AIC(m_slx_red),
    AIC(m_sem),
    AIC(m_sar),
    AIC(m_sdm),
    AIC(m_sdem),
    AIC(m_sac),
    AIC(m_gns)
  ),
  Moran_I_residuos = c(
    moran_ols$statistic,
    moran_slx$statistic,
    moran_slx_red$statistic,
    moran_sem$statistic,
    moran_sar$statistic,
    moran_sdm$statistic,
    moran_sdem$statistic,
    moran_sac$statistic,
    moran_gns$statistic
  ),
  p_Moran_residuos = c(
    moran_ols$p.value,
    moran_slx$p.value,
    moran_slx_red$p.value,
    moran_sem$p.value,
    moran_sar$p.value,
    moran_sdm$p.value,
    moran_sdem$p.value,
    moran_sac$p.value,
    moran_gns$p.value
  )
)

print(tabla_final)

write.csv(
  tabla_modelos,
  "resultados_comparacion_modelos_sar.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Guardar predicciones y residuos
# ------------------------------------------------------------

gdf$pred_OLS <- fitted(m_ols)
gdf$res_OLS <- residuals(m_ols)

gdf$pred_SLX <- fitted(m_slx)
gdf$res_SLX <- residuals(m_slx)

gdf$pred_SEM <- fitted(m_sem)
gdf$res_SEM <- residuals(m_sem)

gdf$pred_SAR <- fitted(m_sar)
gdf$res_SAR <- residuals(m_sar)

gdf$pred_SDM <- fitted(m_sdm)
gdf$res_SDM <- residuals(m_sdm)

gdf$pred_SDEM <- fitted(m_sdem)
gdf$res_SDEM <- residuals(m_sdem)

gdf$pred_SAC <- fitted(m_sac)
gdf$res_SAC <- residuals(m_sac)

gdf$pred_GNS <- fitted(m_gns)
gdf$res_GNS <- residuals(m_gns)

st_write(
  gdf,
  "cuencas_modelos_sar_resultados.gpkg",
  layer = "cuencas_modelos_sar_resultados",
  delete_layer = TRUE
)