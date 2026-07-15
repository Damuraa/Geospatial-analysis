# ============================================================
# MODELOS JERÁRQUICOS PARA FRA - CUENCAS
# ============================================================

# ------------------------------------------------------------
# 1. Paquetes
# ------------------------------------------------------------

paquetes <- c(
  "sf",
  "dplyr",
  "lme4",
  "lmerTest",
  "performance",
  "spdep"
)

instalar <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]

if (length(instalar) > 0) {
  install.packages(instalar, dependencies = TRUE)
}

library(sf)
library(dplyr)
library(lme4)
library(lmerTest)
library(performance)
library(spdep)

# Para evitar problemas de topología con geometrías geográficas
sf::sf_use_s2(FALSE)

# ------------------------------------------------------------
# 2. Rutas
# ------------------------------------------------------------

rutas_posibles <- c(
  "carto/CUENCAS_FINAL.gpkg",
  "CUENCAS_FINAL.gpkg"
)

ruta <- rutas_posibles[file.exists(rutas_posibles)][1]

if (is.na(ruta)) {
  stop("No se encontró CUENCAS_FINAL.gpkg. Revisa si está en la carpeta carto/ o en el directorio actual.")
}

capas <- st_layers(ruta)$name

if ("cuencas" %in% capas) {
  capa <- "cuencas"
} else {
  capa <- capas[1]
}

cat("\nArchivo usado:\n")
print(ruta)

cat("\nCapa usada:\n")
print(capa)

# ------------------------------------------------------------
# 3. Cargar datos
# ------------------------------------------------------------

gdf <- st_read(
  ruta,
  layer = capa
)

gdf <- st_make_valid(gdf)

cat("\nCRS:\n")
print(st_crs(gdf))

cat("\nColumnas:\n")
print(names(gdf))

cat("\nNúmero inicial de cuencas:\n")
print(nrow(gdf))

# ------------------------------------------------------------
# 4. Preparar variables
# ------------------------------------------------------------

# Si la capa está en WGS84, la reproyectamos para análisis espacial posterior
# Esto no cambia las variables, solo mejora las operaciones espaciales.
if (st_is_longlat(gdf)) {
  gdf <- st_transform(gdf, 9377)
}

# Crear area_km2 de forma robusta
# Tu archivo trae la columna como area.km2.
posibles_area <- c(
  "area[km2]",
  "area_km2",
  "area.km2.",
  "area.km2",
  "Area_km2",
  "AREA_KM2"
)

col_area <- posibles_area[posibles_area %in% names(gdf)][1]

if (is.na(col_area)) {
  stop(
    paste(
      "No encuentro columna de área. Las columnas disponibles son:",
      paste(names(gdf), collapse = ", ")
    )
  )
}

cat("\nColumna de área usada:\n")
print(col_area)

gdf$area_km2 <- as.numeric(gdf[[col_area]])

# Verificar columnas necesarias
columnas_necesarias <- c(
  "FRA",
  "area_km2",
  "Beta_media",
  "Beta_std",
  "nom_zh",
  "nom_szh"
)

faltantes <- columnas_necesarias[!(columnas_necesarias %in% names(gdf))]

if (length(faltantes) > 0) {
  stop(
    paste(
      "Faltan estas columnas:",
      paste(faltantes, collapse = ", ")
    )
  )
}

gdf <- gdf %>%
  mutate(
    FRA = as.numeric(FRA),
    area_km2 = as.numeric(area_km2),
    Beta_media = as.numeric(Beta_media),
    Beta_std = as.numeric(Beta_std),
    nom_zh = as.factor(nom_zh),
    nom_szh = as.factor(nom_szh)
  ) %>%
  filter(
    !is.na(FRA),
    !is.na(area_km2),
    !is.na(Beta_media),
    !is.na(Beta_std),
    !is.na(nom_zh),
    !is.na(nom_szh),
    FRA >= 0,
    FRA <= 1
  )

# Estandarizar predictoras
gdf <- gdf %>%
  mutate(
    area_km2_z = as.numeric(scale(area_km2)),
    Beta_media_z = as.numeric(scale(Beta_media)),
    Beta_std_z = as.numeric(scale(Beta_std))
  )

# Data frame sin geometría para modelos lmer
datos <- st_drop_geometry(gdf)

cat("\nNúmero de cuencas usadas:\n")
print(nrow(datos))

cat("\nFrecuencia por zona hidrográfica nom_zh:\n")
print(table(datos$nom_zh))

cat("\nFrecuencia por subzona hidrográfica nom_szh:\n")
print(table(datos$nom_szh))

# ------------------------------------------------------------
# 5. Fórmula base
# ------------------------------------------------------------

formula_base <- FRA ~ area_km2_z + Beta_media_z + Beta_std_z

# ------------------------------------------------------------
# 6. Modelo OLS base
# ------------------------------------------------------------

cat("\n================ OLS BASE ================\n")

m_ols <- lm(
  formula_base,
  data = datos
)

print(summary(m_ols))

# ------------------------------------------------------------
# 7. Función segura para ajustar modelos mixtos
# ------------------------------------------------------------

ajustar_lmer <- function(nombre, formula_modelo, datos) {
  
  cat("\n================", nombre, "================\n")
  
  modelo <- tryCatch(
    {
      lmer(
        formula_modelo,
        data = datos,
        REML = FALSE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 100000)
        )
      )
    },
    error = function(e) {
      cat("\nERROR en", nombre, ":\n")
      print(e$message)
      return(NULL)
    }
  )
  
  if (!is.null(modelo)) {
    print(summary(modelo))
    
    if (isSingular(modelo, tol = 1e-4)) {
      cat("\nADVERTENCIA: El modelo", nombre, "es singular. Puede estar sobreajustado o tener pocos datos por grupo.\n")
    }
  }
  
  return(modelo)
}

# ------------------------------------------------------------
# 8. Modelos jerárquicos
# ------------------------------------------------------------

# J1: intercepto aleatorio por zona hidrográfica
m_j1 <- ajustar_lmer(
  "J1: intercepto aleatorio por nom_zh",
  FRA ~ area_km2_z + Beta_media_z + Beta_std_z + 
    (1 | nom_zh),
  datos
)

# J2: intercepto y pendiente aleatoria de Beta_media por zona
# Se usa || para evitar estimar correlación intercepto-pendiente y mejorar estabilidad.
m_j2 <- ajustar_lmer(
  "J2: intercepto + pendiente Beta_media por nom_zh",
  FRA ~ area_km2_z + Beta_media_z + Beta_std_z + 
    (1 + Beta_media_z || nom_zh),
  datos
)

# J3: intercepto aleatorio por subzona
m_j3 <- ajustar_lmer(
  "J3: intercepto aleatorio por nom_szh",
  FRA ~ area_km2_z + Beta_media_z + Beta_std_z + 
    (1 | nom_szh),
  datos
)

# J4: estructura anidada zona/subzona
m_j4 <- ajustar_lmer(
  "J4: nom_zh / nom_szh anidado",
  FRA ~ area_km2_z + Beta_media_z + Beta_std_z + 
    (1 | nom_zh / nom_szh),
  datos
)

# ------------------------------------------------------------
# 9. Guardar modelos en lista
# ------------------------------------------------------------

modelos <- list(
  OLS = m_ols,
  J1_ZH = m_j1,
  J2_ZH_BETA = m_j2,
  J3_SZH = m_j3,
  J4_ZH_SZH = m_j4
)

modelos <- modelos[!sapply(modelos, is.null)]

# ------------------------------------------------------------
# 10. Función para R2
# ------------------------------------------------------------

extraer_r2 <- function(modelo) {
  
  if (inherits(modelo, "lm") & !inherits(modelo, "merMod")) {
    
    r2_marginal <- summary(modelo)$r.squared
    r2_condicional <- summary(modelo)$r.squared
    
  } else {
    
    r2_obj <- tryCatch(
      {
        performance::r2_nakagawa(modelo)
      },
      error = function(e) {
        return(NULL)
      }
    )
    
    if (is.null(r2_obj)) {
      r2_marginal <- NA
      r2_condicional <- NA
    } else {
      r2_marginal <- as.numeric(r2_obj$R2_marginal)
      r2_condicional <- as.numeric(r2_obj$R2_conditional)
    }
  }
  
  return(
    c(
      R2_marginal = r2_marginal,
      R2_condicional = r2_condicional
    )
  )
}

# ------------------------------------------------------------
# 11. Función para métricas
# ------------------------------------------------------------

calcular_metricas <- function(nombre, modelo, datos) {
  
  y_obs <- datos$FRA
  
  if (inherits(modelo, "lm") & !inherits(modelo, "merMod")) {
    y_pred <- as.numeric(predict(modelo, newdata = datos))
    singular <- FALSE
  } else {
    y_pred <- as.numeric(predict(modelo, newdata = datos, re.form = NULL))
    singular <- isSingular(modelo, tol = 1e-4)
  }
  
  residuos <- y_obs - y_pred
  
  sse <- sum(residuos^2, na.rm = TRUE)
  sst <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
  
  pseudo_r2_pred <- 1 - sse / sst
  mae <- mean(abs(residuos), na.rm = TRUE)
  rmse <- sqrt(mean(residuos^2, na.rm = TRUE))
  
  r2s <- extraer_r2(modelo)
  
  data.frame(
    Modelo = nombre,
    AIC = as.numeric(AIC(modelo)),
    BIC = as.numeric(BIC(modelo)),
    logLik = as.numeric(logLik(modelo)),
    n_param = attr(logLik(modelo), "df"),
    R2_marginal = r2s["R2_marginal"],
    R2_condicional = r2s["R2_condicional"],
    Pseudo_R2_pred = pseudo_r2_pred,
    MAE = mae,
    RMSE = rmse,
    Singular = singular
  )
}

tabla_metricas <- bind_rows(
  lapply(
    names(modelos),
    function(nm) calcular_metricas(nm, modelos[[nm]], datos)
  )
)

cat("\n================ MÉTRICAS MODELOS JERÁRQUICOS ================\n")
print(tabla_metricas)

write.csv(
  tabla_metricas,
  "metricas_modelos_jerarquicos_R.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Moran de residuos
# ------------------------------------------------------------

cat("\n================ MATRIZ QUEEN PARA RESIDUOS ================\n")

queen_nb <- poly2nb(
  gdf,
  queen = TRUE
)

queen_lw <- nb2listw(
  queen_nb,
  style = "W",
  zero.policy = TRUE
)

print(summary(queen_nb))

calcular_moran_residuos <- function(nombre, modelo) {
  
  res <- residuals(modelo)
  
  moran <- moran.mc(
    res,
    listw = queen_lw,
    nsim = 999,
    zero.policy = TRUE
  )
  
  data.frame(
    Modelo = nombre,
    Moran_I_residuos = as.numeric(moran$statistic),
    p_Moran_residuos = as.numeric(moran$p.value)
  )
}

tabla_moran <- bind_rows(
  lapply(
    names(modelos),
    function(nm) calcular_moran_residuos(nm, modelos[[nm]])
  )
)

cat("\n================ MORAN RESIDUOS MODELOS JERÁRQUICOS ================\n")
print(tabla_moran)

write.csv(
  tabla_moran,
  "moran_residuos_modelos_jerarquicos_R.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Agregar predicciones y residuos a capa espacial
# ------------------------------------------------------------

for (nm in names(modelos)) {
  
  modelo <- modelos[[nm]]
  
  pred_col <- paste0("pred_", nm)
  res_col <- paste0("res_", nm)
  abs_col <- paste0("res_abs_", nm)
  
  if (inherits(modelo, "lm") & !inherits(modelo, "merMod")) {
    pred <- as.numeric(predict(modelo, newdata = datos))
  } else {
    pred <- as.numeric(predict(modelo, newdata = datos, re.form = NULL))
  }
  
  gdf[[pred_col]] <- pred
  gdf[[res_col]] <- gdf$FRA - pred
  gdf[[abs_col]] <- abs(gdf[[res_col]])
}

# ------------------------------------------------------------
# 14. Guardar capa final
# ------------------------------------------------------------

if (!dir.exists("carto")) {
  dir.create("carto", recursive = TRUE)
}

st_write(
  gdf,
  "carto/cuencas_modelos_jerarquicos_R.gpkg",
  layer = "cuencas_modelos_jerarquicos_R",
  delete_layer = TRUE
)

# ------------------------------------------------------------
# 15. Mensaje final
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ANÁLISIS JERÁRQUICO FINALIZADO\n")
cat("Archivos generados:\n")
cat("- metricas_modelos_jerarquicos_R.csv\n")
cat("- moran_residuos_modelos_jerarquicos_R.csv\n")
cat("- carto/cuencas_modelos_jerarquicos_R.gpkg\n")
cat("============================================================\n")