# Datos de paper de guerra
epid <- c(43, 102, 122, 142, 162, 238, 333, 421, 500, 585, 658, 743, 816, 872, 923, 1016, 1264, 1466, 1705, 1920, 2168, 2343, 2418, 2430)
# Numero de dias observados
n = as.numeric(length(epid))
# 1. Calculo de R0 usando paqueteria de R y EG
# ---------------------------------------------------------------------------------------------
library(R0)
# Exponential growth
plot(epid)
GT.flu <- generation.time("weibull", c(7, 1))
sensitivity.analysis(epid, GT.flu, begin = 1, end = n, est.method = "EG", sa.type = "time")
EG <- est.R0.EG(epid, GT.flu, begin = 1, end = n, nsim=1000)
print(EG)
plotfit(EG)
# ---------------------------------------------------------------------------------------------
help("generation.time")

# 2. Ahora calcularemos R0 usando el metodo lineal
# ---------------------------------------------------------------------------------------------
plot(epid, xlab='Semana', ylab='Infectados', type='l', col='red', main='Datos de Ebola reportados')
points(epid, xlab='Semana', ylab='Infectados', type='p')
# Ajuste logaritmico a numero de personas infectadas personas infectadas
epid_log <- log(epid)
# Visualizacion de datos logaritmocos
plot(epid_log, xlab='Semana', ylab='Infectados', type='l', col='red', main='Ajuste lineal')
points(epid_log, xlab='Semana', ylab='Infectados', type='p')
# Ajuste lineal
fit <- lm(epid_log~c(1:n))
summary(fit)
coeficientes <- coef(fit)
anova(fit)
# Error estandar de la pendiente
slope_se <- coef(summary(fit))[2,2]
# Ajuste lineal
abline(coeficientes[1], coeficientes[2], col='red')
# Bandas de confianza (creo)
abline(coeficientes[1] - 2.5 * coef(summary(fit))[1,2], coeficientes[2]- 2.5 * coef(summary(fit))[2,2], col='red', lty='dashed')
abline(coeficientes[1] + 2.5 * coef(summary(fit))[1,2], coeficientes[2] + 2.5 * coef(summary(fit))[2,2], col='red', lty='dashed')
# Ahora determianamos R0
# Tomamos 1/mu como 63
# Tomamos 1/gamma = 9
# Recordemos que m = (R0-1) * (mu + gamma)
R0 <- coeficientes[2]/(1/63 + 1/9) + 1
cat('El valor de Ro es', R0)
# Intervalo de confianza de R0
# 2.5 = 1/gamma
IC = c(R0 - 2.5 * slope_se, R0 + 2.5 * slope_se)
cat('Intervalo de confianza', IC)
# ---------------------------------------------------------------------------------------------


# 3. Ahora calcularemos R0 usando el metodo lineal buscando el R0 que minimiza la lognitud del IC
# ---------------------------------------------------------------------------------------------
R0_por_dia <- c()
R0_por_dia_cota_inferior <- c()
R0_por_dia_cota_superior <- c()
dia <- 2
for (dia in c(2:(n-3)))
{
  # epid_log[c(1:dia)]
  # c(1:dia)
  fit <- lm(epid_log[c(1:dia)]~c(1:dia))
  coeficientes <- coef(fit)
  slope_se <- coef(summary(fit))[2,2]
  R0 <- slope_se/(1/63 + 1/9) + 1
  # Intervalo de confianza de R0
  IC = c(R0 - 2.5 * slope_se, R0 + 2.5 * slope_se)
  R0_por_dia <- c(R0_por_dia, R0)
  R0_por_dia_cota_superior <- c(R0_por_dia_cota_superior, IC[2])
  R0_por_dia_cota_inferior <- c(R0_por_dia_cota_inferior, IC[1])
  cat('\n-----------','\nR0 en el dia', dia, ':', R0, '\nIC: [', IC[1], ',', IC[2], "]", '\nLongitud: ', IC[2]-IC[1])
}
plot(R0_por_dia, lty='solid', type='l', col='red', ylab='R0', xlab = 'Días considerados', main='Cálculo de R0 por distinto número de días')
points(R0_por_dia_cota_inferior, lty='dashed', type='l', col='blue')
points(R0_por_dia_cota_superior, lty='dashed', type='l', col='blue')
# plot(R0_por_dia_cota_superior - R0_por_dia_cota_inferior, main='Tamaño del Intervalo de confianza', ylab='Longitud de IC', xlab = 'Días', type='l', col='blue')
# ---------------------------------------------------------------------------------------------