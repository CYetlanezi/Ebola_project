# Simulador del modelo epidemiológico

El notebook interactivo modela el siguiente sistema de ecuaciones que corresponden a un modelo epidemiológico de ébola:

$\frac{dS}{dt} = \pi - (\lambda + v + \mu + q) S$

$\frac{dQ}{dt} = qS - (1-w) \delta_Q Q - \mu Q$

$\frac{dI}{dt} = \lambda S - ((1-w) \delta_I + \gamma_I + \sigma_I + \mu) I$

$\frac{dH}{dt} = (1-w) (\delta_Q Q + \delta_I I) - (\sigma_H + \gamma_H + \mu) H$

$\frac{dR}{dt} = v S + \gamma_I I + \gamma_H H - \mu R$

$\frac{dD}{dt} = \sigma_I I + \sigma_H H - \rho D$

con

$\lambda = \frac{\beta (I + \alpha_H H + \alpha_M D)}{N}$


## Requerimientos

- Anaconda 2.0 o superior.
- Ambiente de anaconda cuya configuración está conetenida en `conf/ebola.yml`. Para crear el ambiente con este archivo se puede utilizar el comando `conda env create -f conf/ebola.yml`. Después se activa con el comando `conda activate ebola`. Si hay problemas con paquetes específicos pueden intentar resolverlos borrando la versión de cada paquete en el archivo `ebola.yml`.
- Instalación de orca. Para instalarlo pueden usar el comando: `conda install -c conda-forge orca`.

