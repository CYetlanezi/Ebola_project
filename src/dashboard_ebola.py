#!/usr/bin/env python
# coding: utf-8


# Preparacion de ambiente
# -----------------------------------------------------------------------------------
# Importacion de librerias
import pandas as pd # Libreria para manejo de datos en forma de DataFrame
import numpy as np # Libreria para operaciones matemaicas
# from sqlalchemy import create_engine # Libreria para conectar con MySQL
import plotly.express as px # Libreria para crear graficos
from dash import Dash # Libreria para desplegar el dashboard
import dash_core_components as dcc # Libreria para agregar componentes al dashboard
import dash_html_components as html # Libreria para usar elementos de html en python 
from dash.dependencies import Input, Output # Libreria para 
from plotly.subplots import make_subplots # Libreria para hacer subgraficas
import plotly.graph_objects as go # Libreria para agregar graficas al dashboard 
from plotly.colors import sequential # Libreria para agregar paletas de colores
import time # Libreria para medir y saber tiempo
import datetime # Libreria para conocer la fecha
import os # Libreria para manipular el sistema operativo
import argparse
import json
Oranges = sequential.Oranges
Greys = sequential.Greys


import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint
import plotly.graph_objects as go
import plotly.io as pio


# Definicion de variables

# S - Poblacion susceptible
S = 86000000
# Q - Poblacion en cuarentena
Q = 0.0
# I - Poblacion Infectada
I = 2619
# H - Poblacion en hospitalizacion
H = 0
# R - Poblacion recuperada
R = 0
# D - Poblacion que murio
D = 1729
# N - Poblacion total
N = S + Q + I + H + R + D

w = 0.8
Gamma = 0.005
Omega = 0.005
epsilon = 0.00001
omega = 0.00001
phi = 0.0001
q = 0.1

# Dias - dias de simulacion
dias = 10


# ------------------------------------------------------------------------------------------


def lambda_x(Beta, I, H, D, N, alpha_h=0.5, alpha_m=0.9):
    """
    Esta función calcula la tasa de contagio de acuerdo a:
    Beta - La tasa de contacto
    I - El número de infectados
    H - El número de hospitalizados
    D - El número de muertos
    N - La población total
    alpha_h - Tasa de contagio de hospitalizados
    alpha_m - Tasa de contaagio de muertos
    """

    return Beta * (I + alpha_h * H + alpha_m * D) / N


def beta_w(w=0.8, beta_max=400, k=10, A=100, tipo=None):

    if tipo == 's':

        resp = beta_max / (1 + A * np.exp(-k * w))

    elif tipo == 'a':

        resp = beta_max / (1 + A * np.exp(-k * w))

    else:

        print('\n\tTipo no válido. Debe ser "a" o "s"\n')

    return resp


def modelo(pob_compartimiento, Gamma, Omega, epsilon, omega, phi, q, pi=1.7, delta_i=0.4, delta_q=0.1, gamma_h=1/4.4, gamma_i=1/10.8, mu=1/756, rho=1/30, sigma_h=1/1.7, sigma_i=1/3.8, v=0.15):
    """
    Sistema de ecuaciones diferenciales con los siguientes parámetros:
    Gamma - Tasa se supesión de la cuarentena
    Omega - Tasa de contagio de los individuos en cuarentena
    delta_i - Tasa de personas hospitalizadas provieniente del compartimiento de infectados
    delta_q - Tasa de personas hospitalizadas provieniente del compartimiento de cuarentena
    epsilon - Tasa de pérdida de inmunidad en vacunados
    gamma_h - Tasa de recuparación para los hispitalizados
    gamma_i - Tasa de recuparación para los infectados
    mu - Tasa de muerte natural
    omega - Nivel de estabilidad en la región
    phi - Tasa de recaida
    pi - Tasa de reclutamiento
    q - Tasa de ingreso a cuarentena
    rho - Tasa de desecho de los cuerpos
    sigma_h - Tasa de mortalidad de los individuos que recibieron atención médica
    sigma_i - Tasa de mortalidad de los individuos infectados no hospitalizados
    v - Tasa de vacunación
    """

    # Cantidades actuales de cada compartimento
    S, Q, I, H, R, D = pob_compartimiento

    # Total de la poblacion
    N = S + Q + I + R + H + D
    
    Beta_s = beta_w(w, tipo='s')               # Tasa de contacto efectivo de subpoblación asintomática a susceptible
    Beta_a = beta_w(w, tipo='a')
    
    lambda_s = lambda_x(Beta_s, I, H, D, N)
    lambda_a = lambda_x(Beta_a, I, H, D, N)
    
    dSdt = pi - (lambda_s + lambda_a + (1 - omega) * v + mu + q + epsilon) * S + Omega * Q

    dQdt = q * S - ((1 - omega) * delta_q - mu - Gamma - Omega) * Q

    dIdt = lambda_s * S - ((1 - omega) * delta_i + gamma_i + sigma_i + mu) * I + phi * R + Gamma * Q

    dHdt = (1 - omega) * (delta_q * Q + delta_i * I) - (sigma_h + gamma_h + mu) * H

    dRdt = (1 - omega) * v * S + gamma_i * I + gamma_h * H - mu * R - phi * R + lambda_a * S - epsilon * R

    dDdt = sigma_i * I + sigma_h * H - rho * D

    return [dSdt, dQdt, dIdt, dHdt, dRdt, dDdt]

# -----------------------------------------------------------------------------------


def ode_solver(t, initial_conditions, params):
    S, Q, I, H, R, D = initial_conditions

    Gamma, Omega, epsilon, omega, phi, q = params

    res = odeint(modelo, [S, Q, I, H, R, D], t, args=(Gamma, Omega, epsilon, omega, phi, q))

    return res


# ---------------------------------------------------------------------------------
initial_conditions = [S, Q, I, H, R, D]
params = (Gamma, Omega, epsilon, omega, phi, q)
# -----------------------------------------------------------------------------------


# Dashboard
# -----------------------------------------------------------------------------------
app = Dash(__name__, routes_pathname_prefix='/modelo_ebola/') # Ruta en la que se expondra el dashboard
server = app.server

# Definimos el layout del dashboard
# Se compone de 4 componentes independientes
app.layout = html.Div([
    html.Div([
        html.Div(className='menus-desplegables3',
            children=[
                html.Div(className='histograma-proceso',
                    children=[
                        dcc.Dropdown(
                            id='dias',
                            options=[{'label':'10', 'value':'10'}] + [{'label':p, 'value':p} for p in list(range(10, 1000, 50))],
                            value='10',
                            clearable=False)
                    ],
                    style=dict(width='25%')),
                ],
            style=dict(display='flex')
        ),

        dcc.Graph('enfermedad-en-tiempo', config={'displayModeBar': False}),

        dcc.Interval(id='interval-component', interval=1*1000)])])
# # -------------------------------------------------------------------------------------------

@app.callback(
    Output('enfermedad-en-tiempo', 'figure'), 
    [
        Input('dias', 'value')
    ])
def update_graph(dias):

    periodo = np.arange(0, int(dias), 1)
    sol = ode_solver(periodo, initial_conditions, params)
    S, Q, I, H, R, D = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], sol[:, 4], sol[:, 5]

    # Iniciamos dos graficos de barras. Uno con la informacion de usuarios conectados por alcaldia y otro con las paginas web que registraron mas visitas
    fig = make_subplots(rows=1, cols=1,
                        specs=[[{"type":"bar"}]])

    fig.add_trace(
        go.Scatter(
            y=S
            , x=list(range(len(S)))
            , marker=dict(color=Greys[0])
            , mode='lines+markers'
            , name='S'
            , hovertext='S'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=Q
            , x=list(range(len(S)))
            , marker=dict(color=Greys[1])
            , mode='lines+markers'
            , name='Q'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=I
            , x=list(range(len(S)))
            , marker=dict(color=Greys[2])
            , mode='lines+markers'
            , name='I'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=H
            , x=list(range(len(S)))
            , marker=dict(color=Greys[3])
            , mode='lines+markers'
            , name='H'
            )
        , row=1
        , col=1)


    fig.add_trace(
        go.Scatter(
            y=R
            , x=list(range(len(S)))
            , marker=dict(color=Greys[4])
            , mode='lines+markers'
            , name='R'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=D
            , x=list(range(len(S)))
            , marker=dict(color=Greys[5])
            , mode='lines+markers'
            , name='D'
            )
        , row=1
        , col=1)

    # fig.add_trace(
    #     go.Scatter(
    #         x=tiempo_dask.num
    #         , y=tiempo_dask.duration
    #         , marker=dict(color=Oranges[4])
    #         , mode='lines+markers'
    #         , name='Dask'
    #         )
    #     , row=1
    #     , col=1)


    fig.update_xaxes(title_text='Número de casos', title_font={'size':12}, showgrid=False, row=1, col=1)
    fig.update_yaxes(title_text='Días', title_font={'size':12}, showgrid=False, row=1, col=1)
    # Actualizamos el tamano y titulo del grafico y ubicacion de los titulos y posicion del bloque
    fig.update_layout(title_text="Modelo Epidemiológico de Ébola"
        , title_font={'size':25}
        , title_x=0.5
        , height=500
        , autosize=True
        # , width=1400
        , template='plotly_white'
        , legend=dict(orientation="h"
            , yanchor="bottom"
            , y=-0.3
            , xanchor="left"
            , x=0.45))

    # Regresamos el bloque de graficos 
    return fig

# # -------------------------------------------------------------------------------------------


if __name__ == '__main__':
        app.run_server(debug=True, port=9092)
