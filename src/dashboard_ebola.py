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
Blues = sequential.Blues
Greens = sequential.Greens
Reds = sequential.Reds
Purples = sequential.Purples



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
# Definicion de variables

# w - Parametro que mide el impacto de la guerra (afecta al numero de infectados y la disminucion en hospitalizacion)
# w esta entre cero y uno
w = 0.8
# beta - Es una funcion de w que representa impacto de la guerra y la tasa de contacto
def beta_w(w, beta_max = 0.4, k =10, A = 100):
    return beta_max / (1 + A * np.exp(-k * w))
# S - Poblacion susceptible
S = 86000000
#Q- Poblacion en cuarentena
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
# delta_i - tasa de hopitalización de los infecciosos
delta_i = 0.008 # 0.4 # 0.1 / 12
#delta_q - tasa de hospitalización de los de cuarentena
delta_q= 0.1
# gamma_i - tasa de recuperacion infecciosos (Tarda 18 dias en recuperarse y lo hace en 60% de los casos)
gamma_i = 0.0003
# gamma_h - tasa de recuperacion de los hospitalizados
gamma_h = 0.0008
# sigma_i - tasa de muerte de los infecciosos
sigma_i = 0.007 # 1/3.8
# sigma_h - tasa de muerte hospitalizados
sigma_h = 0.005 # 1/1.7
#  pi - Tasa de natalidad
# pi = 1/(63*12*30) # 1.7
pi = 1.7
# mu - tasa de muerte natural
# mu = 1/(63*12*30)
mu = 1 / 756
# rho - tasa de deshecho de cuerpos
rho = 1/30
# alpha_h - tasa de infeccion en hospitales 
alpha_h = 0.5
# alpha_m - tasa de infeccion de muertos 
alpha_m = 0.9
# q - tasa de cuarentena
q =  0.3 # 0.0476
# v - tasa de vacunacion
v = 0.15
# Dias - dias se simulacion
dias = 1000
# lambda
# lam = 0.1
# *************************************************************************
# ************************ NUEVOS PARAMETROS ******************************
lam_s = 0.1 # 1.8787826654817788e-05 # Parametro place holder / Tasa de infección de sinomáticos
lam_a = 0.1 # 1.8787826654817788e-05 # Parametro place holder / Tasa de infección de asintomáticos
epsilon = 0.0001 # Parametro place holder / Pérdida de inmunidad de vacunados
Omega = 0.003 # Parametro place holder / Son los que salen de cuarentena
Gamma = 0.0002 # Parametro place holder  / Se infectan estando en cuarentena
psi = 0.0000033 # Parametro place holder / Tasa de recaída (Como Herpes-Zoster)
# *************************************************************************
# ************************ NUEVOS PARAMETROS ******************************

# print('\n\nNUEVO\n\n')

# ------------------------------------------------------------------------------------------
def diff_eqs(INP, t, w, delta_i, delta_q, gamma_i, gamma_h, sigma_h, sigma_i, mu, rho, alpha_h, alpha_m, q, v, pi, lam_s, lam_a, epsilon, Omega, Gamma, psi):
    """
    Sistema de ecuaciones diferenciales.
    """
    # Cantidades actuales de cada compartimento
    S, Q, I, H, R, D = INP
    # Total de la poblacion
    N = S + Q + I + R + H + D
    # print(N)
    beta = beta_w(w)               # Tasa de contacto efectivo de subpoblación asintomática a susceptible

    lam = beta * (I + alpha_h * H + alpha_m * D)/N

    # print(lam)
    
    # Nuevas; 

    dSdt = pi - (lam_s + lam_a + (1 - w) * v + mu + q) * S + epsilon * R + Omega * Q # pi - (lam + v + mu + q) * S
    dQdt = q * S - (mu + Gamma + Omega) * Q # q * S - (1-w) * delta_q * Q - mu * Q
    dIdt = lam_s * S - ((1 - w) * delta_i + gamma_i + sigma_i + mu) * I + (1 - psi) * R + Gamma * Q
    dHdt = (1 - w) * (delta_i * I) - (sigma_h + gamma_h + mu) * H
    dRdt = (1 - w) * v * S + gamma_i * I + gamma_h * H - mu * R - (1 - psi) * R + lam_a * S - epsilon * R
    dDdt = sigma_i * I + sigma_h * H - rho * D

    return [dSdt, dQdt, dIdt, dHdt, dRdt, dDdt]

# -----------------------------------------------------------------------------------


def ode_solver(t, initial_conditions, params):
    S, Q, I, H, R, D = initial_conditions
    w, delta_i, delta_q, gamma_i, gamma_h, sigma_h, sigma_i, mu, rho, alpha_h, alpha_m, q, v, pi, lam_s, lam_a, epsilon, Omega, Gamma, psi = params
    res = odeint(diff_eqs, [S, Q, I, H, R, D], t, args=(w, delta_i, delta_q, gamma_i, gamma_h, sigma_h, sigma_i, mu, rho, alpha_h, alpha_m, q, v, pi, lam_s, lam_a, epsilon, Omega, Gamma, psi))
    return res


# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------


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

                html.Div(className='opciones-dias',
                    children=[
                        html.Label(['Duración:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-dias',
                            options=[{'label':'10 días', 'value':'10'}] + [{'label': str(p) + ' días', 'value':p} for p in list(range(0, 1000, 50))],
                            value='100',
                            clearable=False)
                    ],
                    style=dict(width='10%')),

                html.Div(className='opciones-susceptibles',
                    children=[
                        html.Label(['Susceptibles:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-s',
                            options=[{'label':'10 personas', 'value':'10'}] + [{'label': str(p) + ' personas', 'value':p} for p in list(range(0, int(N), 5000000))],
                            value=S,
                            clearable=False)
                    ],
                    style=dict(width='10%')),

                html.Div(className='opciones-cuarentena',
                    children=[
                        html.Label(['Cuarentena:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-q',
                            options=[{'label':'10 personas', 'value':'10'}] + [{'label': str(p) + ' personas', 'value':p} for p in list(range(0, int(N), 5000000))],
                            value=Q,
                            clearable=False)
                    ],
                    style=dict(width='10%')),

                html.Div(className='opciones-infectados',
                    children=[
                        html.Label(['Infectados:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-i',
                            options=[{'label':'10 personas', 'value':'10'}] + [{'label': str(p) + ' personas', 'value':p} for p in list(range(0, int(N), 5000000))],
                            value=I,
                            clearable=False)
                    ],
                    style=dict(width='10%')),

                html.Div(className='opciones-hospitalizados',
                    children=[
                        html.Label(['Hospitalizados:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-h',
                            options=[{'label':'10 personas', 'value':'10'}] + [{'label': str(p) + ' personas', 'value':p} for p in list(range(0, int(N), 5000000))],
                            value=H,
                            clearable=False)
                    ],
                    style=dict(width='10%')),


                html.Div(className='opciones-recuperados',
                    children=[
                        html.Label(['Recuperados:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-r',
                            options=[{'label':'10 personas', 'value':'10'}] + [{'label': str(p) + ' personas', 'value':p} for p in list(range(0, int(N), 5000000))],
                            value=R,
                            clearable=False)
                    ],
                    style=dict(width='10%')),

                html.Div(className='opciones-fallecidos',
                    children=[
                        html.Label(['Muertes:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-d',
                            options=[{'label':'10 personas', 'value':'10'}] + [{'label': str(p) + ' personas', 'value':p} for p in list(range(0, int(N), 5000000))],
                            value=D,
                            clearable=False)
                    ],
                    style=dict(width='10%')),

                html.Div(className='opciones-inestabilidad',
                    children=[
                        html.Label(['Nivel de inestabilidad:'], style={'font-weight': 'bold', "text-align": "center"}),
                        dcc.Dropdown(
                            id='dropdown-w',
                            options=[{'label': str(p / 10.0), 'value' : p / 10.0 } for p in list(range(0, 11, 1))],
                            value=w,
                            clearable=False)
                    ],
                    style=dict(width='10%')),

                ],
            style=dict(display='flex')
        ),

        dcc.Graph('enfermedad-en-tiempo', config={'displayModeBar': False}),

        dcc.Interval(id='interval-component', interval=1*1000)])])
# # -------------------------------------------------------------------------------------------

@app.callback(
    Output('enfermedad-en-tiempo', 'figure'), 
    [
        Input('dropdown-dias', 'value')
        , Input('dropdown-s', 'value')
        , Input('dropdown-q', 'value')
        , Input('dropdown-i', 'value')
        , Input('dropdown-h', 'value')
        , Input('dropdown-r', 'value')
        , Input('dropdown-d', 'value')
        , Input('dropdown-w', 'value')
    ])
def update_graph(dias, S, Q, I, H, R, D, w):

    initial_conditions = [S, Q, I, H, R, D]
    params = (w, delta_i, delta_q, gamma_i, gamma_h, sigma_h, sigma_i, mu, rho, alpha_h, alpha_m, q, v, pi, lam_s, lam_a, epsilon, Omega, Gamma, psi)
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
            , marker=dict(color=Greens[5])
            , mode='lines+markers'
            , name='S'
            , hovertext='Susceptibles'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=Q
            , x=list(range(len(S)))
            , marker=dict(color=Purples[5])
            , mode='lines+markers'
            , name='Q'
            , hovertext='Cuarentena'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=I
            , x=list(range(len(S)))
            , marker=dict(color=Reds[5])
            , mode='lines+markers'
            , name='I'
            , hovertext='Infectados'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=H
            , x=list(range(len(S)))
            , marker=dict(color=Oranges[5])
            , mode='lines+markers'
            , name='H'
            , hovertext='Hospitalizados'
            )
        , row=1
        , col=1)


    fig.add_trace(
        go.Scatter(
            y=R
            , x=list(range(len(S)))
            , marker=dict(color=Blues[5])
            , mode='lines+markers'
            , name='R'
            , hovertext='Recuperados'
            )
        , row=1
        , col=1)

    fig.add_trace(
        go.Scatter(
            y=D
            , x=list(range(len(S)))
            , marker=dict(color=Greys[6])
            , mode='lines+markers'
            , name='D'
            , hovertext='Muertes'
            )
        , row=1
        , col=1)


    fig.update_xaxes(title_text='Días', title_font={'size':12}, showgrid=False, row=1, col=1)
    fig.update_yaxes(title_text='Número de casos', title_font={'size':12}, showgrid=False, row=1, col=1)
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
