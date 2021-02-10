# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_daq as daq
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import math
import numpy as np
import re
import os
pd.options.plotting.backend = "plotly"

# Create the dashboard
app = dash.Dash(__name__)


def get_cond(rho_,T_):
    A1=-105.161
    A2=+0.9007
    A3=+0.0007
    A4=+3.50e-15
    A5=+3.76e-10
    A6=+0.7500
    A7=+0.0017
    # Evaluate conductivity in mW/(m.K)
    cond=(A1+A2*rho_+A3*np.power(rho_,2)+A4*np.power(rho_,3)*np.power(T_,3)+A5*np.power(rho_,4)+A6*T_+A7*np.power(T_,2))/np.power(T_,0.5)
    k_=0.001*cond
    return k_


# Define here some parameters
Rcst=8.314           # J/(mol.K)
Wmlr=44.01e-3        # kg/mol
Cp=40.0/Wmlr         # J/(kg.K)
Cv=Cp-Rcst/Wmlr      # J/(kg.K)
Gamma=Cp/Cv          # -
Tinlet=430           # K
Tinit=300            # K
Tw=300               # K
MFR=0.2              # kg/m^3
Minit=9.87217        # kg

# Vessel geometry
radius=0.4           # m
length=2.5           # m
area=2.0*math.pi*radius*length+2.0*math.pi*radius*radius
Vtotal=math.pi*radius*radius*length # m^3

# Heat transfer model
dx=2.6/80            # m
k=0.01               # W/(K.m)

AA=Minit
BB=MFR
CC=MFR+k*area/(Cv*dx)
DD=-MFR*Gamma*Tinlet-Tw*k*area/(Cv*dx)

# Temperature graph
def create_Tfig():
    
    # Read the data with Twall=300
    #df=pd.read_csv('monitor_Twall300/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5],names=['Time','Temp','Mass','Pres'])
    #df['Tadia']=(df['Mass'].iloc[0]*df['Temp'].iloc[0]+Gamma*Tinlet*(df['Mass']-df['Mass'].iloc[0]))/df['Mass']
    
    # Create fig
    Tfig=go.Figure()
    
    # Read and plot the Farther Farms data
    df=pd.read_csv('FFdata/temperature.txt',delim_whitespace=True,header=None,skiprows=0,usecols=[0,1],names=['Time','Temp'])
    Tfig.add_trace(go.Scatter(name='Farther Farms experiment',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(width=2)))
    
    #Tfig.add_trace(go.Scatter(name='0D adiabatic model',x=df['Time']/60,y=df['Tadia'],mode='lines',showlegend=True,line=dict(color='firebrick',width=2,dash='dot')))
    #Tfig.add_trace(go.Scatter(name='NGA2 with Twall=300K',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(color='blue',width=2)))
    Tfig.update_layout(width=1200,height=800)
    Tfig.update_xaxes(title_text='Time (min)',title_font_size=24,tickfont_size=24,range=[0,20])
    Tfig.update_yaxes(title_text='Temperature (K)',title_font_size=24,tickfont_size=24)#,range=[280,450])
    Tfig.add_shape(type='line',x0=0,y0=Tinit,x1=df['Time'].iloc[-1]/60,y1=Tinit,line_color='black')
    Tfig.add_annotation(x=7.5,y=Tinit-7,text='Tinit',showarrow=False,font_size=16,font_color='black')
    Tfig.add_shape(type='line',x0=0,y0=Tinlet,x1=df['Time'].iloc[-1]/60,y1=Tinlet,line_color='green')
    Tfig.add_annotation(x=7.5,y=Tinlet+7,text='Tinlet',showarrow=False,font_size=16,font_color='green')
    Tfig.update_layout(legend=dict(font=dict(size=14)))
    
    # Add adiabatic temperature data
    #df=pd.read_csv('monitor_adiabatic/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5],names=['Time','Temp','Mass','Pres'])
    #Tfig.add_trace(go.Scatter(name='NGA2 with adbiatatic walls',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(color='firebrick',width=2)))
    
    # Add temperature data with Twall=350
    #df=pd.read_csv('monitor_Twall350/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5],names=['Time','Temp','Mass','Pres'])
    #Tfig.add_trace(go.Scatter(name='NGA2 with Twall=350K',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(color='navy',width=2)))
    
    # Add wall-modeled temperature data
    #df=pd.read_csv('monitor_wallmodel/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5,6],names=['Time','Temp','Mass','Pres','Twall'])
    #Tfig.add_trace(go.Scatter(name='NGA2 with heat losses - vessel',x=df['Time']/60,y=df['Temp' ],mode='lines',showlegend=True,line=dict(width=2)))
    #Tfig.add_trace(go.Scatter(name='NGA2 with heat losses - wall',  x=df['Time']/60,y=df['Twall'],mode='lines',showlegend=True,line=dict(width=2)))
    
    # Add temperature running now
    df=pd.read_csv('monitor/conservation',delim_whitespace=True,header=None,skiprows=2,usecols=[1,3,4,5,6],names=['Time','Temp','Mass','Pres','Twall'])
    Tfig.add_trace(go.Scatter(name='NGA2 simulation - Twall=300K',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(width=2)))
    #Tfig.add_trace(go.Scatter(name='runnig now - wall',  x=df['Time']/60,y=df['Twall'],mode='lines',showlegend=True,line=dict(color='navy',width=2)))
    
    
    
    # Add various models - analytical
    # df['Mass']=Minit+df['Time']*MFR
    #
    # AA=Minit
    # BB=MFR
    # CC=MFR
    # DD=-MFR*Gamma*Tinlet
    # df['Tadia']=np.power(AA,CC/BB)*(Tinit+DD/CC)*np.power(AA+BB*df['Time'],-CC/BB)-DD/CC
    # Tfig.add_trace(go.Scatter(name='Adia model',x=df['Time']/60,y=df['Tadia'],mode='lines',showlegend=True,line=dict(width=2)))
    #
    # Tw=300
    # k=0.02
    # Ra=radius*Cp*MFR/k
    # coeff=3
    # h=coeff*k/radius
    # AA=Minit
    # BB=MFR
    # CC=MFR+h*area/Cv
    # DD=-MFR*Gamma*Tinlet-Tw*h*area/Cv
    # df['T300']=np.power(AA,CC/BB)*(Tinit+DD/CC)*np.power(AA+BB*df['Time'],-CC/BB)-DD/CC
    # Tfig.add_trace(go.Scatter(name='T=300K model',x=df['Time']/60,y=df['T300'],mode='lines',showlegend=True,line=dict(width=2)))
    
    # Numerical instead - That works well
    # Tw=300
    # coeff=1
    # L_steel=0.0762
    # rho_steel=8050.0
    # Cp_steel=500.0
    # myTime=np.linspace(0,2000,20000)
    # myMass=np.zeros(len(myTime))
    # myMass[0]=Minit
    # myTemp=np.zeros(len(myTime))
    # myTemp[0]=Tinit
    # myMassTemp=myMass*myTemp
    # myTwall=np.zeros(len(myTime))
    # myTwall[0]=Tw
    #
    # timescale_fluid=15  # 15 sec?
    # timescale_wall =100
    #
    # for n in range(0,len(myTime)-1):
    #     k=2#get_cond(myMass[n]/Vtotal,myTemp[n])
    #     h=coeff*k/dx
    #
    #     myTwall[n+1]=myTwall[n]#+(myTime[n+1]-myTime[n])*10*h/(Cp_steel*rho_steel*L_steel)*(myTemp[n]-myTwall[n])
    #     #myTwall[n+1]=myTwall[n]+(myTime[n+1]-myTime[n])*(myTemp[n]-myTwall[n])/timescale_wall
    #
    #     myMass[n+1]=myMass[n]+(myTime[n+1]-myTime[n])*MFR
    #     myMassTemp[n+1]=myMassTemp[n]+(myTime[n+1]-myTime[n])*(MFR*Gamma*Tinlet-h*area/Cv*(myTemp[n]-myTwall[n]))
    #     #myMassTemp[n+1]=myMassTemp[n]+(myTime[n+1]-myTime[n])*(MFR*Gamma*Tinlet-myMass[n+1]*(myTemp[n]-myTwall[n])/timescale_fluid)
    #     myTemp[n+1]=myMassTemp[n+1]/myMass[n+1]
    #
    # df=pd.DataFrame(list(zip(myTime,myTemp,myTwall)),columns=['Time','Temp','Twall'])
    # Tfig.add_trace(go.Scatter(name='Python T=300K - vessel',x=df['Time']/60,y=df['Temp'],mode='lines',showlegend=True,line=dict(width=2)))
    # Tfig.add_trace(go.Scatter(name='Python T=300K - wall',x=df['Time']/60,y=df['Twall'],mode='lines',showlegend=True,line=dict(width=2)))
    
    # Numerical instead
    #Tw=300
    
    #df=pd.read_csv('FFdata/inlet_temp.txt',delim_whitespace=True,header=None,skiprows=0,usecols=[0,1],names=['Time','Tin'])
    
    #myTime=df['Time'].tolist()
    ##myTime=np.linspace(0,2000,20000)
    #myTin=df['Tin'].tolist()
    #myMass=np.zeros(len(myTime))
    #myMass[0]=Minit
    #myTemp=np.zeros(len(myTime))
    #myTemp[0]=Tinit
    #myMassTemp=myMass*myTemp
    #myTwall=np.zeros(len(myTime))
    #myTwall[0]=Tw
    

    #k=45 # W/(K.m) Steel reduced due to teflon
    #L_steel=0.0762
    #coeff=0.5
    #timescale_wall=30
    #timescale_fluid=0.1
    #for n in range(0,len(myTime)-1):
    #    myMass[n+1]    =myMass[n]    +(myTime[n+1]-myTime[n])*MFR
    #    myTwall[n+1]   =myTwall[n]   +(myTime[n+1]-myTime[n])*(myTemp[n]-myTwall[n])/timescale_wall
    
        #myMassTemp[n+1]=myMassTemp[n]+(myTime[n+1]-myTime[n])*(MFR*Gamma*myTin[n]-coeff*k*area/(Cv*L_steel)*(myTemp[n]-myTwall[n]))
    #    myMassTemp[n+1]=myMassTemp[n]+(myTime[n+1]-myTime[n])*(MFR*Gamma*myTin[n]-(myTemp[n]-myTwall[n])/timescale_fluid)
    #    myTemp[n+1]    =myMassTemp[n+1]/myMass[n+1]
    
    #df=pd.DataFrame(list(zip(myTime,myTemp,myTwall,myTin)),columns=['Time','Temp','Twall','Tin'])
    #Tfig.add_trace(go.Scatter(name='Python T=300K',x=df['Time']/60,y=df['Tin']  ,mode='lines',showlegend=True,line=dict(width=2)))
    #Tfig.add_trace(go.Scatter(name='Python T=300K',x=df['Time']/60,y=df['Temp'] ,mode='lines',showlegend=True,line=dict(width=2)))
    #Tfig.add_trace(go.Scatter(name='Python T=300K - wall'  ,x=df['Time']/60,y=df['Twall'],mode='lines',showlegend=True,line=dict(width=2)))
    
    
    return Tfig




# This is where we define the dashboard layout
def serve_layout():
    return html.Div(style={"margin-left": "15px"},children=[
    # Title of doc
    dcc.Markdown('''# Farther Farms Project'''),
    dcc.Markdown('''*NGA2 Dashboard written by O. Desjardins, last updated 02/06/2021*'''),
    # Intro
    dcc.Markdown('''
    ## Overview
    In this dashboard, we post-process the raw data generated by NGA2's pvessel
    case. This simulation is based on an experiment done by Farther Farms where
    a pressure vessel is filled with heated CO2.
    '''),
    # Imbibed volume over time
    dcc.Markdown(f'''
    ---
    ## Average temperature in the vessel
    The graph below shows the evolution of the average temperature inside the pressurized vessel.
    '''),
    #html.Div(create_Tfig(),style={'display':'none'}),
    dcc.Graph(id='Tgraph',figure=create_Tfig()),
])


# This is where we set the layout and run the server
app.layout = serve_layout
if __name__ == '__main__':
    app.run_server(debug=True)