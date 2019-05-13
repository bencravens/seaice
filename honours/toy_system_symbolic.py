import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy
from sympy import init_printing
init_printing(use_latex=True)
"""
python script I am using to calculate the momentum flux tau
of an idealized ice floe with properties taken from the literature.
I will calculate the momentum flux as a function of the form drag for
sails and keels, floe edges, and the atmospheric and oceanic skin drag.
"""

#Defining constants 
#Ratio of keel depth and sail height (Worby (2008))
R_h = 4.4
#Ratio of average distance between keels and average distance between sails
R_d = 1.0
#Weight variables of sails and keels (how many sails vs how many keels)
#alpha in Tsamados Here naming W_s for weighting of sails
W_s = 0.5 
#beta in Tsamados. Here naming W_k for weighting of keels
W_k = 0.5 
#Slope of sails(rad) (Worby (2008))
alpha_r = 0.45
#Slope of keels(rad) (Worby (2008)) 
alpha_k = 0.45 
#Attenuation parameter in sheltering function (given by Tsamados)
s_l = 0.18 
#roughness length of level ice (given by CICE)
z_oi = 5e-4 
#Ice concentration (fraction of 1.0)
A = 0.5 
#Sail height(m) (Worby (2008)) 
H_s = 0.57 
#Ratio of aice to ardg, Ridged ice area fraction
R_f = 1.0/0.11
#dimensionless scaling coefficent for sails
c_ra = 0.2
#dimensionless scaling coefficent for keels
c_kw = 0.2
#dimensionless scaling coefficent for skin drag(atmosphere)
c_sa = 0.0005
#dimensionless scaling coefficent for skin drag(ocean)
c_sw = 0.002
#atmospheric skin drag tunable parameter (from sail height) (Tsamados)
m_a = 20.0
#oceanic skin drag tunable parameter (from keel depth) (Tsamados)
m_w = 10.0
#density of air (kg/m^3) (Wikipedia) (STP)
rho_a = 1.225
#density of ocean water (kg/m^3) (Wikipedia) (STP)
rho_w = 1025
#von karman constant
0.40

"""
Distance between sails (m) 
(Taken from equation in Tsamados) aice/ardg = 0.11 ,
as given in Worby (2008)
"""
D_s, H_s, 

 D_s(H_s):
    return 2.0*H_s*R_f*(W_s/math.tan(alpha_r) + (W_k/math.tan(alpha_k))*(R_h/R_d))

#Keel height(m)
def H_k(H_s):
    return R_h*H_s

#Distance between keels(m)
def D_k(D_s):
    return R_d*D_s

#Now defining individual drag components
#first sheltering function
def S_c(D,H):
    return pow((1.0 - np.exp(-s_l*D/H)),0.5)

#Now the form drag coefficient from sails
def C_dar(H_s,D_s):
    return 0.5*c_ra*pow(S_c(D_s,H_s),2.0)*(H_s/D_s)*A*pow((np.log(H_s/z_oi)/np.log(10/z_oi)),2.0)

#Form drag coefficient from keels
def C_dwr(H_k,D_k):
    return 0.5*c_kw*pow(S_c(D_k,H_k),2.0)*(H_k/D_k)*A*pow((np.log(H_k/z_oi)/np.log(10/z_oi)),2.0)

#Skin drag coefficient for atmosphere
def C_das(H_s,D_s):
    return A*(1-m_a*(H_s/D_s))*c_sa

#Skin drag coefficient for ocean
def C_dws(H_k,D_k):
    return A*(1-m_w*(H_k/D_k))*c_sw

#momentum flux from ice to atmosphere
def tau(z,H_s):
    #function to generate momentum flux at height z from sails to atmosphere
    return 

def plotdrag(H_s,form=False,skin=False):
    #function to plot dependence of drag coefficents on sail height... 2d. Works for both form and skin drag(skin drag kinda)
    assert form or skin, "you need to select either form or skin to plot"
    assert np.logical_xor(form,skin), "you can only plot form or skin... not both."
    D_s_temp = [D_s(i) for i in H_s]
    H_k_temp = [H_k(i) for i in H_s]
    D_k_temp = [D_k(j) for j in D_s_temp]
    #zip stores two arrays ["a","b"] and [1,2] as [["a",1],["b",2]]
    S_zp = zip(H_s,D_s_temp)
    K_zp = zip(H_k_temp,D_k_temp)
    if form:
        totalsail = [C_dar(arr[0],arr[1]) for arr in S_zp]
        totalkeel = [C_dwr(arr[0],arr[1]) for arr in K_zp]
    elif skin:
        totalsail = [C_das(arr[0],arr[1]) for arr in S_zp]
        totalkeel = [C_dws(arr[0],arr[1]) for arr in K_zp]
    totaldrag = np.asarray(totalsail) + np.asarray(totalkeel) 
    fig = plt.figure(dpi=500)
    plt.plot(H_s,totalsail,label="Sails")
    plt.plot(H_s,totalkeel,label="Keels")
    plt.plot(H_s,totaldrag,label="Total")
    plt.legend(loc="center right")
    if form:
        plt.title("How form drag changes with sail height")
    elif skin:
        plt.title("How skin drag changes with sail height")
    plt.xlabel("Height of Sails")
    plt.ylabel("Drag coefficient value")
    plt.xlim([min(H_s),max(H_s)])
    if skin:
        fig.savefig("/home/ben/Desktop/thesis/2D_plot_skin.png")
    elif form:
        fig.savefig("/home/ben/Desktop/thesis/2D_plot_form.png")

plotdrag(np.linspace(0.3,4.0,50),skin=True)
plotdrag(np.linspace(0.3,4.0,50),form=True)
