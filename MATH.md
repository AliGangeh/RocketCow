# Introduction
## Purpose
The purpose of this code is as follows:
- Find equilibrium canditions via cantera
- Calculate Specific heat (without condensed species)
- Calculate thermo derivatives (for above calculations)

## Next Steps
- Calculate CEA rocket output properties: Isp, C*, etc.
- Calculate thermal conductivity, viscocity
- Calculate Bartz heat transfer coefficient

## Relevant Expressions & Variables
### Matrices
$a_{ij}$ is the stoichiometric coefficients matrix. Defined as the number of kilogram-moles of element $i$ per kilogram-mole species $j$ 

### Thermodynamic Properties
#### Specific Heat
Specific heat can be said to have a frozen and Equilibrium component. If $C_p = dH/dT$ where $H = nh$ where $n$ is moles and $h$ is molar enthalpy. Often $n$ is unchanging but if we have a reacting mixture then we get: $C_p = n\frac{dh}{dT} + \frac{dn}{dT}h$ 
This the first term is said to be the frozen specific heat as the $n$ is constant, while the second term is the reacting specific heat. As shown bellow
$$\begin{align}
c_{p,e}=c_{p,f} + c_{p,r}= \sum_{j=1}^{NS}n_jC_{p,j}^o+\sum_{j=1}^{NG} n_j\frac{H_j^o}{T}\left(\frac{\partial\ln n_j}{\partial \ln T}\right)_P+\sum_{j=NG+1}^{NS} \frac{H_j^o}{T}\left(\frac{\partial n_j}{\partial \ln T}\right)_P
\end{align}$$

$$\begin{align}
c_v \equiv c_p + \frac{\frac{PV}{T}\left(\frac{\partial \ln V}{\partial \ln T}\right)_p^2}{\left(\frac{\partial \ln V}{\partial \ln P}\right)_T}\tag{2.70}
\end{align}$$

#### Specific heat ratio
$$\begin{align}
\gamma \equiv \frac{c_p}{c_v}\tag{2.72}
\end{align}$$

$$\begin{align}
\gamma_s = \frac{\gamma}{\left(\frac{\partial \ln V}{\partial \ln P}\right)_T}\tag{2.73}
\end{align}$$

$$\begin{align}
a= \sqrt{nRT\gamma_s} \tag{2.74}
\end{align}$$

### Thermodyanamic Derivatives
#### Derivatives with Respect to Temperature
$$\begin{align}
\sum_{i=1}^\ell \sum_{j=1}^{NG} a_{kj} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \ln T} \right)_P 
&+ \sum_{j=NG+1}^{NS} a_{ij} \left( \frac{\partial n_j}{\partial \ln T} \right)_P \notag
\\ &+ \sum_{j=1}^{NG} a_{kj} n_j \left( \frac{\partial \ln n}{\partial \ln T} \right)_P 
= -\sum_{j=1}^{NG} a_{kj} n_j \frac{H_j^o}{RT} \qquad k = 1, \ldots, \ell
\end{align}$$

$$\begin{align}
\sum_{i=1}^\ell a_{ij} \left( \frac{\partial \pi_i}{\partial \ln T} \right)_P =-\frac{H_j^o}{RT}\qquad j=NG+1,\ldots,NS\tag{2}
\end{align}$$

$$\begin{align}
\sum_{i=1}^\ell \sum_{j=1}^{NG} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \ln T} \right)_P
=-\sum_{j=1}^{NG}\frac{n_jH_j^o}{RT}\tag{3}
\end{align}$$

#### Derivatives with Respect to Pressure
$$\begin{align}
\sum_{i=1}^\ell \sum_{j=1}^{NG} a_{kj} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \ln P} \right)_T
&+ \sum_{j=NG+1}^{NS} a_{ij} \left( \frac{\partial n_j}{\partial \ln P} \right)_T \notag
\\ &+ \sum_{j=1}^{NG} a_{kj} n_j \left( \frac{\partial \ln n}{\partial \ln P} \right)_T
= -\sum_{j=1}^{NG} a_{kj} n_j  \quad k = 1, \ldots, \ell\tag{4}
\end{align}$$

$$\begin{align}
\sum_{i=1}^\ell a_{ij} \left( \frac{\partial  \pi_i}{\partial \ln P} \right)_T 
=0\qquad j=NG+1,\ldots,NS\tag{5}
\end{align}$$

$$\begin{align}
\sum_{i=1}^\ell \sum_{j=1}^{NG} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \ln T} \right)_P
=\sum_{j=1}^{NG}n_j\tag{6}
\end{align}$$

#### Derivatives of volume
$$\begin{align}
\left( \frac{\partial  V}{\partial \ln T} \right)_P = 1 + \left( \frac{\partial  \ln n}{\partial \ln T} \right)_P\tag{2.50}
\end{align}$$

$$\begin{align}
\left( \frac{\partial  V}{\partial \ln P} \right)_T = -1 + \left( \frac{\partial  \ln n}{\partial \ln P} \right)_T\tag{2.51}
\end{align}$$

### Engine conditions
Initial estimate at throat pressure is given as
$$\begin{align}
\frac{P_{inf}}{P_t}=\left(\frac{\gamma_s+1}{2}\right)^{\gamma_s/(\gamma_s-1)}\tag{6.15}
\end{align}$$

### Characteristic velocity
from Sutton
$$\begin{align}
c^*=\frac{\sqrt{\gamma \frac{\bar R}{M}T}}{\gamma \sqrt{\left(\frac2{K+1}\right)^{\frac{\gamma+1}{\gamma-1}}}}\tag{3-32}
\end{align}$$
From which we can derive 
$$\begin{align}
c^*=\sqrt{\frac{RT}{\gamma M}}\left(\frac{2}{\gamma+1}\right)^{-\frac{\gamma+1}{2(\gamma-1)}}
\end{align}$$

### lagrangian multiplier 
$$\pi_i=-\gamma_i/RT$$ 
where $\gamma_i$ are the lagrangian multipliers such that that $\mu_j+\sum_{i=1}^\ell \lambda_ia_{ij}=0$

### Isentropic flow equation
$$\begin{align}
P_t = \frac{P_0}{\frac{\gamma+1}{2}^{\frac{\gamma}{\gamma-1}}} \tag{6.15}
\end{align}$$

https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html

### sub/supersonic nozzle conditions
$$\begin{align}
\ln \frac{P_c}{P_e}=\frac{\ln\left(\frac{P_c}{P_t}\right)}{\frac{A_e}{A_t}+10.587\left(\ln\frac{A_e}{A_t}\right)^3+9.454\ln\frac{A_e}{A_t}}\qquad 1.0001 < \frac{A_e}{A_t} < 1.09\tag{6.20}
\end{align}$$