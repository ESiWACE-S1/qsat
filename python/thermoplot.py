#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import thermo

# plot esl, esi - saturation vapor pressure
# test various iterative solvers of ql from thl
# Fredrik Jansson 2021

T = np.linspace(240, 310, 100)

esl1 = thermo.esatl(T)
esl2 = thermo.esatl_magnus(T)
esl3 = thermo.esatl_huang(T)

re_l1 = (esl1-esl3)/esl3 # Relative error Dales-Huang
re_l2 = (esl2-esl3)/esl3 # Relative error Magnus-Huang

Ti = np.linspace(240, 274, 100)
esi1 = thermo.esati(Ti)
esi2 = thermo.esati_magnus(Ti)
esi3 = thermo.esati_huang(Ti)

re_ice1 = (esi1-esi3)/esi3
re_ice2 = (esi2-esi3)/esi3

plt.semilogy(T,esl1, '-', label='liquid, Dales')
plt.semilogy(T,esl2, '--', label='liquid, Magnus')
plt.semilogy(T,esl3, ':', label='liquid, Huang')
plt.semilogy(Ti,esi1, '-', label='ice, Dales')
plt.semilogy(Ti,esi2, '--', label='ice, Magnus')
plt.semilogy(Ti,esi3, ':', label='ice, Huang')

plt.legend()

# relative error Magnus vs DALES
plt.figure()
plt.title ('Esat relative error compared to Huang')
plt.plot(T, re_l1, '-', label='Dales relative error liquid')
plt.plot(T, re_ice1, '--', label='Dales relative error ice')
plt.plot(T, re_l2, '-', label='Magnus relative error liquid')
plt.plot(T, re_ice2, '--', label='Magnus relative error ice')

plt.legend()

# qsatur(T)   Magnus vs Dales
plt.figure()
T = np.linspace(240, 310, 100)
pres = 101300
qsat = thermo.qsatur(T, pres)
qsat_magnus = thermo.qsatur_magnus(T, pres)
plt.semilogy(T, qsat, '-', label='qsat Dales 1013 mB')
plt.semilogy(T, qsat_magnus, '--', label='qsat Magnus 1013 mB')

pres = 70000
qsat = thermo.qsatur(T, pres)
qsat_magnus = thermo.qsatur_magnus(T, pres)
plt.semilogy(T, qsat, '-', label='qsat Dales 700 mB')
plt.semilogy(T, qsat_magnus, '--', label='qsat Magnus 700 mB')

plt.xlabel('T (K)')
plt.legend()

plt.figure()
re = (qsat_magnus - qsat) / qsat
plt.plot(T, re, '-', label='rel error qsat 700 mB')

T = np.linspace(240, 310, 100)
style=['-', '-', ':']
for p in (101300, 70000):
    plt.figure()
    plt.title(f'{p/100} mB')
    for lw, ql in ((1,0.0001), (2.5,0.001)):
        for niter in (1, 2):
            for col,name,qsatur_fun in (('red', 'Dales', thermo.qsatur), 
	                                ('purple', 'Magnus',thermo.qsatur_magnus),
                                        ('blue', 'Huang', thermo.qsatur_huang)):
                pres = np.ones_like(T) * p
                qsat = thermo.qsatur(T, pres)
                qt = qsat + ql
                thl = T/thermo.exnf(pres) - (thermo.rlv/(thermo.cp*thermo.exnf(pres))) * ql
                qsat_m = thermo.qsat_thl_qt(thl, qt, pres, qsatur_fun=qsatur_fun, niter=niter)
                ql_m = qt - qsat_m
                plt.plot(T, ql_m/ql, style[niter], color=col, lw=lw, 
			label=f'{name} {p/100} mB {ql*1000} g/kg n={niter}')
               

    plt.legend()

plt.show()



