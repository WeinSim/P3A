from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
import random as rd
import math
from LinRegUncertainty import *
from Expressions import *

wavelengths6 = [578.966, 576.960, 546.074, 491.604, 435.405, 404.656]
wavelengths7 = [632.44, 578.966, 576.960, 546.074, 491.604, 435.405, 404.656]

angleUncertainty = 1.0 / 60.0

phi1Deg1 = [34, 34, 34, 33, 32, 31]
phi1Min1 = [31, 29, 31, 29, 25, 30]
phi2Deg1 = [131, 131, 131, 132, 133, 134]
phi2Min1 = [15, 17, 34, 19, 26, 21]

phi1Deg2 = [59, 59, 59, 58, 58, 58]
phi1Min2 = [12, 12, 6, 56, 40, 30]
phi2Deg2 = [106, 106, 106, 106, 106, 106, 106]
phi2Min2 = [10, 10, 19, 30, 41, 53]

phi3Deg = [60, 61, 61, 63, 65, 67, 68]
phi3Min = [30, 41, 47, 26, 27, 29, 34]
phi4Deg = [104, 103, 103, 102, 100, 98, 97]
phi4Min = [59, 21, 17, 7, 8, 10, 0]

def tv2():
    print()
    print("--- Teilversuch 2 ---")

    evaluateTV2("Flintglas", phi1Deg1, phi1Min1, phi2Deg1, phi2Min1)

    evaluateTV2("Methanol", phi1Deg2, phi1Min2, phi2Deg2, phi2Min2)

def evaluateTV2(name, phi1Deg, phi1Min, phi2Deg, phi2Min):
    ior = getIORefractions(phi1Deg, phi1Min, phi2Deg, phi2Min)

    print(f"{name}:")
    fmt = "lambda = %.3f nm: n = %.5f, ∆n = %.5f"
    for i in range(6):
        print(fmt % (wavelengths6[i], ior[0][i], ior[1][i]))

    pp = PdfPages(f"GraphTV2_{name}.pdf")

    plt.figure()
    plt.clf()

    plt.errorbar(wavelengths6, ior[0], ior[1], fmt="-o", label=f"Messwerte {name}")
    # plt.plot(wavelengths6, ior[0], "-o", label=f"Messwerte {name}")

    plt.title(f"Brechungsindex vs. Wellenlänge")
    plt.xlabel('Wellenlänge (nm)')
    plt.ylabel('Brechungsindex')
    plt.legend()

    pp.savefig()
    pp.close()

def getIORefractions(phi1Deg, phi1Min, phi2Deg, phi2Min):
    ioRefractions = []
    uncertainties = []
    theta = Const(60 * math.pi / 180.0)
    for i in range(len(phi1Deg)):
        phi1 = phi1Deg[i] + phi1Min[i] / 60.0
        varphi1 = Var(phi1, angleUncertainty, "phi1")
        phi2 = phi2Deg[i] + phi2Min[i] / 60.0
        varphi2 = Var(phi2, angleUncertainty, "phi2")
        params = [varphi1, varphi2]

        deltaDeg = Div(Sub(varphi2, varphi1) , Const(2))
        delta = Mult(deltaDeg, Const(math.pi / 180.0))
        n = Div(Sin(Div(Add(theta, delta), Const(2.0))), Sin(Div(theta, Const(2.0))))
        ioRefractions.append(n.eval())
        uncertainties.append(gaussian(n, params))
    return (ioRefractions, uncertainties)

def tv3():
    print()
    print("--- Teilversuch 3 ---")

    for i in range(7):
        phi3 = phi3Deg[i] + phi3Min[i] / 60.0
        phi4 = phi4Deg[i] + phi4Min[i] / 60.0

        varphi3 = Var(phi3, angleUncertainty, "phi3")
        varphi4 = Var(phi4, angleUncertainty, "phi4")
        phiDeg = Div(Sub(varphi4, varphi3) , Const(2))
        phi = Mult(phiDeg, Const(math.pi / 180.0))
        g = Const(600 * 1000)
        wl = Div(Sin(phi), g)

        wlNM = wl.eval() * 1e9
        print(f"lambda_{i} = {wlNM} nm")
        params = [varphi3, varphi4]
        unc = gaussian(wl, params) * 1e9
        print(f"∆lambda_{i} = {unc} nm")

def tv4():
    print()
    print("--- Teilversuch 4 ---")

    print("Flintglas:")
    ior = getIORefractions(phi1Deg1, phi1Min1, phi2Deg1, phi2Min1)
    b = 29.8 * 1e-3
    evaluateTV4(ior[0], ior[1], 1, b)
    evaluateTV4(ior[0], ior[1], 4, b)

    print("Methanol:")
    b = 68.4 * 1e-3
    ior = getIORefractions(phi1Deg2, phi1Min2, phi2Deg2, phi2Min2)
    evaluateTV4(ior[0], ior[1], 1, b)
    evaluateTV4(ior[0], ior[1], 4, b)

def evaluateTV4(iorValues, iorUncertainties, wlIndex, b):
    wl = wavelengths6[wlIndex]
    print(f"Auflösungsvermögen bei {wl} nm:")

    ns = []
    wls = []
    for i in range(2):
        ns.append(iorValues[wlIndex + i])
        wls.append(wavelengths6[wlIndex + i] * 1e-9)
    baseLen = Var(b, 0.00005, "b")
    dndl = abs((ns[1] - ns[0]) / (wls[1] - wls[0]))

    r = Mult(baseLen, Const(dndl))

    print(f"R = {r.eval()}")
    unc = gaussian(r, [baseLen])
    print(f"∆R = {unc}")

def tv5():
    print()
    print("--- Teilversuch 5 ---")

    phi5Deg = [84, 85, 87]
    phi5Min = [40, 17, 31]
    phi6Deg = [81, 79, 77]
    phi6Min = [10, 30, 51]
    evaluateTV5(50, "Gelb", phi5Deg, phi5Min, phi6Deg, phi6Min)

    phi5Deg = [84, 85, 86]
    phi5Min = [8, 28, 45]
    phi6Deg = [81, 80, 79]
    phi6Min = [35, 21, 3]
    evaluateTV5(50, "Blau", phi5Deg, phi5Min, phi6Deg, phi6Min)

    phi5Deg = [83, 83]
    phi5Min = [12, 30]
    phi6Deg = [82, 82]
    phi6Min = [30, 12]
    evaluateTV5(10, "Gelb", phi5Deg, phi5Min, phi6Deg, phi6Min)

    phi5Deg = [83, 83]
    phi5Min = [8, 23]
    phi6Deg = [82, 82]
    phi6Min = [38, 24]
    evaluateTV5(10, "Blau", phi5Deg, phi5Min, phi6Deg, phi6Min)

    d1 = [2.95, 2.7, 3.2]
    d2 = [1.55, 1.35, 1.5]
    d = (d1, d2)
    evaluateTV5_2(10, d)

    d1 = [0.7, 0.75, 0.8]
    d2 = [0.4, 0.5, 0.3]
    d3 = [0.25, 0.3, 0.3]
    d = (d1, d2, d3)
    evaluateTV5_2(50, d)


def evaluateTV5(gitterKonst, color, phi5Deg, phi5Min, phi6Deg, phi6Min):
    print(f"{gitterKonst}er-Gitter {color}:")
    mmax = len(phi5Deg)
    g = Const(1.0 / (gitterKonst * 1e3))
    for i in range(mmax):
        phi5 = phi5Deg[i] + phi5Min[i] / 60.0
        varphi5 = Var(phi5, angleUncertainty, "phi5")
        phi6 = phi6Deg[i] + phi6Min[i] / 60.0
        varphi6 = Var(phi6, angleUncertainty, "phi6")
        params = [varphi5, varphi6]
        m = Const(i + 1)

        betaDeg = Div(Sub(varphi5, varphi6), Const(2.0))
        beta = Mult(betaDeg, Const(math.pi / 180.0))
        dpdl = Div(m, Mult(g, Cos(beta)))

        print(f"{i + 1}. Ordnung:")
        print(f"d phi / d lambda = {dpdl.eval()}")
        unc = gaussian(dpdl, params)
        print(f"∆ d phi / d lambda = {unc}")

    print()

def evaluateTV5_2(gitterKonst, d):
    print(f"{gitterKonst}er-Gitter Auflösungsvermögen:")
    mmax = len(d)
    g = Const(gitterKonst * 1e3)
    for i in range(mmax):
        total = None
        params = []
        for j in range(3):
            vard = Var(d[i][j] * 1e-3, 0.00005, "d")
            params.append(vard)
            a = Mult(Mult(Const(i + 1), vard), g)
            if total == None:
                total = a
            else:
                total = Add(total, a)
        total = Div(total, Const(3))
        print(f"A = {total.eval()}")
        unc = gaussian(total, params)
        print(f"∆A = {unc}")
    print()




tv2()
tv3()
tv4()
tv5()
