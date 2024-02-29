from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
import random as rd
import math
from LinRegUncertainty import *
from Expressions import *

# markers - linestyle - color

def tv1():
    print("--- Teilversuch 1 ---")

    xg = 0.1525
    deltaXG = 0.0005
    xb = [1.4700, 1.4200, 1.3700, 1.3200, 1.2700]
    deltaXB = 0.0005
    xsGross = [0.3960, 0.3985, 0.4015, 0.4050, 0.4100]
    xsKlein = [1.2245, 1.1705, 1.1170, 1.0635, 1.0100]
    deltaXS = 0.0005

    print("Kleiner Max:")
    print("g\t∆g\tb\t∆b\tf\t∆f")
    for i in range(len(xb)):
        varXG = Var(xg, deltaXG, "x_g")
        varXB = Var(xb[i], deltaXB, "x_b")
        varXSK = Var(xsKlein[i], deltaXS, "x_sK")
        params = [varXG, varXB, varXSK]

        exprG = Sub(varXSK, varXG)
        exprB = Sub(varXB, varXSK)
        exprF = Pow(Add(Pow(exprG, -1), Pow(exprB, -1)), -1)

        valG = exprG.eval()
        dg = gaussian(exprG, params)

        valB = exprB.eval()
        db = gaussian(exprB, params)

        valF = exprF.eval()
        df = gaussian(exprF, params)

        print("%.4f\t%.5f\t%.4f\t%.5f\t%.4f\t%.5f" % (valG, dg, valB, db, valF, df))

    print()

    g = []
    dg = []
    b = []
    db = []
    a = []
    da = []
    e = []
    de = []

    print("Großer Max:")
    print("g\t∆g\tb\t∆b\tf\t∆f")
    fsum = 0
    for i in range(len(xb)):
        varXG = Var(xg, deltaXG, "x_g")
        varXB = Var(xb[i], deltaXB, "x_b")
        varXS = Var(xsGross[i], deltaXS, "x_s")
        varXSK = Var(xsKlein[i], deltaXS, "x_sK")
        params = [varXG, varXB, varXS, varXSK]

        exprG = Sub(varXS, varXG)
        exprB = Sub(varXB, varXS)
        exprA = Sub(varXB, varXG)
        exprE = Sub(varXSK, varXS)
        exprF = Pow(Add(Pow(exprG, -1), Pow(exprB, -1)), -1)

        valG = exprG.eval()
        g.append(valG)
        dg.append(gaussian(exprG, params))

        valB = exprB.eval()
        b.append(valB)
        db.append(gaussian(exprB, params))

        valA = exprA.eval()
        a.append(valA)
        da.append(gaussian(exprA, params))

        valE = exprE.eval()
        e.append(valE)
        de.append(gaussian(exprE, params))

        valF = exprF.eval()
        df = gaussian(exprF, params)

        f = 1/(1/valG + 1/valB)
        fsum += f
        print("%.4f\t%.5f\t%.4f\t%.5f\t%.4f\t%.5f" % (valG, dg[i], valB, db[i], valF, df))

    print()
    favg = fsum/5
    print(f"f = {favg}")
    print(f"x_S = {xg + 2 * favg}")
    print(f"x_B' = {xg + 4 * favg}")
    print()

    createGraph(g, dg, b, db, True)
    createGraph(g, dg, b, db, False)

    print("Bessel-Verfahren:")
    h = Const(0)
    fParams = []
    fSum = None
    for i in range(len(a)):
        varA = Var(a[i], da[i], "a")
        varE = Var(e[i], de[i], "e")
        params = [varA, varE]

        ah = Sub(varA, h)
        num = Sub(Pow(ah, 2), Pow(varE, 2))
        den = Mult(Const(4), ah)
        f = Div(num, den)

        val = f.eval()
        print(f"f_{i+1}  = {val}")
        unc = gaussian(f, params)
        print(f"∆f_{i+1} = {unc}")

        newF = Var(val, unc, "f")
        fParams.append(newF)
        if fSum == None:
            fSum = newF
        else:
            fSum = Add(newF, fSum)

    f = Div(fSum, Const(len(a)))
    print("Durschnitt:")
    printResult("f", f, fParams)

    print()
    print("Matrizenverfahren:")

    gG = 0.0360
    deltaGG = 0.0005
    bB = -0.150
    deltaBB = 0.001

    print(f"b[0] = {b[0]}")

    varGG = Var(gG, deltaGG, "G")
    varBB = Var(bB, deltaBB, "B")
    varB = Var(b[0], db[0], "b")
    params = [varGG, varB, varBB]

    f = Div(varB, Sub(Const(1), Div(varBB, varGG)))

    printResult("f", f, params)

def createGraph(g, deltaG, b, deltaB, zoom):
    suffix = "Zoom" if zoom else "Normal"
    pp = PdfPages(f"GraphTV1_{suffix}.pdf")

    '''
    plt.figure()
    plt.clf()
    '''
    fig, ax = plt.subplots(figsize=(5, 5), tight_layout=True)

    for i in range(len(g)):
        label = "g = %.1fcm, b = %.1fcm"
        label = label %(g[i] * 100, b[i] * 100)
        ax.plot([0, g[i]], [b[i], 0], "-", label=label)

    if zoom:
        mid = 0.20
        dx = 0.006
        dy = 0.03

        xmin = mid - dx
        xmax = mid + dx
        ymin = mid - dy
        ymax = mid + dy
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

    # plt.plot(xs, ys, "-", label=f"Messwerte")

    plt.title("Bildweite vs. Gegenstandsweite")
    plt.xlabel('Gegenstandsweite (m)')
    plt.ylabel("Bildweite (m)")
    plt.legend()

    pp.savefig()
    pp.close()

# Berechnen der Brennweite aus g und b
def eval_1(g, dg, b, db):
    varG = Var(g, dg, "g")
    varB = Var(b, db, "b")
    params = [varG, varB]
    f = Pow(Add(Pow(varG, -1), Pow(varB, -1)), -1)

    print("f aus g (%.3f) und b (%.3f):" % (g, b))

    print("f = {f.eval()}")
    unc = gaussian(f, params)
    print("∆f = {unc}")

def tv2():
    print()
    print("--- Teilversuch 2 ---")

    # direkt gemessen (Positionen, keine Abstände!)
    xs = Var(0.5490, 0.0005, "x_S")
    xz = Var(0.6695, 0.0005, "x_Z")
    xg = Var(0.1525, 0.0005, "g")
    xb = Var(1.3190, 0.0005, "b")
    gG = Var(0.0300, 0.001, "G")
    bB = Var(0.069, 0.001, "B")
    
    # aus Konstruktion (Positionen, keine Abstände!)
    fs1 = Var(0.35, 0.004, "f_S")
    fs2 = Var(0.748, 0.004, "f_S'")
    fz1 = Var(1.154, 0.004, "f_Z") # Reihenfolge! (negative Brennweite)
    fz2 = Var(0.166, 0.004, "f_Z'")
    f1 = Var(0.256, 0.004, "f")
    f2 = Var(0.764, 0.004, "f'")
    h = Var(0.494, 0.004, "H")
    h_ = Var(0.524, 0.004, "H'")

    params = [xs, xz, xg, xb, gG, bB, fs1, fs2, fz1, fz2, f1, f2, h, h_]

    fs = Div(Sub(fs2, fs1), Const(2))
    fz = Div(Sub(fz2, fz1), Const(2))
    f = Div(Sub(f2, f1), Const(2))

    g = Sub(h, xg)
    b = Sub(xb, h_)
    hAbstand = Sub(h_, h)

    printResult("f_S", fs, params)
    printResult("f_Z", fz, params)
    printResult("f", f, params)
    printResult("g", g, params)
    printResult("b", b, params)
    printResult("h", hAbstand, params)

    left = Add(Pow(g, -1), Pow(b, -1))
    right = Pow(f, -1)
    evalEquation(1, "1/g + 1/b = 1/f", left, right, params)

    left = Div(bB, gG)
    right = Div(b, g)
    evalEquation(2, "B/G = b/g", left, right, params)

    left = Sub(h_, h)
    d = Sub(xz, xs)
    num = Mult(Const(-1), Pow(d, 2))
    den = Sub(Add(fs, fz), d)
    right = Div(num, den)
    evalEquation(9, "h = -d^2 / (f1 + f2 - d)", left, right, params)

    left = Pow(f, -1)
    right = Sub(Add(Pow(fs, -1), Pow(fz, -1)), Div(d, Mult(fs, fz)))
    evalEquation(10, "1/f = 1/f1 + 1/f2 - d/(f1f2)", left, right, params)

def printResult(name, expr, params):
    print(f"{name}  = {expr.eval()}")
    unc = gaussian(expr, params)
    print(f"∆{name} = {unc}")
    print()

def evalEquation(index, equation, left, right, params):
    print(f"Gleichung ({index}): {equation}")
    leftVal = left.eval()
    leftUnc = gaussian(left, params)
    print(f"Links  = {leftVal}")
    print(f"∆Links = {leftUnc}")
    rightVal = right.eval()
    rightUnc = gaussian(right, params)
    print(f"Rechts  = {rightVal}")
    print(f"∆Rechts = {rightUnc}")
    difference = abs(leftVal - rightVal)
    unc = leftUnc + rightUnc
    print(f"Differenz = {difference}")
    factor = difference / unc
    print(f"          = {factor} * (∆Link + ∆Rechts)")
    print()

def tv3():
    print()
    print("--- Teilversuch 3 ---")

    print("Sammellinse:")
    xg = Var(0.1525, 0.0005, "x_G")
    xs = Var(0.3520, 0.0005, "x_S")
    params = [xg, xs]
    f = Sub(xs, xg)
    printResult("f_S", f, params)

    print("Zerstreulinse #1")
    xschirm = Var(1.4205, 0.0005, "x_Schirm")
    xz = Var(0.9420, 0.0005, "x_Z")
    params = [xschirm, xz]
    f = Sub(xz, xschirm)
    printResult("f_Z", f, params)

    print("Zerstreulinse #2")
    xschirm = Var(1.0165, 0.0005, "x_Schirm")
    xz = Var(0.5250, 0.0005, "x_Z")
    params = [xz, xschirm]
    f = Sub(xz, xschirm)
    printResult("f_Z", f, params)


tv1()
tv2()
tv3()
