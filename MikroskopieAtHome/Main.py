from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
import random as rd
import math
from LinRegUncertainty import *
from Expressions import *

# markers - linestyle - color

def tv3():
    print("--- Teilversuch 3 ---")
    s0 = Var(0.25, 0, "s_0")
    n = Var(1.5168, 0, "n")
    d = Var(2.38e-3, 0, "D")

    f = Div(Mult(n, d), Mult(Const(4), Sub(n, Const(1))))
    print(f"f = {f.eval()}")

    gamma = Div(s0, f)

    print("Theoretisch:")
    theo = gamma.eval()
    print(f"gamma = {theo}")

    b = Var(0.025, 0.003, "B")
    g = Div(b, gamma)
    params = [b]
    print("Durchmesser Haar (um):")
    printExpr("G", g.eval() * 1e6, gaussian(g, params) * 1e6)

    g = Var(75e-6, 50e-6, "G")
    params = [b, g]
    gamma = Div(b, g)
    print("Experimentell:")
    val = gamma.eval()
    unc = gaussian(gamma, params)
    printExpr("gamma", val, unc)
    printDiff(val, theo, unc)



def printExpr(name, value, unc):
    print(f"{name}  = {value}")
    print(f"∆{name} = {unc}")

def printVar(var):
    printExpr(var.name, var.value, var.uncertainty)

def printDiff(val, theo, unc):
    u = abs(val - theo) / unc
    print(f"Abweichung vom theoretischen Wert: {u} * Unsicherheit")

tv3()
