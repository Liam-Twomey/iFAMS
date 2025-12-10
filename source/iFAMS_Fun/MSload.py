import numpy as np
import os
from PyQt5 import QtWidgets, QtGui, QtCore

def load_file_norm(self):
    """
    Initiates the load dialog, and names two objects, x and y,
    which are arrays that hold m/z (x) and abundance(y) values.
    A third object, "namebase" is file naming purposes in later
    save def's
    :param self:
    :return: x,y,namebase
    """
    name = QtWidgets.QFileDialog.getOpenFileName(self,'Open File')
    namestr = str(name[0])
    if namestr == '':
        print('No file loaded')
        return
    namebase = os.path.splitext(namestr)[0]
    print(namebase)
    try:
        if namestr.endswith('.csv'):
            try:
                x, y = np.loadtxt(namestr,unpack=True,delimiter=',')
            except ValueError:
                try:
                    index, x, y = np.loadtxt(namestr,unpack=True,delimiter=',')
                except ValueError:
                    data = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
                    if data[0][0] == 'SPECTRUM - MS':
                        x = np.array(data[0][8:len(data[0])],dtype=float)
                        y = np.array(data[1][8:len(data[1])],dtype=float)
                    else:
                        x = np.array(data[0][1:len(data[0])],dtype=float)
                        y = np.array(data[1][1:len(data[1])],dtype=float)
    
        elif namestr.endswith('.txt'):
            x, y = np.loadtxt(namestr,unpack=True)
        return x, y,namebase
    except UnboundLocalError:
        return 'wrong type'

def datapro(x,y):
    """
    A simple function that prepares the x and y for Fourier transform
    1. interpolates that data
    2. zero pads the data
    :param x: m/z values
    :param y: abundance values
    :return: xint,yint,ypadd
    """
    xint = np.linspace(min(x),max(x),len(x))
    yint = np.interp(xint,x,y)
    for i in range(len(yint)):
        po2 = 2**i
        if po2 >len(y):
            break
    zeros = np.zeros(po2-len(y))
    ypadd = np.append(yint,zeros)
    spacing = (xint[1]-xint[0])
    xend = po2-len(xint)
    xpadd = np.linspace(xint[0],xint[-1]+spacing*xend,po2)
    return xint,yint,ypadd,po2,xpadd

def batch_load(name):
    namestr = str(name)
    namebase = os.path.splitext(namestr)[0]

    if namestr.endswith('.csv'):
        try:
            x, y = np.loadtxt(namestr,unpack=True,delimiter=',')
        except ValueError:
            try:
                index, x, y = np.loadtxt(namestr,unpack=True,delimiter=',')
            except ValueError:
                try:
                    data = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
                    x = np.array(data[0][1:len(data[0])],dtype=float)
                    y = np.array(data[1][1:len(data[1])],dtype=float)
                except ValueError:
                    print('Unable to load. Data format needs to be x vs y.')
    elif namestr.endswith('.txt'):
        x, y = np.loadtxt(namestr,unpack=True)
    for i in range(0, len(x)):
        x[i] = round(x[i],7)
        y[i] = round(y[i],7)

    return x,y,namebase,namestr

def load_decon(name):
    #try:

    namestr = str(name)
    namebase = os.path.splitext(namestr)[0]

    if namestr.endswith('.csv'):
        data = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
        if data[0][0] == '0' and data[0][1] == '1':
            x = np.array(data[1][0:len(data[0])],dtype=float)
            y = np.array(data[2][0:len(data[1])],dtype=float)
            bounds = 0
            based = False
            cs = []
            cszero = []
            domain = 'Mass (Da)'
        else:
            x = np.array(data[0][1:len(data[0])],dtype=float)
            y = np.array(data[1][1:len(data[1])],dtype=float)
            domain = data[0][0]
            try:
                boundstr = np.array(data[2][1:len(data[2])],dtype=float)
                bounds = []
                boundstr = boundstr.tolist()
                if int(boundstr[0]) == 1:
                    based = True
                    boundstr.pop(0)
                elif int(boundstr[0]) == 0:
                    based = False
                    boundstr.pop(0)
                else:
                    based = False
                for i in range(0,len(boundstr)):
                    if boundstr[i] > 0:
                        bounds.append(round(boundstr[i],4))
                    else:
                        break
            except IndexError:
                bounds = 0
                based = False
            try:
                cs = []
                cszero = []
                for i in range(3,len(data)):
                    cs.append(int(float(data[i][0])))
                    cszero.append(np.array(data[i][1:len(data[i])],dtype=float))
            except IndexError:
                pass
        
    else:
        x, y = np.loadtxt(namestr,unpack=True)
        for i in range(0, len(x)):
            x[i] = round(x[i],7)
            y[i] = round(y[i],7)
            bounds = 0
            based = False
            cs = []
            cszero = []
            domain = 'Mass (Da)'
            
    return x, y,namebase,bounds,based,cs,cszero,domain

def rec(rlist,cs):
    rectangles = []
    for i in range(0,len(cs)):
        rectemp = []
        for j in range(len(rlist)):
            if rlist[j][4] == i:
                rectemp.append(rlist[j])
            else:
                continue
        rectangles.append(rectemp)
    return rectangles

def agilent(self, pathAgilent):
    namestr = pathAgilent
    print(namestr)
    namebase = pathAgilent.rsplit('.', 1)[0]
    print(namebase)
    try:
        if namestr.endswith('.csv'):
            try:
                x, y = np.loadtxt(namestr,unpack=True,delimiter=',')
            except ValueError:
                index, x, y = np.loadtxt(namestr,unpack=True,delimiter=',')

        elif namestr.endswith('.txt'):
            x, y = np.loadtxt(namestr,unpack=True)
        return x, y,namebase

    except UnboundLocalError:
        print('an exception was found')

