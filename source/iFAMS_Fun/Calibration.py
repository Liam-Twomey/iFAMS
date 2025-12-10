#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jul 21, 2022

@author: Kayd Meldrum
"""

import numpy as np
from scipy import stats as stats
import os
from PyQt5 import QtWidgets, QtGui, QtCore
import scipy.optimize as opt

def load_peaklist(name):
    """
    Initiates the load dialog, and names two objects, x and y,
    which are arrays that hold mass (x) and area (y) values.
    A third object, "namebase" is file naming purposes in later
    save def's
    :param self:
    :return: x,y,namebase
    """
    namebase = os.path.splitext(name)[0]
    try:
        x, y, z, w = np.loadtxt(name,unpack=True,delimiter=',',dtype=str)
        if x[0].lower() == 'centroid':
            x = np.delete(x,0)
            y = np.delete(y,0)
            z = np.delete(z,0)
            w = np.delete(w,0)
        if x[-1].lower() == 'noise:':
            x[-1] = 0.0
        x = np.round(np.array(x,dtype=float),7)
        y = np.round(np.array(y,dtype=float),7)
        z = np.round(np.array(z,dtype=float),7)
        w = np.round(np.array(w,dtype=float),7)

    except ValueError:
        x, y, z = np.loadtxt(name,unpack=True,delimiter=',',dtype=str)
        if x[0].lower() == 'centroid':
            x = np.delete(x,0)
            y = np.delete(y,0)
            z = np.delete(z,0)
        if x[-1].lower() == 'noise:':
            x[-1] = 0.0
        x = np.round(np.array(x,dtype=float),7)
        y = np.round(np.array(y,dtype=float),7)
        z = np.round(np.array(z,dtype=float),7)
        w = (np.array(x,dtype=float))*0+1
        
    except UnboundLocalError:
        print('an exception was found while loading peaklists')
        return
    
    return x, y, z, w, namebase        
        
def sample(x,y,w,sample,tolerance,tunit,num):
    """
    searches through the .txt or .csv file to locate the inputted sample value
    and returns the associated integration value. 
    """
    if tunit=='Da':
        for i in range(len(x)):
            if (sample - tolerance <= x[i]) and (x[i] <= sample + tolerance):
                ysample = y[i]
                peakwidth = int(w[i])
                xcent = x[i]
                try:
                    if abs(sample-x[i+1]) < abs(sample-x[i]):
                        ysample = y[i+1]
                        peakwidth = int(w[i+1])
                        xcent = x[i+1]
                        continue
                except IndexError:
                    pass
                print ('selecting peak ' + str(i+1) + ' as calibrant for concentration ' + str(num+1))
                break
            else:
                ysample = str(0)
                peakwidth = int(0)
                xcent = 0
                #print ('No data in range')
    else:
        tolppm = float(tolerance*sample/(10**6))
        for i in range(len(x)):
            if (sample - tolppm <= x[i]) and (x[i] <= sample + tolppm):
                ysample = y[i]
                peakwidth = int(w[i])
                xcent = x[i]
                try:
                    if abs(sample-x[i+1]) < abs(sample-x[i]):
                        ysample = y[i+1]
                        peakwidth = int(w[i+1])
                        xcent = x[i+1]
                        continue
                except IndexError:
                    pass
                print ('selecting peak ' + str(i+1) + ' as calibrant for concentration ' + str(num+1))
                break
            else:
                ysample = str(0)
                peakwidth = int(0)
                xcent = 0
                #print ('No data in range')

    return ysample, peakwidth, xcent
                
def standard(x,y,standard,tolerance,tunit,totalstds,num):
    """
    searches through the .txt or .csv file to locate the inputted standard value
    and returns the associated integration value. 
    """
    if tunit=='Da':
        for i in range(len(x)):
            if (standard - tolerance <= x[i]) and (x[i] <= standard + tolerance):
                ystandard = y[i]
                xcent = x[i]
                try:
                    if abs(standard-x[i+1]) < abs(standard-x[i]):
                        ystandard = y[i+1]
                        xcent = x[i+1]
                        continue
                except IndexError:
                    pass
                print ('selecting peak ' + str(i+1) + ' as standard for concentration ' + str(num+1))
                break
            else:
                ystandard = str(1/totalstds)
                xcent = 0
                #print ('No data in range')
    else:
        tolppm = float(tolerance*standard/(10**6))
        for i in range(len(x)):
            if (standard - tolppm <= x[i]) and (x[i] <= standard + tolppm):
                ystandard = y[i]
                xcent = x[i]
                try:
                    if abs(standard-x[i+1]) < abs(standard-x[i]):
                        ystandard = y[i+1]
                        xcent = x[i+1]
                        continue
                except IndexError:
                    pass
                print ('selecting peak ' + str(i+1) + ' as standard for concentration ' + str(num+1))
                break
            else:
                ystandard = str(1/totalstds)
                xcent = 0
                #print ('No data in range')
    return ystandard,xcent
                

def calibrate(x,xrep,y,weight,curve,stddev):
    """
    Generates calibration curve and returns calibration information for plotting

    Parameters
    ----------
    x : array
        x-data to be calibrated (concentrations).
    y : array
        y-data to be calibrated (signal heights/areas).
    weight : string
        type of weighting to be used.
    curve : string
        type of curve to generate for fit.

    Returns
    -------
    fitx : array
        calibration x-data for plotting.
    fity : array
        calibration y-data for plotting.
    fitparam : array
        parameters for calibration curve.
    calinfo : string
        curve equation and R^2.
    std_err : float
        standard error about the regression.
    confidence : 2D array
        information for plotting 95% confidence bands about regression.
    uinfo : array
        information to calculate uncertainties.

    """
    x = np.array(x)
    y = np.array(y)
    
    #linear fit
    if curve == 'linear':
        if weight == 'no weighting':
            weights = None
            std_errs = 0
        if weight == '1/x':
            weights = 1/x
        if weight == '1/x^2':
            weights = (1/x)**2
        if weight == 'x':
            weights = x
        if weight == 'x^2':
            weights = x**2
        
        slope, intercept = np.polyfit(x,y,1,w=weights)
        
        yhat = slope*(x)+intercept
        ybar = np.sum(y)/len(y)
        ssres = np.sum((y-yhat)**2)
        sstot = np.sum((y-ybar)**2)
        r_value2 = 1 - ssres/sstot
        
        fitx = np.linspace(-1,max(x)+max(x)/10,100,endpoint=True)
        fity = slope*fitx + intercept
        fitparam = np.array([slope,intercept])
        calinfo = "Calibration Curve \n y = " + str(np.round(slope,7)) + "x + " + str(np.round(intercept,7)) + " \n R^2 = " + str(np.round((r_value2),7))
        print("Calibration Curve \n y = " + str(slope) + "x + " + str(intercept) + " \n R^2 = " + str(r_value2))
        
        std_err = 0
        calmean = np.average(x)
        caldem = 0
        for i in range(0,len(x)):
            std_err += (y[i]-(x[i]*slope+intercept))**2
            caldem += (x[i]-calmean)**2
        std_err = np.sqrt(std_err/(len(x)-2))
        cal_err = np.zeros(len(x))
        try:
            std_errs = std_err*weights
            for i in range(0,len(x)):
                cal_err[i] = std_errs[i]/slope * np.sqrt(1/xrep[i] + 1/len(x) + (y[i]-np.average(y))**2/(slope**2 * caldem))
                ytemp = (x[i] - cal_err[i])*slope+intercept
                cal_err[i] = abs(y[i]-ytemp)
        except TypeError:
            for i in range(0,len(x)):
                cal_err[i] = std_err/slope * np.sqrt(1/xrep[i] + 1/len(x) + (y[i]-np.average(y))**2/(slope**2 * caldem))
                ytemp = (x[i] - cal_err[i])*slope+intercept
                cal_err[i] = abs(y[i]-ytemp)
    
        Avgy = sum(y)/len(y)
        Avgx = sum(x)/len(x)
        xspread = 0
        for i in range(0, len(x)):
            xspread += (x[i]-Avgx)**2
        uinfo = np.array([Avgy,xspread])
        confy = np.linspace(fity[0]-10*std_err,max(fity)+10*std_err,100,endpoint=True)
        confx = (confy-intercept)/slope
        conf = np.zeros(len(confy))
        for i in range(0,len(confy)):
            conf[i] = 1.96*abs(std_err/slope * np.sqrt(1 + 1/len(x) + (confy[i]-np.average(y))**2/(slope**2 * caldem)))

        confR = confx + conf
        confL = confx - conf
        confx = np.linspace(fitx[0],fitx[-1],100,endpoint=True)
        confupper = np.interp(confx,confL,confy)
        conflower = np.interp(confx,confR,confy)
        confidence = np.vstack((confx,conflower,confupper))
    
    #quadratic fit
    if curve == 'quadratic':
        # if weight == 'no weighting':
        #     weights = None
        # if weight == '1/x':
        #     weights = 1/x
        # if weight == '1/x^2':
        #     weights = (1/x)**2
        # if weight == 'x':
        #     weights = x
        # if weight == 'x^2':
        #     weights = x**2
            
        #(coeA, coeB, coeC), cov = np.polyfit(x,y,2,cov='unscaled',w=weights)
        
        if weight == 'no weighting':
            weights = None
        if weight == '1/x':
            weights = x
        if weight == '1/x^2':
            weights = (x)**2
        if weight == 'x':
            weights = 1/x
        if weight == 'x^2':
            weights = 1/x**2
        
        def f(x, a, b, c):
            return a * x**2 + b * x + c
    
        #initial guesses for parameters
        pf = np.array([0,np.average(y/x),0])
    
        (coeA, coeB, coeC), cov = opt.curve_fit(f, x, y,p0=pf,sigma=weights)
        
        yhat = np.poly1d((coeA,coeB,coeC))(x)
        ybar = np.sum(y)/len(y)
        ssres = np.sum((y-yhat)**2)
        sstot = np.sum((y-ybar)**2)
        r_value2 = 1 - ssres/sstot
        std_err = np.sqrt(np.sum((y-yhat)**2)/(len(x)-3))
        
        fitx = np.linspace(min(x)-min(x)/2,max(x)+min(x)/2,200,endpoint=True)
        fity = coeA*fitx**2+coeB*fitx+coeC
        fitparam = np.array([coeA,coeB,coeC])
        print(fitparam)
        calinfo = "Calibration Curve \n y = " +str(np.round(coeA,5))+ "x^2 + " +str(np.round(coeB,5))+ "x + " +str(np.round(coeC,5))+ " \n R^2 = " +str(np.round(r_value2,7))
        print("Calibration Curve \n y = " +str(coeA)+ "x^2 + " +str(coeB)+ "x + " +str(coeC)+ " \n R^2 = " +str(r_value2))
        
        confy = np.linspace(min(fity),max(fity),200,endpoint=True)
        m = np.zeros(len(confy))+1
        confx,confu = unkcalc2(coeA,coeB,coeC,confy,std_err,cov,m)
        uinfo = cov
        confupper = confy + confu
        conflower = confy - confu
        stddev = np.interp(confx,x,stddev)
        for i in range(len(confx)):
            scov = np.average((confupper[i]-confy[i],confy[i]-conflower[i]))
            srep = stddev[i]
            confupper[i] = confy[i]+2*np.sqrt(scov**2+srep**2)
            conflower[i] = confy[i]-2*np.sqrt(scov**2+srep**2)
        confidence = np.vstack((confx,confupper,conflower))

        # confL = confx - confu
        # confR = confx + confu
        # confx = np.linspace(fitx[0],fitx[-1],200,endpoint=True)
        # confupper = np.interp(confx,confL,confy)
        # conflower = np.interp(confx,confR,confy)
        # stddev = np.interp(confx,x,stddev)
        # for i in range(len(confx)):
        #     scov = np.average((confupper[i]-confy[i],confy[i]-conflower[i]))
        #     srep = stddev[i]
        #     confupper[i] = confy[i]+2*np.sqrt(scov**2+srep**2)
        #     conflower[i] = confy[i]-2*np.sqrt(scov**2+srep**2)
        # confidence = np.vstack((confx,confupper,conflower))
        
    #logistic fit
    if curve == 'logistic':
        if weight == 'no weighting':
            weights = None
        if weight == '1/x':
            weights = x
        if weight == '1/x^2':
            weights = (x)**2
        if weight == 'x':
            weights = 1/x
        if weight == 'x^2':
            weights = 1/x**2
    
        def f(x, a, b, c):
            return a / (1. + np.exp(-b * (x))) + c
    
        try:
            #initial guesses for parameters
            pf = np.array([2*y[-1],(2/y[-1]*y[int(np.floor(len(y)/2))]/x[int(np.floor(len(y)/2))]),0])
        
            (a_, b_, c_), cov = opt.curve_fit(f, x, y,p0=pf,sigma=weights)
        
            yhat = f(x, a_, b_, c_)
            ybar = np.sum(y)/len(y)
            ssres = np.sum((y-yhat)**2)
            sstot = np.sum((y-ybar)**2)
            r_value2 = 1 - ssres/sstot
            std_err = np.sqrt(np.sum((y-yhat)**2)/(len(x)-3))
            #std_err = np.sqrt(np.diag(cov))
        
            fitx = np.linspace(0,max(x)+1,100,endpoint=True)
            fity = f(fitx, a_, b_, c_)
            fitparam = np.array([a_,b_,c_])
            calinfo = "Calibration Curve \n y = " +str(np.round(a_,5))+ " / (1 + exp(" +str(-1*np.round(b_,5))+ "x)) + " +str(np.round(c_,5))+ " \n R^2 = " +str(np.round(r_value2,7))
            print("Calibration Curve \n y = " +str(a_)+ " / (1 + exp(" +str(-1*b_)+ "x)) + " +str(c_)+ " \n R^2 = " +str(r_value2))
            
            confy = np.linspace(fity[0]-10*std_err,max(fity)+10*std_err,200,endpoint=True)
            m = np.zeros(len(confy))+1
            confx,confu = unkcalc3(a_,b_,c_,confy,std_err,cov,m)            
            uinfo = cov
            confupper = confy + confu
            conflower = confy - confu
            stddev = np.interp(confx,x,stddev)
            for i in range(len(confx)):
                scov = np.average((confupper[i]-confy[i],confy[i]-conflower[i]))
                srep = stddev[i]
                confupper[i] = confy[i]+2*np.sqrt(scov**2+srep**2)
                conflower[i] = confy[i]-2*np.sqrt(scov**2+srep**2)
            confidence = np.vstack((confx,confupper,conflower))
        except RuntimeError:
            try:        
                #initial guesses for parameters
                pf = np.array([2*y[-1],0.001,0])
                bounds = np.array(([-np.inf,0.001,-np.inf],[np.inf,0.0011,np.inf]))
                print("Optimization Error: constraining parameter 'b' to 0.001")
            
                (a_, b_, c_), cov = opt.curve_fit(f, x, y,p0=pf,sigma=weights,bounds=bounds)
            
                yhat = f(x, a_, b_, c_)
                ybar = np.sum(y)/len(y)
                ssreg = np.sum((yhat-ybar)**2)
                sstot = np.sum((y-ybar)**2)
                r_value2 = ssreg/sstot
                std_err = np.sqrt(np.sum((y-yhat)**2)/(len(x)-3))
                #std_err = np.sqrt(np.diag(cov))
            
                fitx = np.linspace(0,max(x)+1,100,endpoint=True)
                fity = f(fitx, a_, b_, c_)
                fitparam = np.array([a_,b_,c_])
                calinfo = "Calibration Curve \n y = " +str(np.round(a_,5))+ " / (1 + exp(" +str(-1*np.round(b_,5))+ "x)) + " +str(np.round(c_,5))+ " \n R^2 = " +str(np.round(r_value2,7))
               
                confy = np.linspace(fity[0]-10*std_err,max(fity)+10*std_err,200,endpoint=True)
                m = np.zeros(len(confy))+1
                confx,confu = unkcalc3(a_,b_,c_,confy,std_err,cov,m)
                uinfo = cov
                confupper = confy + confu
                conflower = confy - confu
                stddev = np.interp(confx,x,stddev)
                for i in range(len(confx)):
                    scov = np.average((confupper[i]-confy[i],confy[i]-conflower[i]))
                    srep = stddev[i]
                    confupper[i] = confy[i]+2*np.sqrt(scov**2+srep**2)
                    conflower[i] = confy[i]-2*np.sqrt(scov**2+srep**2)
                confidence = np.vstack((confx,confupper,conflower))
            except RuntimeError:
                print('unable to find optimal parameters')
                return

    return fitx, fity, fitparam, calinfo, std_err, confidence, uinfo

        
def unkcalc(slope,intercept,unky,std_err,Avgsignal,concspread,n,m):
    """
    takes the calculated slope and intercept, and the area value calculated to find
    the concentration of the unknown sample
    """
    unkconc = (unky-intercept) / slope
    
    # m = number of replicates averaged for the unknown signal
    # n = number of calibrants used for Avgsignal (average across calibration stds)
    errbarsU = 1.96*(std_err/slope * np.sqrt(1/m + 1/n + (unky-Avgsignal)**2/(slope**2 * concspread)))
    
    return unkconc,errbarsU

def unkcalc2(a,b,c,unky,std_err,cov,m):
    """
    calculates unknown concentration for a quadratic fit
    """
    unkconc2 = (-b+np.sqrt((b)**2-4*a*(c-unky)))/(2*a)
    
    sA = np.sqrt(cov[0][0])
    sB = np.sqrt(cov[1][1])
    sC = np.sqrt(cov[2][2])
    sAB = cov[0][1]
    sAC = cov[0][2]
    sBC = cov[1][2]
    
    # dXdY = 1/np.sqrt(coeB**2-4*coeA*(coeC-unky))
    # dXdC = -1/np.sqrt(coeB**2-4*coeA*(coeC-unky))
    # dXdB = (-1 + coeB/np.sqrt(coeB**2-4*coeA*(coeC-unky)))/(2*coeA)
    # dXdA = (-coeC+unky)/(coeA*np.sqrt(coeB**2-4*coeA*(coeC-unky)))-(-coeB+np.sqrt(coeB**2-4*coeA*(coeC-unky)))/(2*coeA**2)
    # u2 = dXdY**2*std_err**2 + dXdA**2*sA**2+dXdB**2*sB**2 + dXdC**2*sC**2
    # co_u2 = u2 + 2*dXdA*dXdB*sAB + 2*dXdA*dXdC*sAC + 2*dXdB*dXdC*sBC
    
    # errbarsU = 1.96*np.sqrt(abs(co_u2)/m)
       
    def J(x, a, b, c):
        Ja = x**2
        Jb = x
        Jc = 1
        return np.array([[Ja,Jb,Jc]])

    sigma = np.zeros(len(unkconc2))
    for i in range(len(unkconc2)):
        JM = J(unkconc2[i],a,b,c)
        P1 = np.matmul(JM,cov)*3/m[i] 
        # 3 because of the number of parameters, 
        # m is the number of replicates averaged for the unknown signal
        sigma[i] = np.sqrt(np.absolute(np.matmul(P1,JM.T)))
    
    errbarsU = 1.96*sigma
    
    return unkconc2, errbarsU

def unkcalc3(a,b,c,unky,std_err,cov,m):
    """
    calculates unknown concentration for a logistic fit
    """
    unkconc = -np.log(a/(unky-c)-1)/b
    
    def J(x, a, b, c):
        Ja = 1/(1+np.exp(-b*x))
        Jb = a*x*np.exp(-b*x)/(1+np.exp(-b*x))**2
        Jc = 1
        return np.array([[Ja,Jb,Jc]])
    
    sigma = np.zeros(len(unkconc))
    for i in range(len(unkconc)):
        JM = J(unkconc[i],a,b,c)
        P1 = np.matmul(JM,cov)*3/m[i] 
        # 3 because of the number of parameters, 
        # m is the number of replicates averaged for the unknown signal
        sigma[i] = np.sqrt(np.absolute(np.matmul(P1,JM.T)))
    
    errbarsU = 1.96*sigma
    
    
    return unkconc, errbarsU

def load_file_cal(self):
    """
    Initiates the load dialog, and names two objects, x and y,
    which are arrays that hold concentration (x) and signal (y) values.
    A third object, "namebase" is file naming purposes in later
    save def's. SN column added as well as str column with 
    information on signal type, curvetype, and concentration units. 
    """
    name = QtWidgets.QFileDialog.getOpenFileName(self,'Open File')
    namestr = str(name[0])

    calibrants = []
    std = []
    try:
        concstr,concydatastr,SNstr = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
        conc = []
        concydata = []
        SN = []
        for i in range(2,len(concstr)):
            if concstr[i].lower() == 'calibrant masses:' or concstr[i].lower() == 'calibrant masses':
                calIndex = i+1
                break
            else:
                conc.append(round(float(concstr[i]),4))
                concydata.append(round(float(concydatastr[i]),7))
                try:
                    SN.append(round(float(SNstr[i]),2))
                except ValueError:
                    SN.append(0)
        for i in range(calIndex,len(concstr)):
            if int(float(concstr[i])) == 0:
                continue
            else:
                calibrants.append(float(concstr[i]))
            if int(float(concydatastr[i])) == 0:
                continue
            else:
                std.append(float(concydatastr[i]))
        sigtype = concydatastr[1].lower()
        curtype = SNstr[calIndex+3].lower()
        try:
            weight = SNstr[calIndex+5].lower()
        except IndexError:
            weight = 'no weighting'
        except AttributeError:
            weight = 'no weighting'
        try:
            SNstr[calIndex+6].lower()
            Avg = True
        except IndexError:
            Avg = False
        except AttributeError:
            Avg = False
        Calunits = concstr[1]
        caltol = []
        caltol.append('tolerance')
        caltol.append(SNstr[calIndex+1])
        caltol.append(SNstr[calIndex])
        
    except ValueError:
        weight = 'no weighting'
        Avg = False
        try:    
            concstr,concydatastr,SNstr,caltype = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
            caltol = [0,'Da',20]
        except ValueError:
            concstr,concydatastr,SNstr,caltype,caltol,calmasses,stdmasses = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
            concstr = np.delete(concstr,0)
            concydatastr = np.delete(concydatastr,0)
            SNstr = np.delete(SNstr,0)
            for i in range(1,len(calmasses)):
                if float(calmasses[i]) == float(0.0):
                    break
                calibrants.append(float(calmasses[i]))
            for i in range(1,len(stdmasses)):
                if float(stdmasses[i]) == float(0.0):
                    break
                std.append(float(stdmasses[i]))
        conc = []
        concydata = np.zeros(len(concydatastr))
        SN = np.zeros(len(concydatastr)).tolist()
        for i in range(len(concstr)):
            conc.append(round(float(concstr[i]),4))
            concydata[i] = round(float(concydatastr[i]),7)
            try:
                SN[i] = str(round(float(SNstr[i]),2))
            except ValueError:
                SN[i] = 0
        
        sigtype = caltype[0].lower()
        curtype = caltype[1].lower()
        try:
            if round(float(caltype[2]),0) == 0:
                Calunits = '-'
            if round(float(caltype[0]),0) == 1:
                sigtype = 'integration'
            if round(float(caltype[1]),0) == 1:
                curtype = 'linear'
        except ValueError or TypeError:
            Calunits = caltype[2]
         
    return conc, concydata, namestr, SN, sigtype, curtype, Calunits, calibrants, std, caltol, weight, Avg

def noise_calc(xFT,yFT,x1,x2):
    """
    calculates average noise, standard deviation of noise, and the LOQ 
    when provided a box of gabor data
    """
    xFTspacing = xFT[2]-xFT[1]
    x1i = int(np.round(x1/xFTspacing))
    x2i = int(np.round(x2/xFTspacing))
    noiselist = np.array(yFT[x1i:x2i+1])
    Frms = float(np.sqrt(np.mean((np.absolute(noiselist))**2)))
    
    return Frms

def sampleNoise(x,y,w,num):
    """
    searches through the .txt or .csv file to locate the inputted sample noise value
    at mass 0 and returns the associated integration value corrected for total 
    integration width. 
    """
    noise = int(1)
    for i in range(len(x)):
        if x[i] == 0:
            noise = float(y[i]/np.sqrt(w[i]))
            break
        else:
            continue

    if noise == 1 or noise == 0:
        print ('No noise value available for calibrant/unknown ' +str(num+1))
        noise = 1

    return noise

def sampleNoiseH(x,z,num):
    """
    searches through the .txt or .csv file to locate the inputted sample noise value
    at mass 0 and returns the associated height value. 
    """
    noise = int(1)
    for i in range(len(x)):
        if x[i] == 0:
            noise = z[i]
            break
        else:
            continue

    if noise == 1 or noise == 0:
        print ('No noise value available for calibrant/unknown ' +str(num+1))
        noise = 1

    return noise



""""
for noise calculation

will take range of x values and append their y values to a list

def noisefile(x3,y3,minmz,maxmz):
    noiseabun = []
    for i in range(0, len(x3)):
        if (x3[i] >= minmz) and (x3[i] <= maxmz):
            ###print (i)
            noiseabun.append(y3[i])
        
    return noiseabun
"""
