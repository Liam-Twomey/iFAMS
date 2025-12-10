# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:39:35 2024

@author: Kayd Meldrum
"""

import numpy as np

def integrate(y,yc,x,delta,minY,minX, maxX, ntol):
    delta = int(delta*(len(x)/(x[-1]-x[0])))
    if delta < 1:
        delta = 1
    xlist = []
    ylist = []
    xint = []
    xref = []
    yref = []
    sumlist = []
    minL = []
    minR = []
    centx = []
    height = []
    for i in range(0 + delta, len(y) - delta):
        if x[i] < minX:
            continue
        if x[i] > maxX:
            break
        else:
            large = True
            baseline = True
            if y[i] < minY:
                baseline = False
            for j in range(i - delta, i + delta):
                if y[i] < y[j]:
                    large = False
            if large == True and baseline == True:
                xint.append(i)
                xlist.append(x[i])
                ylist.append(y[i])

    for i in range(0,len(xint)):
        sumy = 0
        refnum = xint[i]
        leftnum=0
        rightnum=0
        mincount=0
        maxcount=0
        localminL=[0]
        localmaxL=[refnum]
        localminR=[0]
        localmaxR=[refnum]
        j=0
        while j < len(x):
            if y[localmaxL[maxcount]-j] < y[localmaxL[maxcount]-(j+1)]:
                localminL.append(localmaxL[maxcount]-j)
                mincount+=1
                for k in range(j,len(x)):
                    if (y[localmaxL[maxcount]-k] > y[localmaxL[maxcount]-(k+1)]):
                        localmaxL.append(localmaxL[maxcount]-k)
                        maxcount+=1
                        break
                if maxcount == 0:
                    leftnum = localminL[1]
                else:
                    if (((y[localmaxL[maxcount]]-y[localminL[mincount]])/max(y.tolist()) > ntol) or (y[localmaxL[maxcount]] < minY)):
                            leftnum = localminL[mincount]
                    else:
                        j=0
            try:
                if x[localmaxL[maxcount]-(j)] == minX or y[localmaxL[maxcount]-(j+1)] <= 0:
                    leftnum = localmaxL[maxcount]-j
                    localminL.append(localmaxL[maxcount]-j)
            except IndexError:
                leftnum = localmaxL[maxcount]-j
                localminL.append(localmaxL[maxcount]-j)
            if leftnum > 0:
                mincount=0
                maxcount=0
                j=0
                break
            else:
                j+=1
        while j < len(x):
            if y[localmaxR[maxcount]+j] < y[localmaxR[maxcount]+(j+1)]:
                localminR.append(localmaxR[maxcount]+j)
                mincount+=1
                for k in range(j,len(x)):
                    if (y[localmaxR[maxcount]+k] > y[localmaxR[maxcount]+(k+1)]):
                        localmaxR.append(localmaxR[maxcount]+k)
                        maxcount+=1
                        break
                if maxcount == 0:
                    rightnum = localminR[1]
                else:
                    if (((y[localmaxR[maxcount]]-y[localminR[mincount]])/max(y.tolist()) > ntol) or (y[localmaxR[maxcount]] < minY)):
                            rightnum = localminR[mincount]
                    else:
                        j=0
            try:
                if x[localmaxR[maxcount]+(j)] == maxX or y[localmaxR[maxcount]+(j+1)] <= 0:
                    rightnum = localmaxR[maxcount]+j
                    localminR.append(localmaxR[maxcount]+j)
            except IndexError:
                rightnum = localmaxR[maxcount]+j
                localminR.append(localmaxR[maxcount]+j)
            if rightnum > 0:
                mincount=0
                maxcount=0
                del(j)
                break
            else:
                j+=1

        if ((refnum - leftnum)< (rightnum - refnum)):
            rightnum = refnum + (refnum-leftnum) + 1
        else:
            leftnum = refnum - (rightnum - refnum)
        if np.shape(yc[leftnum:rightnum])[0] == 0:
            leftnum = int(leftnum-1)
            rightnum = int(rightnum+1)
        if (y[refnum]-y[leftnum])/max(y.tolist()) < ntol and (y[refnum]-y[rightnum])/max(y.tolist()) < ntol:
            continue
        xref.append(x[leftnum:rightnum])
        yref.append(yc[leftnum:rightnum])
        height.append(float(y[refnum]))
        weightx = 0
        weight = 0
        for j in range(leftnum,rightnum):
            sumy += y[j]
            if y[j]-min(y[leftnum:rightnum]) >= 0.75*(y[refnum]-min(y[leftnum:rightnum])):
                weightx += x[j]*y[j]
                weight += y[j]
        centx.append(float(weightx/weight))
        sumy = sumy-np.average([y[leftnum],y[rightnum-1]])
        sumlist.append(round(sumy,4))
        minL.append(localminL[-1])
        minR.append(localminR[-1])

    return xlist,ylist,xref,yref,sumlist,xint,minL,minR,centx,height

def boundsmatch(xzero,zerofull,bounds):
    boundtemp = []
    for i in range(0,len(bounds),2):
        if bounds[i] < xzero[0]:
            continue
        elif bounds[i+1] > xzero[-1]:
            break
        else:
            boundtemp.append(bounds[i])
            boundtemp.append(bounds[i+1])
    bounds = boundtemp
    sumlist = np.zeros(int(len(bounds)/2))
    height = np.zeros(int(len(bounds)/2))
    xcent = np.zeros(int(len(bounds)/2))
    xref = []
    yref = []
    peakindex = []
    
    poW = xzero[2]-xzero[1]
    for i in range(0, len(bounds),2):
        peakindex.append(int(i/2))
        for j in range(0,len(xzero)):
            if xzero[j] >= bounds[i]:
                i1 = j
                break
        for j in range(i1,len(xzero)):
            if xzero[j] > bounds[i+1]:
                i2 = j
                break
        sumlist[int(i/2)] = np.sum(zerofull[i1:i2])-np.average([zerofull[i1],zerofull[i2-1]])
        height[int(i/2)] = np.max(zerofull[i1:i2])
        heightI = np.argmax(zerofull[i1:i2])+i1
        xref.append(xzero[i1:i2])
        yref.append(zerofull[i1:i2])
        weightx = 0
        weight = 0
        for j in range(0,len(xref[int(i/2)])):
            if yref[int(i/2)][j]-np.min(yref[int(i/2)]) >= 0.75*(height[int(i/2)]-np.min(yref[int(i/2)])):
                weightx += xref[int(i/2)][j] * yref[int(i/2)][j]
                weight += yref[int(i/2)][j]
        try:
            xcent[int(i/2)] = float(weightx / weight)
        except ZeroDivisionError:
            xcent[int(i/2)] = xzero[heightI]
    sumlist = sumlist.tolist()
    height = height.tolist()
    xcent = xcent.tolist()
    
    return xref, yref, sumlist, height, xcent, peakindex


def smooth(x,y,width):
    """
    This performs a Gaussian smoothing using FT and preserves area under the curve
    :param x: x data of spectrum to be smoothed
    :param y: y data of spectrum to be smoothed
    :param width: width of Gaussian window (FWHM). Values of 1.4-1.8 are typically efficient to smooth over isotope-resolution.
    :return: y_smooth, array of smoothed y data
    """

    x_gauss = np.round(np.linspace(min(x),max(x),len(x)),4)
    FFT = np.fft.fft(y)
    sig = width/(2*np.sqrt(2*np.log(2)))
    gaussian = np.exp(-np.power(x_gauss-x_gauss[0],2) / (2*np.power(sig,2))) + np.exp(-np.power(x_gauss-x_gauss[-1]-((max(x_gauss)-min(x_gauss))/len(x_gauss)),2) / (2*np.power(sig,2)))
    FFT = FFT*np.fft.fft(gaussian)
    area = np.sum(gaussian)
    y_smooth = np.real(np.fft.ifft(FFT))/area

    return y_smooth

def baseline_sub(x,y,xref,yref,peakindex,minL,minR):
    # Old segmented baseline correction
    interminy = []
    interminx = []
    ytemp = []
    xtemp = []
    for i in range(minL[0],peakindex[0]):
        ytemp.append(y[i])
        xtemp.append(x[i])
    interminy.append(min(ytemp))
    minindex = ytemp.index(min(ytemp))
    interminx.append(xtemp[minindex])
    for i in range(0,int(len(peakindex)-1)):
        ytemp = []
        xtemp = []
        for j in range(peakindex[i],peakindex[i+1]):
            ytemp.append(y[j])
            xtemp.append(x[j])
        interminy.append(min(ytemp))
        minindex = ytemp.index(min(ytemp))
        interminx.append(xtemp[minindex])
    ytemp = []
    xtemp = []
    for i in range(peakindex[-1],minR[-1]):
        ytemp.append(y[i])
        xtemp.append(x[i])
    interminy.append(min(ytemp))
    minindex = ytemp.index(min(ytemp))
    interminx.append(xtemp[minindex])
    
    base_m = np.zeros(len(peakindex))
    base_b = np.zeros(len(peakindex))
    newy = np.zeros(len(y))
    yref2 = []
    sumlist = []
      
    for i in range(0,len(peakindex)):
        yreftemp = []
        dy = interminy[i]-interminy[i+1]
        dx = interminx[i]-interminx[i+1]
        base_m[i] = dy/dx
        base_b[i] = interminy[i]-base_m[i]*interminx[i]
        
        Lindex=x.tolist().index(interminx[i])
        Rindex=x.tolist().index(interminx[i+1])
        
        for j in range(Lindex,Rindex):
            newy[j] = y[j] - (base_m[i]*x[j]+base_b[i])
            if newy[j] < 0:
                newy[j] = 0

        for j in range(0,len(yref[i])):
            peaky = yref[i][j] - (base_m[i]*xref[i][j]+base_b[i])
            if peaky < 0:
                peaky = 0
            yreftemp.append(peaky)
        yref2.append(yreftemp)
        sumlist.append(sum(yreftemp))      
                
    return newy, yref2, sumlist
 
def baseline_sub2(x,y,xref,yref,peakindex): 
    # New segmented baseline correction
    base_m = np.zeros(len(peakindex))
    base_b = np.zeros(len(peakindex))
    yref2 = []
    sumlist = []
    height = []
      
    for i in range(0,len(peakindex)):
        try:
            yreftemp = []
            yreftemp2 = []
            dy = yref[i][-1]-yref[i][0]
            dx = xref[i][-1]-xref[i][0]
            base_m[i] = dy/dx
            base_b[i] = yref[i][0]-base_m[i]*xref[i][0]
            
            for j in range(0,len(yref[i])):
                peaky = yref[i][j] - (base_m[i]*xref[i][j]+base_b[i])
                if peaky < 0:
                    peaky = 0
                    yreftemp2.append(yref[i][j])
                else:
                    yreftemp2.append(base_m[i]*xref[i][j]+base_b[i])
                yreftemp.append(peaky)
            yref2.append(yreftemp2)
            sumlist.append(sum(yreftemp))  
            height.append(max(yreftemp))
        except IndexError:
            continue
                
    return yref2, sumlist, height