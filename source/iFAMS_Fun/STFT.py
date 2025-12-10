# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 09:11:06 2021
@author: meghandaniels

Updated on Jan 28 2022
    by Kayd Meldrum
"""

import numpy as np
import scipy
import scipy.signal as sp
import scipy.ndimage as spi

def stft(data,freqchan,window,oversampling):
    """
    This performs the STFT
    :param data: The y data that will be STFT
    :param winnum: The number of windows
    :param window: The built window used for the STFT
    :return: X, a matrix the STFT coefficients
    """
    timehop=freqchan
    if window == 'gaussian':
        stdev = int(freqchan/10)
        window = sp.get_window((window, stdev), freqchan) #don't change to 10*stdev or else it will have rounding errors
    else:
        window = sp.get_window(window, freqchan)
    overlap = int(timehop/oversampling) #this is how far the window is slid over for each new m/z pixel in the GT
    X = scipy.array([scipy.fft(window*data[i:i+timehop])
                     for i in range(0, len(data)-timehop, overlap)])
    
    return X

def re_stft(X,winnum,ypadd,yint,oversampling):
    """
    This function can replot the STFT back into the time domain
    It is good to check that the replot looks like your original spectrum
    If not, then some parameters need to be changed

    :param X: the matrix with the STFT coefficients
    :param winnum: the number of windows
    :param freqnum: The number of frequency points = the number of original data points
    :return:
    """
    inspec = np.fft.fft(ypadd)*0
    freqnum = len(inspec)
    timehop = len(X[1])
    overlap = int(winnum/oversampling)
    for n,i in enumerate(range(0, freqnum-timehop, overlap)):
        try:
            inspec[i:i+timehop] += (scipy.ifft(X[n]))
        except IndexError:
            continue

    # data_remove = np.linspace(len(yint), len(ypadd), len(ypadd) - len(yint) + 1)
    # inspec = np.delete(inspec,data_remove)
    inspec = inspec[0:len(yint)]
    inspec[0] = 0
    inspec[-1] = 0
    return inspec

def line_baseline(y,x,delta):
    xlist = []
    ylist = []
    ilist = []
    baseline = []
    #print (len(y))
    for i in range(len(y)):
        if y[i] > 0 or y[i] < 0:
            minX = x[i]
            break
    for i in range(0 + delta, len(y) - delta):
        if x[i] < minX:
            continue
        else:
            large = True
            negative = True
            if y[i] > 0:
                negative = False
            for j in range(i - delta, i + delta):
                if y[i] > y[j]:
                    large = False
            if large == True and negative == True:
                xlist.append(x[i])
                ylist.append(y[i])
                ilist.append(i)
            
    for i in range(ilist[0]-1):
        baseline.append(0)
    for i in range(0,len(xlist)-1):
        slope = (ylist[i+1]-ylist[i])/(xlist[i+1]-xlist[i])
        b = ylist[i] - slope*xlist[i]
        for j in range(ilist[i],ilist[i+1]):
            baseline.append(slope*x[j] + b)
    for i in range(int(len(x)-ilist[-1])+1):
        baseline.append(0)
    return baseline


def auto_peaks(x,y,X,xG,yG,spec,mzstep,fstep,thresh):
    realF = int(len(yG)/2)
    smoothed = spi.gaussian_filter(spec,2)
    dilated = smoothed.copy()
    for i in range(int(mzstep/2-1),len(smoothed)-int(mzstep/2),mzstep):
        for j in range(int(fstep/2-1),len(smoothed[0]-int(fstep/2+1)),fstep):
            dilated[i-int(mzstep/2-1):i+int(mzstep/2+1),j-int(fstep/2-1):j+int(fstep/2+1)] = np.max(smoothed[i-int(mzstep/2-1):i+int(mzstep/2+1),j-int(fstep/2-1):j+int(fstep/2+1)])
    threshold = np.average(smoothed)*thresh
    localmax = []
    lmAmp = []
    testF = []
    maxyGindex = 0
    minxGindex = len(xG)
    for i in range(len(xG)):
        for j in range(realF):
            if yG[j] < 0.1:
                continue
            if smoothed[i][j] == dilated[i][j] and smoothed[i][j] > threshold and yG[j] > 0.1:
                try:
                    if smoothed[i][j] == np.max(smoothed[i-int(mzstep/2-1):i+int(mzstep/2+1),j-int(fstep/2-1):j+int(fstep/2+1)]):
                        localmax.append([i,j])
                        lmAmp.append(spec[i][j])
                        testF.append(yG[j])
                        if j >= maxyGindex:
                            maxyGindex = j 
                        if i <= minxGindex:
                            minxGindex = i
                    else:
                        continue
                except ValueError:
                    continue
    maxlistx = []
    maxlisty = []
    capx = 1/(xG[minxGindex])
    capy = yG[maxyGindex]
    theta = np.zeros(len(localmax))
    radius = np.zeros(len(localmax))
    for i in range(len(localmax)):
        maxlistx.append(1/xG[localmax[i][0]])
        maxlisty.append(yG[localmax[i][1]])
        theta[i] = (np.arctan((maxlisty[i]/capy)/(maxlistx[i]/capx)))
        radius[i] = (np.sqrt((maxlisty[i]/capy)**2+(maxlistx[i]/capx)**2)/np.sqrt(2))
    Rot = 2*np.pi/max(theta)
    theta = theta*Rot
    Bins = np.linspace(0,2*np.pi,int(1700*capy / Rot**2))
    Hits = np.zeros(len(Bins))
    zGaps = np.zeros(len(Bins)).tolist()
    for i in range(len(Bins)):
        zGaps[i] = []
    binspace = 2*np.pi/(len(Bins)-1)
    for i in range(len(localmax)):
        Hits[int(np.round(theta[i]/binspace,0))] += lmAmp[i]
        localmax[i].append(int(np.round(theta[i]/binspace,0)))
        zGaps[int(np.round(theta[i]/binspace,0))].append(radius[i])
    
    for i in range(len(Bins)):
        radlist = np.sort(zGaps[i])
        gaplist = []
        for j in range(len(radlist)-1):
            gaplist.append(radlist[j+1] - radlist[j])
        count = 0
        for j in range(len(gaplist)-1):
            for k in range(j+1,len(gaplist)-1):
               if np.absolute(gaplist[j]-gaplist[k])/gaplist[j] < 0.15:
                   count += 1
        
        if count > 0:
            Hits[i] = Hits[i]*count
    
    ######
    
    if np.average(testF) > 5:
        iso = True
    else:
        iso = False
    if iso == True:
        Dm = 1.00
    else:
        Dm = 162.1
    norm = np.max(Hits)
    zeroY = Hits.copy()/norm
    zeroX = np.tan(Bins/Rot)*capy/capx * Dm
       
    slopes = []
    avgAngle = []
    testx = []
    testy = []
    testAmp = []
    testTheta = []
    deltaMs = []
    Bthresh = 0.18    
    for i in range(len(Hits)):
        HitIndex = np.argmax(Hits)
        testx.append([])
        testy.append([])
        testAmp.append([])
        testTheta.append([])
        for j in range(len(localmax)):
            try:
                if localmax[j][2] == HitIndex:
                    testx[i].append(maxlistx[j])
                    testy[i].append(maxlisty[j])
                    testAmp[i].append(lmAmp[j])
                    testTheta[i].append(theta[j])
            except IndexError:
                continue
        num = 0
        for j in range(len(testx[i])):
            num += testAmp[i][j]*testTheta[i][j]
        avgTheta = num/np.sum(testAmp[i])
        avgAngle.append(avgTheta)
        slopes.append(np.tan(avgTheta/Rot)*capy/capx)
        deltaMs.append(Dm)
                
        left = 1
        right = 2
        for j in range(1,10):
            try:
                if zeroY[HitIndex-j] < zeroY[HitIndex]*0.2 and zeroY[HitIndex-j-1] > (zeroY[HitIndex-j]*1.2):
                    left = j
                    break
                if zeroY[HitIndex-j] == 0:
                    left = j
                    break
            except IndexError:
                break
        for j in range(1,10):
            try:
                if zeroY[HitIndex+j] < zeroY[HitIndex]*0.2 and zeroY[HitIndex+j+1] > (zeroY[HitIndex+j]*1.2):
                    right = j + 1
                    break
                if zeroY[HitIndex+j] == 0:
                    right = j + 1
                    break
            except IndexError:
                break
        Hits[HitIndex-left:HitIndex+right] = 0
        if zeroY[np.argmax(Hits)] < 0.05:
            break
        # if np.max(Hits) < 2:
        #     break
        
    slicex = []
    slicey = [] 
    slicei = []
    slicem = slopes[0]
    for i in range(len(xG)):
        if xG[i] < 1/max(testx[0])-300:
            continue
        if xG[i] > 1/min(testx[0])+300:
            break
        slicex.append(xG[i])
        for j in range(len(spec[0])):
            if abs(yG[j] - slicem/xG[i]) < (yG[2]-yG[1]):
                slicey.append(spec[i][j])
                slicei.append([i,j])
                break
    stop = False
    ytemp = slicey.copy()
    while stop == False:    
        index1 = np.argmax(ytemp)
        Gi1 = slicei[index1][0]
        Fi1 = slicei[index1][1]
        box1,boxI1 = boxer(Gi1,Fi1,spec,xG,yG,slicey[index1],Bthresh)
        for i in range(1,index1-1):
            if slicex[index1-i] > box1[0]:
                continue
            if slicey[index1-i] < slicey[index1]*0.6:
                continue
            if slicey[index1-i] >= slicey[index1]*0.6 and slicey[index1-i-1] <= slicey[index1-i]:
                index2 = index1-i
                break
        try:
            Gi2 = slicei[index2][0]
            Fi2 = slicei[index2][1]
            box2,boxI2 = boxer(Gi2,Fi2,spec,xG,yG,slicey[index2],Bthresh)
            stop = True
        except UnboundLocalError:
            try:
                ytemp[index1-1:index1+2] = 0
            except TypeError:
                ytemp[index1] = 0
    
    Xfun1 = spec*0
    Xfun2 = spec*0
    xcorners = [box1[0],box1[1],box2[0],box2[1]]
    ycorners = [box1[2],box1[3],box2[2],box2[3]]
    for j in range(0, len(spec)):
        if xG[j] < min(xcorners):
            continue
        if xG[j] > max(xcorners):
            break
        for k in range(0, len(spec[0])):
            if yG[k] < min(ycorners):
                continue
            if yG[k] > max(ycorners):
                break
            if xG[j] >= box1[0] and xG[j] <= box1[1] and yG[k] >= box1[2] and yG[k] <= box1[3]:
                Xfun1[j][k] = spec[j][k]
            if xG[j] >= box2[0] and xG[j] <= box2[1] and yG[k] >= box2[2] and yG[k] <= box2[3]:
                Xfun2[j][k] = spec[j][k]
    
    Xfun1max = np.max(Xfun1)
    Gi1 = int(np.argmax(Xfun1)/len(spec[0]))
    Fi1 = np.argmax(Xfun1)%len(spec[0])
    Xfun2max = np.max(Xfun2)
    Gi2 = int(np.argmax(Xfun2)/len(spec[0]))
    Fi2 = np.argmax(Xfun2)%len(spec[0])
    
    mz1 = []
    abun1 = []
    mz2 = []
    abun2 = []
    for i in range(len(x)):
        if x[i] < min(xcorners):
            continue
        if x[i] > max(xcorners):
            break
        if box1[0] < x[i] < box1[1]:
            mz1.append(x[i])
            abun1.append(y[i])
        if box2[0] < x[i] < box2[1]:
            mz2.append(x[i])
            abun2.append(y[i])
    mz1max = mz1[np.argmax(abun1)]
    mz2max = mz2[np.argmax(abun2)]
    
    z1 = int(abs(np.round(float(1/((mz1max/mz2max)-1)),0)))
    if 1/((mz1max/mz2max)-1) < 0:
        z2 = int(z1 - 1)
        box2.append(0)
        box1.append(1)
        zcountH = z1+1
        zcountL = z2-1
    else:
        z2 = int(z1 + 1)
        box1.append(0)
        box2.append(1)
        zcountH = z2+1
        zcountL = z1-1
    Dm = 1/float(((yG[Fi1]/z1) + (yG[Fi2]/z2))/2)
    m = (mz1max*z1-z1*1.00784 + mz2max*z2-z2*1.00784)/2
    print('Selecting protein with rough mass of '+str(np.round(m,2))+' Da')
    print('with a regular change in mass of '+str(np.round(Dm,1))+' Da')
      
    slicex2 = []
    slicey2 = [] 
    for i in range(len(spec[0])):
        slicex2.append(yG[i])
        slicey2.append(spec[Gi1][i])
            
    for i in range(2,10):
        try:
            if spec[Gi1][i*Fi1] < spec[Gi1][Fi1]*0.05:
                HarmRec = i
                break     
        except IndexError:
            HarmRec = i-1
            break

    rlist = [box1,box2]

    return rlist,maxlistx,maxlisty,slopes,avgAngle,theta,radius,deltaMs,zeroX,zeroY,HarmRec,slicex,slicey,slicex2,slicey2

def auto(x,y,X,xG,yG,spec,mzstep,fstep,thresh,cycles):
    #Bthresh = 0.15
    Bthresh = np.average(spec)
    zthresh = 0.08
    Hthresh = 0.01
    realF = int(len(yG)/2)
    dilated = spec.copy()
    for i in range(int(mzstep/2-1),len(spec)-int(mzstep/2),mzstep):
        for j in range(int(fstep/2-1),len(spec[0]-int(fstep/2+1)),fstep):
            dilated[i-int(mzstep/2-1):i+int(mzstep/2+1),j-int(fstep/2-1):j+int(fstep/2+1)] = np.max(spec[i-int(mzstep/2-1):i+int(mzstep/2+1),j-int(fstep/2-1):j+int(fstep/2+1)])
    threshold = np.average(spec)*thresh
    localmax = []
    lmAmp = []
    testF = []
    maxyGindex = 0
    minxGindex = len(xG)
    for i in range(len(xG)):
        for j in range(realF):
            if yG[j] < 0.1:
                continue
            if spec[i][j] == dilated[i][j] and spec[i][j] > threshold and yG[j] > 0.1:
                try:
                    if spec[i][j] == np.max(spec[i-int(mzstep/2-1):i+int(mzstep/2+1),j-int(fstep/2-1):j+int(fstep/2+1)]):
                        localmax.append([i,j])
                        lmAmp.append(spec[i][j])
                        testF.append(yG[j])
                        if j >= maxyGindex:
                            maxyGindex = j 
                        if i <= minxGindex:
                            minxGindex = i
                    else:
                        continue
                except ValueError:
                    continue
    maxlistx = []
    maxlisty = []
    capx = 1/(xG[minxGindex])
    capy = yG[maxyGindex]
    theta = np.zeros(len(localmax))
    radius = np.zeros(len(localmax))
    for i in range(len(localmax)):
        maxlistx.append(1/xG[localmax[i][0]])
        maxlisty.append(yG[localmax[i][1]])
        theta[i] = (np.arctan((maxlisty[i]/capy)/(maxlistx[i]/capx)))
        radius[i] = (np.sqrt((maxlisty[i]/capy)**2+(maxlistx[i]/capx)**2)/np.sqrt(2))
    Rot = 2*np.pi/max(theta)
    theta = theta*Rot
    Bins = np.linspace(0,2*np.pi,int(1700*capy / Rot**2))
    Hits = np.zeros(len(Bins))
    zGaps = np.zeros(len(Bins)).tolist()
    for i in range(len(Bins)):
        zGaps[i] = []
    binspace = 2*np.pi/(len(Bins)-1)
    for i in range(len(localmax)):
        Hits[int(np.round(theta[i]/binspace,0))] += lmAmp[i]
        localmax[i].append(int(np.round(theta[i]/binspace,0)))
        zGaps[int(np.round(theta[i]/binspace,0))].append(radius[i])
    
    for i in range(len(Bins)):
        radlist = np.sort(zGaps[i])
        gaplist = []
        for j in range(len(radlist)-1):
            gaplist.append(radlist[j+1] - radlist[j])
        count = 0
        for j in range(len(gaplist)-1):
            for k in range(j+1,len(gaplist)-1):
               if np.absolute(gaplist[j]-gaplist[k])/gaplist[j] < 0.15:
                   count += 1
        
        if count > 0:
            Hits[i] = Hits[i]*count
    
    ######
    
    slopes = []
    avgAngle = []
    testx = []
    testy = []
    testAmp = []
    testTheta = []
    HarmRec = []
    zlist = []
    rlist = []
    boxIs = []
    dupcount = 0
    try:
        cycles = int(cycles)
    except ValueError:
        cycles = len(Hits)
    maxHit = np.max(Hits)
        
    for i in range(len(Hits)):
        HitIndex = np.argmax(Hits)
        testx.append([])
        testy.append([])
        testAmp.append([])
        testTheta.append([])
        L = i-dupcount
        for j in range(len(localmax)):
            try:
                if localmax[j][2] == HitIndex:
                    testx[L].append(maxlistx[j])
                    testy[L].append(maxlisty[j])
                    testAmp[L].append(lmAmp[j])
                    testTheta[L].append(theta[j])
            except IndexError:
                continue
        num = 0
        for j in range(len(testx[L])):
            num += testAmp[L][j]*testTheta[L][j]
        avgTheta = num/np.sum(testAmp[L])
        avgAngle.append(avgTheta)
        slopes.append(np.tan(avgTheta/Rot)*capy/capx)
        
        slicex = []
        slicey = [] 
        slicei = []
        slicem = slopes[-1]
        for k in range(len(xG)):
            if xG[k] < 1/max(testx[-1])-300:
                continue
            if xG[k] > 1/min(testx[-1])+300:
                break
            slicex.append(xG[k])
            for j in range(len(spec[0])):
                if abs(yG[j] - slicem/xG[k]) < (yG[2]-yG[1]):
                    slicey.append(spec[k][j])
                    slicei.append([k,j])
                    break
        stop = False
        ytemp = slicey.copy()
        while stop == False:    
            index1 = np.argmax(ytemp)
            Gi1 = slicei[index1][0]
            Fi1 = slicei[index1][1]
            box1,boxI1 = boxer(Gi1,Fi1,spec,xG,yG,slicey[index1],Bthresh)
            for k in range(1,index1-1):
                if slicex[index1-k] > box1[0]:
                    continue
                if slicey[index1-k] < slicey[index1]*0.6:
                    continue
                if slicey[index1-k] >= slicey[index1]*0.6 and slicey[index1-k-1] <= slicey[index1-k]:
                    index2 = index1-k
                    break
            try:
                Gi2 = slicei[index2][0]
                Fi2 = slicei[index2][1]
                box2,boxI2 = boxer(Gi2,Fi2,spec,xG,yG,slicey[index2],Bthresh)
                stop = True
            except UnboundLocalError:
                try:
                    ytemp[index1-1:index1+2] = 0
                except TypeError:
                    ytemp[index1] = 0
        
        Xfun1 = spec*0
        Xfun2 = spec*0
        xcorners = [box1[0],box1[1],box2[0],box2[1]]
        ycorners = [box1[2],box1[3],box2[2],box2[3]]
        for j in range(0, len(spec)):
            if xG[j] < min(xcorners):
                continue
            if xG[j] > max(xcorners):
                break
            for k in range(0, len(spec[0])):
                if yG[k] < min(ycorners):
                    continue
                if yG[k] > max(ycorners):
                    break
                if xG[j] >= box1[0] and xG[j] <= box1[1] and yG[k] >= box1[2] and yG[k] <= box1[3]:
                    Xfun1[j][k] = spec[j][k]
                if xG[j] >= box2[0] and xG[j] <= box2[1] and yG[k] >= box2[2] and yG[k] <= box2[3]:
                    Xfun2[j][k] = spec[j][k]
        
        Xfun1max = np.max(Xfun1)
        Gi1 = int(np.argmax(Xfun1)/len(spec[0]))
        Fi1 = np.argmax(Xfun1)%len(spec[0])
        Xfun2max = np.max(Xfun2)
        Gi2 = int(np.argmax(Xfun2)/len(spec[0]))
        Fi2 = np.argmax(Xfun2)%len(spec[0])
        
        #checking against previous selections
        dup1 = False
        dup2 = False
        for j in range(0,len(rlist)):
            if rlist[j][0] < xG[Gi1] < rlist[j][1] and rlist[j][2] < yG[Fi1] < rlist[j][3]:
                dup1 = True
                break
        for j in range(0,len(rlist)):
            if rlist[j][0] < xG[Gi2] < rlist[j][1] and rlist[j][2] < yG[Fi2] < rlist[j][3]:
                dup2 = True
                break
        if dup1 == True and dup2 == True:
            testx.pop()
            testy.pop()
            testAmp.pop()
            testTheta.pop()
            avgAngle.pop()
            slopes.pop()
            Hits[HitIndex] = 0
            dupcount += 1
            continue
            
        
        mz1 = []
        abun1 = []
        mz2 = []
        abun2 = []
        for k in range(len(x)):
            if x[k] < min(xcorners):
                continue
            if x[k] > max(xcorners):
                break
            if box1[0] < x[k] < box1[1]:
                mz1.append(x[k])
                abun1.append(y[k])
            if box2[0] < x[k] < box2[1]:
                mz2.append(x[k])
                abun2.append(y[k])
        mz1max = mz1[np.argmax(abun1)]
        mz2max = mz2[np.argmax(abun2)]
        
        z1 = int(abs(np.round(float(1/((mz1max/mz2max)-1)),0)))
        if 1/((mz1max/mz2max)-1) < 0:
            z2 = int(z1 - 1)
            box2.append(0)
            box1.append(1)
            zcountH = z1+1
            zcountL = z2-1
        else:
            z2 = int(z1 + 1)
            box1.append(0)
            box2.append(1)
            zcountH = z2+1
            zcountL = z1-1
        Dm = 1/float(((yG[Fi1]/z1) + (yG[Fi2]/z2))/2)
        m = (mz1max*z1-z1*1.00784 + mz2max*z2-z2*1.00784)/2
        print('Selecting protein with rough mass of '+str(np.round(m,2))+' Da')
        print('with a regular change in mass of '+str(np.round(Dm,1))+' Da')
        
        
        BoxL = ((box1[0]*z1+z1*1.00784) + (box2[0]*z2+z2*1.00784))/2
        BoxR = ((box1[1]*z1+z1*1.00784) + (box2[1]*z2+z2*1.00784))/2
        BoxB = ((box1[2]/z1) + (box2[2]/z2))/2
        BoxT = ((box1[3]/z1) + (box2[3]/z2))/2
        
        Zthresh = (Xfun1max+Xfun2max)/2*zthresh
        for k in range(1,index1):
            if abs(slicex[index1-k]-(m+zcountH*1.00784)/zcountH) > (slicex[2]-slicex[1])/2:
                continue
            elif slicey[index1-k] > Zthresh:
                zcountH += 1
                continue
            else:
                break
        for k in range(1,len(slicex)-index1):
            if abs(slicex[index1+k]-(m+zcountL*1.00784)/zcountL) > (slicex[2]-slicex[1])/2:
                continue
            elif slicey[index1+k] > Zthresh:
                zcountL -= 1
                continue
            else:
                zcountL += 1
                break
        if (BoxR+zcountL*1.00784)/zcountL > x[-1]:
            zcountL += 1
        if (BoxL+zcountH*1.00784)/zcountH < x[0]:
            zcountH -= 1
            
        zlisttemp = np.array(range(zcountL,zcountH))
        print('using charge states:')
        print(zlisttemp)
        zlist.append(zlisttemp)
        
        HarmRec.append(2)    
        for k in range(2,20):
            try:
                if spec[Gi1][k*Fi1] < spec[Gi1][Fi1]*Hthresh:
                    HarmRec[L] = k
                    break     
            except IndexError:
                HarmRec[L] = k-1
                break

        Wpixel = xG[2]-xG[1]
        Hpixel = yG[2]-yG[1]
        zeroF = 0
        for k in range(len(zlist[L])):
            for j in range(-zeroF,HarmRec[L]):
                left = (BoxL+zlist[L][k]*1.00784)/zlist[L][k] + Wpixel/(len(zlist[L])+1)*abs(zlist[L][k]-z1)
                right = (BoxR+zlist[L][k]*1.00784)/zlist[L][k] - Wpixel/(len(zlist[L])+1)*abs(zlist[L][k]-z1)
                bottom = (j-1) * ((BoxB*zlist[L][k]+BoxT*zlist[L][k])/2)+BoxB*zlist[L][k] + Hpixel/(len(zlist[L])+1)*abs(zlist[L][k]-z1)
                top = (j-1) * ((BoxB*zlist[L][k]+BoxT*zlist[L][k])/2)+BoxT*zlist[L][k] - Hpixel/(len(zlist[L])+1)*abs(zlist[L][k]-z1)
                if i == 0:
                    r4 = k
                else:
                    r4 = len(zlist[L-1]) + k
                rlist.append([left,right,bottom,top,r4])
                # Wbox = (right-left)
                # Hbox = (top-bottom)
                # zerorlist.append([left+Wbox*0.0,right-Wbox*0.0,bottom+Hbox*0,top-Hbox*0,i])
        
        boxIs.append(boxI1)
        boxIs.append(boxI2)
        
        if len(zlist) == cycles:
            break
        
        Hits[HitIndex] = 0
        if np.max(Hits)/maxHit < 0.05 and cycles == len(Hits):
            break
    
    return rlist,zlist,HarmRec


def masscalc(mzx0,mzx1):
    z0 = int(abs(np.round(float(1/((mzx0/mzx1)-1)),0)))
    if 1/((mzx0/mzx1)-1) < 0:
        z1 = int(z0 - 1)
    else:
        z1 = int(z0 + 1)
    m = float(((mzx0*z0-z0)+(mzx1*z1-z1))/2)

    return m, z0, z1

def boxer(Gi,Fi,X,xG,yG,Gmax,thresh0):
    L = 0
    Li = 0
    #thresh = thresh0*Gmax
    thresh = thresh0
    if yG[Fi] < 1:
        fudgeT = 0.88
        fudgeB = 1.72
    else:
        fudgeT = 0.7
        fudgeB = 0.7
    while L == 0:
        for i in range(1,Fi):
            if (X[Gi-i][Fi]-thresh)/(Gmax-thresh) <= 0.05:
                L = xG[Gi-i-1]
                Li = Gi-i-1
                break
        thresh += 0.05
    R = 0
    Ri = 0
    while R == 0:
        for i in range(1,Fi):
            if (X[Gi+i][Fi]-thresh)/(Gmax-thresh) <= 0.05:
                R = xG[Gi+i+1]
                Ri = Gi+i+1
                break
        thresh += 0.05
    B = 0
    Bi = 0
    while B == 0:
        for i in range(1,Fi):
            if (X[Gi][Fi-i]-thresh)/(Gmax-thresh) <= 0.05:
                B = yG[Fi-i-1]
                Bi = Fi-i-1
                break
        thresh += 0.07
    T = 0
    Ti = 0
    while T == 0:
        for i in range(1,Fi):
            if (X[Gi][Fi+i]-thresh)/(Gmax-thresh) <= 0.05:
                T = yG[Fi+i+1]
                Ti = Fi+i+1
                break
        thresh += 0.03
    """
    try:
        T
    except UnboundLocalError:
        for i in range(1,Fi):
            if X[Gi-1][Fi+i] <= 0.25*Gmax:
                T = yG[Fi+i]
                Ti = Fi+i
                break
    """
    
    return [L,R,B,T],[Li,Ri,Bi,Ti]

def peak_finder(rlist,glist,xgabor,ygabor,X,xMZ,xFT,NFrms,mzchan):
    #completes fundamental series from first two peaks
    try:
        freqchan = len(ygabor)
        mzx0 = xMZ[np.argmax(glist[0])]
        mzx1 = xMZ[np.argmax(glist[1])]
    except TypeError:
        mzx0 = glist
        mzx1 = NFrms * 100
    except IndexError:
        mzx0 = glist
        mzx1 = NFrms * 100
    z0 = int(abs(np.round(float(1/((mzx0/mzx1)-1)),0)))
    if 1/((mzx0/mzx1)-1) < 0:
        z1 = int(z0 - 1)
    else:
        z1 = int(z0 + 1)
    m = float(((mzx0*z0-z0)+(mzx1*z1-z1))/2)
    print('estimating charge states of entered fundamentals to be ' +str(z0)+ ' and ' +str(z1))
    print('with a mass of ' + str(np.round(m,3)) + ' Da')
    print('optimizing fundamental selections for search')
    
    Xfun0 = abs(X*0)
    Xfun1 = abs(X*0)
    for j in range(0, len(X)):
        for k in range(0, len(X[0])):
            if xgabor[j] >= rlist[0][0] and xgabor[j] <= rlist[0][1] and ygabor[k] >= rlist[0][2] and ygabor[k] <= rlist[0][3] and abs(X[j][k]) > float((4*NFrms)/np.sqrt(mzchan)):
                Xfun0[j][k] = abs(X[j][k])
            if xgabor[j] >= rlist[1][0] and xgabor[j] <= rlist[1][1] and ygabor[k] >= rlist[1][2] and ygabor[k] <= rlist[1][3] and abs(X[j][k]) > float((4*NFrms)/np.sqrt(mzchan)):
                Xfun1[j][k] = abs(X[j][k])
   
    Xfun0max = np.max(Xfun0)
    x0max = int(np.argmax(Xfun0)/freqchan)
    y0max = np.argmax(Xfun0)%freqchan
    Xfun1max = np.max(Xfun1)
    x1max = int(np.argmax(Xfun1)/freqchan)
    y1max = np.argmax(Xfun1)%freqchan
    gaborm = float(((xgabor[x0max]*z0-z0) + (xgabor[x1max]*z1-z1))/2)
    print('using a mass of ' + str(np.round(gaborm,3)) + ' Da to find other fundamentals')
    Dm = float(((ygabor[y0max]/z0) + (ygabor[y1max]/z1))/2)
    print('using frequency from a change in mass of '+ str(np.round(float(1/Dm),2)) +' Da')
    """
    try:
        thresh = 0.14
        [mzleft0,mzright0,Flow0,Fhigh0],[Li0,Ri0,Bi0,Ti0] = boxer(x0max,y0max,abs(X),xgabor,ygabor,Xfun0max,thresh)
        [mzleft1,mzright1,Flow1,Fhigh1],[Li1,Ri1,Bi1,Ti1] = boxer(x1max,y1max,abs(X),xgabor,ygabor,Xfun1max,thresh)
        
        print('new boxes')
    except UnboundLocalError:
    """
    thresh = 0.01 #threshold for shrinking boxes to signal
    for i in range(1,x0max):
        if Xfun0[x0max-i][y0max] == 0 or Xfun0[x0max-i][y0max] < thresh*Xfun0max:
            mzleft0 = xgabor[x0max-i]
            break
        else:
            mzleft0 = rlist[0][0]
    for i in range(1,x0max):
        if Xfun0[x0max+i][y0max] == 0 or Xfun0[x0max+i][y0max] < thresh*Xfun0max:
            mzright0 = xgabor[x0max+i]
            break
        else:
            mzright0 = rlist[0][1]
    for i in range(1,y0max):
        if Xfun0[x0max][y0max-i] == 0 or Xfun0[x0max][y0max-i] < thresh*Xfun0max:
            Flow0 = ygabor[y0max-i]
            break
        else:
            Flow0 = rlist[0][2]
    for i in range(1,y0max):
        if Xfun0[x0max][y0max+i] == 0 or Xfun0[x0max][y0max+i] < thresh*Xfun0max:
            Fhigh0 = ygabor[y0max+i]
            break
        else:
            Fhigh0 = rlist[0][3]
            
    for i in range(1,x1max):
        if Xfun1[x1max-i][y1max] == 0 or Xfun1[x1max-i][y1max] < thresh*Xfun1max:
            mzleft1 = xgabor[x1max-i]
            break
        else:
            mzleft1 = rlist[1][0]
    for i in range(1,x1max):
        if Xfun1[x1max+i][y1max] == 0 or Xfun1[x1max+i][y1max] < thresh*Xfun1max:
            mzright1 = xgabor[x1max+i]
            break
        else:
            mzright1 = rlist[1][1]
    for i in range(1,y1max):
        if Xfun1[x1max][y1max-i] == 0 or Xfun1[x1max][y1max-i] < thresh*Xfun1max:
            Flow1 = ygabor[y1max-i]
            break
        else:
            Flow1 = rlist[1][2]
    for i in range(1,y1max):
        if Xfun1[x1max][y1max+i] == 0 or Xfun1[x1max][y1max+i] < thresh*Xfun1max:
            Fhigh1 = ygabor[y1max+i]
            break
        else:
            Fhigh1 = rlist[1][3]

    mleft = float(((mzleft0*z0-z0) + (mzleft1*z1-z1))/2)
    mright = float(((mzright0*z0-z0) + (mzright1*z1-z1))/2)
    Dmlow = float((Flow0/z0 + Flow1/z1)/2)
    Dmhigh = float((Fhigh0/z0 + Fhigh1/z1)/2)
    checkmagnitude = 0.1*min(Xfun0max,Xfun1max) #for determining whether to add fundamental
    
    print('searching...')
    zlist = []
    zcount = z0-30
    while zcount < 1 or float((mright+zcount)/zcount) > xMZ[-1] or float(Dmlow*zcount) < 0.001:
        zcount += 1
    for i in range(1,len(xgabor)):
        if xgabor[-i] > float((gaborm+zcount)/zcount+1.5*(xMZ[-1]-xgabor[0])/mzchan):
            continue
        elif xgabor[-i] < float((gaborm+zcount)/zcount-1.5*(xMZ[-1]-xgabor[0])/mzchan):
            if len(zlist) > 2:
                break
            else:
                zcount += 1
                continue
        for j in range(0,len(ygabor)):
            if ygabor[j] <= float(Dm*zcount-1.5*xFT[-1]/freqchan):
                continue
            elif ygabor[j] >= float(Dm*zcount+1.5*xFT[-1]/freqchan):
                break
            elif abs(X[-i][j]) > checkmagnitude:
                zlist.append(int(zcount))
                zcount += 1
                break
                
    return zlist,mleft,mright,Dmlow,Dmhigh

