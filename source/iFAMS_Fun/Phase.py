"""
Created on Tue Jan 21 15:53:07 2020
Extracts Fourier Phase information from FT/GT data and plots mass defect profiles
@author: Andrew Swansiger
UPDATES:
    added if/else statement (line 38) to fix baseline correction with one or fewer identified minima
    added def cyclic_int and def defect_library (don't incorporate until v 6.4)
    CHANGED DEFINITION of specmax in both def phase and def g_phase from the IFT global maximum to the IFT centroid
"""
import numpy as np

def line_baseline(y,x,delta):
    xlist = []
    ylist = []
    baseline = []
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
            if y[i] >= 0:
                negative = False
            for j in range(i - delta, i + delta):
                if y[i] > y[j]:
                    large = False
            if large == True and negative == True:
                xlist.append(x[i])
                ylist.append(y[i])

    gap = (x[-1]-xlist[-1]) + (xlist[0] - x[0])
    endslope = (ylist[0] - ylist[-1])/gap
    endb = ylist[0] - endslope*xlist[0]
    if len(xlist)<=1:
        baseline = min(y)
    else:
        for i in range(int(xlist[0]-x[0])):
            baseline.append(endslope*x[i]+endb)
        
        for i in range(0,len(xlist)-1):
            slope = (ylist[(i+1)]-ylist[i])/(xlist[(i+1)]-xlist[i])
            b = ylist[i] - slope*xlist[i]
            for j in range(int(xlist[i]-min(x)),int(xlist[(i+1)]-min(x))):
                baseline.append(slope*x[j] + b)
                
        base = []        
        for i in range(int(x[-1] - xlist[-1]+1)):
            base.append(endslope*(-x[i])+endb)
        base = np.flip(base)
        for i in range(len(base)):
            baseline.append(base[i])
    
    return baseline

def cyclic_int(reconmzint,recontot,submass):
    ###for mass cyclic data sets
    y = recontot
    x = reconmzint
    delta = 3
    indlist=[]
    for i in range (-delta,len(y)-delta):
        minima = True    
        for j in range(i-delta, i+delta):
            if i == 0:
                print(j)
                print(y[j])
            if y[j] < y[i]:
                minima = False
        if minima == True:
            indlist.append(i)
    print(indlist)
    print(len(indlist))
    
    xdelta = (max(x)-min(x))/(len(x))
    sumspec = sum(y)*xdelta
    print(xdelta)
    print(sumspec)
    centroid = []
    stddev = []
    FWHM = []
    relabund = []
    spectot = 0
    for i in range(0,len(indlist)):
        rsum = 0
        ysum = 0
        rdiff = 0
        if indlist[i]==indlist[-1]:
            N = len(x) - indlist[-1] + (indlist[0])
            for j in range(indlist[-1],len(x)):
                rsum+=x[j]*y[j]
                ysum+=y[j]
            for j in range(0,indlist[0]):
                rsum+=(x[j]+submass)*y[j]
                ysum+=y[j]
            c = rsum/ysum
            for j in range(indlist[-1],len(x)):
                rdiff+=y[j]*(x[j]-c)**2
            for j in range(0,indlist[0]):
                rdiff+=y[j]*(x[j]+submass-c)**2
            s = np.sqrt(rdiff/((N-1)*ysum/N))
        else:
            N = indlist[i+1]-indlist[i]
            for j in range(indlist[i],indlist[i+1]):
                rsum+=x[j]*y[j]
                ysum+=y[j]
            c = rsum/ysum
            if x[indlist[i]] > x[indlist[i+1]]:
                for j in range(indlist[i],0):
                    rdiff+=y[j]*(x[j]-submass-c)**2
                for j in range(0,indlist[i+1]):
                    rdiff+=y[j]*(x[j]-c)**2
            else:
                for j in range(indlist[i],indlist[i+1]):
                    rdiff+=y[j]*(x[j]-c)**2
            s = np.sqrt(rdiff/((N-1)*ysum/N))
        fwhm = 2*s*np.sqrt(2*np.log(2))
        ###psum is a feature's partial integral relative to the total integral
        psum = ysum*xdelta/sumspec
        if c > submass:
            c = c-submass
        centroid.append(np.real(c))
        stddev.append(s)
        FWHM.append(fwhm)
        relabund.append(psum)
        spectot += relabund[i]
    relabund.append(spectot)
    print(centroid)
    print(len(centroid))
    print(stddev)
    print(len(stddev))
    print(relabund)
    print(len(relabund))
    print(spectot)
    return centroid,stddev,FWHM

def defect_library(submass):
    
    combin = []
    prob = []
    return combin, prob

def phase(mz, abundance, chgabund, lipidmass, cs, numpoints, numharm, mzspacing, mzmin, proton):
    """
    Extracts local phase information from FT spectra after accounting for global phase
    by shifting the mass spectrum such that the single charge state reconstructions
    are centered at mzmin (performed for each charge state)

    Parameters
    ----------
    mz : Evenly-spaced list of mass spectrum x-values used for interpolation
    abundance : Interpolated list of mass spectrum y-values
    chgabund : Matrix of charge-state-specific mass spectrum reconstruntion y-values
    lipidmass: Mass of the polydisperse species used for the FT analysis
    cs : List of charge states determined from FT analysis (see FT.py)
    numpoints : Number of datapoints in the mass distribution
    numharm: Number of FT harmonics considered for phase analysis
    mzspacing: Average differnce between adjacent x-values of the mass spectrum
    mzmin: Minimum x-value (typically m/z) of loaded spectrum
    proton: Contributing proton mass, dependent on ionizing potential mode (+/-)

    Returns
    -------
    reconmzint: Evenly-spaced MMD x-list for interpolation of charge-state-specific profiles
    recontot : Summation of interpolated charge-state-specific MMD profiles
    reconsub : Matrix of charge-state-specific MMD y-lists
    reconmzsub : Matrix of charge-state-specific MMD xlists
    """
    numfreqint = 9
    ###defines delta for peak integration
    recontot = np.zeros(int(lipidmass),dtype=complex)
    reconsub = []
    reconmzsub = []
    reconmzint = np.linspace(0, lipidmass, int(lipidmass), endpoint = False)
    for i in range(len(cs)):
        chargestate = cs[i]
        print(chargestate)
        abund = 0
        for j in range(len(mz)):
            abund+=chgabund[i][j]*mz[j]
        specmax = abund/sum(chgabund[i])
        print(specmax)
        rotmz = specmax - mzmin
        freq = np.linspace(0,2*np.pi/mzspacing,num=numpoints,endpoint=False)
        ftspec = np.fft.fft(abundance)
        rotwave = np.exp(freq*1j*rotmz)
        counterrotftspec = ftspec*rotwave
        """
        Because abundance is real and positive, ftspec has Hermitian symmetry about the 0 of frequency
        """
        reconft = np.zeros(2*numharm+1,dtype=complex)
    
        harmfreqindex = 0
        prevharmfreqindex = 0
        for harm in range(1,numharm+1,1):
            harmfreq = 2*np.pi*chargestate/float(lipidmass)*float(harm)
            print(harmfreq)
            for i in range(prevharmfreqindex,len(freq),1):
                while (freq[i] <= harmfreq):
                    harmfreqindex = i
                    i += 1
    
            for k in range(int(harmfreqindex-(numfreqint-1)/2),int(harmfreqindex+(numfreqint-1)/2),1):
                reconft[harm] += counterrotftspec[k]
            prevharmfreqindex = harmfreqindex
        
        for harm in range(numharm+1,2*numharm+1,1):
            reconft[harm] = np.conj(reconft[2*numharm-harm+1])
        reconmz = np.linspace((0+(specmax)*chargestate-proton*chargestate),(lipidmass+(specmax)*chargestate-proton*chargestate),num=(2*numharm+1),endpoint=False)
        reconmz = np.remainder(reconmz,lipidmass)
        reconmz = reconmz.tolist()
        recon = np.fft.ifft(reconft)
        
        s = reconmz.index(min(reconmz))
        reconmz = [reconmz[s:-1], [reconmz[-1]], reconmz[0:s]]
        reconmz = [val for sublist in reconmz for val in sublist]
        recon = [recon[s:-1], [recon[-1]], recon[0:s]]
        recon = [val for sublist in recon for val in sublist]
        
        reconmzsub.append(reconmz)
        reconsub.append(recon)
        reconint = np.interp(reconmzint, reconmz, recon, period=lipidmass)
        recontot += reconint
        
    return reconmzint, recontot, reconsub, reconmzsub

def g_phase(mz, abundance, chgabund, lipidmass, cs, numpoints, numharm, mzspacing, mzmin, proton):
    """
    Extracts local phase information from STFT spectra after accounting for global phase
    by shifting the mass spectrum such that the single charge state reconstructions
    are centered at mzmin (performed for each charge state)

    Parameters
    ----------
    mz : Evenly-spaced list of mass spectrum x-values used for interpolation
    abundance : Interpolated list of mass spectrum y-values
    chgabund : Matrix of charge-state-specific mass spectrum reconstruntion y-values
    lipidmass: Mass of the polydisperse species used for the FT analysis
    cs : List of charge states determined from STFT analysis (see STFT.py)
         Selected in order from lowest to highest charge state (right to left)
    numpoints : Number of datapoints in the mass distribution
    numharm: Number of FT harmonics considered for phase analysis
             Selected in ascending order starting from fundamentals
    mzspacing: Average differnce between adjacent x-values of the mass spectrum
    mzmin: Minimum x-value (typically m/z) of loaded spectrum
    proton: Contributing proton mass, dependent on ionizing potential mode (+/-)

    Returns
    -------
    reconmzint: Evenly-spaced MMD x-list for interpolation of charge-state-specific profiles
    recontot : Summation of interpolated charge-state-specific MMD profiles
    reconsub : Matrix of charge-state-specific MMD y-lists
    reconmzsub : Matrix of charge-state-specific MMD xlists
    """
    numfreqint = 9
    ###defines delta for peak integration
    recontot = np.zeros(int(lipidmass),dtype=complex)
    reconsub = []
    reconmzsub = []
    reconmzint = np.linspace(0, lipidmass, int(lipidmass), endpoint = False)
    for i in range(len(cs)):
        chargestate = cs[i]
        print(chargestate)
        abund = 0
        for j in range(len(mz)):
            abund+=chgabund[i][j]*mz[j]
        specmax = abund/sum(chgabund[i])
        print(specmax)
        rotmz = specmax - mzmin
        freq = np.linspace(0,2*np.pi/mzspacing,num=numpoints,endpoint=False)
        ftspec = np.fft.fft(chgabund[i])
        rotwave = np.exp(freq*1j*rotmz)
        counterrotftspec = ftspec*rotwave
        """
        Because abundance is real and positive, ftspec has Hermitian symmetry about the 0 of frequency
        """
        reconft = np.zeros(2*numharm+1,dtype=complex)
    
        harmfreqindex = 0
        prevharmfreqindex = 0
        for harm in range(1,numharm+1,1):
            harmfreq = 2*np.pi*chargestate/float(lipidmass)*float(harm)
            print(harmfreq)
            for i in range(prevharmfreqindex,len(freq),1):
                while (freq[i] <= harmfreq):
                    harmfreqindex = i
                    i += 1
    
            for k in range(int(harmfreqindex-(numfreqint-1)/2),int(harmfreqindex+(numfreqint-1)/2),1):
                reconft[harm] += counterrotftspec[k]
            prevharmfreqindex = harmfreqindex
        
        for harm in range(numharm+1,2*numharm+1,1):
            reconft[harm] = np.conj(reconft[2*numharm-harm+1])
        reconmz = np.linspace((0+(specmax)*chargestate-proton*chargestate),(lipidmass+(specmax)*chargestate-proton*chargestate),num=(2*numharm+1),endpoint=False)
        reconmz = np.remainder(reconmz,lipidmass)
        reconmz = reconmz.tolist()
        recon = np.fft.ifft(reconft)
        
        s = reconmz.index(min(reconmz))
        reconmz = [reconmz[s:-1], [reconmz[-1]], reconmz[0:s]]
        reconmz = [val for sublist in reconmz for val in sublist]
        recon = [recon[s:-1], [recon[-1]], recon[0:s]]
        recon = [val for sublist in recon for val in sublist]
        
        reconmzsub.append(reconmz)
        reconsub.append(recon)
        reconint = np.interp(reconmzint, reconmz, recon, period=lipidmass)
        recontot += reconint
    return reconmzint, recontot, reconsub, reconmzsub