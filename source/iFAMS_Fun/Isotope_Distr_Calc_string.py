# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 15:53:07 2020
Computes the Mass Distribution of a Protein and Oligonucleotide sequences
Calls atom numbers for each letter from a separate dictionary script
@author: Andrew Swansiger
UPDATES:
    added def ssDNA_seq and def ssRNA_seq
    adjusted FFT calculation to allow for non-integer amounts of atoms from 'averaged' constituents
"""

import numpy as np
from . import Isotope_Library_string as IsoLib

def datawindow(sum_elem,elem_list,multimer):
    """
    Defines the list size used for constructing the mass distribution,
    such that it includes all significant peak contributions (>0.001%) at a minimum

    Parameters
    ----------
    sum_elem : list of total atoms of each element in a monomeric unit
    elem_list : ordered list of atomic symbols in the element dictionary
    multimer : number of monomeric units that make up the system 

    Returns
    -------
    iso_min : Smallest possible mass of the system, used to define start of the xlist
    iso_max :  Largest possible mass of the system, used to determine overall xlist size
    iso_avg :  Expected average mass of the system based on isotopic probabilities
    distr_max : Span of the mass distribution
    distr_spc : Number of datapoints in the mass distribution,
                defined as 1000 times the span of the distribution
    mass_corr : mass correction applied for mass shifts due to rounding error 
    """
    iso_min = 0
    iso_max = 0
    iso_avg = 0
    mass_corr = 0
    
    for i in range(len(elem_list)):
        if sum_elem[i] > 0:
            print(elem_list[i],sum_elem[i],len(IsoLib.mass_drift[elem_list[i]]),len(IsoLib.elem[elem_list[i]][1]))
            iso_min+=int(sum_elem[i])*(IsoLib.elem[elem_list[i]][1][0][0])
            iso_max+=sum_elem[i]*(IsoLib.elem[elem_list[i]][1][-1][0])
            iso_err = 0
            for j in range(len(IsoLib.elem[elem_list[i]][1])):
                iso_avg+=(sum_elem[i])*(IsoLib.elem[elem_list[i]][1][j][1])*(IsoLib.elem[elem_list[i]][1][j][0])
                iso_err+=(sum_elem[i])*(IsoLib.elem[elem_list[i]][1][j][1])*(IsoLib.mass_drift[elem_list[i]][j])
            mass_corr+=iso_err
    iso_min = iso_min*multimer
    iso_max = iso_max*multimer
    iso_avg = iso_avg*multimer
    mass_corr = np.round(mass_corr*multimer, decimals=3)
    distr_max = np.round(2*(iso_avg-iso_min), decimals=3)
    distr_spc = int(distr_max*1000)
    
    for i in range(100):
        q = 2**i
        if q>=distr_spc:
            distr_spc = q
            distr_max = np.round(distr_spc/1000, decimals=3)
            break
    
    if distr_max<=300:
        distr_max = 300
        distr_spc = distr_max*1000
    '''    
    if distr_spc>=10**6:
        distr_max = int(2**(np.floor(np.log2(iso_max-iso_min + 2))))
        distr_spc = int(distr_max*1000)
    '''    
    return iso_min, iso_max, iso_avg, distr_max, distr_spc, mass_corr

def FT_tot_string(x,y,k,sum_elem,elem_list,multimer):
    """
    Generates the system mass distribution 
    by multiplying atom distributions in the Fourier domain
    
    read an empirical formula
    
    Parameters
    ----------
    x : range of mass values
    y : list of zeros 
    sum_elem : list of total atoms of each element in a monomeric unit
    elem_list : ordered list of atomic symbols in the element dictionary
    multimer : number of monomeric units that make up the system 

    Returns
    -------
    FTsystem : total mass distribution Fourier transform
               for all atoms in the system
    """
    FTsystem = IsoLib.null(x,y)
    
    for i in range(len(elem_list)):
        if sum_elem[i] > 0:
            y_prime = y
            for j in range(len(IsoLib.elem[elem_list[i]][1])):
                y_prime = np.delete(y_prime,-1)
                y_prime = np.insert(y_prime,
                                    int(np.argwhere(x==IsoLib.elem[elem_list[i]][1][j][0])),
                                    IsoLib.elem[elem_list[i]][1][j][1])
            FT = (np.fft.fft(y_prime))**np.floor((sum_elem[i]))
            if np.remainder(sum_elem[i],1) != 0:
                #fracFT = np.fft.fft(y_prime)*np.exp(1j*k*(1-np.remainder(sum_elem[i],1)))
                y_prime = y
                for j in range(len(IsoLib.elem[elem_list[i]][1])):
                    y_prime = np.delete(y_prime,-1)
                    y_prime = np.insert(y_prime,
                                        int(np.argwhere(x==np.round(IsoLib.elem[elem_list[i]][1][j][0]*np.remainder(sum_elem[i],1),decimals=3))),
                                        IsoLib.elem[elem_list[i]][1][j][1])
                fracFT = np.fft.fft(y_prime)
                FT = FT*fracFT
            FTsystem = FTsystem*FT
            
    FTsystem = FTsystem**multimer
    
    return FTsystem

def shftdata(xP, FFT, iso_min, distr_max, distr_spc):
    """
    Adjusts the data such that the minimum mass of the system is near the start of the list,
    then takes the Inverse Fourier Transform of this counterrotated list

    Parameters
    ----------
    xP : Mass value frame shifted to keep minimum mass labeled with the correct mass (Da)
    FFT : Fourier transformed mass distribution
    iso_min : Smallest possible mass of the system, used to define start of the xlist
    distr_max : Maximum value of the mass distribution
    distr_spc : Number of datapoints in the mass distribution,
                defined as 1000 times the span of the distribution

    Returns
    -------
    xrot : Mass value frame shifted to accomodate for the counterrotated Fourier data
    index_min : index of the floor of the minimum mass value in the original mass frame, 
                used to define the counterrotation
    shiftdatafft : Counterrotated Fourier data
    IFFT : Inverse Fourier Transform of the counterrotated data
    """
    for i in range(len(xP)):
        if xP[i]==np.floor(iso_min):
            index_min = i
            break
    xrot = np.round(np.linspace(xP[index_min], xP[index_min]+distr_max-0.001, distr_spc),decimals=3)
    
    shiftfft = []
    for i in range(len(FFT)):
        fft = np.exp(1j*2*np.pi*((index_min)*i/len(xP)))
        shiftfft.append(fft)
    shiftdatafft = FFT*shiftfft

    IFFT = np.absolute(np.fft.ifft(shiftdatafft))
    return xrot, index_min, shiftdatafft, IFFT

def distr_stats(xP,IFFT):
    """
    Calculates the average mass and mass variance of the system
    accounting for masses down to a threshold probability of occurrance

    Parameters
    ----------
    xP : modified xlist, starting close to the lowest possible system mass
    IFFT : inverse Fourier transform of the multiplied atom distributions
    
    Returns
    -------
    MWavg : average molecular weight of the system
    MWvar : mass variance of the system
    """
    MWtot = 0
    num = 0
    for i in range(len(IFFT)):
        MWtot += IFFT[i]*xP[i]
        num += IFFT[i]
    MWavg = MWtot/num
    
    MWvar = 0
    for i in range(len(IFFT)):
        MWvar += ((xP[i]**2)*IFFT[i])
    MWvar = MWvar/num-MWavg**2
    return MWavg, MWvar

def protein_seq(sum_elem,protein):
###atom count adjustments are for 5' and 3' end groups formed from hydrolysis
    sumAA = [0, 0, 0, 0, 0]    
    for i in range(len(protein)):
        amino = IsoLib.AA[protein[i]]
        sumAA = np.add(sumAA,amino)
    sum_elem[0] += sumAA[0]+2
    sum_elem[5] += sumAA[1]
    sum_elem[6] += sumAA[2]
    sum_elem[7] += sumAA[3]+1
    sum_elem[15] += sumAA[4]
    return sum_elem

def ssRNA_seq(sum_elem,protein):
###atom count adjustments are for 5' oxide protonation and 3' phosphate hydroxyl substitution
    sumNB = [0, 0, 0, 0, 0]    
    for i in range(1,len(protein)):
        ssRNA = IsoLib.RNA[protein[i]]
        sumNB = np.add(sumNB,ssRNA)
    sum_elem[0] += sumNB[0]+1
    sum_elem[5] += sumNB[1]
    sum_elem[6] += sumNB[2]
    sum_elem[7] += sumNB[3]-2
    sum_elem[14] += sumNB[4]-1
    return sum_elem

def ssDNA_seq(sum_elem,protein):
###atom count adjustments are for 5' oxide protonation and 3' phosphate hydroxyl substitution
    sumNB = [0, 0, 0, 0, 0]    
    for i in range(1,len(protein)):
        ssDNA = IsoLib.DNA[protein[i]]
        sumNB = np.add(sumNB,ssDNA)
    sum_elem[0] += sumNB[0]+1
    sum_elem[5] += sumNB[1]
    sum_elem[6] += sumNB[2]
    sum_elem[7] += sumNB[3]-2
    sum_elem[14] += sumNB[4]-1
    return sum_elem

def adduct(add,cluster,sum_elem):
    formula = list(add.split())
    molec = []
    num_atom = []
    for i in range(len(formula)):
        if np.remainder(i,2)==0:
            molec.append(str(formula[i]))
        if np.remainder(i,2)==1:
            num_atom.append(int(formula[i]))
    for i in range(len(molec)):
        sum_elem[IsoLib.elem[molec[i]][0]]+=num_atom[i]*cluster
    return sum_elem
