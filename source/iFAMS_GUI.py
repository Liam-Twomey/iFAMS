"""
Authors: Sean P. Cleary, Andrew K. Swansiger, Meghan M. Daniels, Kayd L. Meldrum, Lily Miller, James S. Prell
Copyright 2016-2023, University of Oregon

iFAMS Software Academic License Agreement

The iFAMS software ("Software") has been developed by the contributing researchers James S. Prell, Sean P. Cleary,
Andrew K. Swansiger, Meghan M. Daniels, Kayd L. Meldrum, and Lily Miller ("Developers") and made available through the University of Oregon ("UO") for your internal,
non-profit research use. The Software was developed with funding support from UO. UO and the Developers allow
researchers at your Institution to run, display, copy and modify Software on the following conditions:

1. The Software remains at your Institution and is not published, distributed, or otherwise transferred or
made available to other than Institution employees and students involved in research under your supervision.

2. You agree to make results generated using Software available to other academic researchers for non-profit
research purposes. If You wish to obtain Software for any commercial purposes, including fee-based service
projects, You will need to execute a separate licensing agreement with the University of Oregon and pay a
fee. In that case please contact: techtran@uoregon.edu. 

3.  You retain in Software and any modifications to Software, the copyright, trademark, or other notices
pertaining to Software as provided by UO and Developers.

4. You provide the Developers with feedback on the use of the Software in your research, and that the
Developers and UO are permitted to use any information You provide in making changes to the Software. All bug
reports and technical questions shall be sent to the email address: jprell@uoregon.edu

5. You acknowledge that the Developers, UO and its licensees may develop modifications to Software that may be
substantially similar to your modifications of Software, and that the Developers, UO and its licensees shall
not be constrained in any way by You in Developer’s, UO’s or its licensees’ use or management of such
modifications. You acknowledge the right of the Developers and UO to prepare and publish modifications
to Software that may be substantially similar or functionally equivalent to your modifications and
improvements, and if You obtain patent protection for any modification or improvement to Software You agree
not to allege or enjoin infringement of your patent by the Developers, UO or by any of UO’s licensees
obtaining modifications or improvements to Software from the UO or the Developers.

6. You agree to acknowledge the contribution Developers and Software make to your research, and cite
appropriate references about the Software in your publications. The current citations for the Software
can be found at:

https://dx.doi.org/10.1021/acs.analchem.6b01088
https://dx.doi.org/10.1007/s13361-018-2018-7
https://dx.doi.org/10.1002/cphc.201900022

7. Any risk associated with using the Software at your institution is with You and your Institution.
Software is experimental in nature and is made available as a research courtesy "AS IS," without obligation
by UO to provide accompanying services or support.

8. UO AND THE DEVELOPERS EXPRESSLY DISCLAIM ANY AND ALL WARRANTIES REGARDING THE SOFTWARE, WHETHER EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES PERTAINING TO NON-INFRINGEMENT, MERCHANTABILITY OR FITNESS
FOR A PARTICULAR PURPOSE.
"""

from .iFAMS_Fun import MSload as load
from .iFAMS_Fun import FT as FT
from .iFAMS_Fun import STFT as STFT
from .iFAMS_Fun import Integrate as Int
from .iFAMS_Fun import Phase as phase
from .iFAMS_Fun import Calibration as cal
from .iFAMS_Fun import Isotope_Distr_Calc_string as IsoCalc
from .iFAMS_Fun import batch_param_template as batch
import os
import sys
import numpy as np
from PyQt5 import QtWidgets, QtGui, QtCore 
from PyQt5.Qt import QApplication
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas  ###
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar  ###
from matplotlib.figure import Figure
import scipy.signal as sp
from matplotlib.widgets import RectangleSelector
from matplotlib.patches import  Rectangle
import matplotlib.pyplot as plt
import datetime
import subprocess
from matplotlib.font_manager import FontProperties
color=np.concatenate((plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),
                      plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),
                      plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20)),plt.cm.tab20(np.linspace(0,1,20))))
font = FontProperties()
font.set_name('Arial')

class hamu(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(hamu, self).__init__(parent)
        self.setWindowTitle("Harmonic Analysis")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        self.harmnum_tag = QtWidgets.QLabel()
        self.harmnum_tag.setText('Highest Harmonic to Include')
        self.harmnum = QtWidgets.QTextEdit()
        self.harmnum.setFixedSize(150,35)
        self.harmnum.setText('1')
        self.harmnum.setToolTip('iFAMS selects harmonics up to the entered value.\n'
                                'The fundamentals are considered the first harmonic. ')

        self.button = QtWidgets.QPushButton('Run iFAMS Analysis')
        self.button.clicked.connect(self.harmav)

        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.harmnum_tag)
        self.layout.addWidget(self.harmnum)
        self.layout.addWidget(self.button)
        self.setLayout(self.layout)
    
    def harmav(self):
        ovnum = np.array(range(1,int(self.harmnum.toPlainText())+1))
        MW.harmnum = int(self.harmnum.toPlainText())
        MW.IFT = FT.higher_harmonic(MW.xFT,MW.yFT,MW.submass,MW.cs,MW.yint,MW.po2,ovnum)
        MW.tIFT = np.absolute(MW.IFT)
        MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(MW.tIFT, MW.xint, MW.cs,MW.proton)

        MW.GaborFraction = len(ovnum)*len(MW.cs)*(1/MW.submass)/(max(MW.xFT)/2)

        MW.fig.clf()
        MW.ax = MW.fig.add_subplot(121)
        MW.ax2 = MW.fig.add_subplot(122)
        MW.ax.plot(MW.x, MW.y, label='Original Spectrum')
        for i in range(len(MW.tIFT)):
            label = 'Charge State ' + str(int(MW.cs[i]))
            MW.ax.plot(MW.xint, MW.tIFT[i], label=label)
        MW.ax.set_title('Mass Spectrum')
        MW.ax.set_ylabel('Relative Abundance')
        MW.ax.set_xlabel('m/z')
        MW.ax.legend(loc='upper right')

        MW.ax2.plot(MW.xzero, MW.zerofull, label='Full Zero',color='mediumseagreen')
        for i in range(len(MW.cszero)):
            label = 'Charge State ' + str(int(MW.cs[i]))
            MW.ax2.plot(MW.xzero, MW.cszero[i], label=label)
        MW.ax2.set_title('Zero Charge Spectrum')
        MW.ax2.set_ylabel('Relative Abundance')
        MW.ax2.set_xlabel('Mass (Da)')
        MW.ax2.legend(loc='upper right')
        MW.fig.tight_layout()
        MW.canvas.draw()

class defect_analysis(QtWidgets.QDialog):
    """
    Conducts Phase Analysis for 1D Fourier Transform selections
    """
    def __init__(self, parent=None):
        super(defect_analysis, self).__init__(parent)
        self.setWindowTitle("Parameters for Mass Defect Analysis")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.harmnum_tag = QtWidgets.QLabel()
        self.harmnum_tag.setText('Highest Harmonic to Include')
        self.harmnum = QtWidgets.QTextEdit()
        self.harmnum.setText('1')
        self.harmnum.setFixedSize(150,35)
        self.harmnum.setTabChangesFocus(True)
        self.harmnum.setToolTip('iFAMS selects harmonics up to the entered value.\n'
                                'The fundamentals are considered the first harmonic. ')
        
        self.cb = QtWidgets.QCheckBox('Baseline Correction')

        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.harmnum_tag)
        self.layout.addWidget(self.harmnum)
        self.layout.addWidget(self.cb)
        self.setLayout(self.layout)

        self.button = QtWidgets.QPushButton('Run Mass Defect Analysis')
        self.button.clicked.connect(self.run_phase)
        self.layout.addWidget(self.button)
    
    def run_phase(self):
        """
        runs iFAMS molecular mass defect analysis
        """
##########################
        try:
            """
            MW.IFT,FTcentroid,FTstddev = FT.inverseFT(MW.xFT,MW.yFT,MW.submass,MW.cs,MW.y,MW.po2)
            MScentroid = np.dot(MW.xint,MW.yint)/sum(MW.yint)
            rdiff = 0
            for j in range(len(MW.yint)):
                rdiff+=MW.yint[j]*(MW.xint[j]-MScentroid)**2
            MSstddev = np.sqrt(rdiff/((len(MW.yint)-1)*sum(MW.yint)/len(MW.yint)))
            print(FTcentroid, FTstddev, MScentroid, MSstddev)
            for j in range(len(FTstddev)):
                slope = -FTstddev[j]/MSstddev
                print(slope, FTstddev[0]/FTstddev[j], abs((FTcentroid[j])**2/slope))
            """
            ###reestablished MW.IFT calculation to allow for MMD analysis before FT deconvolution
            ###changed MW.y to MW.yint
            MW.IFT = FT.inverseFT(MW.xFT,MW.yFT,MW.submass,MW.cs,MW.yint,MW.po2)
            MW.tIFT = np.absolute(MW.IFT)
            MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(MW.tIFT, MW.xint, MW.cs, MW.proton)
            
            numpoints = len(MW.xint)
            mzmin = min(MW.xint)
            ###mzspacing = np.round((max(MW.xint)-mzmin)/numpoints, decimals=4)
            mzspacing = MW.xint[1] - MW.xint[0]
            numharm = int(str(self.harmnum.toPlainText()))
            MW.reconmzint, MW.recontot, MW.reconsub, MW.reconmzsub = phase.phase(MW.xint, MW.yint, MW.tIFT, MW.submass, MW.cs, numpoints, numharm, mzspacing, mzmin, MW.proton)
            
            if self.cb.isChecked():
                MW.basecor = True
                
        ###baseline correction before charge state summation
                MW.recontot1 = np.zeros(int(MW.submass), dtype=complex)
                for i in range(len(MW.cs)):
                    MW.reconint = np.interp(MW.reconmzint,MW.reconmzsub[i],MW.reconsub[i],period=MW.submass)
                    baseline = phase.line_baseline(MW.reconint,MW.reconmzint,5)
                    MW.reconint = MW.reconint - baseline
                    MW.recontot1+=MW.reconint
                
        ###baseline correction after charge state summation 
                baseline = phase.line_baseline(MW.recontot,MW.reconmzint,5)
                MW.recontot = MW.recontot-baseline
        ###integration and identification of mass defect peaks 
            """
            MW.MMDcent,MW.MMDstddev,MW.MMDfwhm = phase.cyclic_int(MW.reconmzint,MW.recontot,MW.submass)
            MW.MMDcent = np.remainder(MW.MMDcent,MW.submass)*MW.submass
            height = np.interp(MW.MMDcent,MW.reconmz,MW.recontot)
            """
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)
            MW.ax.plot(MW.x, MW.y, label = 'Original Spectrum')
            for i in range(len(MW.tIFT)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax.plot(MW.xint,MW.tIFT[i], label=label)
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')
            MW.ax.legend(loc = 'upper right')

            MW.ax2.plot(MW.reconmzint, MW.recontot, label = 'Full Zero')
            ###MW.ax2.plot.vlines(centroid,0,max(y),colors='k')
            '''
            if self.cb.isChecked():
                MW.ax2.plot(MW.reconmzint, MW.recontot1, label = 'full zero*')
            '''
            for i in range(len(MW.cs)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax2.plot(MW.reconmzsub[i], MW.reconsub[i], label = label)
            
            MW.ax2.set_title('Molecular Mass Defect Profile')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.ax2.legend(loc = 'upper right')
            MW.fig.tight_layout()
            MW.canvas.draw()

        except AttributeError:
            print('No charge states or subunit mass found')
            print('Please either calculate charge states and subunit mass or manually enter them')
            MW.message = 'No charge states or subunit mass found'
            MW.submessage = 'Please either calculate charge states and subunit mass or manually enter them'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

class g_defect_analysis(QtWidgets.QDialog):
    """
    Conducts Phase Analysis for Gabor Transform selections
    currently must run prior to harmonic selection
    """
    def __init__(self, parent=None):
        super(g_defect_analysis, self).__init__(parent)
        self.setWindowTitle("Parameters for Mass Defect Analysis")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.harmnum_tag = QtWidgets.QLabel()
        self.harmnum_tag.setText('Highest Harmonic to Include')
        self.harmnum = QtWidgets.QTextEdit()
        self.harmnum.setText('1')
        self.harmnum.setFixedSize(150,35)
        self.harmnum.setTabChangesFocus(True)
        self.harmnum.setToolTip('iFAMS selects harmonics up to the entered value.\n'
                                'The fundamentals are considered the first harmonic. ')
        
        self.submass_tag = QtWidgets.QLabel()
        self.submass_tag.setText('Repeated Subunit Mass')
        self.submass_tag.setFixedSize(300,35)
        self.submass = QtWidgets.QTextEdit()
        self.submass.setText('0')
        self.submass.setFixedSize(150,35)
        self.submass.setTabChangesFocus(True)
        
        self.newlist =[]
        self.newlist2=[]
        
        for i in range(0,len(MW.glist)):
            self.newlist.append('self.selection'+str(i))
            self.newlist2.append('self.selnum'+str(i))
        for i in range(0,len(self.newlist)):
            label = 'Selection ' + str(i+1) +' Charge State'
            self.newlist[i] = QtWidgets.QLabel()
            self.newlist[i].setText(label)
            self.newlist[i].setFixedSize(300,25)
            self.newlist2[i] = QtWidgets.QTextEdit()
            self.newlist2[i].setFixedSize(150,35)
            self.newlist2[i].setTabChangesFocus(True)
            try:
                self.newlist2[i].setText(str((MW.zlist[i])))
            except AttributeError:
                continue
            except IndexError:
                continue

        self.cb = QtWidgets.QCheckBox('Baseline Correction')
        self.button = QtWidgets.QPushButton('Run Mass Defect Analysis')
        self.button.clicked.connect(self.run_phase)
        
        layout = QtWidgets.QGridLayout()
        count = 0
        for i in range(0,len(self.newlist)):
            if i <= 7+8*count:
                    if count == 0:
                        layout.addWidget(self.newlist[i],i*2,0)
                        layout.addWidget(self.newlist2[i],i*2+1,0)
                    else:
                        layout.addWidget(self.newlist[i],i*2-(16*count),count)
                        layout.addWidget(self.newlist2[i],i*2-(16*count-1),count)
            else:
                count+=1
                layout.addWidget(self.newlist[i],i*2-(16*count),count)
                layout.addWidget(self.newlist2[i],i*2-(16*count-1),count)
                
        layout.addWidget(self.submass_tag,16,0)
        layout.addWidget(self.submass,17,0)
        layout.addWidget(self.harmnum_tag,18,0)
        layout.addWidget(self.harmnum,19,0)
        layout.addWidget(self.cb,20,0)
        layout.addWidget(self.button,21,0)
            
        self.setLayout(layout)
    
    def run_phase(self):
        """
        runs iFAMS molecular mass defect analysis
        """
##########################
        try:    
            MW.cs = []
            for i in range(0,len(self.newlist2)):
                MW.cs.append(int(str(self.newlist2[i].toPlainText())))
            MW.submass = float(str(self.submass.toPlainText()))
            """
            MW.IFT, FTcentroid, FTstddev = FT.inverseFT(MW.xFT,MW.yFT,MW.submass,MW.cs,MW.y,MW.po2)
            MScentroid = np.dot(MW.xint,MW.yint)/sum(MW.yint)
            rdiff = 0
            for j in range(len(MW.yint)):
                rdiff+=MW.yint[j]*(MW.xint[j]-MScentroid)**2
            MSstddev = np.sqrt(rdiff/((len(MW.yint)-1)*sum(MW.yint)/len(MW.yint)))
            print(FTcentroid, FTstddev, MScentroid, MSstddev)
            for j in range(len(FTstddev)):
                slope = FTstddev[j]/MSstddev
                print(abs((FTcentroid[j])**2/slope))
            """
            MW.tIFT = np.absolute(MW.glist)
            MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(MW.tIFT, MW.xint, MW.cs, MW.proton)
            
            numpoints = len(MW.xint)
            mzmin = min(MW.xint)
            ###mzspacing = np.round((max(MW.xint)-mzmin)/numpoints, decimals=4)
            mzspacing = MW.xint[1] - MW.xint[0]
            numharm = int(str(self.harmnum.toPlainText()))
            MW.reconmzint, MW.recontot, MW.reconsub, MW.reconmzsub = phase.phase(MW.xint, MW.yint, MW.tIFT, MW.submass, MW.cs, numpoints, numharm, mzspacing, mzmin, MW.proton)
            if self.cb.isChecked():
                MW.basecor = True
                baseline = phase.line_baseline(MW.recontot,MW.reconmzint,5)
                MW.recontot = MW.recontot-baseline
            """
            MW.MMDcent,MW.MMDstddev,MW.MMDfwhm = phase.cyclic_int(MW.reconmzint,MW.recontot,MW.submass)
            MW.MMDcent = np.remainder(MW.MMDcent,MW.submass)*MW.submass
            height = np.interp(MW.MMDcent,MW.reconmz,MW.recontot)
            """
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)
            MW.ax.plot(MW.x, MW.y, label = 'Original Spectrum')
            for i in range(len(MW.tIFT)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax.plot(MW.xint,MW.tIFT[i], label=label)
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')
            MW.ax.legend(loc = 'upper right')

            MW.ax2.plot(MW.reconmzint, MW.recontot, label = 'Full Zero')
            ###MW.ax2.plot.vlines(MW.MMDcent,0,max(MW.recontot),colors = 'k')

            for i in range(len(MW.cs)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax2.plot(MW.reconmzsub[i], MW.reconsub[i], label = label)
            
            MW.ax2.set_title('Molecular Mass Defect Profile')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.ax2.legend(loc = 'upper right')
            MW.fig.tight_layout()
            MW.canvas.draw()

        except AttributeError:
            print('No charge states or subunit mass found')
            print('Please either calculate charge states and subunit mass or manually enter them')
            MW.message = 'No charge states or subunit mass found'
            MW.submessage = 'Please either calculate charge states and subunit mass or manually enter them'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()


class integrate(QtWidgets.QDialog):
    """
    This builds the iFAMS integration parameters and deconvolution baseline corrections menu
    """

    def __init__(self, parent=None):
        super(integrate, self).__init__(parent)
        self.setWindowTitle("Parameters for Integration")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))

        self.minx_tag = QtWidgets.QLabel()
        self.minx_tag.setText('Minimum x')
        self.minx_tag.setToolTip('The x-value the program will start looking for\n' 
                                'peaks to integrate')
        
        self.maxx_tag = QtWidgets.QLabel()
        self.maxx_tag.setText('Maximum x')
        self.maxx_tag.setToolTip('The x-value the program will stop looking for\n' 
                                'peaks to integrate')

        self.miny_tag = QtWidgets.QLabel()
        self.miny_tag.setText('Minimum y')
        self.miny_tag.setToolTip('This y-value is the smallest peak height\n' 
                                'the program will consider an integratable peak')

        self.delta_tag = QtWidgets.QLabel()
        self.delta_tag.setText('Minimum Distance Between Peaks')
        self.delta_tag.setToolTip('This value is the minimum distance between possible\n' 
                                'peaks needed for the program to consider them\n'
                                'two different peaks, integrating them seperately')

        self.ntol_tag = QtWidgets.QLabel()
        self.ntol_tag.setText('Noise Tolerance for Peaks')
        self.ntol_tag.setToolTip('This value is the minimum relative height for \n' 
                                'a neighboring local maximum to be considered a new \n'
                                'peak. A larger value will tolerate more noise in a peak. \n'
                                'Typical value between 0 and 0.1')

        self.minx = QtWidgets.QTextEdit()
        self.minx.setFixedSize(150,35)
        self.minx.setTabChangesFocus(True)
        self.maxx = QtWidgets.QTextEdit()
        self.maxx.setFixedSize(150,35)
        self.maxx.setTabChangesFocus(True)
        try:
            self.minx.setText(str(np.round(MW.minimumX,0)))
            self.maxx.setText(str(np.round(MW.maximumX,0)))
        except AttributeError:
            pass
        self.miny = QtWidgets.QTextEdit()
        self.miny.setFixedSize(150,35)
        self.miny.setTabChangesFocus(True)
        self.delta = QtWidgets.QTextEdit()
        self.delta.setFixedSize(150,35)
        self.delta.setTabChangesFocus(True)
        self.ntol = QtWidgets.QTextEdit()
        self.ntol.setFixedSize(150,35)
        self.ntol.setTabChangesFocus(True)
        try:
            self.miny.setText(str(MW.minimumY))
            self.delta.setText(str(MW.deltaY))
            self.ntol.setText(str(MW.ntol))
        except AttributeError:
            self.miny.setText(str(0.01))
            self.delta.setText(str(10))
            self.ntol.setText(str(0.005))
    
        self.button1 = QtWidgets.QPushButton('Perform Peak Integration')
        self.button1.clicked.connect(self.integrate)
        
        self.button2 = QtWidgets.QPushButton('Perform Simple Integration')
        self.button2.clicked.connect(self.simpInt)
        self.button2.setToolTip('This integrates over the zero-charge spectrum within\n'
                                'given minimum and maximum x values entered above.\n'
                                'Other parameters are not needed.')
        
        self.separator = QtWidgets.QLabel()
        self.separator.setText('Baseline and Resolution Tools: ')
        
        self.button3 = QtWidgets.QPushButton('Smooth Isotope Resolution')
        self.button3.clicked.connect(self.gaussian_smooth)
        self.button3.setToolTip('Convolves the zero-charge spectrum by a gaussian \n'
                                'with a full-width at half-max of 1.8 Da \n'
                                'to smooth out isotope modulations prior to integration ')
        
        self.buttonMBC = QtWidgets.QPushButton('Minima Baseline Correction')
        self.buttonMBC.clicked.connect(self.Min_baseline_sub)
        self.buttonMBC.setToolTip('This moves the baseline of the spectrum up to \n'
                                'zero by drawing a curve between minima below zero \n'
                                'and subtracting the curve from the spectrum.')
        
        self.button4 = QtWidgets.QPushButton('Fourier Baseline Correction')
        self.button4.clicked.connect(self.NZ_baseline_sub)
        self.button4.setToolTip('This moves the base of each peak to the baseline \n'
                                'by subtracting half of the Near-Zero frequency \n'
                                'from each charge state.')
        
        self.button5 = QtWidgets.QPushButton('Segmented Baseline Correction')
        self.button5.clicked.connect(self.baseline_sub)
        self.button5.setToolTip('This moves the base of each peak to the baseline \n'
                                'by subtracting the lines between adjacent minima \n'
                                'and subtracts the removed area from the integration.\n'
                                'Requires peak integration')

        self.button6 = QtWidgets.QPushButton('Undo Adjustments')
        self.button6.clicked.connect(self.sub_undo)
        
                
        self.minx2_tag = QtWidgets.QLabel()
        self.minx2_tag.setText('Minimum x for 2nd range')
        self.minx2_tag.setToolTip('the x value the program will start integrating. \n'
                                  'Optional 2nd range for simple integration only')
        
        self.maxx2_tag = QtWidgets.QLabel()
        self.maxx2_tag.setText('Maximum x for 2nd range')
        self.maxx2_tag.setToolTip('the x value the program will stop integrating. \n'
                                  'Optional 2nd range for simple integration only')

        self.minx2 = QtWidgets.QTextEdit()
        self.minx2.setFixedSize(150,35)
        self.minx2.setTabChangesFocus(True)
        self.maxx2 = QtWidgets.QTextEdit()
        self.maxx2.setFixedSize(150,35)
        self.maxx2.setTabChangesFocus(True)

        self.button7 = QtWidgets.QPushButton('Toggle Charge-state Spectra')
        self.button7.clicked.connect(self.CS_toggle)
        self.button7.setToolTip('Adds or removes the individual charge-state zero-charge spectra \n'
                                'that are added up for the full zero-charge spectrum (or combined \n'
                                'deconvolution), if possible.')

        self.window_tag = QtWidgets.QLabel()
        self.window_tag.setText('Batch Integration Type ')
        self.window = QtWidgets.QComboBox()
        self.window.addItem('Select type')
        self.window.addItem('Parameter Match')
        self.window.addItem('Bounds Match')
        self.window.activated[str].connect(self.batch_integration)
        self.window.setToolTip('These types are used for batching with peak integration. \n'
                               'Parameter match will use the same parameters entered here \n'
                               'to search for peaks in each spectrum batched. \n'
                               'Bounds match will use the bounds of integration determined \n'
                               'with this spectrum and will integrate each batched spectrum \n'
                               'over the same bounds.')

        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.minx_tag,1,0)
        self.layout.addWidget(self.minx,2,0)
        self.layout.addWidget(self.maxx_tag,3,0)
        self.layout.addWidget(self.maxx,4,0)
        self.layout.addWidget(self.miny_tag,5,0)
        self.layout.addWidget(self.miny,6,0)
        self.layout.addWidget(self.delta_tag,7,0)
        self.layout.addWidget(self.delta,8,0)
        self.layout.addWidget(self.ntol_tag,9,0)
        self.layout.addWidget(self.ntol,10,0)
        self.layout.addWidget(self.button1,11,0,1,2)
        self.layout.addWidget(self.button2,12,0,1,2)
        self.layout.addWidget(self.separator,13,0,1,2)
        self.layout.addWidget(self.button3,14,0,1,2)
        self.layout.addWidget(self.buttonMBC,15,0,1,2)
        self.layout.addWidget(self.button4,16,0,1,2)
        self.layout.addWidget(self.button5,17,0,1,2)
        self.layout.addWidget(self.button6,18,0,1,2)
        self.layout.addWidget(self.minx2_tag,1,1)
        self.layout.addWidget(self.minx2,2,1)
        self.layout.addWidget(self.maxx2_tag,3,1)
        self.layout.addWidget(self.maxx2,4,1)
        self.layout.addWidget(self.button7,8,1)
        self.layout.addWidget(self.window_tag,9,1)
        self.layout.addWidget(self.window,10,1)
        
        self.window.setCurrentText('Bounds Match')
        self.batch_integration('Bounds Match')
        self.CSshown = True
        
        self.setLayout(self.layout)
        
    def integrate(self):
        MW.minimumX = float(str(self.minx.toPlainText()))
        MW.maximumX = float(str(self.maxx.toPlainText()))
        MW.minimumY = float(str(self.miny.toPlainText()))
        MW.deltaY = float(str(self.delta.toPlainText()))
        MW.ntol = float(str(self.ntol.toPlainText()))
        
        MW.ynorm = MW.zerofull/max(MW.zerofull)
        massrange = np.nonzero(MW.zerofull)
        MW.xlist,MW.ylist,MW.xref,MW.yref,MW.sumlist,MW.peakindex,MW.minL,MW.minR,MW.xcent,MW.height = Int.integrate(MW.ynorm[massrange],MW.zerofull[massrange], 
                                                                                            MW.xzero[massrange], MW.deltaY, MW.minimumY,
                                                                                            MW.minimumX, MW.maximumX, MW.ntol)
        max_value = max(MW.zerofull)
        for i in range(0,len(MW.sumlist)):
            MW.sumlist[i] = max_value*MW.sumlist[i]
            MW.height[i] = max_value*MW.height[i]
         
        try:
            MW.integratedx = 0  
            MW.peakwidths = np.zeros(len(MW.xref))
            for i in range(0, len(MW.xref)):
                MW.integratedx += len(MW.xref[i])
                MW.peakwidths[i] = len(MW.xref[i])
            MW.xpoints = MW.integratedx
            MW.IntFraction = MW.integratedx/len(MW.xzero[massrange])
            MW.N = np.sqrt(MW.IntFraction)*np.sqrt(MW.GaborFraction)*MW.NFrms*2*np.pi
            MW.Nh = np.sqrt(MW.GaborFraction)*MW.NFrms/np.sqrt(len(MW.xzero[massrange]))*2*np.pi
            if MW.NFrms > 0:
                print('Frequency Noise RMSD: ' +str(MW.NFrms))
                print('Mass Noise RMSD: ' +str(MW.Nrms))
                print('Fraction of Gabor used: ' + str(MW.GaborFraction))
                print('Fraction of Zero Charge Spectrum integrated: ' + str(MW.IntFraction))
                print('Noise for integrated spectrum: ' + str(MW.N))
        except AttributeError:
            print('Unable to perform noise analysis')
        
        if MW.based == True:
            self.baseline_sub()
        else:
            self.int_show()
        
        if MW.ReInt == False:
            MW.ReInt = True  
        if MW.paramMatch == False:
            MW.bounds = []
            for i in range(0,len(MW.xref)):
                MW.bounds.append(MW.xref[i][0])
                MW.bounds.append(MW.xref[i][-1])
      
    def gaussian_smooth(self):
        try:
            try:
                MW.zerofullstore
            except AttributeError:
                MW.zerofullstore = MW.zerofull.copy()
            MW.smoothed = True
            
            MW.zerofull = Int.smooth(MW.xzero,MW.zerofull,1.8)
            
            MW.ax2.clear()
            MW.zero_norm = MW.zerofull/max(MW.zerofullstore)
            MW.zero_norm2 = MW.zerofullstore/max(MW.zerofullstore)
            MW.ax2.plot(MW.xzero,MW.zero_norm2, label='Original Spectrum')
            MW.ax2.plot(MW.xzero,MW.zero_norm, label='Smoothed Spectrum')
            MW.ax2.set_title('Zero Charge Spectrum')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.ax2.legend(loc = 'upper right')
            MW.fig.tight_layout()
            MW.canvas.draw()
            
        except AttributeError:
            print('Failed to apply Gaussian convolution')
            MW.message = 'Failed to apply Gaussian convolution'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def baseline_sub(self):
        try:
            try:
                MW.zerofullstore
            except AttributeError:
                MW.zerofullstore = MW.zerofull.copy()
            MW.yrefstore = MW.yref.copy()
            MW.sumliststore = MW.sumlist.copy()
            MW.heightstore = MW.height.copy()
            
            MW.yref2, MW.sumlist, MW.height = Int.baseline_sub2(MW.xzero, MW.zerofull, MW.xref, 
                                                                 MW.yref, MW.peakindex)            
             
            MW.based = True
            self.int_show() 
        except IndexError:
            MW.minimumX = float(str(self.minx.toPlainText()))
            MW.maximumX = float(str(self.maxx.toPlainText()))
            MW.minimumY = float(str(self.miny.toPlainText()))
            MW.deltaY = float(str(self.delta.toPlainText()))
            MW.ntol = float(str(self.ntol.toPlainText()))
            
            MW.ynorm = MW.zerofull/max(MW.zerofull)
            massrange = np.nonzero(MW.zerofull)
            MW.xlist,MW.ylist,MW.xref,MW.yref,MW.sumlist,MW.peakindex,MW.minL,MW.minR,MW.xcent,MW.height = Int.integrate(MW.ynorm[massrange],MW.zerofull[massrange], 
                                                                                                MW.xzero[massrange], MW.deltaY, MW.minimumY,
                                                                                                MW.minimumX, MW.maximumX, MW.ntol)
            max_value = max(MW.zerofull)
            for i in range(0,len(MW.sumlist)):
                MW.sumlist[i] = max_value*MW.sumlist[i]
                MW.height[i] = max_value*MW.height[i]
            MW.zerofullstore = MW.zerofull.copy()
            MW.yrefstore = MW.yref.copy()
            MW.sumliststore = MW.sumlist.copy()
            MW.heightstore = MW.height.copy()
            
            MW.yref2, MW.sumlist, MW.height = Int.baseline_sub2(MW.xzero, MW.zerofull, MW.xref, 
                                                                 MW.yref, MW.peakindex)            
             
            MW.based = True
            self.int_show() 
        except AttributeError:
            print('unable to perform Segmented Baseline Correction. Try performing peak integration first')
            MW.message = 'Unable to perform Segmented Baseline Correction. \nTry performing peak integration first'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

    def Min_baseline_sub(self):
        try:
            MW.zerofullstore
        except AttributeError:
            MW.zerofullstore = MW.zerofull.copy()
        MW.baseline = STFT.line_baseline(MW.zerofull, MW.xzero, 5)
        MW.zerofull = MW.zerofull - MW.baseline
        MW.basecor = True

        MW.ax2.clear()
        MW.ax2.plot(MW.xzero,MW.zerofullstore, label='Before Subtraction',color='darkgray')
        MW.ax2.plot(MW.xzero,MW.baseline, label='Subtracted Curve',color='orange')
        MW.ax2.plot(MW.xzero,MW.zerofull, label='After Subtraction',color='mediumseagreen')
        MW.ax2.legend(loc='upper right', title='Charge State')
        MW.ax2.set_title('Zero Charge Spectrum')
        MW.ax2.set_ylabel('Relative Abundance')
        MW.ax2.set_xlabel('Mass (Da)')
        MW.fig.tight_layout()
        MW.canvas.draw()

    def NZ_baseline_sub(self):
        try:
            try:
                MW.zerofullstore
            except AttributeError:
                MW.zerofullstore = MW.zerofull.copy()
            MW.gliststore = MW.glist.copy()
            MW.NZbased = True
            
            print('Selecting near-zero frequencies')
            if len(MW.batchstore) > 0:
                piW = MW.xshort[2]-MW.xshort[1]
                piH = MW.yshort[2]-MW.yshort[1]
                MW.NZrex = []
                nzrex = []
                MW.NZcs = []
                nzcs = []
                zcount = -1
                zcount2 = 0
                for i in range(0,len(MW.rlist)):
                    if zcount > MW.rlist[i][4]:
                        zcount2 += zcount+1
                        zcount = -1
                        MW.NZrex.append(nzrex)
                        MW.NZcs.append(nzcs)
                        nzrex = []
                        nzcs = []
                    if zcount == MW.rlist[i][4]:
                        continue
                    else:
                        nzrex.append([int(np.round((MW.rlist[i][0]-MW.xshort[0])/piW)),int(np.round((MW.rlist[i][1]-MW.xshort[0])/piW)),0, int(np.round((MW.rlist[i][3]-((MW.rlist[i][2]+MW.rlist[i][3])/2))/piH)), MW.rlist[i][4]])
                        zcount += 1
                        nzcs.append(MW.cs[zcount+zcount2])
                MW.NZrex.append(nzrex)
                MW.NZcs.append(nzcs)
                MW.Xtemp = MW.X*0
                MW.XtempBL2 = MW.XBL2*0
                MW.NZrecon = []
                MW.NZreconBL2 = []
                MW.NZzerofull = []
                for i in range(len(MW.NZrex)):
                    MW.NZrecon.append([])
                    MW.NZreconBL2.append([])
                    for j in range(len(MW.NZrex[i])):
                        MW.Xtemp[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1] = MW.X[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1]
                        MW.Xtemp[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])] = MW.X[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])]
                        MW.XtempBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1] = MW.XBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1]
                        MW.XtempBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])] = MW.XBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])]
                        MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                        MW.NZrecon[i].append(MW.yrecon/2)
                        MW.Xtemp = MW.X*0
                        MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                        MW.NZreconBL2[i].append(MW.yreconBL2/2)
                        MW.XtempBL2 = MW.XBL2*0
                
                    MW.NZxzero, MW.NZcszero, NZzerofull = FT.zerocharge(np.real(MW.NZrecon[i]), MW.xint, MW.NZcs[i],MW.proton)
                    MW.NZxzeroBL2, MW.NZcszeroBL2, MW.NZzerofullBL2 = FT.zerocharge(np.real(MW.NZreconBL2[i]), MW.xint, MW.NZcs[i],MW.proton)
                    massrange=np.nonzero(MW.NZzerofullBL2)
                    MW.NZxzero = MW.NZxzero[massrange]
                    NZzerofull = NZzerofull[massrange]
                    MW.NZxzeroBL2 = MW.NZxzeroBL2[massrange]
                    MW.NZzerofullBL2 = MW.NZzerofullBL2[massrange]
                    check = np.max(MW.NZzerofullBL2)/4
                    for j in range(len(MW.NZxzero)):
                        if MW.NZzerofullBL2[j] > check:
                            mid1 = j
                            break
                    for j in range(mid1+1,len(MW.NZxzero)):
                        if MW.NZzerofullBL2[j] <= check or MW.NZzerofullBL2[j] <= MW.NZzerofullBL2[mid1]:
                            mid2 = j
                            break
                    int1 = min(sum(NZzerofull[0:mid1]),sum(NZzerofull[mid2:len(MW.NZxzero)]))
                    int2 = (sum(MW.NZzerofullBL2[0:mid1])+sum(MW.NZzerofullBL2[mid2:len(MW.NZxzero)]))/2
                    MW.NZzerofullBL2 = MW.NZzerofullBL2*int1/int2
                    MW.NZzerofull.append(np.interp(MW.deconxlist[i],MW.NZxzero,NZzerofull - MW.NZzerofullBL2))
                    MW.deconylist[i] = np.asarray(MW.deconylist[i]) - np.asarray(MW.NZzerofull[i])
                y0s = []
                NZ0s = []
                for i in range(len(MW.deconylist)):
                    y0s.append(np.interp(MW.xzero,MW.deconxlist[i],MW.deconylist[i]))
                    NZ0s.append(np.interp(MW.xzero,MW.deconxlist[i],MW.NZzerofull[i]))
                MW.zerofull = np.sum(y0s,axis=0)
                MW.NZzerofull = np.sum(NZ0s,axis=0)
                MW.zero_norm = MW.zerofullstore/max(MW.zerofull)
            
            else:
                piW = MW.xshort[2]-MW.xshort[1]
                piH = MW.yshort[2]-MW.yshort[1]
                MW.NZrex = []
                zcount = -1
                for i in range(0,len(MW.rlist)):
                    if zcount == MW.rlist[i][4]:
                        continue
                    else:
                        MW.NZrex.append([int(np.round((MW.rlist[i][0]-MW.xshort[0])/piW)),int(np.round((MW.rlist[i][1]-MW.xshort[0])/piW)),0, int(np.round((MW.rlist[i][3]-((MW.rlist[i][2]+MW.rlist[i][3])/2))/piH)), MW.rlist[i][4]])
                        zcount += 1
                        if MW.rlist[i][4] < zcount:
                            zcount = 0
                MW.Xtemp = MW.X*0
                MW.XtempBL2 = MW.XBL2*0
                MW.NZrecon = []
                MW.NZreconBL2 = []
                for j in range(len(MW.NZrex)):
                    MW.Xtemp[MW.NZrex[j][0]:MW.NZrex[j][1]+1,0:MW.NZrex[j][3]+1] = MW.X[MW.NZrex[j][0]:MW.NZrex[j][1]+1,0:MW.NZrex[j][3]+1]
                    MW.Xtemp[MW.NZrex[j][0]:MW.NZrex[j][1]+1,-MW.NZrex[j][3]:len(MW.X[0])] = MW.X[MW.NZrex[j][0]:MW.NZrex[j][1]+1,-MW.NZrex[j][3]:len(MW.X[0])]
                    MW.XtempBL2[MW.NZrex[j][0]:MW.NZrex[j][1]+1,0:MW.NZrex[j][3]+1] = MW.XBL2[MW.NZrex[j][0]:MW.NZrex[j][1]+1,0:MW.NZrex[j][3]+1]
                    MW.XtempBL2[MW.NZrex[j][0]:MW.NZrex[j][1]+1,-MW.NZrex[j][3]:len(MW.X[0])] = MW.XBL2[MW.NZrex[j][0]:MW.NZrex[j][1]+1,-MW.NZrex[j][3]:len(MW.X[0])]
                    MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                    MW.NZrecon.append(MW.yrecon/2)
                    MW.Xtemp = MW.X*0
                    MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                    MW.NZreconBL2.append(MW.yreconBL2/2)
                    MW.XtempBL2 = MW.XBL2*0
                
                MW.NZxzero, MW.NZcszero, MW.NZzerofull = FT.zerocharge(np.real(MW.NZrecon), MW.xint, MW.cs,MW.proton)
                MW.NZxzeroBL2, MW.NZcszeroBL2, MW.NZzerofullBL2 = FT.zerocharge(np.real(MW.NZreconBL2), MW.xint, MW.cs,MW.proton)
                NZspec = np.sum(MW.NZrecon,axis=0)
                
                massrange=np.nonzero(MW.NZzerofullBL2)
                MW.NZxzero = MW.NZxzero[massrange]
                MW.NZzerofull = MW.NZzerofull[massrange]
                MW.NZxzeroBL2 = MW.NZxzeroBL2[massrange]
                MW.NZzerofullBL2 = MW.NZzerofullBL2[massrange]
                check = np.max(MW.NZzerofullBL2)/4
                for j in range(len(MW.NZxzero)):
                    if MW.NZzerofullBL2[j] > check:
                        mid1 = j
                        break
                for j in range(mid1+1,len(MW.NZxzero)):
                    if MW.NZzerofullBL2[j] <= check or MW.NZzerofullBL2[j] <= MW.NZzerofullBL2[mid1]:
                        mid2 = j
                        break
                int1 = min(sum(MW.NZzerofull[0:mid1]),sum(MW.NZzerofull[mid2:len(MW.NZxzero)]))
                int2 = (sum(MW.NZzerofullBL2[0:mid1])+sum(MW.NZzerofullBL2[mid2:len(MW.NZxzero)]))/2
                MW.NZzerofullBL2 = MW.NZzerofullBL2*int1/int2
                MW.NZzerofull = MW.NZzerofull - MW.NZzerofullBL2
                MW.zero_norm = MW.zerofull/np.max(MW.zerofull - MW.NZzerofull)
                MW.zerofull = MW.zerofull - MW.NZzerofull
             
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)
            MW.ax.plot(MW.xint,MW.yint,label='Mass Spectrum',color='darkgray')
            if len(MW.batchstore) == 0:
                for i in range(0,len(MW.cs)):
                    label = 'Selection' + str(MW.cs[i])
                    MW.ax.plot(MW.xint,MW.tIFT[i],label=label)
                try:
                    MW.ax.plot(MW.xint,NZspec, label='Halved Near-Zero Frequency',color='r')
                except AttributeError:
                    pass
            else:
                for i in range(0,len(MW.batchstore)):
                    NZspec = np.sum(MW.NZrecon[i],axis=0)
                    label = 'Series ' + str(i+1)
                    MW.ax.plot(MW.xint,MW.tIFT[i],label=label)
                    try:
                        MW.ax.plot(MW.xint,NZspec, label='Halved Near-Zero Frequency',color='r')
                    except AttributeError:
                        pass
            MW.ax.legend(loc='upper right', title='Charge State')
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Relative Abundance')    
    
            MW.ax2.plot(MW.xzero,MW.zero_norm, label='Before Subtraction')
            MW.ax2.plot(MW.xzero,MW.NZzerofull/max(MW.zerofull), label='Subtracted Curve')
            MW.ax2.plot(MW.xzero,MW.zerofull/max(MW.zerofull), label='After Subtraction')
            MW.ax2.legend(loc='upper right', title='Charge State')
            MW.ax2.set_title('Zero Charge Spectrum')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.fig.tight_layout()
            MW.canvas.draw()
        except AttributeError:
            print('Unable to find near-zero frequencies')
            MW.message = 'Unable to find near-zero frequencies to subtract'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            
    def sub_undo(self):
        MW.zerofull = MW.zerofullstore.copy()
        if MW.based == True:
            MW.yref = MW.yrefstore.copy()
            MW.sumlist = MW.sumliststore.copy()
            MW.height = MW.heightstore.copy()
            MW.based = False
            self.int_show()
            
        if MW.basecor == True or MW.NZbased == True:
            if MW.NZbased == True:
                MW.glist = MW.gliststore.copy()
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)
            MW.ax.plot(MW.xint,MW.yint,label='Mass Spectrum',color='darkgray')
            try:
                for i in range(0,len(MW.cs)):
                    label = 'Selection' + str(MW.cs[i])
                    MW.ax.plot(MW.xint,MW.tIFT[i],label=label)
                MW.ax.legend(loc='upper right', title='Charge State')
            except IndexError:
                pass
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Relative Abundance')
    
            MW.zero_norm = MW.zerofull/max(MW.zerofull)
            MW.ax2.plot(MW.xzero,MW.zero_norm, label='Zero Charge Spectrum')
            MW.ax2.set_title('Zero Charge Spectrum')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.fig.tight_layout()
            MW.canvas.draw()
         
            
        if MW.smoothed == True:  
            MW.ax2.clear()
            MW.zero_norm = MW.zerofull/max(MW.zerofull)
            MW.ax2.plot(MW.xzero,MW.zero_norm, label='Zero Charge Spectrum')
            MW.ax2.set_title('Zero Charge Spectrum')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.fig.tight_layout()
            MW.canvas.draw()
            
        MW.NZbased = False
        MW.smoothed = False
        MW.basecor = False
   
    def int_show(self):
            # plottop = 1.1
            # plotbott = -0.05
            # plotxL = MW.xzero[0]
            # plotxR = MW.xzero[-1]

        MW.time_stamp = str(datetime.datetime.now())

        normfact = np.max(MW.zerofull)
        MW.fig.clf()
        MW.ax2 = MW.fig.add_subplot(111)
        MW.graphinfo = ("Zero Charge Spectrum with Integrated Peaks \n"
                        "Total number of points in integration =" + str(MW.xpoints)+"\n"
                        "Y-axis normalized to "+str("{:e}".format(np.round(normfact)))+"\n"
                        +MW.time_stamp)
        MW.ax2.set_title(MW.graphinfo)
        MW.ax2.plot(MW.xzero, MW.zerofull/normfact,color='mediumseagreen')
        try:
            if self.CSshown == True:
                if len(MW.batchstore) == 0:
                    for i in range(len(MW.cszero)):
                        MW.ax2.plot(MW.xzero,MW.cszero[i]/normfact,alpha=0.5)
                        MW.ax2.text(MW.xzero[np.argmax(MW.cszero[i])],np.max(MW.cszero[i])/normfact,str(MW.cs[i]),color=color[i],clip_on=True,ha='center')
                else:
                    for i in range(len(MW.batchstore)):
                        for j in range(len(MW.batchstore[i][5])):
                            MW.ax2.plot(MW.batchstore[i][3],MW.batchstore[i][5][j]/normfact,color=color[j+1],alpha=0.5)
                            MW.ax2.text(MW.batchstore[i][3][np.argmax(MW.batchstore[i][5][j])],np.max(MW.batchstore[i][5][j])/normfact,str(MW.batchstore[i][0][j]),color=color[j+1],clip_on=True,ha='center')
                self.CSshown = True
        except AttributeError:
            pass
        except IndexError:
            pass
        except TypeError:
            pass
        MW.ax2.set_ylabel('Relative Abundance')
        MW.ax2.set_xlabel('Mass (Da)')        
        
        if MW.based == True:
            for i in range(0, len(MW.xref)):
                MW.ax2.fill_between(MW.xref[i], MW.yref[i]/normfact, MW.yref2[i]/normfact)
        else:
            for i in range(0, len(MW.xref)):
                MW.ax2.fill_between(MW.xref[i], MW.yref[i]/normfact)

        try: 
            peaktxt = []
            for i in range(0, len(MW.sumlist)):
                peaktxt.append(i)
            for i in peaktxt:
                if MW.based == True:
                    try:
                        MW.ax2.text(MW.xcent[i],MW.heightstore[i]/normfact,str(np.round(MW.xcent[i],4)),rotation='vertical',ha='center',clip_on=True)
                    except AttributeError:
                        MW.ax2.text(MW.xcent[i],MW.height[i]/normfact,str(np.round(MW.xcent[i],4)),rotation='vertical',ha='center',clip_on=True)
                else:
                    MW.ax2.text(MW.xcent[i],MW.height[i]/normfact,str(np.round(MW.xcent[i],4)),rotation='vertical',ha='center',clip_on=True)
        except IndexError:
            pass

        # MW.ax2.set_ylim(plotbott,plottop)
        # MW.ax2.set_xlim(plotxL,plotxR)
        MW.masslist=[]
        MW.intlist=[]
        MW.heightlist=[]
        for i in range(len(MW.xcent)): #aligns values in table by decimal
            MW.masslist.append("%.4f" % MW.xcent[i])
            MW.intlist.append("%.4f" % MW.sumlist[i])
            MW.heightlist.append("%.4f" % MW.height[i])  
        try:   
            MW.intlist.append("%.4f" % MW.N)
            MW.heightlist.append("%.4f" % MW.Nh)
            MW.masslist.append('Noise')    
        except AttributeError:
            pass
            
        integrate.newtable(self)

        MW.fig.tight_layout()
        MW.canvas.draw() 
        
    def CS_toggle(self):
        try:
            if len(MW.cszero) == 0:
                MW.noname
            self.CSshown
            if self.CSshown == True:
                self.CSshown = False
            else:
                self.CSshown = True
            
            plottop = MW.ax2.get_ylim()[1]
            plotbott = MW.ax2.get_ylim()[0]
            plotxL = MW.ax2.get_xlim()[0]
            plotxR = MW.ax2.get_xlim()[1]
            normfact = np.max(MW.zerofull)
            MW.ax2.clear()
            MW.graphinfo = ("Zero Charge Spectrum with Integrated Peaks \n"
                            "Total number of points in integration =" + str(MW.xpoints)+"\n"
                            "Y-axis normalized to "+str("{:e}".format(np.round(normfact))))
            MW.ax2.set_title(MW.graphinfo)
            MW.ax2.plot(MW.xzero, MW.zerofull/normfact,color='mediumseagreen')
            try:
                if self.CSshown == True:
                    if len(MW.batchstore) == 0:
                        for i in range(len(MW.cszero)):
                            MW.ax2.plot(MW.xzero,MW.cszero[i]/normfact,color=color[i],alpha=0.5)
                            MW.ax2.text(MW.xzero[np.argmax(MW.cszero[i])],np.max(MW.cszero[i])/normfact,str(MW.cs[i]),color=color[i],clip_on=True,ha='center')
                    else:
                        for i in range(len(MW.batchstore)):
                            for j in range(len(MW.batchstore[i][5])):
                                MW.ax2.plot(MW.batchstore[i][3],MW.batchstore[i][5][j]/normfact,color=color[j+1],alpha=0.5)
                                MW.ax2.text(MW.batchstore[i][3][np.argmax(MW.batchstore[i][5][j])],np.max(MW.batchstore[i][5][j])/normfact,str(MW.batchstore[i][0][j]),color=color[j+1],clip_on=True,ha='center')
                    self.CSshown = True
            except AttributeError:
                pass
            except IndexError:
                pass
            except TypeError:
                pass
            MW.ax2.set_ylim(plotbott,plottop)
            MW.ax2.set_xlim(plotxL,plotxR)
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')        
            
            if MW.based == True:
                for i in range(0, len(MW.xref)):
                    MW.ax2.fill_between(MW.xref[i], MW.yref[i]/normfact, MW.yref2[i]/normfact)
            else:
                for i in range(0, len(MW.xref)):
                    MW.ax2.fill_between(MW.xref[i], MW.yref[i]/normfact)
    
            try: 
                peaktxt = []
                for i in range(0, len(MW.sumlist)):
                    peaktxt.append(i)
                for i in peaktxt:
                    if MW.based == True:
                        try:
                            MW.ax2.text(MW.xcent[i],MW.heightstore[i]/normfact,str(np.round(MW.xcent[i],4)),rotation='vertical',ha='center',clip_on=True)
                        except AttributeError:
                            MW.ax2.text(MW.xcent[i],MW.height[i]/normfact,str(np.round(MW.xcent[i],4)),rotation='vertical',ha='center',clip_on=True)
                    else:
                        MW.ax2.text(MW.xcent[i],MW.height[i]/normfact,str(np.round(MW.xcent[i],4)),rotation='vertical',ha='center',clip_on=True)
            except IndexError:
                pass
            except AttributeError:
                pass
            MW.fig.tight_layout()
            MW.canvas.draw() 
            
        except AttributeError:
            print('Unable to plot charge-state data')
            MW.message = 'Unable to plot charge-state data'
            MW.submessage = str('Either no charge-state data has been stored \n'
                                'or an integration needs to be performed first. ')
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def simpInt(self):
        try:
            MW.x0 = []
            MW.y0 = []
            MW.x02 = []
            MW.y02 = []
            MW.zeroWidth = []
            MW.IntWidth = []
            MW.IntWidth2 = []
            MW.minimumX = float(str(self.minx.toPlainText()))
            MW.maximumX = float(str(self.maxx.toPlainText()))
            try:
                MW.minimumX2 = float(str(self.minx2.toPlainText()))
                MW.maximumX2 = float(str(self.maxx2.toPlainText()))
                MW.range2 = True
            except ValueError:
                MW.range2 = False
    
            for i in range(0,len(MW.xzero)):
                if MW.zerofull[i] > 0:
                    MW.zeroWidth.append(MW.zerofull[i])
                if MW.xzero[i] >= MW.minimumX and MW.xzero[i] <= MW.maximumX:
                    MW.x0.append(MW.xzero[i])
                    MW.y0.append(MW.zerofull[i])
                if MW.zerofull[i] > 0 and MW.xzero[i] >= MW.minimumX and MW.xzero[i] <= MW.maximumX:              
                    MW.IntWidth.append(MW.xzero[i])
                
            MW.fullsum = 0
            MW.yref = [[]]
            MW.xref = [[]]
            for i in range(0, len(MW.y0)):
                MW.fullsum += MW.y0[i]
                MW.yref[0].append(MW.y0[i])
                MW.xref[0].append(MW.x0[i])
            xzero = MW.xzero.tolist()
            MW.peakindex = [xzero.index(MW.xref[0][np.argmax(MW.yref[0])])]
    
    
            sumnum = 0
            for i in range(0,len(MW.x0)):
                sumnum += MW.x0[i] * MW.y0[i]
            MW.height = [max(MW.y0)]
            MW.xcent = [float(sumnum / MW.fullsum )]
            MW.sumlist = [MW.fullsum]
    
            if MW.range2 == True:
                for i in range(0,len(MW.xzero)):
                    if MW.xzero[i] >= MW.minimumX2 and MW.xzero[i] <= MW.maximumX2:
                        MW.x02.append(MW.xzero[i])
                        MW.y02.append(MW.zerofull[i])
                    if MW.zerofull[i] > 0 and MW.xzero[i] >= MW.minimumX2 and MW.xzero[i] <= MW.maximumX2:              
                        MW.IntWidth2.append(MW.xzero[i])
            
                MW.fullsum2 = 0
                MW.yref.append([])
                MW.xref.append([])
                for i in range(0, len(MW.y02)):
                    MW.fullsum2 += MW.y02[i]
                    MW.yref[1].append(MW.y0[i])
                    MW.xref[1].append(MW.x0[i])
                MW.peakindex.append(xzero.index(MW.xref[1][np.argmax(MW.yref[1])]))
            
                sumnum = 0
                for i in range(0,len(MW.x02)):
                    sumnum += MW.x02[i] * MW.y02[i]
                MW.height.append(max(MW.y02))
                MW.xcent.append(float(sumnum / MW.fullsum2))
                MW.sumlist.append(MW.fullsum2)
    
            MW.xpoints = int(len(MW.IntWidth) + len(MW.IntWidth2))
            MW.integratedx = MW.xpoints
            if MW.range2 == True:
                MW.peakwidths = np.array([len(MW.IntWidth),len(MW.IntWidth2)])
            else:
                MW.peakwidths = np.array([len(MW.IntWidth)])   
    
            try:
                MW.IntFraction = MW.integratedx/len(MW.zeroWidth)
                MW.N = np.sqrt(MW.IntFraction)*np.sqrt(MW.GaborFraction)*MW.NFrms*2*np.pi
                MW.Nh = np.sqrt(MW.GaborFraction)*MW.NFrms/np.sqrt(len(MW.zeroWidth))*2*np.pi
                if MW.NFrms > 0:
                    print('Frequency Noise RMSD: ' +str(MW.NFrms))
                    print('Mass Noise RMSD: ' +str(MW.Nrms))
                    print('Fraction of Gabor used: ' + str(MW.GaborFraction))
                    print('Fraction of Zero Charge Spectrum integrated: '+str(MW.IntFraction))
                    print('Noise for integrated spectrum: ' + str(MW.N))
            except AttributeError:
                print('Unable to perform noise analysis')
                MW.N = 0
                MW.Nh = 0
            print('Points integrated over: ' + str(MW.integratedx))
    
            MW.paramMatch = False
            MW.bounds = []
            MW.bounds.append(MW.minimumX)
            MW.bounds.append(MW.maximumX)
            if MW.range2 == True:
                MW.bounds.append(MW.minimumX2)
                MW.bounds.append(MW.maximumX2)
    
            MW.based = False
            self.int_show()
        
        except AttributeError:
            print('Unable to perform simple integration. Enter minimum and maximum x values')
            MW.message = 'Unable to perform simple integration.\nPlease enter minimum and maximum x values'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
      
    def batch_integration(self,text):
        self.text = text
        if self.text == 'Parameter Match':
            MW.paramMatch = True
        if self.text == 'Bounds Match':
            MW.paramMatch = False
            try:
                MW.bounds = []
                for i in range(0,len(MW.xref)):
                    MW.bounds.append(MW.xref[i][0])
                    MW.bounds.append(MW.xref[i][-1])
            except AttributeError:
                pass
                
    def newtable(self):
        self.dialog = Int_Table(self)
        self.dialog.show()
        
class Int_Table(QtWidgets.QDialog):
    """
    Creates a pop out table for the integration, with the Mass(Da), 
    Integration, and Height as the columns and a scroll bar option to see all of the entries. 
    """
    def __init__(self, parent=None):
        super(Int_Table, self).__init__(parent)
        self.setWindowTitle("Integration Table")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.setGeometry(100, 100, 800, 800)

        self.time_stamp = QtWidgets.QLabel(MW.time_stamp)
        self.time_stamp.setFixedSize(400,35)
    
        self.createTable()
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.time_stamp)
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout)
        self.show()

    def createTable(self):
        self.tableWidget = QtWidgets.QTableWidget()
        self.tableWidget.setGeometry(QtCore.QRect(310, 10, 311, 161))

        numrows = len(MW.masslist)
        numcols = 3

        self.tableWidget.setRowCount(numrows)
        self.tableWidget.setColumnCount(numcols)
        
        self.tableWidget.setHorizontalHeaderLabels(("Mass (Da)", "Peak Area", "Peak Height"))

        MW.xcent_str = []
        MW.sumlist_str = []
        MW.height_str = []
        
        for i in MW.masslist:
            MW.xcent_str.append(str(i))
            
        for i in MW.intlist:
            MW.sumlist_str.append(str(i))
            
        for i in MW.heightlist:
            MW.height_str.append(str(i))
        
        table_list = [list(i) for i in zip(MW.xcent_str, MW.sumlist_str, MW.height_str)]

        for i in range(numrows):
            for m in range(numcols):
                item = QtWidgets.QTableWidgetItem((table_list[i][m]))
                item.setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
                self.tableWidget.setItem(i, m, item)
        
        self.tableWidget.move(0, 0)
        
class iFAMS_STFT(QtWidgets.QDialog):
    """
    This builds the iFAMS parameters inverse Gabor transform calculation window
    """

    def __init__(self, parent=None):
        super(iFAMS_STFT, self).__init__(parent)
        self.setWindowTitle("Parameters for iFAMS Analysis")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        self.cscb = []
        self.chargestate = []
       
        try:
            if len(MW.inclusion) == len(MW.glist_store) and len(MW.inclusion) == len(MW.zlist_store):
                pass
            else:
                del(MW.inclusion)
        except AttributeError:
            pass
        
        try:
            MW.inclusion
            MW.glist = MW.glist_store.copy()
            try:
                MW.zlist_store[np.nonzero(MW.inclusion)] = MW.zlist
                MW.zlist = MW.zlist_store.copy()
            except ValueError:
                MW.zlist = MW.zlist_store.copy()
            for i in range(0,len(MW.glist)):
                self.cscb.append('self.selection'+str(i))
                self.chargestate.append('self.selnum'+str(i))
            
        except AttributeError:
            MW.inclusion = []
            for i in range(0,len(MW.glist)):
                self.cscb.append('self.selection'+str(i))
                self.chargestate.append('self.selnum'+str(i))
                MW.inclusion.append(True)
            MW.inclusion = np.array(MW.inclusion)
         
        for i in range(0,len(self.cscb)):
            label = 'Selection ' + str(i+1) +' Charge State'
            self.cscb[i] = QtWidgets.QCheckBox(label)
            self.cscb[i].setChecked(MW.inclusion[i])
            self.cscb[i].setFixedSize(300,35)
            self.chargestate[i] = QtWidgets.QTextEdit()
            self.chargestate[i].setFixedSize(150,35)
            self.chargestate[i].setTabChangesFocus(True)
            self.chargestate[i].setToolTip("Enter integer value of charge state for charge correction.\n"
                                        "Enter '0' if no charge correction is needed.")
            try:
                self.chargestate[i].setText(str((MW.zlist[i])))
            except AttributeError:
                continue
            except IndexError:
                continue
        if len(self.cscb) < 9:
            self.note = QtWidgets.QLabel("*Selections near ends of spectrum \n"
                                         "might create windowing artifacts")
        else:
            self.note = QtWidgets.QLabel("*Selections near ends of spectrum might create windowing artifacts")
        self.note.setToolTip('Due to the Gaussian windowing, IFFTs extend further in the m/z direction \n'
                             'than their corresponding pixels in the Gabor spectrogram. If a selection \n'
                             'is too close to the end of the mass spectrum, it might generate windowing \n'
                             'artifacts that can interfere with the automatic baseline correction.\n'
                             'If a deconvolution seems to fail, try reprocessing after deselecting \n'
                             'charge states near the ends of the spectrum.')

        self.button1 = QtWidgets.QPushButton('Run iFAMS Analysis')
        self.button1.clicked.connect(self.replot)

        self.cbN = QtWidgets.QCheckBox('Update Noise Calculation')
        self.cbN.setToolTip('If checked, iFAMS re-estimates the noise included \n'
                            'by calculating from a region in the Fourier spectrum \n'
                            'closer to the Gabor selections')

        self.button2 = QtWidgets.QPushButton('Open Charge Adjuster')
        self.button2.clicked.connect(self.CS_adjust)

        layout = QtWidgets.QGridLayout()
        count = 0
        layout.addWidget(self.note,0,0,1,2)
        for i in range(0,len(self.cscb)):
            if i <= 7+8*count:
                    if count == 0:
                        layout.addWidget(self.cscb[i],i*2+1,0)
                        layout.addWidget(self.chargestate[i],i*2+2,0)
                    else:
                        layout.addWidget(self.cscb[i],i*2-(16*count-1),count)
                        layout.addWidget(self.chargestate[i],i*2-(16*count-2),count)
            else:
                count+=1
                layout.addWidget(self.cscb[i],i*2-(16*count-1),count)
                layout.addWidget(self.chargestate[i],i*2-(16*count-2),count)
                
        layout.addWidget(self.cbN,18,0)
        layout.addWidget(self.button1,19,0)
        try:
            MW.zerofull
            layout.addWidget(self.button2,20,0)
        except AttributeError:
            pass
            
        self.setLayout(layout)

        MW.selecArea = 0
        MW.GaborArea = (MW.xint[-1]-MW.xshort[0])*(MW.xFT[-1]-MW.xFT[0])
        for i in range(0,len(MW.rlist)):
            MW.selecArea += (MW.rlist[i][1]-MW.rlist[i][0])*(MW.rlist[i][3]-MW.rlist[i][2])*2
        MW.GaborFraction = float(MW.selecArea/MW.GaborArea)


    def replot(self):
        MW.cs = []
        rtemp = []
        gtemp = []
        try:
            for i in range(len(self.cscb)):
                if self.cscb[i].isChecked():
                    MW.cs.append(int(str(self.chargestate[i].toPlainText())))
                    gtemp.append(MW.glist[i])
                    for j in range(len(MW.rlist)):
                        if MW.rlist[j][4] == i:
                            rtemp.append(MW.rlist[j])
                            rtemp[-1][4] = len(gtemp)-1
                        elif MW.rlist[j][4] > i:
                            break
        except ValueError:
            MW.message = 'Unable to perform iFAMS analysis.\nPlease enter values for every charge state'
            print(MW.message)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
        del(MW.inclusion)
        MW.fig.clf()
        MW.fig.tight_layout()
        MW.canvas.draw()
        MW.ax = MW.fig.add_subplot(121)
        MW.ax2 = MW.fig.add_subplot(122)
        MW.rlist = rtemp.copy()
        MW.glist = gtemp.copy()
        MW.zlist = np.array(MW.cs)
        MW.ax.plot(MW.xint,MW.yint,label='Mass Spectrum',color='darkgray')
        if MW.realsel.isChecked() == False:
            MW.tIFT = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.tIFT = np.real(MW.glist)
        for i in range(0,len(MW.cs)):
            label = 'Charge ' + str(MW.cs[i])
            MW.ax.plot(MW.xint,MW.tIFT[i],label=label)
        MW.ax.legend(loc='upper right', title='Charge State')
        MW.ax.set_title('Mass Spectrum')
        MW.ax.set_xlabel('m/z')
        MW.ax.set_ylabel('Relative Abundance')

        MW.Ilist = []
        piW = MW.xshort[2]-MW.xshort[1]
        piH = MW.yshort[2]-MW.yshort[1]
        for j in range(len(MW.rlist)):
            lowmzI = int(np.round((MW.rlist[j][0]-MW.xshort[0])/piW))
            himzI = int(np.round((MW.rlist[j][1]-MW.xshort[0])/piW))
            lowFI = int(np.round(MW.rlist[j][2]/piH))
            hiFI = int(np.round(MW.rlist[j][3]/piH))
            MW.Ilist.append([lowmzI,himzI,lowFI,hiFI,MW.rlist[j][4]])
      
        try:
            MW.XBL2
        except AttributeError:
            MW.yBL2 = np.zeros(len(MW.xint))+ min(MW.yint)+ 1+ 3*MW.Nrms
            MW.xintBL2,MW.yintBL2,MW.ypaddBL2,MW.po2BL2,MW.xpaddBL2 = load.datapro(MW.xint,MW.yBL2)
            MW.XBL2 = STFT.stft(MW.ypaddBL2, MW.winnum, MW.window, MW.OF)
        
        MW.glistBL2 = []
        MW.XtempBL2 = MW.XBL2*0
        Icount = 0
        for j in range(len(MW.Ilist)):
            if MW.Ilist[j][4] > Icount:
                MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                MW.glistBL2.append(MW.yreconBL2)
                MW.XtempBL2 = MW.XBL2*0
                Icount += 1
            if MW.Ilist[j][2] <= 1:
                MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1]
                MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:len(MW.X[0])] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:len(MW.X[0])]
            else:
                MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1]
                MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1]
        MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
        MW.glistBL2.append(MW.yreconBL2)
        MW.XtempBL2 = MW.XBL2*0

        MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(MW.tIFT, MW.xint, MW.cs,MW.proton)
        try:
            MW.xzeroBL2, MW.cszeroBL2, MW.zerofullBL2 = FT.zerocharge(np.real(MW.glistBL2), MW.xint, MW.cs,MW.proton)
        except ValueError:
            MW.xzeroBL2 = MW.xzero.copy()
            MW.cszeroBL2 = MW.cszero*0
            MW.zerofullBL2 = MW.zerofull*0
        try:
            massrange=np.nonzero(MW.zerofullBL2)
            MW.xzero = MW.xzero[massrange]
            MW.zerofull = MW.zerofull[massrange]
            MW.xzeroBL2 = MW.xzeroBL2[massrange]
            MW.zerofullBL2 = MW.zerofullBL2[massrange]
            MW.BSint = []
            for j in range(len(MW.cszero)):
                MW.cszero[j] = MW.cszero[j][massrange]
                MW.cszeroBL2[j] = MW.cszeroBL2[j][massrange]
                check = np.max(MW.cszeroBL2[j])/4
                if check == 0:
                    continue
                mid1 = int(len(MW.xzero)/8)+1
                mid2 = int(len(MW.xzero)*7/8)-1
                for k in range(len(MW.cszero[j])):
                    if MW.cszeroBL2[j][k] > check:
                        mid1 = k
                        break
                for k in range(mid1+1,len(MW.cszero[j])):
                    if MW.cszeroBL2[j][k] <= check or MW.cszeroBL2[j][k] <= MW.cszeroBL2[j][mid1]:
                        mid2 = k
                        break
                int1 = min(sum(MW.cszero[j][0:mid1]),sum(MW.cszero[j][mid2:len(MW.xzero)]))
                int2 = (sum(MW.cszeroBL2[j][0:mid1])+sum(MW.cszeroBL2[j][mid2:len(MW.xzero)]))/2
                MW.BSint.append(int1/int2)
            BSint = np.average(MW.BSint)
            MW.zerofullBL2 = MW.zerofullBL2*BSint
            correct_fact = max(MW.zerofull)
            MW.BL_corrected = (MW.zerofull - MW.zerofullBL2)/correct_fact
            MW.zero_norm = MW.zerofull/correct_fact
            MW.BL_norm = MW.zerofullBL2/correct_fact
            MW.zerofull = MW.zerofull - MW.zerofullBL2
        except IndexError:
            pass
        
        if self.cbN.isChecked():
            try:
                self.harm = int(len(MW.rlist)/len(MW.glist))
                if MW.rlist[0][3] < 2:               
                    self.high_funF = MW.rlist[-self.harm][3]
                    MW.NFmin = self.high_funF*(self.harm+1)
                    self.Fgap = MW.rlist[1][2] - MW.rlist[0][3]
                    if self.Fgap > float(30*(MW.xFT[2]-MW.xFT[1])):
                        MW.NFmax = MW.NFmin + self.Fgap
                    else:
                        MW.NFmax = MW.NFmin + float(40*(MW.xFT[2]-MW.xFT[1]))
                else:
                    MW.NFmin = MW.xFT[int((MW.rlist[-self.harm][3]+0.1)*len(MW.xFT)/MW.xFT[-1])]
                    MW.NFmax = MW.xFT[int((MW.rlist[-self.harm][3]+0.3)*len(MW.xFT)/MW.xFT[-1])]
                print('Re-estimating noise between frequencies of '+str(np.round(MW.NFmin,4))+' and '+str(np.round(MW.NFmax,4)))
                MW.NFrms = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin,MW.NFmax),4)
                MW.Nrms = np.round(MW.NFrms/np.sqrt(len(MW.x)),4)
                print('Frequency Noise RMSD: ' + str(MW.NFrms))
                print('Mass Noise RMSD: ' + str(MW.Nrms))
            except AttributeError:
                print('Unable to re-estimate noise. Try manually updating noise in the noise calculator')

        MW.ax2.plot(MW.xzero,MW.zero_norm, label='Full zero-charge')
        try:
            MW.ax2.plot(MW.xzero,MW.BL_norm, label='Modeled Baseline',color='orange')
            MW.ax2.plot(MW.xzero,MW.BL_corrected, label='Baseline-corrected',color='mediumseagreen')
        except AttributeError:
            pass
        except ValueError:
            pass
        MW.cs_norm = []
        for i in range(len(MW.cszero)):
            label = 'Charge '+str(MW.cs[i])
            MW.cs_norm.append(MW.cszero[i]/correct_fact)
            MW.ax2.plot(MW.xzero,MW.cs_norm[i],label=label,alpha=0.5)
        MW.ax2.legend(loc='upper right', title='Charge State')
        MW.ax2.set_title('Zero Charge Spectra')
        MW.ax2.set_ylabel('Relative Abundance')
        MW.ax2.set_xlabel('Mass (Da)')
        
        MW.deltaY = 5
        MW.minimumY = 0.03
        MW.minimumX = MW.xzero[0]
        MW.maximumX = MW.xzero[-1]
        MW.ntol = 0.0
        massrange = np.nonzero(MW.zerofull)
        try:
            MW.xlist,MW.ylist,MW.xref,MW.yref,MW.sumlist,MW.peakindex,MW.minL,MW.minR,MW.xcent,MW.height = Int.integrate(MW.zero_norm[massrange],MW.zerofull[massrange],MW.xzero[massrange], MW.deltaY, MW.minimumY,MW.minimumX, MW.maximumX, MW.ntol)
                
            peaktxt = []
            sumlist = MW.sumlist.copy()
            for i in range(0, min(8,len(sumlist))):
                peaktxt.append(np.argmax(sumlist))
                sumlist[np.argmax(sumlist)] = 0
            for i in peaktxt:
                if len(MW.batchstore) == 0:
                    txty = MW.height[i]*max(MW.BL_corrected)
                else:
                    txty = MW.height[i]*max(MW.zerofull)
                MW.ax2.text(MW.xcent[i],txty,str(np.round(MW.xcent[i],3)),rotation='vertical',ha='center',clip_on=True)
    
            max_value = max(MW.zerofull)
            for i in range(0,len(MW.sumlist)):
                MW.sumlist[i] = max_value*MW.sumlist[i]
                MW.height[i] = max_value*MW.height[i]
        except IndexError:
            pass

        if MW.realsel.isChecked() == False:
            MW.tIFT = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.tIFT = np.real(MW.glist)

        #Noise Calculations
        try:
            MW.selecArea = 0
            MW.GaborArea = (MW.xint[-1]-MW.xshort[0])*(MW.xFT[-1]-MW.xFT[0])
            for i in range(0,len(MW.rlist)):
                MW.selecArea += (MW.rlist[i][1]-MW.rlist[i][0])*(MW.rlist[i][3]-MW.rlist[i][2])*2
            MW.GaborFraction = float(MW.selecArea/MW.GaborArea)
    
            MW.integratedx = 0  
            MW.peakwidths = np.zeros(len(MW.xref))
            for i in range(0, len(MW.xref)):
                MW.integratedx += len(MW.xref[i])
                MW.peakwidths[i] = len(MW.xref[i])
            MW.IntFraction = MW.integratedx/len(MW.xzero[massrange])
            MW.N = np.sqrt(MW.IntFraction)*np.sqrt(MW.GaborFraction)*MW.NFrms*2*np.pi
            MW.Nh = np.sqrt(MW.GaborFraction)*MW.NFrms/np.sqrt(len(MW.xzero[massrange]))*2*np.pi
            if MW.NFrms > 0:
                print('Fraction of Gabor used: ' + str(np.round(MW.GaborFraction,5)))
                print('Fraction of Zero Charge Spectrum integrated: ' + str(np.round(MW.IntFraction,5)))
                print('Noise for integrated spectrum: ' + str(np.round(MW.N,3)))
            
            MW.paramMatch = False
            MW.bounds = []
            for i in range(0,len(MW.xref)):
                MW.bounds.append(MW.xref[i][0])
                MW.bounds.append(MW.xref[i][-1])
        except AttributeError:
            pass
        
        self.close()
        MW.fig.tight_layout()
        MW.canvas.draw()

    def CS_adjust(self):
        MW.cs = []
        gtemp = []
        MW.zlist_store = []
        MW.glist_store = MW.glist.copy()
        try:
            for i in range(len(self.cscb)):
                if self.cscb[i].isChecked():
                    MW.cs.append(int(str(self.chargestate[i].toPlainText())))
                    MW.zlist_store.append(int(str(self.chargestate[i].toPlainText())))
                    gtemp.append(MW.glist[i])
                else:
                    MW.inclusion[i] = False
                    MW.zlist_store.append(int(str(self.chargestate[i].toPlainText())))
        except ValueError:
            MW.message = 'Unable to adjust charge states.\nPlease enter values for every charge state'
            print(MW.message)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
        MW.glist = gtemp.copy()
        MW.zlist = np.array(MW.cs)
        MW.zlist_store = np.array(MW.zlist_store)

        if MW.realsel.isChecked() == False:
            MW.glist = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.glist = np.real(MW.glist)
        self.close()
        self.dialog = Z_Adjust(self)
        self.dialog.show()
        

class STFT_parm(QtWidgets.QDialog):
    """
    This builds the iFAMS parameters recalculation window
    """

    def __init__(self, parent=None):
        super(STFT_parm, self).__init__(parent)
        self.setGeometry(1100, 100, 500, 600)
        self.setWindowTitle("STFT Parameters")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        self.xlist = np.arange(0,5)
        self.window_tag = QtWidgets.QLabel()
        self.window_tag.setText('Window Type')
        self.window_tag.setFixedSize(200,35)
        self.window = QtWidgets.QComboBox()
        self.window.addItem('Select window')
        self.window.addItem('Gaussian')
        self.window.addItem('Hann')
        self.window.addItem('Blackman')
        self.window.addItem('Bartlett')
        self.window.activated[str].connect(self.window_selection)
        self.window.setFixedSize(200,35)
        self.window.setToolTip('The window type used for the Gabor Transform')

        self.winnumTag = QtWidgets.QLabel()
        self.winnumTag.setText('Number of FFT Channels')
        self.winnum2 = QtWidgets.QTextEdit()
        self.winnum2.setFixedSize(150,35)
        self.winnum2.setText(str(MW.winnum))
        self.winnum2.setTabChangesFocus(True)
        self.winnum2.setToolTip('This value defines frequency resolution. Increasing \n'
                                'this number will increase the frequency resolution, but \n '
                                'lower the m/z resolution.  The opposite is true for \n'
                                'a lower number.')

        self.maxtag = QtWidgets.QLabel()
        self.maxtag.setText('Color Scale Relative Maximum')
        self.maxnum = QtWidgets.QTextEdit()
        self.maxnum.setFixedSize(150,35)
        self.maxnum.setText(str(MW.vmax))
        self.maxnum.setTabChangesFocus(True)
        self.maxnum.setToolTip('Adjusting this value will vary the relative maximum \n'
                               'amplitude of the spectrogram.  Less abundant signals \n'
                               'will be seen if this number is smaller.')

        self.OFtag = QtWidgets.QLabel()
        self.OFtag.setText('Oversampling Factor')
        self.OF = QtWidgets.QTextEdit()
        self.OF.setFixedSize(150,35)
        self.OF.setText(str(MW.OF))
        self.OF.setTabChangesFocus(True)
        self.OF.setToolTip('This value determines how many more pixels there are in \n'
                           'the Gabor spectrogram relative to the number of data \n'
                           'data points in the mass spectrum. Higher values result \n'
                           'in a smoother-looking spectrogram but will take longer \n'
                           'to compute. \n'
                           '10x oversampling is recommended for most applications.')

        self.fig2 = Figure()
        self.canvas2 = FigureCanvas(self.fig2)
        self.ax5 = self.fig2.add_subplot(121)
        self.ax6 = self.fig2.add_subplot(122)
        self.canvas2.setFixedSize(480,300)
        self.button2 = QtWidgets.QPushButton('Replot Gabor Spectrogram')
        self.button2.clicked.connect(self.replotgabor)

        self.autotag = QtWidgets.QLabel()
        self.autotag.setText('Auto Adjust: \n'
                             'Enter frequency to resolve')
        self.automass = QtWidgets.QTextEdit()
        self.automass.setFixedSize(150,35)
        self.automass.setText(str(np.round(MW.FoI,2)))
        self.automass.setToolTip('Estimate middle of the range of frequencies containing \n'
                                 'the signal of interest. Auto adjust will attempt to\n'
                                 'resolve the signal around the entered frequency.')
        self.cbMassRes = QtWidgets.QCheckBox('High Mass Resolution')
        self.cbMassRes.setToolTip('If checked, the auto adjust will increase the mass \n'
                                  'resolution at the expense of the frequency resolution. \n'
                                  'This can be useful for isolation ions that are over- \n'
                                  'lapped in the mass spectrum.')
        self.button3 = QtWidgets.QPushButton('Auto Adjust')
        self.button3.clicked.connect(self.replotgaborAuto)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.canvas2,0,0,2,2)
        layout.addWidget(self.window_tag,3,0)
        layout.addWidget(self.window,4,0)
        layout.addWidget(self.winnumTag,5,0)
        layout.addWidget(self.winnum2,6,0)
        layout.addWidget(self.maxtag,7,0)
        layout.addWidget(self.maxnum,8,0)
        layout.addWidget(self.OFtag,9,0)
        layout.addWidget(self.OF,10,0)
        layout.addWidget(self.button2,11,0,1,2)
        layout.addWidget(self.autotag,4,1)
        layout.addWidget(self.automass,5,1)
        layout.addWidget(self.cbMassRes,6,1)
        layout.addWidget(self.button3,7,1)

        self.window.setCurrentText('Gaussian')
        self.window_selection('Gaussian')

        self.setLayout(layout)


    def replotgabor(self):

        try:
            plottop = MW.ax.get_ylim()[1]
            plotbott = MW.ax.get_ylim()[0]
            plotxL = MW.ax.get_xlim()[0]
            plotxR = MW.ax.get_xlim()[1]
            MW.ax.clear()
            MW.window = self.text
            MW.winnum = int(str(self.winnum2.toPlainText()))
            MW.vmax = int(str(self.maxnum.toPlainText()))
            MW.OF = int(str(self.OF.toPlainText()))
            MW.X = STFT.stft(MW.ypadd, MW.winnum, MW.window,MW.OF)
            ytemp = STFT.re_stft(MW.X, MW.winnum, MW.ypadd, MW.yint,MW.OF)
            MW.correction_factor = max(MW.yint)/max(ytemp)
            MW.Xtemp = MW.X*0
            MW.ax.set_title('Gabor spectrogram')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Frequency')
            MW.mzspacing = float((MW.xint[-1] - MW.xint[0])/len(MW.xint))
            MW.xshort = np.linspace((MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing + min(MW.xpadd), max(MW.xpadd) - (MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing, len(MW.X))
            MW.yshort = np.linspace(0 - (max(MW.xFT)/len(MW.X[0]))/2, max(MW.xFT) + (max(MW.xFT)/len(MW.X[0]))/2, len(MW.X[0]))
            MW.mzchan = float((MW.xint[-1]-MW.xshort[0])/((MW.xshort[-1]-MW.xshort[0])/len(MW.xshort)))
            MW.ax.imshow(np.absolute(MW.X.T), origin='lower', aspect='auto',
                     interpolation='nearest',extent= [MW.xshort[0],MW.xshort[-1],MW.yshort[0],MW.yshort[-1]],cmap='jet',vmax=MW.vmax)
            MW.ax.set_ylim(plotbott,plottop)
            MW.ax.set_xlim(plotxL,plotxR)
            MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)
            MW.fig.tight_layout()
            MW.canvas.draw()

        except ValueError:
            print('Not all parameters entered. please enter a value for all parameters')

    def replotgaborAuto(self):
        try:
            MW.FoI = float(str(self.automass.toPlainText()))
            if MW.FoI >= 8:
                MW.winnum = int((MW.maxF)*20)
                MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -5)*len(MW.yFT)/MW.maxF):int((MW.FoI +5)*len(MW.yFT)/MW.maxF)]))/20)
            elif MW.FoI < 8 and MW.FoI > 1:
                MW.winnum = int((MW.maxF)*30)
                MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -1)*len(MW.yFT)/MW.maxF):int((MW.FoI +1)*len(MW.yFT)/MW.maxF)]))/10)
            elif MW.FoI <= 1 and MW.FoI >= 0.1:
                MW.winnum = int((MW.maxF)*80)
                MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -0.05)*len(MW.yFT)/MW.maxF):int((MW.FoI +0.05)*len(MW.yFT)/MW.maxF)]))/10)
            else:
                MW.winnum = int((MW.maxF)*1000)
                MW.vmax = int(max(abs(MW.yFT[int((MW.FoI)*len(MW.yFT)/MW.maxF)-3:int((MW.FoI)*len(MW.yFT)/MW.maxF)+3]))/10)
            if len(MW.xint)/MW.winnum < 8:
                MW.winnum = int(len(MW.xint)/8)
            if self.cbMassRes.isChecked():
                MW.winnum = int(MW.winnum/3)
            plottop = MW.ax.get_ylim()[1]
            plotbott = MW.ax.get_ylim()[0]
            plotxL = MW.ax.get_xlim()[0]
            plotxR = MW.ax.get_xlim()[1]
            MW.ax.clear()
            MW.window = str('gaussian')
            MW.OF = int(10)
            MW.X = STFT.stft(MW.ypadd, MW.winnum,MW.window,MW.OF)
            ytemp = STFT.re_stft(MW.X, MW.winnum, MW.ypadd, MW.yint,MW.OF)
            MW.correction_factor = max(MW.yint)/max(ytemp)
            MW.Xtemp = MW.X*0
            MW.ax.set_title('Gabor Spectrogram')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Frequency')
            MW.mzspacing = float((MW.xint[-1] - MW.xint[0])/len(MW.xint))
            MW.xshort = np.linspace((MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing + min(MW.xpadd), max(MW.xpadd) - (MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing, len(MW.X))
            MW.yshort = np.linspace(0 - (max(MW.xFT)/len(MW.X[0]))/2, max(MW.xFT) + (max(MW.xFT)/len(MW.X[0]))/2, len(MW.X[0]))
            MW.mzchan = float((MW.xint[-1]-MW.xshort[0])/((MW.xshort[-1]-MW.xshort[0])/len(MW.xshort)))
            MW.ax.imshow(np.absolute(MW.X.T), origin='lower', aspect='auto',
                     interpolation='nearest',extent= [MW.xshort[0],MW.xshort[-1],MW.yshort[0],MW.yshort[-1]],cmap='jet',vmax=MW.vmax)
            MW.ax.set_ylim(plotbott,plottop)
            MW.ax.set_xlim(plotxL,plotxR)
            MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)
            self.winnum2.setText(str(MW.winnum))
            self.maxnum.setText(str(MW.vmax))
            self.OF.setText(str(MW.OF))
            
            self.window.setCurrentText('Gaussian')
            self.window_selection('gaussian')
            
            MW.fig.tight_layout()
            MW.canvas.draw()

        except ValueError:
            print('Not all parameters entered. please enter a value for frequency to resolve')


    def window_selection(self,text):
        self.ax5.clear(),self.ax6.clear()
        self.text = text.lower()
        if self.text == 'gaussian':
            self.winDraw = sp.get_window((self.text,15), 100)
        else:
            self.winDraw = sp.get_window(self.text,100)
        self.winDrawFT = abs(np.fft.fft(self.winDraw))
        self.winDrawFT2 = []
        for i in range(0,len(self.winDrawFT)):
            self.winDrawFT2.append(self.winDrawFT[i])
        for i in range(0,50):
            self.winDrawFT2.append(self.winDrawFT2[0])
            del self.winDrawFT2[0]
        self.ax5.plot(self.winDraw)
        self.ax5.set_title('Window')
        self.ax6.plot(self.winDrawFT2)
        self.ax6.set_title('Window FT')
        self.canvas2.draw()

    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y2 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y1 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)


class manu(QtWidgets.QDialog):
    """
    This builds the manual input window
    """

    def __init__(MW_3, parent=None):
        super(manu, MW_3).__init__(parent)
        MW_3.setGeometry(50, 50, 300, 200)
        MW_3.setWindowTitle("Parameters for Finding Maxima")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))

        MW_3.minchar_tag = QtWidgets.QLabel()
        MW_3.minchar_tag.setText('Minimum Charge State')
        MW_3.minchar_tag.setFixedSize(300,35)
        MW_3.minchar = QtWidgets.QTextEdit()
        MW_3.minchar.setFixedSize(150,35)
        MW_3.minchar.setTabChangesFocus(True)

        MW_3.maxchar_tag = QtWidgets.QLabel()
        MW_3.maxchar_tag.setText('Maximum Charge State')
        MW_3.maxchar_tag.setFixedSize(300,35)
        MW_3.maxchar = QtWidgets.QTextEdit()
        MW_3.maxchar.setFixedSize(150,35)
        MW_3.maxchar.setTabChangesFocus(True)

        MW_3.sm_tag = QtWidgets.QLabel()
        MW_3.sm_tag.setText('Subunit Mass')
        MW_3.sm_tag.setFixedSize(300,35)
        MW_3.sm = QtWidgets.QTextEdit()
        MW_3.sm.setFixedSize(150,35)
        MW_3.sm.setTabChangesFocus(True)

        MW_3.button1 = QtWidgets.QPushButton('Update Finder')
        MW_3.button1.clicked.connect(MW_3.recalculate)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(MW_3.minchar_tag)
        layout.addWidget(MW_3.minchar)
        layout.addWidget(MW_3.maxchar_tag)
        layout.addWidget(MW_3.maxchar)
        layout.addWidget(MW_3.sm_tag)
        layout.addWidget(MW_3.sm)
        layout.addWidget(MW_3.button1)
        MW_3.setLayout(layout)

    def recalculate(MW_3):
        try:
            mincharge = float(str(MW_3.minchar.toPlainText()))
            maxcharge = float(str(MW_3.maxchar.toPlainText()))
            MW.submass = float(str(MW_3.sm.toPlainText()))
            MW.cs = np.linspace(mincharge,maxcharge,int(maxcharge-mincharge+1))
            MW.xvalues,MW.yvalues = FT.rel_freq(MW.cs,MW.submass,MW.xFT,MW.yFT)

            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)

            MW.ax.plot(MW.xint, MW.yint)
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')

            MW.ax2.plot(MW.xFT, abs(MW.yFT))
            MW.ax2.scatter(MW.xvalues,MW.yvalues,color='g')
            MW.ax2.set_title('Fourier Spectrum')
            MW.ax2.set_ylabel('Relative Amplitude')
            MW.ax2.set_xlabel('Frequency')
            MW.ax2.set_xlim(0, MW.xFT[-1] / 2)
            MW.fig.tight_layout()
            MW.canvas.draw()

        except ValueError:
            print('Missing a value.  Make sure each box has a value in it')
            MW.message = 'Missing a value. Make sure each box has a value in it'
            MW.Err = True
            MW_3.dialog = iFAMS_message(MW_3)
            MW_3.dialog.show()

        except AttributeError:
            print('No data found for Fourier transform.  Please load in data set')
            MW.message = 'No data found for Fourier transform. Please load in data set'
            MW.Err = True
            MW_3.dialog = iFAMS_message(MW_3)
            MW_3.dialog.show()

class Maxi(QtWidgets.QDialog):
    """
    This builds the iFAMS parameters recalculation window
    """

    def __init__(MW_2, parent=None):
        super(Maxi, MW_2).__init__(parent)
        MW_2.setGeometry(50, 50, 300, 200)
        MW_2.setWindowTitle("Parameters for Finding Maxima")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))


        MW_2.minFreq_tag = QtWidgets.QLabel()
        MW_2.minFreq_tag.setText('Minimum Frequency')
        MW_2.minFreq_tag.setFixedSize(300,35)
        MW_2.minFreq = QtWidgets.QTextEdit()
        MW_2.minFreq.setText('0.001')
        MW_2.minFreq.setFixedSize(150,35)
        MW_2.minFreq.setTabChangesFocus(True)
        MW_2.minFreq.setToolTip('Sets minimum x-value to consider \n'
                                'during search for local maxima.')

        MW_2.minAmp_tag = QtWidgets.QLabel()
        MW_2.minAmp_tag.setText('Minimum Amplitude')
        MW_2.minAmp_tag.setFixedSize(300,35)
        MW_2.minAmp = QtWidgets.QTextEdit()
        MW_2.minAmp.setText('1000')
        MW_2.minAmp.setFixedSize(150,35)
        MW_2.minAmp.setTabChangesFocus(True)
        MW_2.minAmp.setToolTip('Sets minimum y-value to consider \n'
                               'during search for local maxima.')

        MW_2.delta_tag = QtWidgets.QLabel()
        MW_2.delta_tag.setText('Minimum Change in Frequency')
        MW_2.delta_tag.setFixedSize(300,35)
        MW_2.deltatext = QtWidgets.QTextEdit()
        MW_2.deltatext.setText('5')
        MW_2.deltatext.setFixedSize(150,35)
        MW_2.deltatext.setTabChangesFocus(True)
        MW_2.deltatext.setToolTip('Sets search radius around each datapoint \n'
                                  'when searching for local maxima.')

        MW_2.button1 = QtWidgets.QPushButton('Update Finder')
        MW_2.button1.clicked.connect(MW_2.chgnsubmass)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(MW_2.minFreq_tag)
        layout.addWidget(MW_2.minFreq)
        layout.addWidget(MW_2.minAmp_tag)
        layout.addWidget(MW_2.minAmp)
        layout.addWidget(MW_2.delta_tag)
        layout.addWidget(MW_2.deltatext)
        layout.addWidget(MW_2.button1)
        MW_2.setLayout(layout)

    def chgnsubmass(MW_2):
        try:
            minFreq = float(str(MW_2.minFreq.toPlainText()))
            minAmp = float(str(MW_2.minAmp.toPlainText()))
            delta = int(str(MW_2.deltatext.toPlainText()))
            MW.xvalues, MW.yvalues = FT.peakfinder(MW.xFT, np.abs(MW.yFT), minFreq, minAmp, delta)
            MW.cs, MW.submass, stdevmass = FT.chgnsubmass(MW.xvalues, MW.xFT, np.abs(MW.yFT))

            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)

            MW.ax.plot(MW.xint, MW.yint)
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')

            MW.ax2.plot(MW.xFT, abs(MW.yFT))
            MW.ax2.scatter(MW.xvalues, np.abs(MW.yvalues), color='g')
            MW.ax2.set_title('Fourier spectrum')
            MW.ax2.set_ylabel('Relative Amplitude')
            MW.ax2.set_xlabel('Frequency')
            MW.ax2.set_xlim(0, MW.xFT[-1] / 2)
            text = 'The charge states are ' + str(MW.cs) +' \nand the subunit mass is ' + str(round(MW.submass,3)) \
                   + '+/-' + str(round(stdevmass,3))
            MW.ax2.text(0.5, 0.9,text, horizontalalignment='center', verticalalignment='center',
                        transform=MW.ax2.transAxes,color='g',fontsize=12)
            MW.fig.tight_layout()
            MW.canvas.draw()
        except TypeError:
            pass
        except ValueError:
            print('Not enough parameters entered.  Please enter a value in each box.')
            MW.message = 'Missing a value. Make sure each box has a value in it'
            MW.Err = True
            MW_2.dialog = iFAMS_message(MW_2)
            MW_2.dialog.show()
    
        
class testplot(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(testplot, self).__init__(parent)
        self.setGeometry(1100, 100, 800, 900)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot()
        self.ax.plot(MW.px,MW.py)
        self.ax.plot([MW.px[0],MW.px[-1]],[np.average(MW.py),np.average(MW.py)],color='r')
        
        
        self.layout = QtWidgets.QGridLayout()
        self.mpl_toolbar = NavigationToolbar(self.canvas,parent=None)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.mpl_toolbar)
        self.setLayout(self.layout)
        self.canvas.draw()
        MW.fig.tight_layout()
        MW.canvas.draw()

class HarmonicFinder(QtWidgets.QDialog):
    """
    This builds the automatic harmonic finder parameter box. Needs the fundamental peaks
    to be selected first. 
    
    This menu includes the option to select the number of harmonics to add to 
    plot.
    """
    
    def __init__(self, parent=None):
        super(HarmonicFinder, self).__init__(parent)
        self.setWindowTitle("Harmonic Finder Parameters")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        self.setGeometry(1400, 300, 400, 200)
        
        self.harmnum_tag = QtWidgets.QLabel()
        self.harmnum_tag.setText('Highest Harmonic to Include')
        self.harmnum = QtWidgets.QTextEdit()
        self.harmnum.setFixedSize(100,35)
        try:
            self.harmnum.setText(str(MW.HarmRec))
        except AttributeError:
            self.harmnum.setText('2')
        self.harmnum.setToolTip('iFAMS selects harmonics up to the entered value.\n'
                                'The fundamentals are considered the first harmonic. ')
        
        self.cbZeroHarm = QtWidgets.QCheckBox('Include Near-zero Frequencies')
        self.cbZeroHarm.setChecked(True)
        self.cbZeroHarm.setToolTip('Including the near-zero frequencies is recommended \n'
                                   'unless there is a lot of overlap in m/z between \n'
                                   'different species.')

        self.button1 = QtWidgets.QPushButton('Add Harmonics to Plot')
        self.button1.clicked.connect(self.harmadder)

        self.button2 = QtWidgets.QPushButton('Remove Harmonics')
        self.button2.clicked.connect(self.harmremove)
        
        self.button3 = QtWidgets.QPushButton('Start iFAMS Analysis')
        self.button3.clicked.connect(self.harmcs)

        self.layout = QtWidgets.QGridLayout() 
        self.layout.addWidget(self.harmnum_tag)
        self.layout.addWidget(self.harmnum)
        self.layout.addWidget(self.cbZeroHarm)
        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)
        self.layout.addWidget(self.button3)
        self.setLayout(self.layout)


    def harmadder(self):
        """
        This section of codes calculates the harmonic boxes and adds them to the plot. 
        creates a copy of rlist (list containing box points and charge state info) to use for indexing
        creates another copy of rlist for indexing purposes
        """
        print('Finding harmonics')
        try:
            MW.harmonicnum = int(self.harmnum.toPlainText())
        except AttributeError:
            MW.harmonicnum = MW.HarmRec
        self.rliststore = MW.rlist.copy()
        self.gliststore = MW.glist.copy()
        try:
            if self.cbZeroHarm.isChecked() == True:
                start = 0
            else:
                start = 1
        except AttributeError:
            start = 0

        MW.Ilist = []
        piW = MW.xshort[2]-MW.xshort[1]
        piH = MW.yshort[2]-MW.yshort[1]
        MW.Xtemp = MW.X*0
        rtemp = []
        Itemp = []
        gtemp = []
        for i in range(len(MW.rlist)):
            i1 = int(np.round((MW.rlist[i][0]-MW.xshort[0])/piW))
            i2 = int(np.round((MW.rlist[i][1]-MW.xshort[0])/piW))
            i3 = int(np.round(MW.rlist[i][2]/piH))
            i4 = int(np.round(MW.rlist[i][3]/piH))
            MW.Ilist.append([i1,i2,i3,i4,i])
            for j in range(start,MW.harmonicnum+1):
                if j == 0:
                    bottom = 0
                else:
                    bottom = int(np.round(j*MW.Ilist[i][2]+(j-1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                top = int(np.round(j*MW.Ilist[i][2]+(j+1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                rtemp.append([MW.rlist[i][0],MW.rlist[i][1],MW.yshort[bottom],MW.yshort[top],len(gtemp)])
                Itemp.append([MW.Ilist[i][0],MW.Ilist[i][1],bottom,top,len(gtemp)])
                if j == 0:
                    MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1]
                    MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])]
                else:
                    MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1]
                    MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1]
            MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
            gtemp.append(MW.yrecon)
            MW.Xtemp = MW.X*0
        
        MW.Ilist = Itemp
        MW.rlist = rtemp
        MW.glist = gtemp
        
        for i in range(len(MW.rlist)):
            rex = Rectangle((MW.rlist[i][0], MW.rlist[i][2]), (MW.rlist[i][1] - MW.rlist[i][0]), (MW.rlist[i][3] - MW.rlist[i][2]),facecolor='none', edgecolor=color[MW.rlist[i][4]], linewidth=3)
            MW.ax.add_artist(rex)
        
        MW.ax2.clear()
        MW.ax2.plot(MW.x, MW.y, color='darkgrey')
        if MW.realsel.isChecked() == False:
            for i in range(len(MW.glist)):
                MW.ax2.plot(MW.xint,abs(MW.glist[i]),color = color[i])
        if MW.realsel.isChecked() == True:
            for i in range(len(MW.glist)):
                MW.ax2.plot(MW.xint,np.real(MW.glist[i]),color = color[i])
        MW.ax2.set_title('Reconstructed Spectrum')
        MW.ax2.set_xlabel('m/z')
        MW.ax2.set_ylabel('Abundance')

        MW.fig.tight_layout()
        MW.canvas.draw()
        
    def harmremove(self):
        MW.rlist = self.rliststore.copy()
        MW.glist = self.gliststore.copy()
        MW.ax.clear()
        MW.ax2.clear()
        
        MW.ax.imshow(abs(MW.X.T), origin='lower', aspect='auto',
                     interpolation='nearest', extent=(MW.xshort[0], MW.xshort[-1], MW.yshort[0], MW.yshort[-1]),
                     cmap='jet', vmax=MW.vmax)
        MW.ax.set_title('Gabor Spectrogram')
        MW.ax.set_xlabel('m/z')
        MW.ax.set_ylabel('Frequency')
        MW.ax.set_xlim(MW.xint[0], MW.xint[-1])
        MW.ax.set_ylim(MW.xFT[0]-MW.maxF/MW.winnum/2,MW.xFT[int(np.ceil(len(MW.xFT)/2))]+MW.maxF/MW.winnum/2)
        MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)
       
        MW.ax2.plot(MW.x, MW.y, color='darkgrey')
        for i in range(len(MW.rlist)):
            rex = Rectangle((MW.rlist[i][0], MW.rlist[i][2]), (MW.rlist[i][1] - MW.rlist[i][0]), (MW.rlist[i][3] -
                                                                                                  MW.rlist[i][2]),
                            facecolor='none', edgecolor=color[MW.rlist[i][4]], linewidth=3)
            MW.ax.add_artist(rex)
        
        if MW.realsel.isChecked() == False:
            for i in range(len(MW.glist)):
                MW.ax2.plot(MW.xint,abs(MW.glist[i]),color = color[i])
        if MW.realsel.isChecked() == True:
            for i in range(len(MW.glist)):
                MW.ax2.plot(MW.xint,np.real(MW.glist[i]),color = color[i])
        MW.ax2.set_title('Reconstructed Spectrum')
        MW.ax2.set_xlabel('m/z')
        MW.ax2.set_ylabel('Abundance')

        MW.fig.tight_layout()
        MW.canvas.draw()
        
    def harmcs(self):
        self.close()
        self.dialog = iFAMS_STFT(self)
        self.dialog.show()

    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = MW.toggle_selector.RS.y1
            MW.y2 = MW.toggle_selector.RS.y2
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = MW.toggle_selector.RS.y1
            MW.y1 = MW.toggle_selector.RS.y2
     
        
class Calibration(QtWidgets.QDialog):
    """
    This builds the iFAMS calibration curve creation window
    """
    
    def __init__(self, parent=None):
        super(Calibration, self).__init__(parent)
        self.setWindowTitle("Calibration Menu")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))

        self.button1 = QtWidgets.QPushButton('Select Peaklist Files')
        self.button1.clicked.connect(self.load_files)
        self.button1.setFixedSize(400,50)
        self.button1.setToolTip('Select all the peaklist files to use \n'
                                'as calibrants and unknowns/QCs.')

        self.button2 = QtWidgets.QPushButton('Use Batched Files')
        self.button2.clicked.connect(self.batched_files)
        self.button2.setFixedSize(400,50)
        self.button2.setToolTip('Uses the stored files from the \n'
                                'batch just performed.')

        self.layout = QtWidgets.QGridLayout() 
        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)
        self.setLayout(self.layout)

        try:
            MW.batchpeaklists
        except AttributeError:
            self.load_files()

    def load_files(self):
        try:
            try:
                del(MW.batchfolder)
            except AttributeError:
                pass
            MW.Calfiles, extra = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select calibrant peaklist files')
            MW.Calfiles[0]
            MW.cal_dir = os.path.dirname(os.path.dirname(MW.Calfiles[0]))
            self.calparams1()
        except IndexError:
            print('No files selected')
        
    def batched_files(self):
        try:
            MW.Calfiles = MW.batchpeaklists
            MW.cal_dir = MW.batchfolder
            self.calparams1()
        except AttributeError:
            print('No batch to use. Try selecting files')
            MW.message = 'No batch to use. Try selecting files'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def calparams1(self):
        self.button1.hide()
        self.button2.hide()
        self.layout.removeWidget(self.button1)
        self.layout.removeWidget(self.button2)
        
        self.mcalibrant_tag = QtWidgets.QLabel()
        self.mcalibrant_tag.setText('Number of Calibrant Masses')
        self.mcalibrant = QtWidgets.QTextEdit()
        self.mcalibrant.setFixedSize(150,35)
        self.mcalibrant.setText(str(1))
        self.mcalibrant.setTabChangesFocus(True)
        self.mcalibrant.setToolTip('Input number of mass peaks to be combined for a \n' 
                                      'single calibrant concentration. Include all \n'
                                      'peaks of interest. Peaks to include in \n'
                                      'calibration can be toggled after batch.')
        
        self.mstd_tag = QtWidgets.QLabel()
        self.mstd_tag.setText('Number of Internal Standard Masses')
        self.mstd = QtWidgets.QTextEdit()
        self.mstd.setFixedSize(150,35)
        self.mstd.setText(str(0))
        self.mstd.setTabChangesFocus(True)
        self.mstd.setToolTip('Input number of mass peaks to be combined for the \n' 
                                      'internal standard in a single spectrum. \n'
                                      'Enter zero if no internal standard.')

        self.namelist = np.zeros(len(MW.Calfiles)).tolist()
        self.labellist = np.zeros(len(MW.Calfiles)).tolist()
        self.Files_Label = QtWidgets.QLabel()
        self.Files_Label.setText('Add concentration labels and check calibrant data files:')
                
        for i in range(0,len(MW.Calfiles)):
            self.namelist[i] = QtWidgets.QCheckBox(os.path.basename(str(MW.Calfiles[i])))
            self.namelist[i].setFixedSize(320,25)
            self.namelist[i].setChecked(True)
            self.labellist[i] = QtWidgets.QTextEdit()
            self.labellist[i].setFixedSize(150,35)
            self.labellist[i].setTabChangesFocus(True)
         
        self.CalLabel_tag = QtWidgets.QLabel()
        self.CalLabel_tag.setText('Concentration Units')
        self.CalLabel_tag.setFixedSize(300,35) 
        self.CalLabel = QtWidgets.QTextEdit()
        self.CalLabel.setFixedSize(150,35)
        self.CalLabel.setTabChangesFocus(True)
               
        self.button3 = QtWidgets.QPushButton('Next')
        self.button3.clicked.connect(self.calparams2)
        self.button3.setFixedSize(150,35)
        
        self.buttonLL = QtWidgets.QPushButton('Load File Label List')
        self.buttonLL.clicked.connect(self.label_load1)
        self.buttonLL.setFixedSize(300,35)
        self.buttonLL.setToolTip('Option to label file concentrations by \n'
                                'loading list of concentrations as a \n'
                                '.csv or .txt (only numeric concentrations)')
        
        self.setGeometry(50, 50, 300, 200)
        
        self.layout.addWidget(self.Files_Label,0,1)
        self.layout.addWidget(self.buttonLL,1,1)
        count = 0
        for i in range(0,len(MW.Calfiles)):
            if i < 10*(count+1):
                self.layout.addWidget(self.namelist[i],2*i+2-20*count,count+1)
                self.layout.addWidget(self.labellist[i],2*i+3-20*count,count+1)
            else:
                count+=1
                self.layout.addWidget(self.namelist[i],2*i+2-20*count,count+1)
                self.layout.addWidget(self.labellist[i],2*i+3-20*count,count+1)
            
        self.layout.addWidget(self.mcalibrant_tag,0,0)
        self.layout.addWidget(self.mcalibrant,1,0)
        self.layout.addWidget(self.mstd_tag,2,0)
        self.layout.addWidget(self.mstd,3,0)
        self.layout.addWidget(self.CalLabel_tag,4,0)
        self.layout.addWidget(self.CalLabel,5,0)
        self.layout.addWidget(self.button3,6,0)
               
        
    def calparams2(self):
        try:
            MW.fileconc = []
            MW.calconc = []
            MW.calfileIndex = []
            MW.unkfilename = []
            for i in range(0,len(MW.Calfiles)):
                MW.fileconc.append(str(self.labellist[i].toPlainText()))
                if self.namelist[i].isChecked() == True:
                    MW.calconc.append(float(self.labellist[i].toPlainText()))
                    MW.calfileIndex.append(i)
                try:
                    if MW.Calfiles == MW.batchpeaklists and self.namelist[i].isChecked() == False:
                        MW.unkfilename.append(str(self.labellist[i].toPlainText()))
                except AttributeError:
                    pass
                
            self.calMnum = int(self.mcalibrant.toPlainText())
            MW.standards = int(self.mstd.toPlainText())
            MW.Calunits = str(self.CalLabel.toPlainText())
            
            self.button3.hide()
            self.buttonLL.hide()
            self.layout.removeWidget(self.button3)
            self.layout.removeWidget(self.buttonLL)
            self.layout.removeWidget(self.Files_Label)
            self.layout.removeWidget(self.mcalibrant_tag)
            self.layout.removeWidget(self.mcalibrant)
            self.layout.removeWidget(self.mstd_tag)
            self.layout.removeWidget(self.mstd)
            self.layout.removeWidget(self.CalLabel_tag)
            self.layout.removeWidget(self.CalLabel)
            self.Files_Label.setParent(None)
            self.mcalibrant_tag.setParent(None)
            self.mcalibrant.setParent(None)
            self.mstd_tag.setParent(None)
            self.mstd.setParent(None)
            self.CalLabel_tag.setParent(None)
            self.CalLabel.setParent(None)
    
            for i in range(0,len(MW.Calfiles)):
                self.layout.removeWidget(self.namelist[i])
                self.namelist[i].setParent(None)
                self.layout.removeWidget(self.labellist[i])
                self.labellist[i].setParent(None)
            
            self.resize(100,100)
            
            self.tolerance_tag = QtWidgets.QLabel()
            self.tolerance_tag.setText('Tolerance')
            self.tolerance = QtWidgets.QTextEdit()
            self.tolerance.setFixedSize(200,35)
            self.tolerance.setTabChangesFocus(True)
            self.tolerance.setToolTip('How far from expected sample or standard mass \n' 
                                          'is acceptable to use for the calibration \n'
                                          'curve?')
    
            self.tunit = QtWidgets.QComboBox()
            self.tunit.addItem('Select tolerance units')
            self.tunit.addItem('Da')
            self.tunit.addItem('ppm')
            self.tunit.setCurrentText('Da')
            
            self.sigtype = QtWidgets.QComboBox()
            self.sigtype.addItem('Select signal type')
            self.sigtype.addItem('Integration')
            self.sigtype.addItem('Height')
            self.sigtype.setCurrentText('Integration')
            
            self.curtype = QtWidgets.QComboBox()
            self.curtype.addItem('Select curve type')
            self.curtype.addItem('Linear')
            self.curtype.addItem('Quadratic')
            self.curtype.addItem('Logistic')
            self.curtype.setCurrentText('Linear')
            
            self.weight = QtWidgets.QComboBox()
            self.weight.addItem('Select fit weighting')
            self.weight.addItem('No weighting')
            self.weight.addItem('1/x')
            self.weight.addItem('1/x^2')
            self.weight.addItem('x')
            self.weight.addItem('x^2')
            self.weight.setCurrentText('No weighting')
            
            self.buttonLLcal = QtWidgets.QPushButton('Load Calibrant Mass List')
            self.buttonLLcal.clicked.connect(self.label_load2)
            self.buttonLLcal.setFixedSize(300,35)
            self.buttonLLcal.setToolTip('Option to enter calibrant masses by \n'
                                    'loading list of masses as a .csv or .txt')
            
            self.buttonLLstd = QtWidgets.QPushButton('Load Standard Mass List')
            self.buttonLLstd.clicked.connect(self.label_load3)
            self.buttonLLstd.setFixedSize(300,35)
            self.buttonLLstd.setToolTip('Option to enter internal standard masses by \n'
                                    'loading list of masses as a .csv or .txt')
            
            self.callist1 = np.zeros(self.calMnum).tolist()
            self.callist2 = np.zeros(self.calMnum).tolist()
            self.layout.addWidget(self.buttonLLcal,0,1)
            count = 0
            for i in range(0,self.calMnum):
                self.callist1[i] = QtWidgets.QLabel()
                self.callist1[i].setText('Calibrant Mass '+str(i+1))
                self.callist1[i].setFixedSize(200,25)
                self.callist2[i] = QtWidgets.QTextEdit()
                self.callist2[i].setFixedSize(150,35)
                self.callist2[i].setTabChangesFocus(True)
                if i <= 9+10*count:
                    self.layout.addWidget(self.callist1[i],i*2-(20*count-1),count+1)
                    self.layout.addWidget(self.callist2[i],i*2-(20*count-2),count+1)
                else:
                    count+=1
                    self.layout.addWidget(self.callist1[i],i*2-(20*count-1),count+1)
                    self.layout.addWidget(self.callist2[i],i*2-(20*count-2),count+1)
    
            self.stdlist1 = np.zeros(MW.standards).tolist()
            self.stdlist2 = np.zeros(MW.standards).tolist()
            if MW.standards > 0:
                self.layout.addWidget(self.buttonLLstd,0,count+2)
            start_count = count
            for i in range(0,MW.standards):
                self.stdlist1[i] = QtWidgets.QLabel()
                self.stdlist1[i].setText('Standard Mass '+str(i+1))
                self.stdlist1[i].setFixedSize(200,25)
                self.stdlist2[i] = QtWidgets.QTextEdit()
                self.stdlist2[i].setFixedSize(150,35)
                self.stdlist2[i].setTabChangesFocus(True)
                if i <= 9+10*(count-start_count):
                    self.layout.addWidget(self.stdlist1[i],2*i-(20*(count-start_count)-1),count+2)
                    self.layout.addWidget(self.stdlist2[i],2*i-(20*(count-start_count)-2),count+2)
                else:
                    count+=1
                    self.layout.addWidget(self.stdlist1[i],2*i-(20*(count-start_count)-1),count+2)
                    self.layout.addWidget(self.stdlist2[i],2*i-(20*(count-start_count)-2),count+2)

            self.setGeometry(50, 50, 300, 200)
    
            self.button4 = QtWidgets.QPushButton('Create Calibration Curve')
            self.button4.clicked.connect(self.calbuild)
    
            self.layout.addWidget(self.tolerance_tag,0,0)
            self.layout.addWidget(self.tolerance,1,0)
            self.layout.addWidget(self.tunit,2,0)
            self.layout.addWidget(self.sigtype,3,0)
            self.layout.addWidget(self.curtype,4,0)
            self.layout.addWidget(self.weight,5,0)
            self.layout.addWidget(self.button4,6,0)

            MW.canvas.draw()
        except ValueError:
            print('Missing values. Every checked file needs a concentration label')
            MW.message = 'Missing values.\nEvery checked file needs a concentration label'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def label_load1(self):
        try:
            name = QtWidgets.QFileDialog.getOpenFileName(self,'Open Concentration List File')
            namestr = str(name[0])
            if namestr == '':
                print('No file loaded')
                return
            namebase = os.path.splitext(namestr)[0]
            print(namebase)
            if namestr.endswith('.csv'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
            elif namestr.endswith('.txt'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,dtype=str)
            else:
                MW.message = 'Unsupported data format. Please load .csv or .txt'
                print(MW.message)
                MW.Err = True
                self.dialog = iFAMS_message(self)
                self.dialog.show()
            for i in range(0,len(self.labellist)):
                self.labellist[i].setText(self.loaded_label[i])
        except ValueError:
            MW.message = 'Unable to load data. Please check data type'
            print(MW.message)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def label_load2(self):
        try:
            name = QtWidgets.QFileDialog.getOpenFileName(self,'Open Calibrant Mass List File')
            namestr = str(name[0])
            if namestr == '':
                print('No file loaded')
                return
            namebase = os.path.splitext(namestr)[0]
            print(namebase)
            if namestr.endswith('.csv'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
            elif namestr.endswith('.txt'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,dtype=str)
            else:
                MW.message = 'Unsupported data format. Please load .csv or .txt'
                print(MW.message)
                MW.Err = True
                self.dialog = iFAMS_message(self)
                self.dialog.show()
            for i in range(0,len(self.callist2)):
                self.callist2[i].setText(self.loaded_label[i])
        except ValueError:
            MW.message = 'Unable to load data. Please check data type'
            print(MW.message)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def label_load3(self):
        try:
            name = QtWidgets.QFileDialog.getOpenFileName(self,'Open Standard Mass List File')
            namestr = str(name[0])
            if namestr == '':
                print('No file loaded')
                return
            namebase = os.path.splitext(namestr)[0]
            print(namebase)
            if namestr.endswith('.csv'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
            elif namestr.endswith('.txt'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,dtype=str)
            else:
                MW.message = 'Unsupported data format. Please load .csv or .txt'
                print(MW.message)
                MW.Err = True
                self.dialog = iFAMS_message(self)
                self.dialog.show()
            for i in range(0,len(self.stdlist2)):
                self.stdlist2[i].setText(self.loaded_label[i])
        except ValueError:
            MW.message = 'Unable to load data. Please check data type'
            print(MW.message)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def calbuild(self):
        MW.tolerance = float(self.tolerance.toPlainText())
        MW.tunit = str(self.tunit.currentText())
        MW.sigtype = str(self.sigtype.currentText()).lower()
        MW.curtype = str(self.curtype.currentText()).lower()
        MW.weight = str(self.weight.currentText()).lower()

        MW.calmasses = []
        for i in range(0,self.calMnum):
            MW.calmasses.append(float(self.callist2[i].toPlainText()))
        MW.stdmasses = []
        for i in range(0,MW.standards):
            MW.stdmasses.append(float(self.stdlist2[i].toPlainText()))
            #print(self.stdlist2[i].toPlainText())
        
        if MW.sigtype == 'select signal type' or MW.curtype == 'select curve type' or MW.weight == 'select regression weighting' or MW.tunit == 'Select tolerance units':
            print('Missing parameters. Please select an option for every drop-down box')
            MW.message = 'Missing parameters.\nPlease select an option for every drop-down box'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
            
        self.close()
        MW.init_curve = True
        Cal_adjuster(self).show()
                
        
class Cal_adjuster(QtWidgets.QDialog):
    """
    This builds the calibration adjustment window
    """
    def __init__(self, parent = None):
        super(Cal_adjuster, self).__init__(parent)
        self.setWindowTitle("Calibration Curve Menu")
        self.setGeometry(50, 50, 300, 200)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.calHeader = QtWidgets.QLabel()
        self.calHeader.setText('Calibrant Peaks:')
        self.cbcal = np.zeros(len(MW.calmasses)).tolist()
        for i in range(0,len(MW.calmasses)):
            self.cbcal[i] = QtWidgets.QCheckBox(str(MW.calmasses[i]))
            self.cbcal[i].setToolTip('Check to use peak in calibration')
            self.cbcal[i].setChecked(True)
         
        if MW.standards > 0:
            self.stdHeader = QtWidgets.QLabel()
            self.stdHeader.setText('Standard Peaks:')
            self.stdHeader.setFixedSize(200,35)
            self.cbstd = np.zeros(len(MW.stdmasses)).tolist()
            for i in range(0,len(MW.stdmasses)):
                self.cbstd[i] = QtWidgets.QCheckBox(str(MW.stdmasses[i]))
                self.cbstd[i].setToolTip('Check to use peak in calibration')
                self.cbstd[i].setChecked(True)
                
        self.concHeader = QtWidgets.QLabel()
        self.concHeader.setText('Calibrant Concentrations:')
        self.cbconc = np.zeros(len(MW.calconc)).tolist()
        for i in range(0,len(MW.calconc)):
            self.cbconc[i] = QtWidgets.QCheckBox(str(MW.calconc[i]))
            self.cbconc[i].setToolTip('Check to use concentration data point in calibration.')
            self.cbconc[i].setChecked(True)
        
        self.button1 = QtWidgets.QPushButton('Update Calibration')
        self.button1.clicked.connect(self.calAdjust)
        
        self.button2 = QtWidgets.QPushButton('Calculate Concentrations')
        self.button2.clicked.connect(self.calcUnknowns_Menu)
        self.button2.setToolTip('Load peaklist files to calculate sample concentrations \n'
                                'from calibration curve.')
        
        self.tolerance_tag = QtWidgets.QLabel()
        self.tolerance_tag.setText('Tolerance')
        self.tolerance_tag.setFixedSize(200,35) 
        self.tolerance = QtWidgets.QTextEdit()
        self.tolerance.setFixedSize(200,35)
        self.tolerance.setText(str(MW.tolerance))
        self.tolerance.setTabChangesFocus(True)
        self.tolerance.setToolTip('How far from expected sample or standard mass \n' 
                                      'is acceptable to use for the calibration \n'
                                      'curve?')

        self.tunit = QtWidgets.QComboBox()
        self.tunit.addItem('Select tolerance units')
        self.tunit.addItem('Da')
        self.tunit.addItem('ppm')
        self.tunit.setCurrentText(MW.tunit)
        self.tunit.setFixedSize(200,35)
        
        self.sigtype = QtWidgets.QComboBox()
        self.sigtype.addItem('Select signal type')
        self.sigtype.addItem('Integration')
        self.sigtype.addItem('Height')
        self.sigtype.setCurrentText(MW.sigtype.capitalize())
        self.sigtype.setFixedSize(200,35)
        
        self.curtype = QtWidgets.QComboBox()
        self.curtype.addItem('Select curve type')
        self.curtype.addItem('Linear')
        self.curtype.addItem('Quadratic')
        self.curtype.addItem('Logistic')
        self.curtype.setCurrentText(MW.curtype.capitalize())
        self.curtype.setFixedSize(200,35)
        
        self.weight = QtWidgets.QComboBox()
        self.weight.addItem('Select fit weighting')
        self.weight.addItem('No weighting')
        self.weight.addItem('1/x')
        self.weight.addItem('1/x^2')
        self.weight.addItem('x')
        self.weight.addItem('x^2')
        self.weight.setCurrentText(MW.weight.capitalize())
        self.weight.setFixedSize(200,35)

        self.cbCurveSave = QtWidgets.QCheckBox('Overwrite previous curve file')
        self.cbCurveSave.setChecked(True)
        self.cbCurveSave.setToolTip('If checked, the calibration curve will be saved with \n'
                                    'the same name as the previously generated curve.\n'
                                    'If unchecked, the calibration curve will be saved \n'
                                    'to an iterable name to ensure no file is overwritten. ')
        
        self.cbFitAvg = QtWidgets.QCheckBox('Fit averaged data')
        self.cbFitAvg.setChecked(True)
        self.cbFitAvg.setToolTip('Averages calibrants of same concentration before fitting.')
        
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.calHeader,0,0)
        for i in range(len(self.cbcal)):
            self.layout.addWidget(self.cbcal[i],i+1,0)
        n = len(self.cbcal)+2
        self.layout.addWidget(self.tolerance_tag,n,0)
        self.layout.addWidget(self.tolerance,n+1,0)
        self.layout.addWidget(self.tunit,n+2,0)
        self.layout.addWidget(self.sigtype,n+3,0)
        self.layout.addWidget(self.curtype,n+4,0)
        self.layout.addWidget(self.weight,n+5,0)
        self.layout.addWidget(self.cbCurveSave,n+6,0)
        self.layout.addWidget(self.cbFitAvg,n+7,0)
        self.layout.addWidget(self.button1,n+8,0)
        self.layout.addWidget(self.button2,n+9,0)
        if MW.standards > 0:
            self.layout.addWidget(self.stdHeader,0,1)
            for i in range(len(self.cbstd)):
                self.layout.addWidget(self.cbstd[i],i+1,1)
            self.layout.addWidget(self.concHeader,0,2)
            for i in range(len(self.cbconc)):
                self.layout.addWidget(self.cbconc[i],i+1-18*int(i/18),2+int(i/18))
        else:
            self.layout.addWidget(self.concHeader,0,1)
            for i in range(len(self.cbconc)):
                self.layout.addWidget(self.cbconc[i],i+1-18*int(i/18),1+int(i/18))
        self.setLayout(self.layout)
        self.cbFitAvg.hide()
        
        if MW.init_curve == True:
            self.calAdjust()
            MW.init_curve = False
    
    def calAdjust(self):
        try:
            del(Unknown_calculator.unkconc)
        except AttributeError:
            pass
        
        MW.tolerance = float(str(self.tolerance.toPlainText()))
        MW.tunit = str(self.tunit.currentText())
        MW.sigtype = str(self.sigtype.currentText()).lower()
        MW.curtype = str(self.curtype.currentText()).lower()
        MW.weight = str(self.weight.currentText()).lower()
        if MW.sigtype == 'select signal type' or MW.curtype == 'select curve type' or MW.weight == 'select regression weighting' or MW.tunit == 'Select tolerance units':
            print('Missing parameters. Please select an option for every drop-down box')
            MW.message = 'Missing parameters.\nPlease select an option for every drop-down box'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return

        if self.cbCurveSave.isChecked() == False or MW.init_curve == True:
            MW.OverwriteCurve = False
        else:
            MW.OverwriteCurve = True
        
        self.calmasses = []
        for i in range(len(self.cbcal)):
            if self.cbcal[i].isChecked():
                self.calmasses.append(MW.calmasses[i])

        if MW.standards > 0:
            self.stdmasses = []
            for i in range(len(self.cbstd)):
                if self.cbstd[i].isChecked():
                    self.stdmasses.append(MW.stdmasses[i])
          
        self.calconc = []
        MW.concIndex = []
        for i in range(len(self.cbconc)):
            if self.cbconc[i].isChecked():
                if float(MW.calconc[i]) == 0.0 and MW.weight == '1/x' or float(MW.calconc[i]) == 0.0 and MW.weight == '1/x^2':
                    self.cbconc[i].setChecked(False)
                    print('Cannot include concentration of 0 with selected weighting.')
                    MW.message = 'Cannot include concentration of 0 with selected weighting.'
                    MW.Err = False
                    self.dialog = iFAMS_message(self)
                    self.dialog.show()
                    continue
                self.calconc.append(MW.calconc[i])
                MW.concIndex.append(MW.calfileIndex[i])
                
        if len(self.calconc) <= 1 or len(self.calmasses) < 1:
            print('Not enough calibrant information. Please check more calibrant masses or files')
            MW.message = 'Not enough calibrant information.\nPlease check more calibrant masses or files'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
        
        MW.x = []
        MW.y = []
        MW.z = []
        MW.w = []
        MW.ysample = []
        MW.xcent1 = []
        MW.ylist1 = []
        MW.ystandard = []
        MW.xcent2 = []
        MW.ylist2 = []
        MW.areapoint = []
        MW.ysamples = np.zeros(len(self.calmasses))
        MW.xcentsamples = np.zeros(len(self.calmasses))
        if MW.standards > 0:
            MW.ystandards = np.zeros(len(self.stdmasses))
            MW.xcentstandards = np.zeros(len(self.stdmasses))
        MW.sampleNoise = []
        MW.SN = []
        MW.samplewidths = []
        for i in range(0,len(self.calconc)):
            MW.x.append('self.xdata'+str(i))
            MW.y.append('self.ydata'+str(i))
            MW.z.append('self.zdata'+str(i))
            MW.w.append('self.widthdata' +str(i))
            MW.ysample.append('self.ysamp'+str(i))
            MW.xcent1.append('centroid'+str(i))
            MW.ylist1.append('abundance'+str(i))
            MW.ystandard.append('self.ystnd'+str(i))
            MW.xcent2.append('centroid'+str(i))
            MW.ylist2.append('abundance'+str(i))
            MW.areapoint.append('self.areadata'+str(i))
            MW.sampleNoise.append(0)
            MW.SN.append(0)
            MW.samplewidths.append(0)
              
        for i in range(0,len(self.calconc)):
            MW.x[i], MW.y[i], MW.z[i], MW.w[i],namebase = cal.load_peaklist(MW.Calfiles[MW.concIndex[i]])
            
            if MW.sigtype == 'height':
                MW.sampleNoise[i] = cal.sampleNoiseH(MW.x[i],MW.z[i],i)
                MW.sampleNoise[i] = float(MW.sampleNoise[i]*np.sqrt(len(self.calmasses)))

                for j in range(0,len(self.calmasses)):
                    MW.ysamples[j],MW.pwidth,MW.xcentsamples[j] = cal.sample(MW.x[i],MW.z[i],MW.w[i],self.calmasses[j],MW.tolerance,MW.tunit,i)
                MW.ysample[i] = sum(MW.ysamples.copy())
                MW.xcent1[i] = MW.xcentsamples.copy()
                MW.ylist1[i] = MW.ysamples.copy()
                if MW.ysample[i] == 0:
                    print('Warning: No calibrant data within range of tolerance for concentration ' + str(self.calconc[i]) + '. Setting value equal to zero.')
                if MW.sampleNoise[i] == 0:
                    MW.SN[i] == '-'
                else:
                    MW.SN[i] = MW.ysample[i]/MW.sampleNoise[i]

                if MW.standards > 0:
                    if len(self.stdmasses) > 0:
                        for j in range(0,len(self.stdmasses)):
                            MW.ystandards[j],MW.xcentstandards[j] = cal.standard(MW.x[i],MW.z[i],self.stdmasses[j],MW.tolerance,MW.tunit,len(self.stdmasses),i)
                        MW.ystandard[i] = sum(MW.ystandards.copy())
                        MW.xcent2[i] = MW.xcentstandards.copy()
                        MW.ylist2[i] = MW.ystandards.copy()
                        MW.areapoint[i] = MW.ysample[i]/MW.ystandard[i]
                        if MW.ystandard[i] == 1:
                            print('Warning: No internal standard data within range of tolerance for concentration ' + str(self.calconc[i]) + '. Setting value equal to one.')
    
            else:
                MW.sampleNoise[i] = cal.sampleNoise(MW.x[i],MW.y[i],MW.w[i],i)

                for j in range(0,len(self.calmasses)):
                    MW.ysamples[j],pwidth,MW.xcentsamples[j] = cal.sample(MW.x[i],MW.y[i],MW.w[i],self.calmasses[j],MW.tolerance,MW.tunit,i)
                    MW.samplewidths[i] += pwidth
                MW.ysample[i] = sum(MW.ysamples.copy())
                MW.xcent1[i] = MW.xcentsamples.copy()
                MW.ylist1[i] = MW.ysamples.copy()
                if MW.ysample[i] == 0:
                    print('Warning: No calibrant data within range of tolerance for concentration ' + str(self.calconc[i]) + '. Setting value equal to zero.')
                MW.sampleNoise[i] = float(MW.sampleNoise[i]*np.sqrt(MW.samplewidths[i]))
                if MW.sampleNoise[i] == 0:
                    MW.SN[i] == '-'
                else:
                    MW.SN[i] = MW.ysample[i]/MW.sampleNoise[i]
                
                if MW.standards > 0:
                    if len(self.stdmasses) > 0:
                        for j in range(0,len(self.stdmasses)):
                            MW.ystandards[j],MW.xcentstandards[j] = cal.standard(MW.x[i],MW.y[i],self.stdmasses[j],MW.tolerance,MW.tunit,len(self.stdmasses),i)
                        MW.ystandard[i] = sum(MW.ystandards.copy())
                        MW.xcent2[i] = MW.xcentstandards.copy()
                        MW.ylist2[i] = MW.ystandards.copy()
                        MW.areapoint[i] = MW.ysample[i]/MW.ystandard[i]
                        if MW.ystandard[i] == 1:
                            print('Warning: No internal standard data within range of tolerance for concentration ' + str(self.calconc[i]) + '. Setting value equal to one.')

            if MW.sampleNoise[i]==np.sqrt(len(self.calmasses)) or MW.sampleNoise[i]==np.sqrt(MW.samplewidths[i]):
                MW.sampleNoise[i] = str('-')
                MW.SN[i] = str('-')
            else:
                MW.sampleNoise[i] = str("%.2f" % MW.sampleNoise[i])
                MW.SN[i] = str("%.2f" % MW.SN[i])  

        MW.calibrants_used = self.calmasses 
        if MW.standards > 0:
            MW.std_used = self.stdmasses
        else:
            MW.std_used = []

        if MW.standards > 0 and len(MW.std_used) > 0:
            MW.concydata = MW.areapoint
        else:
            MW.concydata = MW.ysample
   
        print('Abundances: ')
        print(np.round(MW.concydata,2))
        MW.conc = self.calconc
        self.calcurve()

    def calcurve(self):
        """
        plots the calibration curve and stores the slope and intercept
        """
        try:
            MW.fig.tight_layout()
            conc = [MW.conc[0]]
            concytemp = [MW.concydata[0]]
            concy = []
            xrep = []
            for i in range(1,len(MW.conc)):
                if MW.conc[i] == conc[-1]:
                    concytemp.append(MW.concydata[i])
                else:
                    conc.append(MW.conc[i])
                    xrep.append(len(concytemp))
                    concy.append(concytemp)
                    concytemp = [MW.concydata[i]]
            concy.append(concytemp)
            xrep.append(len(concytemp))
            MW.concAvg = conc
            MW.concydataAvg = np.zeros(len(conc))
            ystddev = np.zeros(len(conc))
            for i in range(len(MW.concydataAvg)):
                MW.concydataAvg[i] = np.mean(concy[i])
                ystddev[i] = np.std(concy[i])
            if self.cbFitAvg.isChecked() == False:
                MW.concAvg = MW.conc
                MW.concydataAvg = MW.concydata
                ystddev = np.zeros(len(MW.concAvg))
                xrep = 1+np.zeros(len(MW.concAvg))
            try:
                MW.calx,MW.caly,MW.calparam,MW.calinfo,MW.std_err,MW.confidence,MW.uinfo = cal.calibrate(MW.concAvg,xrep,MW.concydataAvg,MW.weight,MW.curtype,ystddev)
                MW.ystddev = ystddev
            except TypeError:
                MW.message = 'Unable to optimize curve. Try different curve type'
                print(MW.message)
                MW.Err = True
                self.dialog = iFAMS_message(self)
                self.dialog.show()
                return
   
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(1,4,(1,3))
            MW.ax.fill_between(MW.confidence[0],MW.confidence[1],MW.confidence[2],color = color[3],alpha=0.3)
            print (MW.calinfo)
            if MW.conc == MW.concAvg:
                if MW.weight == 'no weighting':
                    bars = MW.std_err
                if MW.weight == '1/x':
                    bars = 1/np.array(MW.conc)*MW.std_err
                if MW.weight == '1/x^2':
                    bars = 1/np.array(MW.conc)**2*MW.std_err
                if MW.weight == 'x':
                    bars = MW.std_err*np.array(MW.conc)
                if MW.weight == 'x^2':
                    bars = MW.std_err*np.array(MW.conc)**2
    
                MW.ax.errorbar(MW.conc, MW.concydata, yerr = bars,linestyle='None',marker='o',label='calibrant')
                print(' Standard Deviation About the Regression = ' +str(np.round(MW.std_err,4))+'\n'
                      ' Error bars indicate one standard deviation.')
            else:
                MW.ax.errorbar(MW.conc, MW.concydata, yerr = 0,linestyle='None',marker='o',label='calibrant')
                MW.ax.errorbar(MW.concAvg, MW.concydataAvg, yerr = ystddev,linestyle='None',marker='o',color='r',label='averaged')
                print(' Error bars on averaged data indicate one standard deviation at that concentration.')
                print(' Standard Deviation About the Regression = ' +str(np.round(MW.std_err,4)))
                self.cbFitAvg.show()
            MW.ax.plot(MW.calx, MW.caly)
            MW.ax.set_title(MW.calinfo)
            print(' Band indicates 95% confidence interval.')
            MW.ax.set_xlabel('Concentration ('+MW.Calunits+')')
            if MW.sigtype == 'integration':
                siglabel = 'Peak-area Relative Abundance'
            else:
                siglabel = 'Peak-height Relative Abundance'
            MW.ax.set_ylabel(siglabel)
            MW.ax.legend(loc='lower right')
            try:
                err = max(MW.std_err)
            except TypeError:
                err = MW.std_err
            MW.ax.set_ylim(min(MW.concydata)-3*err,max(MW.concydata)+3*err)
            MW.ax.set_xlim(min(MW.conc)-max(MW.conc)/100,max(MW.conc)+max(MW.conc)/100)
            
            MW.axN = MW.fig.add_subplot(1,4,4)
            columnsN = ('Calibrant \nConcentration ('+MW.Calunits+')','Response')
            MW.axN.axis('off')
            response = []
            for i in range(len(MW.concydata)):
                response.append(str(np.format_float_scientific(MW.concydata[i],precision=2,unique=False)))
            MW.Ntable = MW.axN.table(cellText=np.c_[np.concatenate((MW.conc,[str('Regression SD')])),np.concatenate((response,[np.format_float_scientific(err,precision=2,unique=False)]))],colLabels=columnsN,loc='center')
            MW.Ntable.auto_set_font_size(False)
            MW.Ntable.set_fontsize(12)
            MW.axN.axis('tight')
            MW.Ntable.scale(1.0,2.2)
                       
            MW.file_display.setText('File: '+str(MW.cal_dir)) 
            MW.canvas.draw()

            MW.save_cali(self)
            MW.fig.tight_layout()

        except AttributeError:
            print('No calibration data loaded. Please load integration files first')
            MW.message = 'No calibration data loaded.\nPlease load integration files first'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        except ValueError:
            print('Not enough calibrant information. Please check more calibrants')
            MW.message = 'Not enough calibrant information.\nPlease check more calibrants'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

    def calcUnknowns_Menu(self):
        Unknown_calculator(self).show()
        
        
class Unknown_calculator(QtWidgets.QDialog):
    """
    This builds the unknown calculator
    """
    def __init__(self, parent = None):
        super(Unknown_calculator, self).__init__(parent)
        self.setWindowTitle("Concentration Calculator")
        self.setGeometry(50, 50, 300, 200)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.button1 = QtWidgets.QPushButton('Calculate Concentrations')
        self.button1.clicked.connect(self.calcUnknowns)
        self.button1.setFixedSize(300,35)

        self.button2 = QtWidgets.QPushButton('Load New Peaklist Files')
        self.button2.clicked.connect(self.fileAdjust)
        self.button2.setFixedSize(300,35)

        self.button3 = QtWidgets.QPushButton('Export Concentrations')
        self.button3.clicked.connect(self.export)
        self.button3.setFixedSize(300,35)
        self.button3.setToolTip('Exports a CSV of Unknown/QC calculated concentrations /n'
                                'to the same location as the calibration curve.')

        self.QCcb = QtWidgets.QCheckBox('Plot Expected Concentrations')
        self.QCcb.setToolTip('If concentrations are entered for the file labels, \n'
                             'this option will add vertical dashed lines to \n'
                             'indicate where on the calibration curve that \n'
                             'concentration should be.')
        
        self.AVGcb = QtWidgets.QCheckBox('Average Signal')
        self.AVGcb.setChecked(True)
        self.AVGcb.setToolTip('Averages the signal across samples of identical \n'
                              'labels before calculating the concentration.')
        
        self.layout = QtWidgets.QGridLayout()
        
        self.unkHeader = QtWidgets.QLabel()
        self.unkHeader.setText('Files to Calculate:')
        self.unkHeader.setToolTip('Add labels to the files for concentration calculation \n'
                                  '(any alphanumeric label is acceptable). \n'
                                  'If a number is entered for the label, a line \n'
                                  'indicating the concentration on the x-axis can be \n'
                                  'plotted using the checkbox below.\n'
                                  'If no label is entered, the file name will be used.')
        
        self.buttonLLunk = QtWidgets.QPushButton('Load Label List')
        self.buttonLLunk.clicked.connect(self.label_load4)
        self.buttonLLunk.setFixedSize(300,35)
        self.buttonLLunk.setToolTip('Option to enter unknown labels or QC concentrations \n'
                                    'by loading .csv or .txt list.')

        self.cbUnk = []
        self.unkLabel = []
        try:    
            try:
                if MW.Calfiles == MW.batchpeaklists:
                    MW.unkfiles = []
                    for i in range(len(MW.Calfiles)):
                        try:
                            MW.calfileIndex.index(i)
                            continue
                        except ValueError:
                            MW.unkfiles.append(MW.Calfiles[i])
                else:
                    MW.unkfiles, extra = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select Peaklist Files')
                    MW.unkfiles[0]
            except AttributeError:
                MW.unkfiles, extra = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select Peaklist Files')
                MW.unkfiles[0]  
                    
            for i in range(len(MW.unkfiles)):
                self.cbUnk.append(QtWidgets.QCheckBox(os.path.basename(str(MW.unkfiles[i]))))
                self.unkLabel.append(QtWidgets.QTextEdit())
                self.unkLabel[i].setFixedSize(150,35)
            try:
                if MW.Calfiles == MW.batchpeaklists:
                    for i in range(len(self.unkLabel)):
                        self.unkLabel[i].setText(MW.unkfilename[i])
            except AttributeError:
                pass
                
            self.layout.addWidget(self.unkHeader,0,0)
            self.layout.addWidget(self.buttonLLunk,1,0)
            
            count = 0
            for i in range(len(self.cbUnk)):
                self.cbUnk[i].setChecked(True)
                if i <= 9+10*count:
                    self.layout.addWidget(self.cbUnk[i],i*2-(20*count-2),count)
                    self.layout.addWidget(self.unkLabel[i],i*2-(20*count-3),count)
                else:
                    count+=1
                    self.layout.addWidget(self.cbUnk[i],i*2-(20*count-2),count)
                    self.layout.addWidget(self.unkLabel[i],i*2-(20*count-3),count)

        except IndexError:
            print('No files selected')
            self.err = QtWidgets.QLabel()
            self.err.setText('No files selected')
            self.layout.addWidget(self.err,0,0)

        self.layout.addWidget(self.QCcb,len(self.cbUnk)*2+2,0)
        self.layout.addWidget(self.AVGcb,len(self.cbUnk)*2+3,0)
        self.layout.addWidget(self.button1,len(self.cbUnk)*2+4,0)
        self.layout.addWidget(self.button2,len(self.cbUnk)*2+5,0)
        self.layout.addWidget(self.button3,len(self.cbUnk)*2+6,0)
        self.setLayout(self.layout)
            
    def calcUnknowns(self):
        self.unklab = []
        self.unkfileIndex = []
        for i in range(len(self.cbUnk)):
            if self.cbUnk[i].isChecked():
                self.unklab.append(str(self.unkLabel[i].toPlainText()))
                self.unkfileIndex.append(i)
        
        if self.QCcb.isChecked():
            MW.QClines = True
        else:
            MW.QClines = False
        if self.AVGcb.isChecked():
            self.AvgQC = True
        else:
            self.AvgQC = False
        
        MW.x2 = []
        MW.y2 = []
        MW.z2 = []
        MW.w2 = []
        MW.ysample2 = []
        MW.ystandard2 = []
        MW.areapoint2 = []
        MW.ysamples2 = np.zeros(len(MW.calibrants_used))
        if MW.standards > 0:
            MW.ystandards2 = np.zeros(len(MW.std_used))
        MW.unkNoise = []
        MW.unkSN = []
        MW.unkwidths = []
        MW.conc2 = MW.conc.copy()
        MW.calSN = MW.SN.copy()
        for i in range(0,len(self.unklab)):
            MW.x2.append('self.xdata'+str(i))
            MW.y2.append('self.ydata'+str(i))
            MW.z2.append('self.zdata'+str(i))
            MW.w2.append('self.widthdata'+str(i))
            MW.ysample2.append('self.ysamp'+str(i))
            MW.ystandard2.append('self.ystnd'+str(i))
            MW.areapoint2.append('self.areadata'+str(i))
            MW.unkNoise.append(0)
            MW.unkSN.append(0)
            MW.unkwidths.append(0)
      
        self.unkfilename = []
        for i in range(0,len(self.unklab)):
            MW.x2[i], MW.y2[i], MW.z2[i], MW.w2[i],namebase = cal.load_peaklist(MW.unkfiles[self.unkfileIndex[i]])
            self.unkfilename.append(str(os.path.basename(str(namebase))))

            if MW.sigtype == 'height':
                MW.unkNoise[i] = cal.sampleNoiseH(MW.x2[i],MW.z2[i],i)
                MW.unkNoise[i] = float(MW.unkNoise[i]*np.sqrt(len(MW.calibrants_used)))
                
                for j in range(0,len(MW.calibrants_used)):
                    MW.ysamples2[j],pwidth,xcent = cal.sample(MW.x2[i],MW.z2[i],MW.w2[i],MW.calibrants_used[j],MW.tolerance,MW.tunit,i)
                MW.ysample2[i] = sum(MW.ysamples2)
                if MW.ysample2[i] == 0:
                    print('Warning: No sample data within range of tolerance for unknown ' + str(i+1) + '. Setting value equal to zero.')
                if MW.unkNoise[i] == 0:
                    MW.unkSN[i] == '-'
                else:
                    MW.unkSN[i] = MW.ysample2[i]/MW.unkNoise[i]
                
                if MW.standards > 0:
                    for j in range(0,len(MW.std_used)):
                        MW.ystandards2[j],xcent = cal.standard(MW.x2[i],MW.z2[i],MW.std_used[j],MW.tolerance,MW.tunit,MW.standards,i)
                    MW.ystandard2[i] = sum(MW.ystandards2)
                    MW.areapoint2[i] = MW.ysample2[i]/MW.ystandard2[i]
                    if MW.ystandard2[i] == 1:
                        print('Warning: No internal standard data within range of tolerance for unknown ' + str(i+1) + '. Setting value equal to one.')
            
            else:
                MW.unkNoise[i] = cal.sampleNoise(MW.x2[i],MW.y2[i],MW.w2[i],i)
                
                for j in range(0,len(MW.calibrants_used)):
                    MW.ysamples2[j],pwidth,xcent = cal.sample(MW.x2[i],MW.y2[i],MW.w2[i],MW.calibrants_used[j],MW.tolerance,MW.tunit,i)
                    MW.unkwidths[i] += pwidth
                MW.ysample2[i] = sum(MW.ysamples2)
                if MW.ysample2[i] == 0:
                    print('Warning: No sample data within range of tolerance for unknown ' + str(i+1) + '. Setting value equal to zero.')
                MW.unkNoise[i] = float(MW.unkNoise[i]*np.sqrt(MW.unkwidths[i]))
                if MW.unkNoise[i] == 0:
                    MW.unkSN[i] == '-'
                else:
                    MW.unkSN[i] = MW.ysample2[i]/MW.unkNoise[i]
                
                if MW.standards > 0:
                    for j in range(0,len(MW.std_used)):
                        MW.ystandards2[j],xcent = cal.standard(MW.x2[i],MW.y2[i], MW.std_used[j],MW.tolerance,MW.tunit,MW.standards,i)
                    MW.ystandard2[i] = sum(MW.ystandards2)
                    MW.areapoint2[i] = MW.ysample2[i]/MW.ystandard2[i]
                    if MW.ystandard2[i] == 1:
                        print('Warning: No internal standard data within range of tolerance for unknown ' + str(i+1) + '. Setting value equal to one.')

            if MW.unkNoise[i]==np.sqrt(len(MW.calibrants_used)) or MW.unkNoise[i]==np.sqrt(MW.unkwidths[i]):
                MW.unkNoise[i] = str('-')
                MW.unkSN[i] = str('-')
            else:
                MW.unkNoise[i] = str("%.2f" % MW.unkNoise[i])
                MW.unkSN[i] = str("%.2f" % MW.unkSN[i])  
            if self.unklab[i] == '':
                self.unklab[i] = os.path.basename(str(MW.unkfiles[self.unkfileIndex[i]]))
            MW.conc2.append(self.unklab[i]) 
            MW.calSN.append(MW.unkSN[i])
            
        if MW.standards > 0:
            MW.concydata2 = MW.areapoint2            
        else:
            MW.concydata2 = MW.ysample2

        #average
        if self.AvgQC == True:
            conc = [self.unklab[0]]
            concytemp = [MW.concydata2[0]]
            concy = []
            self.m = []
            for i in range(1,len(MW.concydata2)):
                if self.unklab[i] == conc[-1]:
                    concytemp.append(MW.concydata2[i])
                else:
                    conc.append(self.unklab[i])
                    self.m.append(len(concytemp))
                    concy.append(concytemp)
                    concytemp = [MW.concydata2[i]]
            concy.append(concytemp)
            self.m.append(len(concytemp))
            self.unklab = conc
            MW.concydata2 = np.zeros(len(self.unklab))
            self.ystddev = np.zeros(len(self.unklab))
            for i in range(len(MW.concydata2)):
                MW.concydata2[i] = np.mean(concy[i])
                self.ystddev[i] = np.std(concy[i])
            self.m = np.array(self.m)
        else:
            self.m = np.zeros(len(MW.concydata2))+1
            self.ystddev = np.zeros(len(MW.concydata2))

        self.unkerr = []
        if MW.curtype == 'linear':
            self.unkconc,MW.errbarsU = cal.unkcalc(MW.calparam[0],MW.calparam[1],MW.concydata2,MW.std_err,MW.uinfo[0],MW.uinfo[1],len(MW.conc),self.m)
            for i in range(len(self.unkconc)):
                self.unkerr.append(MW.errbarsU[i])
        if MW.curtype == 'quadratic':
            self.unkconc,MW.errbarsU = cal.unkcalc2(MW.calparam[0],MW.calparam[1],MW.calparam[2],MW.concydata2,MW.std_err,MW.uinfo,self.m)
            for i in range(len(self.unkconc)):
                self.unkerr.append(MW.errbarsU[i])
        if MW.curtype == 'logistic':
            self.unkconc,MW.errbarsU = cal.unkcalc3(MW.calparam[0],MW.calparam[1],MW.calparam[2],MW.concydata2,MW.std_err,MW.uinfo,self.m)
            for i in range(len(self.unkconc)):
                yerr1 = MW.concydata2[i]+MW.errbarsU[i]
                yerr2 = MW.concydata2[i]-MW.errbarsU[i]
                xerr1 = -np.log(MW.calparam[0]/(yerr1-MW.calparam[2])-1)/MW.calparam[1]-self.unkconc[i]
                xerr2 = self.unkconc[i]+np.log(MW.calparam[0]/(yerr2-MW.calparam[2])-1)/MW.calparam[1]
                self.unkerr.append(np.average((xerr1,xerr2)))
        print('Abundances: ')
        print (np.round(MW.concydata2,2))
        print('Calculated concentrations: ')
        for i in range(len(self.unkconc)):
            print(np.round(self.unkconc[i],4))
        self.QCliney = []
        self.QClinex = []

        MW.fig.clf()
        MW.ax2 = MW.fig.add_subplot(1,4,(1,3))
        MW.ax2.fill_between(MW.confidence[0],MW.confidence[1],MW.confidence[2],color = color[3],alpha=0.3)
        MW.ax2.plot(MW.calx, MW.caly,color=color[2])
        if MW.conc == MW.concAvg:
            if MW.weight == 'no weighting':
                bars = MW.std_err
            if MW.weight == '1/x':
                bars = 1/np.array(MW.conc)*MW.std_err
            if MW.weight == '1/x^2':
                bars = 1/np.array(MW.conc)**2*MW.std_err
            if MW.weight == 'x':
                bars = MW.std_err*np.array(MW.conc)
            if MW.weight == 'x^2':
                bars = MW.std_err*np.array(MW.conc)**2

            MW.ax2.errorbar(MW.conc, MW.concydata, yerr = bars,linestyle='None',marker='o',label='calibrant')
        else:
            MW.ax2.errorbar(MW.conc, MW.concydata, yerr = 0,linestyle='None',marker='o',label='calibrant')
            MW.ax2.errorbar(MW.concAvg, MW.concydataAvg, yerr = MW.ystddev,linestyle='None',marker='o',color='r',label='averaged')
        
        count = 0
        colorlist = []
        for i in range(0, len(self.unkconc)):
            if i > 0 and self.unklab[i] == self.unklab[i-1]:
                count += 2
                MW.ax2.errorbar(self.unkconc[i], MW.concydata2[i], yerr=self.ystddev[i], linestyle = 'None', marker='s',color = color[i*2+2-count])
            else:
                count = 0
                MW.ax2.errorbar(self.unkconc[i], MW.concydata2[i], yerr=self.ystddev[i], linestyle = 'None', marker='s',color = color[i*2+2-count],label=str(self.unklab[i]))
            colorlist.append(i*2+2-count)
            try:
                self.QClinex.append([float(self.unklab[i]),float(self.unklab[i])])
                if MW.curtype == 'linear':
                    QCy = float(MW.calparam[0]*float(self.unklab[i])+MW.calparam[1])
                    self.QCliney.append([MW.caly[0],QCy])
                if MW.curtype == 'quadratic':
                    QCy = float(MW.calparam[0]*float(self.unklab[i])**2+MW.calparam[1]*float(self.unklab[i])+MW.calparam[2])
                    self.QCliney.append([MW.caly[0],QCy])
                if MW.curtype == 'logistic':
                    QCy = float(MW.calparam[0]/(1+np.exp(-MW.calparam[1]*float(self.unklab[i])))+MW.calparam[2])
                    self.QCliney.append([MW.caly[0],QCy])
            except ValueError:
                self.QClinex.append(0)
                self.QCliney.append(0)
        MW.ax2.set_title(MW.calinfo)
        MW.ax2.legend(loc='lower right')
        MW.ax2.set_xlabel('Concentration ('+MW.Calunits+')')
        try:
            err = max(MW.std_err)
        except TypeError:
            err = MW.std_err
        MW.ax2.set_ylim(min(MW.concydata)-3*err,max(MW.concydata)+3*err)
        MW.ax2.set_xlim(min(MW.conc)-max(MW.conc)/100,max(MW.conc)+max(MW.conc)/100)
        if MW.sigtype == 'integration':
            siglabel = 'Peak-area Relative Abundance'
        else:
            siglabel = 'Peak-height Relative Abundance'
        MW.ax2.set_ylabel(siglabel)
        
        if MW.QClines == True:
            for i in range(0,len(self.unkconc)):
                if self.QCliney[i] == 0:
                    continue
                else:
                    MW.ax2.plot(self.QClinex[i],self.QCliney[i],linestyle='dashed',marker='None',color = color[colorlist[i]])
                    MW.ax2.text(self.QClinex[i][1],float(self.QCliney[i][1]-count),str(str(self.QClinex[i][1])+' QC'),color = color[colorlist[i]],horizontalalignment = 'right',clip_on=True)
        
        MW.axU = MW.fig.add_subplot(1,4,4)
        columnsN = ('Calibrant \nConcentration ('+MW.Calunits+')','Response')
        MW.axU.axis('off')
        response = []
        for i in range(len(MW.concydata)):
            response.append(str(np.format_float_scientific(MW.concydata[i],precision=2,unique=False)))
        MW.Utable = MW.axU.table(cellText=np.c_[np.concatenate((MW.conc,[str('Regression SD')])),np.concatenate((response,[np.format_float_scientific(err,precision=2,unique=False)]))],colLabels=columnsN,loc='center')
        MW.Utable.auto_set_font_size(False)
        MW.Utable.set_fontsize(12)
        MW.axU.axis('tight')
        MW.Utable.scale(1.0,2.2)
        
        self.unklabT,self.unkconcT,self.unkerrT,self.bias,self.cv = [],[],[],[],[]
        unkconctemp,unkerrtemp = [],[]
        for i in range(len(self.unklab)):
            if i > 0 and self.unklab[i] == self.unklab[i-1]:
                pass
            elif i == 0:
                self.unklabT.append(self.unklab[i])
            else:
                self.unkconcT.append(np.round(np.average(unkconctemp),4))
                self.unkerrT.append(np.round(np.average(unkerrtemp)*np.sqrt(1/len(unkerrtemp)),4))
                try:
                    labeltemp = float(self.unklabT[-1])
                    self.bias.append(np.round((self.unkconcT[-1]-labeltemp)/labeltemp*100,1))
                    self.cv.append(np.round((np.std(unkconctemp)/labeltemp*100),1))
                except ValueError:
                    self.bias.append('N/A')
                    self.cv.append('N/A')
                unkconctemp,unkerrtemp = [],[]
                self.unklabT.append(self.unklab[i])
            unkconctemp.append(self.unkconc[i])
            unkerrtemp.append(self.unkerr[i])
        self.unkconcT.append(np.round(np.average(unkconctemp),4))
        self.unkerrT.append(np.round(np.average(unkerrtemp)*np.sqrt(1/len(unkerrtemp)),4))
        try:
            labeltemp = float(self.unklabT[-1])
            self.bias.append(np.round((self.unkconcT[-1]-labeltemp)/labeltemp*100,1))
            self.cv.append(np.round((np.std(unkconctemp)/labeltemp*100),1))
        except ValueError:
            self.bias.append('N/A')
            self.cv.append('N/A')
        MW.UNK_len = len(self.unklabT)
        MW.UNK_TABLE = np.c_[self.unklabT,self.unkconcT,self.unkerrT,self.bias,self.cv]
        Unknown_calculator.newtable(self)

        if self.AvgQC == True:
            print('Error bars on square points indicate one standard deviation ' 
                  'in the averaged signal')

        MW.canvas.draw() 

        
    def fileAdjust(self):
        try:
            MW.unkfilesstore = MW.unkfiles.copy()
            MW.unkfiles, extra = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select Peaklist Files')
            MW.unkfiles[0]
            
            try:
                self.layout.removeWidget(self.err)
                self.err.setParent(None)
            except AttributeError:
                self.layout.removeWidget(self.unkHeader)
                self.unkHeader.setParent(None)
                for i in range(len(self.cbUnk)):
                    self.layout.removeWidget(self.cbUnk[i])
                    self.layout.removeWidget(self.unkLabel[i])
                    self.cbUnk[i].setParent(None)
                    self.unkLabel[i].setParent(None)
           
            self.layout.removeWidget(self.QCcb)
            self.layout.removeWidget(self.AVGcb)
            self.button1.hide()
            self.button2.hide()
            self.button3.hide()
            self.buttonLLunk.hide()
            self.QCcb.setParent(None)
            self.AVGcb.setParent(None)
            self.button1.setParent(None)
            self.button2.setParent(None)
            self.button3.setParent(None)
            self.buttonLLunk.setParent(None)

            self.cbUnk = []
            self.unkLabel = []
            for i in range(len(MW.unkfiles)):
                self.cbUnk.append(QtWidgets.QCheckBox(os.path.basename(str(MW.unkfiles[i]))))
                self.unkLabel.append(QtWidgets.QTextEdit())
                self.unkLabel[i].setFixedSize(150,35)
                
            self.layout.addWidget(self.unkHeader,0,0)
            self.layout.addWidget(self.buttonLLunk,1,0)
            
            count = 0
            for i in range(len(self.cbUnk)):
                self.cbUnk[i].setChecked(True)
                if i <= 9+10*count:
                    self.layout.addWidget(self.cbUnk[i],i*2-(20*count-2),count)
                    self.layout.addWidget(self.unkLabel[i],i*2-(20*count-3),count)
                else:
                    count+=1
                    self.layout.addWidget(self.cbUnk[i],i*2-(20*count-2),count)
                    self.layout.addWidget(self.unkLabel[i],i*2-(20*count-3),count)
            
            self.layout.addWidget(self.QCcb,len(self.cbUnk)*2+2,0)
            self.layout.addWidget(self.AVGcb,len(self.cbUnk)*2+3,0)
            self.layout.addWidget(self.button1,len(self.cbUnk)*2+4,0)
            self.layout.addWidget(self.button2,len(self.cbUnk)*2+5,0)
            self.layout.addWidget(self.button3,len(self.cbUnk)*2+6,0)
            self.button1.show()
            self.button2.show()
            self.button3.show()
            self.buttonLLunk.show()

        except IndexError:
            print('No files selected')
            MW.unkfiles = MW.unkfilesstore.copy()
        
    def label_load4(self):
        try:
            name = QtWidgets.QFileDialog.getOpenFileName(self,'Open Label List File')
            namestr = str(name[0])
            if namestr == '':
                print('No file loaded')
                return
            namebase = os.path.splitext(namestr)[0]
            print(namebase)
            if namestr.endswith('.csv'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,delimiter=',',dtype=str)
            elif namestr.endswith('.txt'):
                self.loaded_label = np.loadtxt(namestr,unpack=True,dtype=str)
            else:
                MW.message = 'Unsupported data format. Please load .csv or .txt'
                print(MW.message)
                MW.Err = True
                self.dialog = iFAMS_message(self)
                self.dialog.show()
            for i in range(0,len(self.unkLabel)):
                self.unkLabel[i].setText(self.loaded_label[i])
        except ValueError:
            MW.message = 'Unable to load data. Please check data type'
            print(MW.message)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def export(self):
        try:
            unklab = np.concatenate((['Unknown/QC Label'],self.unklab))
            if MW.sigtype == 'height':
                unkresponse = np.concatenate((['Response (height)'],MW.concydata2))
            else:
                unkresponse = np.concatenate((['Response (area)'],MW.concydata2))
            unkconc = np.concatenate((['Calc. Concentration ('+str(MW.Calunits)+')'],self.unkconc))
            temp_path = os.path.join(os.path.dirname(MW.cal_dir), str('Calc_Concs_from_'+str(os.path.basename(MW.cal_dir))+'.csv'))
            np.savetxt(temp_path, np.c_[unklab, unkresponse, unkconc], delimiter=',',fmt='%s')
            MW.message = str('Exported calculated concentration data to:\n' + temp_path)
            print(MW.message)
            MW.Err = False
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        except AttributeError:
            print('No calculation to save.\nPlease recalculate concentrations of unknown/QC data first')
            MW.message = 'No calculation to save.\nPlease recalculate concentrations of unknown/QC data first'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
               
    def newtable(self):
        self.dialog = Unk_Table(self)
        self.dialog.show()
            
class Unk_Table(QtWidgets.QDialog):
    """
    Creates a pop out table for the unknown/QC concentrations
    """
    def __init__(self, parent=None):
        super(Unk_Table, self).__init__(parent)
        self.setWindowTitle("Calculated Concentrations Table")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.setGeometry(100, 100, 1200, 800)

        self.cal_stamp = QtWidgets.QLabel(MW.calinfo)
        self.cal_stamp.setFixedSize(400,100)
    
        self.createTable()
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.cal_stamp)
        self.layout.addWidget(self.tableWidget)
        self.setLayout(self.layout)
        self.show()

    def createTable(self):
        self.tableWidget = QtWidgets.QTableWidget()
        self.tableWidget.setGeometry(QtCore.QRect(310, 10, 311, 161))

        numrows = MW.UNK_len
        numcols = 5

        self.tableWidget.setRowCount(numrows)
        self.tableWidget.setColumnCount(numcols)
        
        self.tableWidget.setHorizontalHeaderLabels(('Label','Concentration\n('+MW.Calunits+')','Uncertainty\n(95% conf.)','% Bias', '% CV'))
        table_list = MW.UNK_TABLE

        for i in range(numrows):
            for m in range(numcols):
                item = QtWidgets.QTableWidgetItem((table_list[i][m]))
                item.setTextAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
                self.tableWidget.setItem(i, m, item)
        
        self.tableWidget.move(0, 0)


class Spec_domain(QtWidgets.QDialog):
    """
    This builds the window for adjusting the mass spectrum domain, replotting, and saving a truncated data file
    """
    def __init__(self, parent = None):
        try:
            super(Spec_domain, self).__init__(parent)
            self.setGeometry(50, 50, 400, 400)
            self.setWindowTitle("Mass Spectrum Domain")
            QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
            
            self.minx_tag = QtWidgets.QLabel()
            self.minx_tag.setText('Minimum x-value')
            self.minx_tag.setFixedSize(400,45)
            self.minx = QtWidgets.QTextEdit()
            self.minx.setFixedSize(400,45)
            self.minx.setText(str(np.round(MW.x[0],4)))
            self.minx.setTabChangesFocus(True)
            
            self.maxx_tag = QtWidgets.QLabel()
            self.maxx_tag.setText('Maximum x-value')
            self.maxx_tag.setFixedSize(400,45)
            self.maxx = QtWidgets.QTextEdit()
            self.maxx.setFixedSize(400,45)
            self.maxx.setText(str(np.round(MW.x[-1],4)))
            self.maxx.setTabChangesFocus(True)
            
            self.mzstep_tag = QtWidgets.QLabel()
            self.mzstep_tag.setText('Interpolated m/z sampling')
            self.mzstep_tag.setFixedSize(400,45)
            self.mzstep = QtWidgets.QTextEdit()
            self.mzstep.setFixedSize(400,45)
            self.mzstep.setText(str(np.round(1/MW.maxF,4)))
            self.mzstep.setTabChangesFocus(True)
            
            self.button1 = QtWidgets.QPushButton('Replot')
            self.button1.clicked.connect(self.replot)
            self.button1.setToolTip('Re-interpolates data over entered domain with entered m/z sampling and replots spectra ')
            
            self.button2 = QtWidgets.QPushButton('Save Raw Truncated Data')
            self.button2.clicked.connect(self.trunc_save1)
            self.button2.setToolTip('Saves data in adjusted domain as a .csv')
            
            self.button3 = QtWidgets.QPushButton('Save Interpolated Data')
            self.button3.clicked.connect(self.trunc_save2)
            self.button3.setToolTip('Saves data in adjusted domain with specified interpolation as a .csv')
            
            self.button4 = QtWidgets.QPushButton('Restore Spectra')
            self.button4.clicked.connect(self.restore_plot)
            self.button4.setToolTip('Replots spectra over original domain')
            
            layout = QtWidgets.QGridLayout()
            layout.addWidget(self.minx_tag)
            layout.addWidget(self.minx)
            layout.addWidget(self.maxx_tag)
            layout.addWidget(self.maxx)
            layout.addWidget(self.mzstep_tag)
            layout.addWidget(self.mzstep)
            layout.addWidget(self.button1)
            layout.addWidget(self.button2)
            layout.addWidget(self.button3)
            layout.addWidget(self.button4)
            self.setLayout(layout)
            
            MW.xstore = MW.x.copy()
            MW.ystore = MW.y.copy()
            if MW.stft_pro == True:
                MW.winnumstore = MW.winnum
                MW.vmaxstore = MW.vmax
        except TypeError:
            MW.message = 'Please load mass spectrum as FT or STFT first'
            print(MW.message)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def replot(self):
        try:
            MW.x = []
            MW.y = []
            
            minx = float(str(self.minx.toPlainText()))
            maxx = float(str(self.maxx.toPlainText()))
            mzstep = float(str(self.mzstep.toPlainText()))
            lenx = int(np.round((maxx-minx)/mzstep))

            for i in range(len(MW.xstore)):
                if (minx <= MW.xstore[i]) and (MW.xstore[i] <= maxx):
                     MW.x.append(MW.xstore[i])
                     MW.y.append(MW.ystore[i])
            
            MW.xint = np.linspace(minx,maxx,lenx)
            MW.yint = np.interp(MW.xint,MW.x,MW.y)
            for i in range(len(MW.yint)):
                MW.po2 = 2**i
                if MW.po2 >lenx:
                    break
            zeros = np.zeros(MW.po2-lenx)
            MW.ypadd = np.append(MW.yint,zeros)
            spacing = mzstep
            xend = MW.po2-len(MW.xint)
            MW.xpadd = np.linspace(MW.xint[0],MW.xint[-1]+spacing*xend,MW.po2)

            MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)
            MW.maxF = max(MW.xFT)

            if MW.ft_pro == True:
                MW.fig.clf()
                MW.ax = MW.fig.add_subplot(121)#put in upper left hand corner
                MW.ax2 = MW.fig.add_subplot(122)
    
                if MW.MS_interp.isChecked():
                    MW.ax.plot(MW.xint,MW.yint)
                    MW.ax.set_title('Interpolated Mass Spectrum')
                else:
                    MW.ax.plot(MW.x,MW.y)
                    MW.ax.set_title('Raw Mass Spectrum')
                MW.ax.set_ylabel('Relative Abundance')
                MW.ax.set_xlabel('m/z')
    
                MW.ax2.plot(MW.xFT,abs(MW.yFT))
                MW.ax2.set_title('Fourier Spectrum')
                MW.ax2.set_ylabel('Relative Amplitude')
                MW.ax2.set_xlabel('Frequency')
                MW.ax2.set_xlim(0, MW.xFT[-1] / 2)
                MW.fig.tight_layout()
                MW.canvas.draw() #shows the graphs

            if MW.stft_pro == True:
                #Estimating Gabor parameters
                MW.maxF = float(MW.xFT[-1])
                MW.maxab = float(max(MW.yint))
                if MW.maxF < 20:
                    MW.winnum = int((MW.maxF)**(4/3)*20)
                    MW.vmax = int(MW.maxab*2)
                elif max(MW.yFT[10*int(len(MW.yFT)/MW.maxF):int(len(MW.yFT)/2)]) > 10*MW.NFrms:
                    MW.winnum = int((MW.maxF)**(2/3)*80)
                    MW.vmax = int(MW.maxab*2)
                elif max(MW.yFT[1*int(len(MW.yFT)/MW.maxF):(10*int(len(MW.yFT)/MW.maxF))]) > 10*MW.NFrms:
                    MW.winnum = int((MW.maxF)**(4/3)*20)
                    MW.vmax = int(MW.maxab*10)
                else:
                    MW.winnum = int((MW.maxF)**(4/3)*80)
                    MW.vmax = int(MW.maxab*20)
                if MW.winnum < 10:
                    MW.winnum = 10
    
                try:
                    MW.plot_gabor(self)
                except IndexError:
                    MW.winnum = 10
                    MW.OF = 10
                    MW.plot_gabor(self)

            print('Domain: {0} to {1}'.format(np.round(MW.x[0],4),np.round(MW.x[-1],4)))

            MW.domaintrunc = True
            MW.domainminx = float(str(self.minx.toPlainText()))
            MW.domainmaxx = float(str(self.maxx.toPlainText()))


        except ValueError:
            print('Not all parameters entered. please enter a value for all parameters')
            MW.message = 'Not all parameters entered.\nPlease enter a value for all parameters'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = MW.toggle_selector.RS.y1
            MW.y2 = MW.toggle_selector.RS.y2
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = MW.toggle_selector.RS.y1
            MW.y1 = MW.toggle_selector.RS.y2

    def trunc_save1(self):
        try:
            if MW.namebase == 'From Clipboard':
                MW.ReturnFunc = Spec_domain.trunc_save
                Folder_Select(self).show()
                return
            print('Saving files')
            label = os.path.basename(MW.namebase) + " trunc.csv"
            temp_path = os.path.join(os.path.dirname(MW.namebase), label)
            np.savetxt(temp_path, np.c_[np.round(MW.x, 7), np.round(MW.y,7)], delimiter=',')
            print('Saved truncated data as ' + label)
            MW.message = str('Saved truncated data as: \n' + label)
            MW.Err = False
            self.dialog = iFAMS_message(self)
            self.dialog.show()

        except AttributeError:
            print('No data to save')
            MW.message = 'No data to save'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            
    def trunc_save2(self):
        try:
            if MW.namebase == 'From Clipboard':
                MW.ReturnFunc = Spec_domain.trunc_save
                Folder_Select(self).show()
                return
            print('Saving files')
            label = os.path.basename(MW.namebase) + " interp_trunc.csv"
            temp_path = os.path.join(os.path.dirname(MW.namebase), label)
            np.savetxt(temp_path, np.c_[np.round(MW.xint, 7), np.round(MW.yint,7)], delimiter=',')
            print('Saved truncated data as ' + label)
            MW.message = str('Saved truncated data as: \n' + label)
            MW.Err = False
            self.dialog = iFAMS_message(self)
            self.dialog.show()

        except AttributeError:
            print('No data to save')
            MW.message = 'No data to save'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            
    def restore_plot(self):
        MW.x = MW.xstore.copy()
        MW.y = MW.ystore.copy()
        
        MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
        MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)

        MW.maxF = float(MW.xFT[-1])
        MW.maxab = float(max(MW.yint))

        if MW.ft_pro == True:
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)#put in upper left hand corner
            MW.ax2 = MW.fig.add_subplot(122)

            if MW.MS_interp.isChecked():
                MW.ax.plot(MW.xint,MW.yint)
                MW.ax.set_title('Interpolated Mass Spectrum')
            else:
                MW.ax.plot(MW.x,MW.y)
                MW.ax.set_title('Raw Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')

            MW.ax2.plot(MW.xFT,abs(MW.yFT))
            MW.ax2.set_title('Fourier Spectrum')
            MW.ax2.set_ylabel('Relative Amplitude')
            MW.ax2.set_xlabel('Frequency')
            MW.ax2.set_xlim(0, MW.xFT[-1] / 2)
            MW.fig.tight_layout()
            MW.canvas.draw() #shows the graphs

        if MW.stft_pro == True:
            MW.winnum = MW.winnumstore
            MW.vmax = MW.vmaxstore    
            MW.plot_gabor(self)
       
        self.minx.setText(str(np.round(MW.x[0],4)))
        self.maxx.setText(str(np.round(MW.x[-1],4)))
        self.mzstep.setText(str(np.round(1/MW.maxF,4)))
        

class Noise(QtWidgets.QDialog):
    """
    This builds the window for adjusting the frequency domain over which the noise is determined
    """
    def __init__(self, parent=None):
        super(Noise, self).__init__(parent)
        self.setGeometry(1100, 100, 500, 600)
        self.setWindowTitle("Noise Calculator")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))

        self.fig = Figure()
        self.canvas3 = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(MW.xFT,abs(MW.yFT))
        self.ax.set_ylim(0,max(abs(MW.yFT)))
        self.ax.set_title('Fourier Spectrum')
        self.ax.set_ylabel('Amplitude')
        self.ax.set_xlabel('Frequency')
        self.canvas3.setFixedSize(500,400)

        self.Noise_tag = QtWidgets.QLabel()
        self.Noise_tag.setText('Frequency range for Noise determination:')
        self.Noise_tag.setToolTip('Enter frequency range that only contains white noise \n'
                                  'i.e. appears free of signal')

        self.minx_tag = QtWidgets.QLabel()
        self.minx_tag.setText('Minimum Frequency')
        self.minx_tag.setFixedSize(200,35)
        self.minx = QtWidgets.QTextEdit()
        self.minx.setFixedSize(150,35)
        self.minx.setText(str(np.round(MW.NFmin,4)))
        self.minx.setTabChangesFocus(True)
        
        self.maxx_tag = QtWidgets.QLabel()
        self.maxx_tag.setText('Maximum Frequency')
        self.maxx_tag.setFixedSize(200,35)
        self.maxx = QtWidgets.QTextEdit()
        self.maxx.setFixedSize(150,35)
        self.maxx.setText(str(np.round(MW.NFmax,4)))
        self.maxx.setTabChangesFocus(True)

        self.NFrms_tag = QtWidgets.QLabel()
        self.NFrms_tag.setText('Frequency Noise RMSD: ' + str(MW.NFrms))
        
        self.Nrms_tag = QtWidgets.QLabel()
        self.Nrms_tag.setText('Mass Noise RMSD: ' + str(MW.Nrms))

        self.button = QtWidgets.QPushButton('Replot Fourier')
        self.button.clicked.connect(self.replot)
        self.button.setToolTip('Displays Fourier spectrum in entered domain')

        self.button2 = QtWidgets.QPushButton('Update Noise Calculation')
        self.button2.clicked.connect(self.noise_calc)

        self.separator_tag = QtWidgets.QLabel()
        self.separator_tag.setText('Baseline Corrections to Mass Spectrum: ')

        self.button3 = QtWidgets.QPushButton('Subtract Average Noise')
        self.button3.clicked.connect(self.noise_sub)
        self.button3.setToolTip('Subtracts the input multiples of the mass noise RMSD from the \n'
                                'mass spectrum and recalculates the Fourier and Gabor spectra')

        self.button4 = QtWidgets.QPushButton('Undo Subtraction')
        self.button4.clicked.connect(self.sub_undo)
        
        self.button5 = QtWidgets.QPushButton('Linear Baseline Subtraction')
        self.button5.clicked.connect(self.linear_sub)
        self.button5.setToolTip('Subtracts the line connecting the minima within \n' 
                                'the first and last 20th of the mass spectrum \n' 
                                'and recalculates the Fourier and Gabor spectra')
        
        self.button6 = QtWidgets.QPushButton('Segmented Baseline Subtraction')
        self.button6.clicked.connect(self.segmented_sub)
        self.button6.setToolTip('Subtracts the lines connecting the minima within \n' 
                                '30 evenly spaced intervals in the mass spectrum and \n' 
                                'recalculates the Fourier and Gabor spectra. \n'
                                'SUGGESTION: update noise before adjusting the baseline')

        self.factor_tag = QtWidgets.QLabel()
        self.factor_tag.setText('   Multiples of Noise RMSD: ')
        self.factor_tag.setFixedSize(200,35)
        self.factor = QtWidgets.QTextEdit()
        self.factor.setFixedSize(100,35)
        self.factor.setText(str(3))
        self.factor.setTabChangesFocus(True)
        self.factor.setToolTip('Enter the number of RMSDs the noise extends \n'
                               'beyond the average. This is used to subtract \n'
                               'the average noise from the mass spectrum')

        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.canvas3,0,0,2,2)
        self.layout.addWidget(self.Noise_tag,3,0,1,2)
        self.layout.addWidget(self.minx_tag,4,0)
        self.layout.addWidget(self.minx,5,0)
        self.layout.addWidget(self.maxx_tag,4,1)
        self.layout.addWidget(self.maxx,5,1)
        self.layout.addWidget(self.NFrms_tag,6,0,1,2)
        self.layout.addWidget(self.Nrms_tag,7,0,1,2)
        self.layout.addWidget(self.button,8,0,1,2)
        self.layout.addWidget(self.button2,9,0,1,2)
        self.layout.addWidget(self.separator_tag,10,0,1,2)
        self.layout.addWidget(self.factor_tag,11,0)
        self.layout.addWidget(self.factor,11,1)
        self.layout.addWidget(self.button5,12,0)
        self.layout.addWidget(self.button6,12,1)
        self.layout.addWidget(self.button3,13,0)
        self.layout.addWidget(self.button4,13,1)

        MW.ycopy = MW.y.copy()
        MW.xcopy = MW.x.copy()

        self.setLayout(self.layout)
        self.canvas3.draw()
        
    def replot(self):
        self.xFT = []
        self.yFT = []
        for i in range(0,len(MW.xFT)):
            if MW.xFT[i] >= float(str(self.minx.toPlainText())) and MW.xFT[i] <= float(str(self.maxx.toPlainText())):
                self.xFT.append(MW.xFT[i])
                self.yFT.append(MW.yFT[i])
        self.yFT = np.array(self.yFT)
        
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(self.xFT,abs(self.yFT))
        self.ax.set_title('Fourier Spectrum')
        self.ax.set_ylabel('Amplitude')
        self.ax.set_xlabel('Frequency')
        self.canvas3.draw()
        
    def noise_calc(self):
        MW.NFmin = float(str(self.minx.toPlainText()))
        MW.NFmax = float(str(self.maxx.toPlainText()))
        MW.NFrms = cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin,MW.NFmax)
        MW.Nrms = MW.NFrms/np.sqrt(len(MW.x))
        self.NFrms_tag.setText('Frequency Noise RMSD: ' + str(np.round(MW.NFrms,4)))
        self.NFrms_tag.show()
        self.Nrms_tag.setText('Mass Noise RMSD: ' + str(np.round(MW.Nrms,4)))
        self.Nrms_tag.show()

        print('Frequency Noise RMSD: ' + str(MW.NFrms))
        print('Mass Noise RMSD: ' + str(MW.Nrms))

    def noise_sub(self):
        try:
            MW.RMSDmult = float(str(self.factor.toPlainText()))
            MW.NoiseAvg = MW.RMSDmult*MW.Nrms
            
            MW.y = MW.y - MW.NoiseAvg
            MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
            MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)
                        
            MW.plot_gabor(self)
    
            MW.Nsubtraction = True
        except TypeError:
            print('Unable to subtract noise. Update noise calculation')
            MW.message = 'Unable to subtract noise. Update noise calculation'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

    def sub_undo(self):        
        try:
            MW.y = MW.ycopy.copy()
            MW.x = MW.xcopy.copy()
            MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
            MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)
                        
            MW.plot_gabor(self)
    
            if MW.Nsubtraction == True:
                MW.Nsubtraction = False
            if MW.Nlinear_sub == True:
                MW.Nlinear_sub = False
            if MW.Nseg_sub == True:
                MW.Nseg_sub = False

        except AttributeError:
            print('No subtraction to undo')
            MW.message = 'No subtraction to undo'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

    def linear_sub(self):
        try:
            MW.globalmin1 = min(MW.y[0:int(len(MW.y)/20)])
            MW.globalmin2 = min(MW.y[-int(len(MW.y)/20):-1])
            dy = MW.globalmin2-MW.globalmin1
            dx = MW.x[-1]-MW.x[0]
            base_m = dy/dx
            base_b = MW.globalmin1-base_m*MW.x[0]
            for i in range(0,len(MW.y)):
                MW.y[i] = MW.y[i] - (base_m*MW.x[i]+base_b)
                
            MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
            MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)
                        
            MW.plot_gabor(self)

            MW.Nlinear_sub = True       
        except AttributeError:
            print('Unable to perform linear baseline correction')

    def segmented_sub(self):
        try:
            print('Starting segmented baseline subtraction')
            MW.mzminy = []
            MW.mzminx = []
            MW.mzslopes = []
            MW.mzintercepts = []
            
            MW.mzminy.append(min(MW.yint[0:int(len(MW.yint)/50)]))
            MW.mzminx.append(MW.xint[0])
            for i in range(1,30):
                MW.mzminy.append(min(MW.yint[int(len(MW.yint)*i/30-len(MW.yint)/100):int(len(MW.yint)*i/30+len(MW.yint)/100)]))
                MW.mzminx.append(MW.xint[int(len(MW.yint)*i/30)])
            MW.mzminy.append(min(MW.yint[-int(len(MW.yint)/50):-1]))
            MW.mzminx.append(MW.xint[-1])
            for i in range(0,30):
                dy = MW.mzminy[i+1]-MW.mzminy[i]
                dx = MW.mzminx[i+1]-MW.mzminx[i]
                MW.mzslopes.append(float(dy/dx))
                MW.mzintercepts.append(float(MW.mzminy[i]-(dy/dx)*MW.mzminx[i]))
            j = 0
            for i in range(0,len(MW.yint)):
                if MW.xint[i] > MW.mzminx[j+1]:
                    j += 1
                    if j == 30:
                        j -= 1
                MW.yint[i] = MW.yint[i] - (MW.mzslopes[j]*MW.xint[i]+MW.mzintercepts[j])
            
            MW.y = MW.yint
            MW.x = MW.xint
            
            MW.Nseg_sub = True
            
            MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
            MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)
                        
            MW.plot_gabor(self)

        except AttributeError:
            print('Unable to perform segmented baseline correction')

    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = MW.toggle_selector.RS.y1
            MW.y2 = MW.toggle_selector.RS.y2
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = MW.toggle_selector.RS.y1
            MW.y1 = MW.toggle_selector.RS.y2

class iso(QtWidgets.QDialog):
    """
    This builds the isotope distribution calculator window 
    UPDATES: 
        modified MW_4.protein_tag.setText, MW.4.protein.setToolTip, and MW_4.adductstip1_tag.setText
        new elif statements in def isotope_calculator for DNA/RNA
        added _temp objects to store initial FFT, IFFT and IFFTnorm
        added dirac delta spectrum for plotting of isotope fine structure
        adjusted MW.ax.plot lines to plot gaussian-smoothed distributions along with unsmoothed plot, after MW.fig.clf()
        see other updates in Isotope_Library_string.py and Isotope_Distr_Calc_string.py
    """
    def __init__(MW_4, parent=None):
        super(iso, MW_4).__init__(parent)
        MW_4.setGeometry(50, 50, 800, 800)
        MW_4.setWindowTitle("Molecular System Parameters")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        MW_4.protein_tag = QtWidgets.QLabel()
        MW_4.protein_tag.setText('Protein/Nucleotide Sequence \n'
                                 '[amino acid sequence (UPPERCASE)]\n'
                                 '[RNA/DNA sequence (r/d + lowercase)]')
        MW_4.protein_tag.setFixedSize(800,110)
        MW_4.protein = QtWidgets.QTextEdit()
        MW_4.protein.setText('')
        MW_4.protein.setFixedSize(800,200)
        MW_4.protein.setTabChangesFocus(True)
        MW_4.protein.setToolTip('Input amino acid sequence (UpperCase) or single strand RNA/DNA sequence (r-/d- lowercase) \n'
                                'of your protein or oligonucleotide monomeric unit, without spaces. \n'
                                'Use X or x for an unknown amino acid or nucleic base, respectively. \n'
                                'Masses determined at neutral pH with typical end caps (hydroxyls or DMT, depending). \n'
                                'All other modifications must be described in the "Adducts" input. \n'
                                'If your analyte is not a protein or oligonucleotide, leave this box blank.')
        
        MW_4.oligo_tag = QtWidgets.QLabel()
        MW_4.oligo_tag.setText('Oligomeric State')
        MW_4.oligo_tag.setFixedSize(800,35)
        MW_4.oligo = QtWidgets.QTextEdit()
        MW_4.oligo.setText('1')
        MW_4.oligo.setFixedSize(150,35)
        MW_4.oligo.setTabChangesFocus(True)
        MW_4.oligo.setToolTip('Input number of repeated monomeric units in your assembly. \n'
                              'If your analyte is not a protein or oligonucleotide, enter 1 in this box.')
        
        MW_4.adducts_tag = QtWidgets.QLabel()
        MW_4.adducts_tag.setText('Adduct Composition')
        MW_4.adducts_tag.setFixedSize(800,35)
        MW_4.adductstip1_tag = QtWidgets.QLabel()
        MW_4.adductstip1_tag.setText('(e.g. "Na 1 Cl 1")')
        MW_4.adductstip1_tag.setFixedSize(800,35)
        MW_4.adducts = QtWidgets.QTextEdit()
        MW_4.adducts.setText('')
        MW_4.adducts.setFixedSize(800,100)
        MW_4.adducts.setTabChangesFocus(True)
        MW_4.adducts.setToolTip('Input empirical formula of your adduct. \n'
                                'Leave spaces between each number and atomic symbol \n'
                                'and explicitly include 1s for single atoms. \n'
                                'If you are not including adducts, leave this box blank.')
        
        MW_4.cluster_tag = QtWidgets.QLabel()
        MW_4.cluster_tag.setText('Adduct Density (per monomer)')
        MW_4.cluster_tag.setFixedSize(800,35)
        MW_4.cluster = QtWidgets.QTextEdit()
        MW_4.cluster.setText('1')
        MW_4.cluster.setFixedSize(150,35)
        MW_4.cluster.setTabChangesFocus(True)
        MW_4.cluster.setToolTip('Input number of adducts (specified above) \n'
                                'per monomeric unit of your protein. \n'
                                'If you are not including adducts, enter 1 in this box.')
        
        MW_4.width_tag = QtWidgets.QLabel()
        MW_4.width_tag.setText('Set Resolution Width')
        MW_4.width_tag.setFixedSize(800,35)
        MW_4.width = QtWidgets.QTextEdit()
        MW_4.width.setText('0.1')
        MW_4.width.setFixedSize(150,35)
        MW_4.width.setTabChangesFocus(True)
        MW_4.width.setToolTip('Input the desired full-width at half-max for your isotope peaks. \n'
                                'Input 0 to keep isotope fine structure. \n'
                                'Note: cyclic nature of the data window \n'
                                'means tails near boundaries will bleed over')
        
        MW_4.button1 = QtWidgets.QPushButton('Calculate Mass Distribution')
        MW_4.button1.clicked.connect(MW_4.isotope_calculator)
        
        MW_4.button2 = QtWidgets.QPushButton('Export Mass Distribution')
        MW_4.button2.clicked.connect(MW_4.isotope_exporter)
        
        layout = QtWidgets.QGridLayout()
        layout.addWidget(MW_4.protein_tag)
        layout.addWidget(MW_4.protein)
        layout.addWidget(MW_4.oligo_tag)
        layout.addWidget(MW_4.oligo)
        layout.addWidget(MW_4.adducts_tag)
        layout.addWidget(MW_4.adductstip1_tag)
        ###layout.addWidget(MW_4.adductstip2_tag) REMOVED
        layout.addWidget(MW_4.adducts)
        layout.addWidget(MW_4.cluster_tag)
        layout.addWidget(MW_4.cluster)
        layout.addWidget(MW_4.width_tag)
        layout.addWidget(MW_4.width)
        layout.addWidget(MW_4.button1)
        layout.addWidget(MW_4.button2)
        MW_4.setLayout(layout)
        
    def isotope_calculator(MW_4):
        try:
            protein = str(MW_4.protein.toPlainText())
            protein = [protein[i] for i in range(len(protein))]
            multimer = int(str(MW_4.oligo.toPlainText()))
            add = str(MW_4.adducts.toPlainText())
            cluster = int(str(MW_4.cluster.toPlainText()))
            width = float(str(MW_4.width.toPlainText()))
            
            elem_list = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg',
                         'Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
                         'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
                         'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                         'In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',
                         'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
                         'Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
                         'At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu']
            
            sum_elem = np.zeros(94)
            
            if protein == []:
                multimer = 1
                subunit_type = ' subunits'
            elif protein[0] == 'r':
                sum_elem = IsoCalc.ssRNA_seq(sum_elem,protein)
                subunit_type = ' nucleotides'
            elif protein[0] == 'd':
                sum_elem = IsoCalc.ssDNA_seq(sum_elem,protein)
                subunit_type = ' nucleotides'
            else:
                sum_elem = IsoCalc.protein_seq(sum_elem,protein)
                subunit_type = ' amino acids'
            if add == []:
                cluster = 1
            else:
                sum_elem = IsoCalc.adduct(add,cluster,sum_elem)
            print("Monomer is "+str(len(protein))+subunit_type+"  in length")
            iso_min, iso_max, iso_avg, distr_max, distr_spc, mass_corr = IsoCalc.datawindow(sum_elem,elem_list,multimer)
            #print(distr_max, distr_spc,iso_min,iso_max,iso_avg)
            print("Finished Calculating Data Window")
            x = np.round(np.linspace(0, distr_max-0.001, distr_spc),decimals=3)
            y = np.zeros(len(x))
            mzspan = x[-1]-x[0]
            k = np.linspace(0, 2*np.pi/mzspan, len(x))
            xP = x + np.round(distr_max*np.floor(iso_min/distr_max), decimals=3)

            FFT = IsoCalc.FT_tot_string(x,y,k,sum_elem,elem_list,multimer)
            FFT_temp = FFT
            print("Finished Fourier domain calculations")
            
            if abs(width) > 0:
                sig = width/(2*np.sqrt(2*np.log(2)))
                gaussian = np.exp(-np.power(x,2) / (2*np.power(sig,2))) + np.exp(-np.power(x-x[-1],2) / (2*np.power(sig,2)))
                FFT = FFT*np.fft.fft(gaussian)
            else:
                dirac = np.zeros(len(x))
                dirac[0]+=1
                FFT = FFT*np.fft.fft(dirac)
            
            xrot,index_min,shiftdatafft,iso.IFFT = IsoCalc.shftdata(xP, FFT, iso_min, distr_max, distr_spc)
            iso.xrot1 = xrot + mass_corr
            print("Finished Data Counterrotation and IFFT")
            
            MWavg,MWvar = IsoCalc.distr_stats(iso.xrot1,iso.IFFT)
            print("Average Mass",MWavg)
            print("Standard Deviation",MWvar)
            print("Finished Statistical Calculations")
            
            iso.IFFT_norm1 = iso.IFFT/max(iso.IFFT)
            
            xrot,index_min,shiftdatafft,iso.IFFT_temp = IsoCalc.shftdata(xP, FFT_temp, iso_min, distr_max, distr_spc)
            iso.xrot2 = xrot + mass_corr
            iso.IFFT_norm2 = iso.IFFT_temp/max(iso.IFFT_temp)

            MW.iso_plot(MW_4)

        except ValueError:
            MW.message = 'Not enough parameters entered.  Please check each input value.'
            print(MW.message)
            MW.Err = True
            MW_4.dialog = iFAMS_message(MW_4)
            MW_4.dialog.show()
        except IndexError:
            MW.message = 'Not enough parameters entered.  Please check each input value.'
            print(MW.message)
            MW.Err = True
            MW_4.dialog = iFAMS_message(MW_4)
            MW_4.dialog.show()
        except KeyError:
            MW.message = 'Invalid character entered in sequence.  \nPlease remove invalid character or extra space and try again.'
            print(MW.message)
            MW.Err = True
            MW_4.dialog = iFAMS_message(MW_4)
            MW_4.dialog.show()
            
    def isotope_exporter(MW_4):
        MW.iso_export(MW_4)
            
class batch_menu(QtWidgets.QDialog):
    """
    builds the batch menu window
    """
    def __init__(self, parent = None):
        super(batch_menu, self).__init__(parent)
        self.setWindowTitle("Batch Menu")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.saving = QtWidgets.QLabel()
        self.saving.setText('In the new batch folder, processing parameters will \n'
                            'be saved along with folders for the deconvolved spectra \n'
                            'and integrated peaklists. \n'
                            'Additional data from each file can be saved if specified: ')
        self.cbinterp = QtWidgets.QCheckBox('Match original spectrum density')
        self.cbinterp.setToolTip('USE ONLY if the original spectrum is the sum-average of the batched spectra \n'
                                 '(e.g. sum across MS imaging pixels) \n'
                                 'Increases frequency range of the Fourier transform to match \n'
                                 'the original spectrum, avoiding STFT selection bound errors')
        self.cbFT = QtWidgets.QCheckBox('Save Fourier spectrum')
        self.cbIFFT = QtWidgets.QCheckBox('Save IFFT spectra')
        self.cbIFFT.setToolTip('Saves the inverse fast Fourier transform \n'
                               'for each charge-state of a series in a \n'
                               'single file that can be plotted as overlays \n'
                               'on top of the original mass spectrum.')
        self.cbZ = QtWidgets.QCheckBox('Save charge-state zero-charge spectra')
        self.cbZ.setToolTip('Saves the zero-charge spectrum for each charge state')
        self.cbseries = QtWidgets.QCheckBox('Save individual ion series deconvolution')
        self.cbseries.setToolTip('If multiple ion series selected, this saves the \n'
                                 'zero-charge spectrum for each ion series \n'
                                 'individually in addition to the combined \n'
                                 'zero-charge spectrum.')
                
        self.button1 = QtWidgets.QPushButton('Start Batch')
        self.button1.clicked.connect(self.Batch)
                 
        self.spacer = QtWidgets.QLabel()
        self.spacer.setFixedSize(300,5)
        
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.saving,0,0)
        self.layout.addWidget(self.cbFT,1,0)
        self.layout.addWidget(self.cbIFFT,2,0)
        self.layout.addWidget(self.cbZ,3,0)
        self.layout.addWidget(self.cbseries,4,0)
        #self.layout.addWidget(self.cbinterp,5,0)
        self.layout.addWidget(self.spacer,6,0)
        self.layout.addWidget(self.button1,7,0)
        self.setLayout(self.layout)
    
    def Batch(self):
        if self.cbFT.isChecked() == True:
            MW.saveFT = True
        else:
            MW.saveFT = False
        if self.cbIFFT.isChecked() == True:
            MW.saveIFFT = True
        else:
            MW.saveIFFT = False
        if self.cbZ.isChecked() == True:
            MW.saveZspec = True
        else:
            MW.saveZspec = False
        if self.cbseries.isChecked() == True:
            MW.saveSeries = True
        else:
            MW.saveSeries = False
        if self.cbinterp.isChecked() == True:
            MW.setdatarho = True
        else:
            MW.setdatarho = False
        
        self.close()
        MW.batch_function(self)   

class iFAMS_message(QtWidgets.QDialog):
    """
    creates a dialog box when various tasks are completed or unable to complete
    """
    def __init__(self, parent=None):
        super(iFAMS_message, self).__init__(parent)
        if MW.Err == True:
            self.setWindowTitle("Oops")
        else:
            self.setWindowTitle("iFAMS Message")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.message = QtWidgets.QLabel()
        self.message.setText(MW.message)
        self.message.setFont(QtGui.QFont('Arial',8)) 
        self.message.setAlignment(QtCore.Qt.AlignCenter)
               
        self.spacer = QtWidgets.QLabel()
        self.spacer.setFixedSize(300,5)
        
        self.button1 = QtWidgets.QPushButton('OK')
        self.button1.clicked.connect(self.OK)
                       
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.message)
        try:
            self.submessage = QtWidgets.QLabel()
            self.submessage.setText(MW.submessage)
            self.submessage.setFont(QtGui.QFont('Arial',7)) 
            self.submessage.setAlignment(QtCore.Qt.AlignCenter)
            self.layout.addWidget(self.submessage)
            del(MW.submessage)
        except AttributeError:
            pass
        self.layout.addWidget(self.spacer)
        self.layout.addWidget(self.button1)
        self.setLayout(self.layout)
        MW.fig.tight_layout()
        MW.canvas.draw()
    
    def OK(self):
        self.hide()
        self.close()
      
class Folder_Select(QtWidgets.QDialog):     
    """
    creates a dialog box for selecting a location to save data
    """
    def __init__(self, parent=None):
        super(Folder_Select, self).__init__(parent)
        self.setWindowTitle("Save Menu")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.message = QtWidgets.QLabel()
        self.message.setText("Select save destination")
        self.message.setFont(QtGui.QFont('Arial',10)) 
        self.message.setAlignment(QtCore.Qt.AlignCenter)
        
        self.message2 = QtWidgets.QLabel()
        self.message2.setText("Data will be saved to a new folder in the file's directory.\n"
                             "The folder will be named after the 'File Name' shown below. ")
        self.message2.setFont(QtGui.QFont('Arial',8)) 
        self.message2.setAlignment(QtCore.Qt.AlignCenter)
        try:
            self.namebase = os.path.basename(MW.namebase)
            self.directory = os.path.dirname(MW.namebase)
            if self.directory == '':
                self.directory = 'None'
        except AttributeError:
            self.namebase = 'None'
            self.directory = 'None'
        
        self.submessageD = QtWidgets.QLabel()
        self.submessageD.setText('Directory Name: '+str(self.directory))
        self.submessageD.setFont(QtGui.QFont('Arial',8)) 
        self.submessageD.setAlignment(QtCore.Qt.AlignCenter)
        self.submessageD.setToolTip("Determined by the 'file' location and \n"
                                    "will be the location for the new save folder." )
        
        self.submessageF = QtWidgets.QLabel()
        self.submessageF.setText('File Name: '+str(self.namebase))
        self.submessageF.setFont(QtGui.QFont('Arial',8)) 
        self.submessageF.setAlignment(QtCore.Qt.AlignCenter)
        self.submessageF.setToolTip("Used to name the new save folder. \n"
                                    "There does not need to be an actual file with this name.")
        
        self.spacer = QtWidgets.QLabel()
        self.spacer.setFixedSize(300,5)
        self.savemess = QtWidgets.QLabel()
        self.savemess.setText('The deconvolved spectrum (combined full zero-charge spectrum), \n'
                              'integrated peaklist, and processing parameters are saved to the folder\n'
                              'by default. Additional data can be selected for export below:')
        self.savemess.setFont(QtGui.QFont('Arial',8)) 
        self.savemess.setAlignment(QtCore.Qt.AlignCenter)
        
        self.cbFT = QtWidgets.QCheckBox('Save Fourier spectrum')
        self.cbIFFT = QtWidgets.QCheckBox('Save IFFT spectra')
        self.cbIFFT.setToolTip('Saves the inverse fast Fourier transform \n'
                               'for each charge-state of a series in a \n'
                               'single file that can be plotted as overlays \n'
                               'on top of the original mass spectrum.')
        self.cbZ = QtWidgets.QCheckBox('Save charge-state zero-charge spectra')
        self.cbZ.setToolTip('Saves the zero-charge spectrum for each charge state')
        self.cbseries = QtWidgets.QCheckBox('Save individual ion series deconvolution')
        self.cbseries.setToolTip('If multiple ion series selected, this saves the \n'
                                 'zero-charge spectrum for each ion series \n'
                                 'individually in addition to the combined \n'
                                 'zero-charge spectrum.')
        
        self.button1 = QtWidgets.QPushButton('Select Save Destination')
        self.button1.clicked.connect(self.select)
        
        self.button2 = QtWidgets.QPushButton('Accept')
        self.button2.clicked.connect(self.OK)
                       
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.message)
        self.layout.addWidget(self.message2)
        self.layout.addWidget(self.submessageD)
        self.layout.addWidget(self.submessageF)
        if MW.full_save == True:
            self.layout.addWidget(self.spacer)
            self.layout.addWidget(self.savemess)
            self.layout.addWidget(self.cbFT)
            self.layout.addWidget(self.cbIFFT)
            self.layout.addWidget(self.cbZ)
            self.layout.addWidget(self.cbseries)
        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)
        self.setLayout(self.layout)

    
    def select(self):
        self.namebase, extra = QtWidgets.QFileDialog.getSaveFileName(self, 'Select Save Destination')
        if self.namebase == '':
            pass
        else:
            MW.namebase = self.namebase
            self.namebase = os.path.basename(MW.namebase)
            self.directory = os.path.dirname(MW.namebase)
            self.submessageD.setText('Directory Name: '+str(self.directory))
            self.submessageF.setText('File Name: '+str(self.namebase))
            self.submessageD.show()
            self.submessageF.show()
    
    def OK(self):
        if self.directory == 'None':
            MW.message = 'Please select folder first'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
        if self.cbFT.isChecked() == True:
            MW.saveFT = True
        else:
            MW.saveFT = False
        if self.cbIFFT.isChecked() == True:
            MW.saveIFFT = True
        else:
            MW.saveIFFT = False
        if self.cbZ.isChecked() == True:
            MW.saveZspec = True
        else:
            MW.saveZspec = False
        if self.cbseries.isChecked() == True:
            MW.saveSeries = True
        else:
            MW.saveSeries = False
        self.hide()
        self.close()
        MW.ReturnFunc(self)
    

class guided_dialog(QtWidgets.QDialog):
    """
    creates a dialog box for the guided search tool
    """
    def __init__(self, parent=None):
        super(guided_dialog, self).__init__(parent)
        self.setWindowTitle("iFAMS Guided Search")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.setGeometry(1100, 400, 400, 200)
        
        self.message = QtWidgets.QLabel()
        self.message.setText('Select signal or enter ion mass, then press Go')
        self.message.setFont(QtGui.QFont('Arial',8)) 
        self.message.setAlignment(QtCore.Qt.AlignCenter)

        self.submessage = QtWidgets.QLabel()
        self.submessage.setText('To select signal, add and save boxes in the Gabor spectrogram \n'
                                'of two adjacent signals in the charge state series of interest. \n'
                                'Signal can be a fundamental frequency or near-zero frequencies. ')
        self.submessage.setFont(QtGui.QFont('Arial',7)) 
        self.submessage.setAlignment(QtCore.Qt.AlignCenter)
        self.submessage.setToolTip('For best results, select signal at or near the center \n'
                                   'of the distribution and as far away from neighboring \n'
                                   'distributions as possible.')
               
        self.ion_tag = QtWidgets.QLabel()
        self.ion_tag.setText('Ion Mass(es) (Da):')
        self.ion_tag.setFixedSize(400,25)
        self.ion = QtWidgets.QTextEdit()
        self.ion.setFixedSize(500,100)
        self.ion.setTabChangesFocus(True)
        self.ion.setToolTip('Enter the most abundant mass of the ion of interest in Da. \n'
                            'Multiple ion masses can be entered if separated by a comma, \n'
                            'for example: 8565,12360,148200 \n'
                            'Entering multiple ions will skip the charge state toggle \n'
                            'menu but will automatically use the recommended charge states\n'
                            'and harmonics.')
        
        self.cbBox = QtWidgets.QCheckBox('Automatically Adjust Boxes to Signal')
        self.cbBox.setChecked(True)
        self.cbBox.setToolTip('If checked, each selection box will be adjusted \n'
                              'to fit the signal based on abundances.\n'
                              'If unchecked, an average box size will be used for \n'
                              'every charge state with only slight adjustments for \n'
                              'expected spreading/compressing of signal in the \n'
                              'mass spectrum')
        
        self.cbRange = QtWidgets.QCheckBox('Automatically Adjust Range of Selection')
        self.cbRange.setChecked(True)
        self.cbRange.setToolTip('If checked, the range of charge states selected will \n'
                              'be determined by thresholds relative to the most \n'
                              'abundant charge state. This option may result in less \n'
                              'charge states selected. \n'
                              'If unchecked, every possible charge state in the mass \n'
                              'domain will be selected. Then, the charge states to be \n'
                              'included in the deconvolution can be manually toggled. \n'
                              'It may be necessary to truncate the mass domain before \n'
                              'selecting this option.')
        
        self.button1 = QtWidgets.QPushButton('Go')
        self.button1.clicked.connect(self.Go)
        self.button1.setToolTip('For best results, select signal at or near the center \n'
                                   'of the distribution and as far away from neighboring \n'
                                   'distributions as possible.')
        
        self.button2 = QtWidgets.QPushButton('Finish Search')
        self.button2.clicked.connect(self.Finish)
        self.button2.setToolTip('Finalizes deconvolutions of saved selections without searching \n'
                                'for an additional series.')
                       
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.message)
        self.layout.addWidget(self.submessage)
        self.layout.addWidget(self.ion_tag)
        self.layout.addWidget(self.ion)
        self.layout.addWidget(self.cbBox)
        self.layout.addWidget(self.cbRange)
        self.layout.addWidget(self.button1)
        if len(MW.batchstore) > 0:
            self.layout.addWidget(self.button2)
        self.setLayout(self.layout)
    
    def Go(self):
        try:
            try:
                MW.ionmass = str(self.ion.toPlainText())
            except AttributeError:
                MW.ionmass = False
            if MW.ionmass == '':
                MW.ionmass = False
            else:
                try:
                    MW.ionmass = float(MW.ionmass)
                except ValueError:
                    MW.ionmass = MW.ionmass.split(',')
                except TypeError:
                    MW.message = 'Unexpected value for ion mass'
                    MW.Err = True
                    MW.submessage = 'Please enter a number and try again'
                    self.dialog = iFAMS_message(self)
                    self.dialog.show()
            if MW.ionmass == False and len(MW.rlist) < 2:
                MW.message = 'No mass or signal selection for search'
                MW.Err = True
                MW.submessage = 'Please select signal or enter ion mass and try again'
                self.dialog = iFAMS_message(self)
                self.dialog.show()
            else:
                if self.cbBox.isChecked():
                    MW.boxadjust = True
                else:
                    MW.boxadjust = False
                if self.cbRange.isChecked():
                    MW.rangeadjust = True
                else:
                    MW.rangeadjust = False
                MW.stft_selec_man = False
                MW.stft_box_thresh = MW.boxadjust
                MW.stft_range_thresh = MW.rangeadjust
                self.close()
                try:
                    if len(MW.ionmass) > 1:
                        masslist = MW.ionmass.copy()
                        MW.autosearch = True
                        for i in range(len(masslist)):
                            MW.ionmass = float(masslist[i])
                            MW.rlist = []
                            MW.glist = []
                            MW.zlist = []
                            MW.Go_search(self)
                        MW.auto_check_open(self)
                except TypeError:
                    MW.Go_search(self)
        except AttributeError:
            self.ErrReset()
        except ValueError:
            self.ErrReset()
        except IndexError:
            self.ErrReset()
        except TypeError:
            self.ErrReset()    
                
    def ErrReset(self):
        del(MW.ionmass)
        MW.plot_gabor(self)
        for j in range(len(MW.batchstore)):
            for i in range(len(MW.batchstore[j][1])):    
                rex = Rectangle((MW.batchstore[j][1][i][0],MW.batchstore[j][1][i][2]), (MW.batchstore[j][1][i][1] - MW.batchstore[j][1][i][0]), (MW.batchstore[j][1][i][3] - MW.batchstore[j][1][i][2]),
                            facecolor='none', edgecolor=color[-2*(j+1)], linewidth=1)
                MW.ax.add_artist(rex) 
        self.close()
        MW.guided_search(self)
        MW.fig.tight_layout()
        MW.canvas.draw()
        MW.message = 'Unable to complete search'
        MW.Err = True
        MW.submessage = 'Please adjust search parameters and try again.'
        self.dialog = iFAMS_message(self)
        self.dialog.show()

    def Finish(self):
        MW.fig.clf()
        MW.fig.tight_layout()
        MW.canvas.draw()
        MW.ax = MW.fig.add_subplot(121)
        MW.ax2 = MW.fig.add_subplot(122)
        
        MW.ax.plot(MW.xint,MW.yint,label='Mass Spectrum',color='darkgray')
        for i in range(len(MW.batchstore)):
            label = 'Series '+str(i+1)
            MW.ax.plot(MW.xint,MW.batchstore[i][2],label=label)
        MW.ax.legend(loc='upper right', title='Selections')
        MW.ax.set_title('Mass Spectrum')
        MW.ax.set_xlabel('m/z')
        MW.ax.set_ylabel('Relative Abundance')
        
        min0 = []
        max0 = []
        for i in range(len(MW.batchstore)):
            min0.append(np.min(MW.batchstore[i][3]))
            max0.append(np.max(MW.batchstore[i][3]))
        min0 = np.min(min0)
        max0 = np.max(max0)
        space0 = MW.xzero[2]-MW.xzero[1]
        MW.xzero = np.linspace(min0,max0,int(2*np.round((max0-min0)/space0)),endpoint=True)
        y0s = []
        MW.deconylist = []
        MW.deconxlist = []
        for i in range(len(MW.batchstore)):
            y0s.append(np.interp(MW.xzero,MW.batchstore[i][3],MW.batchstore[i][4]))
            MW.deconxlist.append(MW.batchstore[i][3])
            MW.deconylist.append(MW.batchstore[i][4])
        MW.zerofull = np.sum(y0s,axis=0)
        MW.zero_norm = MW.zerofull/max(MW.zerofull)
        MW.ax2.plot(MW.xzero,MW.zerofull, label='Combined',color='darkgray',linewidth=3)
        for i in range(len(MW.batchstore)):
            label = 'Series ' + str(i+1)
            MW.ax2.plot(MW.batchstore[i][3],MW.batchstore[i][4],label=label)
            for j in range(len(MW.batchstore[i][5])):
                MW.ax2.plot(MW.batchstore[i][3],MW.batchstore[i][5][j],color=color[j+1],alpha=0.5)
        MW.ax2.legend(loc='upper right', title='Deconvolutions')
        MW.ax2.set_title('Zero Charge Spectra')
        MW.ax2.set_ylabel('Relative Abundance')
        MW.ax2.set_xlabel('Mass (Da)')
                
        MW.deltaY = 5
        MW.minimumY = 0.03
        MW.minimumX = MW.xzero[0]
        MW.maximumX = MW.xzero[-1]
        MW.ntol = 0.0
        massrange = np.nonzero(MW.zerofull)
        try:
            MW.xlist,MW.ylist,MW.xref,MW.yref,MW.sumlist,MW.peakindex,MW.minL,MW.minR,MW.xcent,MW.height = Int.integrate(MW.zero_norm[massrange],MW.zerofull[massrange],MW.xzero[massrange], MW.deltaY, MW.minimumY,MW.minimumX, MW.maximumX, MW.ntol)
            
            peaktxt = []
            sumlist = MW.sumlist.copy()
            for i in range(0, min(8,len(sumlist))):
                peaktxt.append(np.argmax(sumlist))
                sumlist[np.argmax(sumlist)] = 0
            for i in peaktxt:
                if len(MW.batchstore) == 0:
                    txty = MW.height[i]*max(MW.BL_corrected)
                else:
                    txty = MW.height[i]*max(MW.zerofull)
                MW.ax2.text(MW.xcent[i],txty,str(np.round(MW.xcent[i],3)),rotation='vertical',ha='center',clip_on=True)
    
            max_value = max(MW.zerofull)
            for i in range(0,len(MW.sumlist)):
                MW.sumlist[i] = max_value*MW.sumlist[i]
                MW.height[i] = max_value*MW.height[i]
        except IndexError:
            pass

        rtemp = []
        cstemp = []
        gtemp = []
        for i in range(len(MW.batchstore)):
            cstemp.append(MW.batchstore[i][0])
            rtemp.append(MW.batchstore[i][1])
            gtemp.append(MW.batchstore[i][2])
        MW.cs = np.concatenate((cstemp))
        MW.rlist = np.concatenate((rtemp))
        MW.glist = gtemp

        if MW.realsel.isChecked() == False:
            MW.tIFT = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.tIFT = np.real(MW.glist)

        #Noise Calculations
        try:
            MW.selecArea = 0
            MW.GaborArea = (MW.xint[-1]-MW.xshort[0])*(MW.xFT[-1]-MW.xFT[0])
            for i in range(0,len(MW.rlist)):
                MW.selecArea += (MW.rlist[i][1]-MW.rlist[i][0])*(MW.rlist[i][3]-MW.rlist[i][2])*2
            MW.GaborFraction = float(MW.selecArea/MW.GaborArea)
    
            MW.integratedx = 0  
            MW.peakwidths = np.zeros(len(MW.xref))
            for i in range(0, len(MW.xref)):
                MW.integratedx += len(MW.xref[i])
                MW.peakwidths[i] = len(MW.xref[i])
            MW.IntFraction = MW.integratedx/len(MW.xzero[massrange])
            MW.N = np.sqrt(MW.IntFraction)*np.sqrt(MW.GaborFraction)*MW.NFrms*2*np.pi
            MW.Nh = np.sqrt(MW.GaborFraction)*MW.NFrms/np.sqrt(len(MW.xzero[massrange]))*2*np.pi
            if MW.NFrms > 0:
                print('Fraction of Gabor used: ' + str(np.round(MW.GaborFraction,5)))
                print('Fraction of Zero Charge Spectrum integrated: ' + str(np.round(MW.IntFraction,5)))
                print('Noise for integrated spectrum: ' + str(np.round(MW.N,3)))
            
            MW.paramMatch = False
            MW.bounds = []
            for i in range(0,len(MW.xref)):
                MW.bounds.append(MW.xref[i][0])
                MW.bounds.append(MW.xref[i][-1])
        except AttributeError:
            pass
        
        self.close()
        MW.fig.tight_layout()
        MW.canvas.draw()

    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y2 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y1 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
    
class guided_decharger(QtWidgets.QDialog):
    """
    guided search tool dialog box for editing charge state selections/harmonics
    """
    def __init__(self, parent=None):
        super(guided_decharger, self).__init__(parent)
        self.setWindowTitle("iFAMS Guided Search")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.setGeometry(1200, 100, 400, 200)
        self.layout = QtWidgets.QGridLayout()
                
        if len(MW.zlist) < 16:
            self.message = QtWidgets.QLabel('Select charge states to \n'
                                            'include in deconvolution: ')
            self.note = QtWidgets.QLabel("*Selections near ends of spectrum \n"
                                         "might create windowing artifacts")   
        else:
            self.message = QtWidgets.QLabel('Select charge states to include in deconvolution: ')
            self.note = QtWidgets.QLabel("*Selections near ends of spectrum might create windowing artifacts")
        self.message.setFont(QtGui.QFont('Arial',8)) 
        self.note.setFont(QtGui.QFont('Arial',8))
        self.note.setToolTip('Due to the Gaussian windowing, IFFTs extend further in the m/z direction \n'
                             'than their corresponding pixels in the Gabor spectrogram. If a selection \n'
                             'is too close to the end of the mass spectrum, it might generate windowing \n'
                             'artifacts that can interfere with the automatic baseline correction.\n'
                             'If a deconvolution seems to fail, try reprocessing after deselecting \n'
                             'charge states near the ends of the spectrum.')

        self.layout.addWidget(self.message,0,0,1,2)
        self.layout.addWidget(self.note,1,0,1,2)

        self.cbCS = np.zeros(len(MW.zlist)).tolist()
        count = 0
        for i in range(len(MW.zlist)):
            if i == 15*(count+1):
                count+=1
            self.cbCS[i] = QtWidgets.QCheckBox('Charge State '+str(MW.zlist[i]))
            self.cbCS[i].setChecked(True)
            self.layout.addWidget(self.cbCS[i],i+2-15*count,count)

        try:
            MW.Harmrec
            self.harm_tag = QtWidgets.QLabel()
            self.harm_tag.setText('Highest Harmonic to Include')
            self.layout.addWidget(self.harm_tag,17,0)
            
            self.harm = QtWidgets.QTextEdit()
            self.harm.setFixedSize(100,35)
            self.harm.setTabChangesFocus(True)
            self.harm.setText(str(MW.Harmrec))
            self.harm.setToolTip('iFAMS selects harmonics up to the entered value.\n'
                                'The fundamentals are considered the first harmonic (1). \n'
                                'Value entered must be an integer.')
            self.layout.addWidget(self.harm,18,0)
        
            self.cbZeroHarm = QtWidgets.QCheckBox('Include Near-zero Frequencies')
            self.cbZeroHarm.setChecked(True)
            self.cbZeroHarm.setToolTip('Including the near-zero frequencies is recommended \n'
                                       'unless there is a lot of overlap in m/z between \n'
                                       'different species.')
            self.layout.addWidget(self.cbZeroHarm,19,0)
            
            self.button1 = QtWidgets.QPushButton('Run iFAMS Deconvolution')
            self.button1.clicked.connect(self.iFAMS_Decon)
            self.layout.addWidget(self.button1,20,0,1,3)
            
            self.button2 = QtWidgets.QPushButton('Deconvolve and Select Additional Series')
            self.button2.clicked.connect(self.Decon_2)
            self.button2.setToolTip('Finishes and saves current selection, then resets spectrogram \n'
                                    'for selecting another charge state series.')
            self.layout.addWidget(self.button2,21,0,1,3)
            
            self.button3 = QtWidgets.QPushButton('Adjust Assigned Charges')
            self.button3.clicked.connect(self.CS_Adjust)
            self.button3.setToolTip('Opens window to help adjust charges from the suggested \n'
                                    'values to the correct values.')
            self.layout.addWidget(self.button3,22,0,1,3)
            
            self.button4 = QtWidgets.QPushButton('Remove Selection')
            self.button4.clicked.connect(self.Remove)
            self.button4.setToolTip('Removes current search results from saved selections and \n'
                                    'restores previous saved selections, if any.')
            self.layout.addWidget(self.button4,23,0,1,3)
            
        except AttributeError:
            self.button1 = QtWidgets.QPushButton('Run iFAMS Deconvolution')
            self.button1.clicked.connect(self.iFAMS_Decon)
            self.layout.addWidget(self.button1,17,0,1,3)
        
            self.button2 = QtWidgets.QPushButton('Deconvolve and Select Additional Series')
            self.button2.clicked.connect(self.Decon_2)
            self.button2.setToolTip('Finishes and saves current selection, then resets spectrogram \n'
                                    'for selecting another charge state series.')
            self.layout.addWidget(self.button2,18,0,1,3)
            
            self.button3 = QtWidgets.QPushButton('Adjust Assigned Charges')
            self.button3.clicked.connect(self.CS_Adjust)
            self.button3.setToolTip('Opens window to help adjust charges from the suggested \n'
                                    'values to the correct values.')
            self.layout.addWidget(self.button3,19,0,1,3)
            
            self.button4 = QtWidgets.QPushButton('Remove Selection')
            self.button4.clicked.connect(self.Remove)
            self.button4.setToolTip('Removes current search results from saved selections and \n'
                                    'restores previous saved selections, if any.')
            self.layout.addWidget(self.button4,20,0,1,3)
                    
        self.setLayout(self.layout)
    
    def iFAMS_Decon(self):
        'Finishes Guided Search and performs the final iFAMS analysis'
        try:
            MW.Harmrec
            try:
                MW.Harmrec = int(self.harm.toPlainText())
            except TypeError:
                MW.message = 'Unexpected value for harmonic number'
                MW.Err = True
                MW.submessage = 'Please enter an integer and try again'
                self.dialog = iFAMS_message(self)
                self.dialog.show()
                return
            if self.cbZeroHarm.isChecked():
                start = 0
            else:
                start = 1
            ztemp = []
            rtemp = []
            Itemp = []
            gtemp = []
            gBL2temp = []
            MW.Xtemp = MW.X*0
            MW.XtempBL2 = MW.XBL2*0
            for i in range(len(self.cbCS)):
                if self.cbCS[i].isChecked():
                    ztemp.append(MW.zlist[i])
                    for j in range(start,MW.Harmrec+1):
                        if j == 0:
                            bottom = 0
                        else:
                            bottom = int(np.round(j*MW.Ilist[i][2]+(j-1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                        top = int(np.round(j*MW.Ilist[i][2]+(j+1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                        rtemp.append([MW.rlist[i][0],MW.rlist[i][1],MW.yshort[bottom],MW.yshort[top],len(gtemp)])
                        Itemp.append([MW.Ilist[i][0],MW.Ilist[i][1],bottom,top,len(gtemp)])
                        if j == 0:
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1]
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])]
                        else:
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1]
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1]
                    MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                    gtemp.append(MW.yrecon)
                    MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                    gBL2temp.append(MW.yreconBL2)
                    MW.Xtemp = MW.X*0
                    MW.XtempBL2 = MW.XBL2*0

        except AttributeError:
            start = 0
            ztemp = []
            rtemp = []
            Itemp = []
            gtemp = []
            gBL2temp = []
            for i in range(len(self.cbCS)):
                if self.cbCS[i].isChecked():
                    ztemp.append(MW.zlist[i])
                    rtemp.append(MW.rlist[i])
                    rtemp[-1][4] = len(rtemp)-1
                    Itemp.append(MW.Ilist[i])
                    Itemp[-1][4] = len(Itemp)-1
                    gtemp.append(MW.glist[i])
                    gBL2temp.append(MW.glistBL2[i])

        MW.zlist = ztemp.copy()
        MW.cs = MW.zlist
        MW.rlist = rtemp.copy()
        MW.Ilist = Itemp.copy()
        MW.glist = gtemp.copy()
        MW.glistBL2 = gBL2temp.copy()
        
        if MW.realsel.isChecked() == False:
            MW.tIFT = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.tIFT = np.real(MW.glist)
        
        MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(MW.tIFT, MW.xint, MW.cs,MW.proton)
        MW.xzeroBL2, MW.cszeroBL2, MW.zerofullBL2 = FT.zerocharge(np.real(MW.glistBL2), MW.xint, MW.cs,MW.proton)
       
        massrange=np.nonzero(MW.zerofullBL2)
        MW.xzero = MW.xzero[massrange]
        MW.zerofull = MW.zerofull[massrange]
        MW.xzeroBL2 = MW.xzeroBL2[massrange]
        MW.zerofullBL2 = MW.zerofullBL2[massrange]
        MW.BSint = []
        for j in range(len(MW.cszero)):
            MW.cszero[j] = MW.cszero[j][massrange]
            MW.cszeroBL2[j] = MW.cszeroBL2[j][massrange]
            check = np.max(MW.cszeroBL2[j])/4
            if check == 0:
                continue
            mid1 = int(len(MW.xzero)/8)+1
            mid2 = int(len(MW.xzero)*7/8)-1
            for k in range(len(MW.cszero[j])):
                if MW.cszeroBL2[j][k] > check:
                    mid1 = k
                    break
            for k in range(mid1+1,len(MW.cszero[j])):
                if MW.cszeroBL2[j][k] <= check or MW.cszeroBL2[j][k] <= MW.cszeroBL2[j][mid1]:
                    mid2 = k
                    break
            int1 = min(sum(MW.cszero[j][0:mid1]),sum(MW.cszero[j][mid2:len(MW.xzero)]))
            int2 = (sum(MW.cszeroBL2[j][0:mid1])+sum(MW.cszeroBL2[j][mid2:len(MW.xzero)]))/2
            MW.BSint.append(int1/int2)
        BSint = np.average(MW.BSint)
        MW.zerofullBL2 = MW.zerofullBL2*BSint
        MW.BL_norm = MW.zerofullBL2/max(MW.zerofull)
        MW.BL_corrected = (MW.zerofull - MW.zerofullBL2)/max(MW.zerofull)
        MW.zerofull0 = MW.zerofull
        MW.zerofull = MW.zerofull - MW.zerofullBL2
        
        MW.fig.clf()
        MW.fig.tight_layout()
        MW.canvas.draw()
        MW.ax = MW.fig.add_subplot(121)
        MW.ax2 = MW.fig.add_subplot(122)
        
        MW.ax.plot(MW.xint,MW.yint,label='Mass Spectrum',color='darkgray')
        if len(MW.batchstore) == 0:
            for i in range(0,len(MW.cs)):
                label = 'Charge ' + str(MW.cs[i])
                MW.ax.plot(MW.xint,MW.tIFT[i],label=label)
        else:
            MW.glist = np.sum(MW.tIFT,axis=0)
            for i in range(len(MW.batchstore)):
                label = 'Series '+str(i+1)
                MW.ax.plot(MW.xint,MW.batchstore[i][2],label=label)
            label = 'Series '+str(len(MW.batchstore)+1)
            MW.ax.plot(MW.xint,MW.glist,label=label)
        MW.ax.legend(loc='upper right', title='Selections')
        MW.ax.set_title('Mass Spectrum')
        MW.ax.set_xlabel('m/z')
        MW.ax.set_ylabel('Relative Abundance')
        
        if len(MW.batchstore) == 0:
            MW.cs_norm = np.zeros(len(MW.cszero)).tolist()
            MW.zero_norm = MW.zerofull0/max(MW.zerofull)
            MW.ax2.plot(MW.xzero,MW.zero_norm, label='Full zero-charge',color='darkgray')
            try:
                MW.ax2.plot(MW.xzero,MW.BL_norm, label='Modeled Baseline',color='orange')
                MW.ax2.plot(MW.xzero,MW.BL_corrected, label='Baseline-corrected',color='mediumseagreen')
            except AttributeError:
                pass
            for i in range(len(MW.cszero)):
                label = 'Charge '+str(MW.cs[i])
                MW.cs_norm[i] = MW.cszero[i]/max(MW.zerofull)
                MW.ax2.plot(MW.xzero,MW.cs_norm[i],label=label,alpha=0.5)
            MW.ax2.legend(loc='upper right', title='Charge State')
            MW.ax2.set_title('Zero Charge Spectra')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
        
        else:
            storetemp = []
            storetemp.append(MW.cs)
            storetemp.append(MW.rlist)
            storetemp.append(MW.glist)
            storetemp.append(MW.xzero)
            storetemp.append(MW.zerofull)
            storetemp.append(MW.cszero)
            storetemp.append(MW.cszeroBL2)
            storetemp.append(len(MW.tIFT))
            storetemp.append(MW.tIFT)
            MW.batchstore.append(storetemp)
            min0 = []
            max0 = []
            for i in range(len(MW.batchstore)):
                min0.append(np.min(MW.batchstore[i][3]))
                max0.append(np.max(MW.batchstore[i][3]))
            min0 = np.min(min0)
            max0 = np.max(max0)
            space0 = MW.xzero[2]-MW.xzero[1]
            MW.xzero = np.linspace(min0,max0,int(2*np.round((max0-min0)/space0)),endpoint=True)
            y0s = []
            MW.deconylist = []
            MW.deconxlist = []
            for i in range(len(MW.batchstore)):
                y0s.append(np.interp(MW.xzero,MW.batchstore[i][3],MW.batchstore[i][4]))
                MW.deconxlist.append(MW.batchstore[i][3])
                MW.deconylist.append(MW.batchstore[i][4])
            MW.zerofull = np.sum(y0s,axis=0)
            MW.zero_norm = MW.zerofull/max(MW.zerofull)
            MW.ax2.plot(MW.xzero,MW.zerofull, label='Combined',color='darkgray',linewidth=3)
            for i in range(len(MW.batchstore)):
                label = 'Series ' + str(i+1)
                MW.ax2.plot(MW.batchstore[i][3],MW.batchstore[i][4],label=label)
                for j in range(len(MW.batchstore[i][5])):
                    MW.ax2.plot(MW.batchstore[i][3],MW.batchstore[i][5][j],color=color[j+1],alpha=0.5)
            MW.ax2.legend(loc='upper right', title='Deconvolutions')
            MW.ax2.set_title('Zero Charge Spectra')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
                
        MW.deltaY = 5
        MW.minimumY = 0.03
        MW.minimumX = MW.xzero[0]
        MW.maximumX = MW.xzero[-1]
        MW.ntol = 0.0
        massrange = np.nonzero(MW.zerofull)
        try:
            MW.xlist,MW.ylist,MW.xref,MW.yref,MW.sumlist,MW.peakindex,MW.minL,MW.minR,MW.xcent,MW.height = Int.integrate(MW.zero_norm[massrange],MW.zerofull[massrange],MW.xzero[massrange], MW.deltaY, MW.minimumY,MW.minimumX, MW.maximumX, MW.ntol)
                
            peaktxt = []
            sumlist = MW.sumlist.copy()
            for i in range(0, min(8,len(sumlist))):
                peaktxt.append(np.argmax(sumlist))
                sumlist[np.argmax(sumlist)] = 0
            for i in peaktxt:
                if len(MW.batchstore) == 0:
                    txty = MW.height[i]*max(MW.BL_corrected)
                else:
                    txty = MW.height[i]*max(MW.zerofull)
                MW.ax2.text(MW.xcent[i],txty,str(np.round(MW.xcent[i],3)),rotation='vertical',ha='center',clip_on=True)
    
            max_value = max(MW.zerofull)
            for i in range(0,len(MW.sumlist)):
                MW.sumlist[i] = max_value*MW.sumlist[i]
                MW.height[i] = max_value*MW.height[i]
        except IndexError:
            pass

        if len(MW.batchstore) > 0:
            rtemp = []
            cstemp = []
            gtemp = []
            for i in range(len(MW.batchstore)):
                cstemp.append(MW.batchstore[i][0])
                rtemp.append(MW.batchstore[i][1])
                gtemp.append(MW.batchstore[i][2])
            MW.cs = np.concatenate((cstemp))
            MW.rlist = np.concatenate((rtemp))
            MW.glist = gtemp

        if MW.realsel.isChecked() == False:
            MW.tIFT = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.tIFT = np.real(MW.glist)

        #Noise Calculations
        try:
            MW.selecArea = 0
            MW.GaborArea = (MW.xint[-1]-MW.xshort[0])*(MW.xFT[-1]-MW.xFT[0])
            for i in range(0,len(MW.rlist)):
                MW.selecArea += (MW.rlist[i][1]-MW.rlist[i][0])*(MW.rlist[i][3]-MW.rlist[i][2])*2
            MW.GaborFraction = float(MW.selecArea/MW.GaborArea)
    
            MW.integratedx = 0  
            MW.peakwidths = np.zeros(len(MW.xref))
            for i in range(0, len(MW.xref)):
                MW.integratedx += len(MW.xref[i])
                MW.peakwidths[i] = len(MW.xref[i])
            MW.IntFraction = MW.integratedx/len(MW.xzero[massrange])
            MW.N = np.sqrt(MW.IntFraction)*np.sqrt(MW.GaborFraction)*MW.NFrms*2*np.pi
            MW.Nh = np.sqrt(MW.GaborFraction)*MW.NFrms/np.sqrt(len(MW.xzero[massrange]))*2*np.pi
            if MW.NFrms > 0:
                print('Fraction of Gabor used: ' + str(np.round(MW.GaborFraction,5)))
                print('Fraction of Zero Charge Spectrum integrated: ' + str(np.round(MW.IntFraction,5)))
                print('Noise for integrated spectrum: ' + str(np.round(MW.N,3)))
            
            MW.paramMatch = False
            MW.bounds = []
            for i in range(0,len(MW.xref)):
                MW.bounds.append(MW.xref[i][0])
                MW.bounds.append(MW.xref[i][-1])
        except AttributeError:
            pass
        
        self.close()
        MW.fig.tight_layout()
        MW.canvas.draw()
        
    def Decon_2(self):
        'Performs iFAMS analysis, saves results, and resets spectrogram for continued searching'
        try:
            MW.Harmrec
            try:
                MW.Harmrec = int(self.harm.toPlainText())
            except TypeError:
                MW.message = 'Unexpected value for harmonic number'
                MW.Err = True
                MW.submessage = 'Please enter an integer and try again'
                self.dialog = iFAMS_message(self)
                self.dialog.show()
                return
            if self.cbZeroHarm.isChecked():
                start = 0
            else:
                start = 1
            ztemp = []
            rtemp = []
            Itemp = []
            gtemp = []
            gBL2temp = []
            MW.Xtemp = MW.X*0
            MW.XtempBL2 = MW.XBL2*0
            for i in range(len(self.cbCS)):
                if self.cbCS[i].isChecked():
                    ztemp.append(MW.zlist[i])
                    for j in range(start,MW.Harmrec+1):
                        if j == 0:
                            bottom = 0
                        else:
                            bottom = int(np.round(j*MW.Ilist[i][2]+(j-1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                        top = int(np.round(j*MW.Ilist[i][2]+(j+1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                        rtemp.append([MW.rlist[i][0],MW.rlist[i][1],MW.yshort[bottom],MW.yshort[top],len(gtemp)])
                        Itemp.append([MW.Ilist[i][0],MW.Ilist[i][1],bottom,top,len(gtemp)])
                        if j == 0:
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1]
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-1]
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,-1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,-1]
                        else:
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1]
                            MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1]
                            MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1]
                    MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                    gtemp.append(MW.yrecon)
                    MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                    gBL2temp.append(MW.yreconBL2)
                    MW.Xtemp = MW.X*0
                    MW.XtempBL2 = MW.XBL2*0

        except AttributeError:
            ztemp = []
            rtemp = []
            Itemp = []
            gtemp = []
            gBL2temp = []
            for i in range(len(self.cbCS)):
                if self.cbCS[i].isChecked():
                    ztemp.append(MW.zlist[i])
                    rtemp.append(MW.rlist[i])
                    rtemp[-1][4] = len(rtemp)-1
                    Itemp.append(MW.Ilist[i])
                    Itemp[-1][4] = len(Itemp)-1
                    gtemp.append(MW.glist[i])
                    gBL2temp.append(MW.glistBL2[i])

        MW.zlist = ztemp.copy()
        MW.cs = MW.zlist
        MW.rlist = rtemp.copy()
        MW.Ilist = Itemp.copy()
        MW.glist = gtemp.copy()
        MW.glistBL2 = gBL2temp.copy()
        
        if MW.realsel.isChecked() == False:
            MW.tIFT = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.tIFT = np.real(MW.glist)
        
        MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(MW.tIFT, MW.xint, MW.cs,MW.proton)
        MW.xzeroBL2, MW.cszeroBL2, MW.zerofullBL2 = FT.zerocharge(np.real(MW.glistBL2), MW.xint, MW.cs,MW.proton)
       
        massrange=np.nonzero(MW.zerofullBL2)
        MW.xzero = MW.xzero[massrange]
        MW.zerofull = MW.zerofull[massrange]
        MW.xzeroBL2 = MW.xzeroBL2[massrange]
        MW.zerofullBL2 = MW.zerofullBL2[massrange]
        MW.BSint = []
        for j in range(len(MW.cszero)):
            MW.cszero[j] = MW.cszero[j][massrange]
            MW.cszeroBL2[j] = MW.cszeroBL2[j][massrange]
            check = np.max(MW.cszeroBL2[j])/4
            if check == 0:
                continue
            mid1 = int(len(MW.xzero)/8)+1
            mid2 = int(len(MW.xzero)*7/8)-1
            for k in range(len(MW.cszero[j])):
                if MW.cszeroBL2[j][k] > check:
                    mid1 = k
                    break
            for k in range(mid1+1,len(MW.cszero[j])):
                if MW.cszeroBL2[j][k] <= check or MW.cszeroBL2[j][k] <= MW.cszeroBL2[j][mid1]:
                    mid2 = k
                    break
            int1 = min(sum(MW.cszero[j][0:mid1]),sum(MW.cszero[j][mid2:len(MW.xzero)]))
            int2 = (sum(MW.cszeroBL2[j][0:mid1])+sum(MW.cszeroBL2[j][mid2:len(MW.xzero)]))/2
            MW.BSint.append(int1/int2)
        BSint = np.average(MW.BSint)
        MW.zerofullBL2 = MW.zerofullBL2*BSint
        MW.BL_corrected = (MW.zerofull - MW.zerofullBL2)/max(MW.zerofull)
        MW.zerofull = MW.zerofull - MW.zerofullBL2

        MW.glist = np.sum(MW.tIFT,axis=0)

        #save selection
        storetemp = []
        storetemp.append(MW.cs)
        storetemp.append(MW.rlist)
        storetemp.append(MW.glist)
        storetemp.append(MW.xzero)
        storetemp.append(MW.zerofull)
        storetemp.append(MW.cszero)
        storetemp.append(MW.cszeroBL2)
        storetemp.append(len(MW.tIFT))
        storetemp.append(MW.tIFT)
        MW.batchstore.append(storetemp)
                
        MW.plot_gabor(self)
        
        for j in range(len(MW.batchstore)):
            for i in range(len(MW.batchstore[j][1])):    
                rex = Rectangle((MW.batchstore[j][1][i][0],MW.batchstore[j][1][i][2]), (MW.batchstore[j][1][i][1] - MW.batchstore[j][1][i][0]), (MW.batchstore[j][1][i][3] - MW.batchstore[j][1][i][2]),
                            facecolor='none', edgecolor=color[-2*(j+1)], linewidth=1)
                MW.ax.add_artist(rex) 
                    
        self.close()
        MW.guided_search(self)
        MW.fig.tight_layout()
        MW.canvas.draw()
        
    def CS_Adjust(self):
        ztemp = []
        rtemp = []
        Itemp = []
        gtemp = []
        gBL2temp = []
        slicetext = []
        for i in range(len(self.cbCS)):
            if self.cbCS[i].isChecked():
                ztemp.append(MW.zlist[i])
                rtemp.append(MW.rlist[i])
                rtemp[-1][4] = len(rtemp)-1
                Itemp.append(MW.Ilist[i])
                Itemp[-1][4] = len(Itemp)-1
                gtemp.append(MW.glist[i])
                gBL2temp.append(MW.glistBL2[i])
                slicetext.append(MW.slicetext[i])

        MW.zlist = ztemp.copy()
        MW.cs = MW.zlist
        MW.rlist = rtemp.copy()
        MW.Ilist = Itemp.copy()
        MW.glist = gtemp.copy()
        MW.glistBL2 = gBL2temp.copy()
        MW.slicetext = slicetext.copy()
        if MW.realsel.isChecked() == False:
            MW.glist = np.absolute(MW.glist)
        if MW.realsel.isChecked() == True:
            MW.glist = np.real(MW.glist)
        self.close()
        self.dialog = Z_Adjust(self)
        self.dialog.show()
        
    def Remove(self):
        MW.plot_gabor(self)
        
        for j in range(len(MW.batchstore)):
            for i in range(len(MW.batchstore[j][1])):    
                rex = Rectangle((MW.batchstore[j][1][i][0],MW.batchstore[j][1][i][2]), (MW.batchstore[j][1][i][1] - MW.batchstore[j][1][i][0]), (MW.batchstore[j][1][i][3] - MW.batchstore[j][1][i][2]),
                            facecolor='none', edgecolor=color[-2*(j+1)], linewidth=1)
                MW.ax.add_artist(rex) 
        
        self.close()
        MW.guided_search(self)
        MW.fig.tight_layout()
        MW.canvas.draw()
        
    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y2 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y1 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)

       
class autosearch_check(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(autosearch_check, self).__init__()
        self.setWindowTitle("iFAMS Guided Search")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.setGeometry(1200, 100, 400, 200)
        self.layout = QtWidgets.QGridLayout()
                
        self.message = QtWidgets.QLabel()
        self.message.setText('Select charge state series to keep: ')
        self.message.setFont(QtGui.QFont('Arial',8)) 
        self.layout.addWidget(self.message,0,0)

        self.cbCS = np.zeros(len(MW.batchstore)).tolist()
        count = 0
        for i in range(len(self.cbCS)):
            if i == 20*(count+1):
                count+=1
            seriesM = str(np.round(MW.batchstore[i][3][np.argmax(MW.batchstore[i][4])]/1000,1))
            self.cbCS[i] = QtWidgets.QCheckBox('Series '+str(i+1)+' mass '+seriesM+' kDa')
            self.cbCS[i].setChecked(True)
            self.layout.addWidget(self.cbCS[i],i+1-20*count,count)

        self.button1 = QtWidgets.QPushButton('Finish Search')
        self.button1.clicked.connect(self.Finish_Auto)
        self.layout.addWidget(self.button1,21,0)
        
        self.button2 = QtWidgets.QPushButton('Select Additional Series')
        self.button2.clicked.connect(self.Setup_Search)
        self.layout.addWidget(self.button2,22,0)

        #Adding series to plots
        MW.plot_gabor(self)
        for j in range(len(MW.batchstore)):
            for i in range(len(MW.batchstore[j][1])):    
                rex = Rectangle((MW.batchstore[j][1][i][0],MW.batchstore[j][1][i][2]), (MW.batchstore[j][1][i][1] - MW.batchstore[j][1][i][0]), (MW.batchstore[j][1][i][3] - MW.batchstore[j][1][i][2]),
                            facecolor='none', edgecolor=color[2*j+1], linewidth=1)
                MW.ax.add_artist(rex) 
        MW.ax2.clear()
        min0=[]
        max0=[]
        for i in range(len(MW.batchstore)):
            min0.append(np.min(MW.batchstore[i][3]))
            max0.append(np.max(MW.batchstore[i][3]))
        min0 = np.min(min0)
        max0 = np.max(max0)
        space0 = MW.xzero[2]-MW.xzero[1]
        MW.xzero = np.linspace(min0,max0,int(2*np.round((max0-min0)/space0)),endpoint=True)
        y0s = []
        MW.deconylist = []
        MW.deconxlist = []
        for i in range(len(MW.batchstore)):
            y0s.append(np.interp(MW.xzero,MW.batchstore[i][3],MW.batchstore[i][4]))
            MW.deconxlist.append(MW.batchstore[i][3])
            MW.deconylist.append(MW.batchstore[i][4])
        MW.zerofull = np.sum(y0s,axis=0)
        MW.zero_norm = MW.zerofull/max(MW.zerofull)
        MW.ax2.plot(MW.xzero,MW.zerofull, label='Combined',color='darkgray',linewidth=3)
        for i in range(len(MW.batchstore)):
            label = 'series ' + str(i+1)
            MW.ax2.plot(MW.batchstore[i][3],MW.batchstore[i][4],label=label,color=color[2*i+1])
        MW.ax2.legend(loc='upper right', title='Deconvolutions')
        MW.ax2.set_title('Zero Charge Spectra')
        MW.ax2.set_ylabel('Relative Abundance')
        MW.ax2.set_xlabel('Mass (Da)')

        self.setLayout(self.layout)
        self.show()
        MW.fig.tight_layout()
        MW.canvas.draw()
    
    def Finish_Auto(self):
        bstemp = []
        for i in range(len(self.cbCS)):
            if self.cbCS[i].isChecked():
                bstemp.append(MW.batchstore[i])
        MW.batchstore = bstemp.copy()

        self.close()
        guided_dialog.Finish(self)
        
    def Setup_Search(self):
        bstemp = []
        for i in range(len(self.cbCS)):
            if self.cbCS[i].isChecked():
                bstemp.append(MW.batchstore[i])
        MW.batchstore = bstemp.copy()
        
        MW.plot_gabor(self)
        
        for j in range(len(MW.batchstore)):
            for i in range(len(MW.batchstore[j][1])):    
                rex = Rectangle((MW.batchstore[j][1][i][0],MW.batchstore[j][1][i][2]), (MW.batchstore[j][1][i][1] - MW.batchstore[j][1][i][0]), (MW.batchstore[j][1][i][3] - MW.batchstore[j][1][i][2]),
                            facecolor='none', edgecolor=color[-2*(j+1)], linewidth=1)
                MW.ax.add_artist(rex) 
                    
        self.close()
        MW.guided_search(self)
        MW.fig.tight_layout()
        MW.canvas.draw()
        
    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y2 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y1 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
 
class Z_Adjust(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(Z_Adjust, self).__init__(parent)
        self.setGeometry(1100,100,1000,1200)
        self.setWindowTitle("Charge Adjuster")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))

        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        self.ax.set_title('Mass Spectrum')
        self.ax.set_xlabel('m/z')
        self.ax.set_ylabel('Relative Abundance')
        self.ax2.set_title('Zero Charge Spectrum')
        self.ax2.set_xlabel('Mass (Da)')
        self.ax2.set_ylabel('Relative Abundance')
        
        self.ax.plot(MW.x, MW.y, color='darkgrey')
        for i in range(len(MW.glist)):
            self.ax.plot(MW.xint,np.real(MW.glist[i]),color = color[i+1])
            self.ax.text(MW.xint[np.argmax(MW.glist[i])],max(MW.glist[i]),str(MW.zlist[i]),color=color[i+1],clip_on=True)
        
        self.xzero, self.cszero, self.zerofull = FT.zerocharge(np.real(MW.glist), MW.xint, MW.zlist,MW.proton)
        massrange=np.nonzero(self.zerofull)
        self.xzero = self.xzero[massrange]
        self.zerofull = self.zerofull[massrange]
        for i in range(len(self.cszero)):
            self.cszero[i] = np.array(self.cszero[i])[massrange]
        self.ax2.plot(self.xzero,self.zerofull)
        for i in range(len(self.cszero)):
            self.ax2.plot(self.xzero,self.cszero[i],color=color[i+1])
        self.score = FT.cszero_score(self.cszero)
            
        self.label1 = QtWidgets.QLabel()
        self.label1.setText('Charge Range: ' +str(MW.zlist[-1])+'-'+str(MW.zlist[0]))
        self.label1.setFixedSize(250,30)

        self.labelscore = QtWidgets.QLabel()
        self.labelscore.setText('Score: ' +str(self.score))
        self.labelscore.setFixedSize(250,30)
        self.labelscore.setToolTip('Score is calculated from multiplying charge state spectra, \n'
                                   'integrating result, and taking the natural log of the sum.')

        self.label2 = QtWidgets.QLabel()
        self.label2.setText('Adjust Range:')
        self.label2.setFixedSize(250,30)
        
        self.button1 = QtWidgets.QPushButton('+1')
        self.button1.clicked.connect(self.Higher_Z)
        self.button2 = QtWidgets.QPushButton('-1')
        self.button2.clicked.connect(self.Lower_Z)
        
        self.button3 = QtWidgets.QPushButton('Finish')
        self.button3.clicked.connect(self.close)
        
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.label1)
        self.layout.addWidget(self.labelscore)
        self.layout.addWidget(self.label2)
        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)
        self.layout.addWidget(self.button3)
        self.mpl_toolbar = NavigationToolbar(self.canvas,parent=None)
        self.layout.addWidget(self.mpl_toolbar)
        self.setLayout(self.layout)
        self.canvas.draw()
        
    def Higher_Z(self):
        MW.zlist = np.array(MW.zlist)+1
        self.label1.setText('Charge range: ' +str(MW.zlist[-1])+'-'+str(MW.zlist[0]))
        
        self.ax.clear()
        self.ax2.clear()
        self.ax.set_title('Mass Spectrum')
        self.ax.set_xlabel('m/z')
        self.ax.set_ylabel('Relative Abundance')
        self.ax2.set_title('Zero Charge Spectrum')
        self.ax2.set_xlabel('Mass (Da)')
        self.ax2.set_ylabel('Relative Abundance')
        
        self.ax.plot(MW.x, MW.y, color='darkgrey')
        for i in range(len(MW.glist)):
            self.ax.plot(MW.xint,np.real(MW.glist[i]),color = color[i+1])
            self.ax.text(MW.xint[np.argmax(MW.glist[i])],max(MW.glist[i]),str(MW.zlist[i]),color=color[i+1],clip_on=True)
        
        self.xzero, self.cszero, self.zerofull = FT.zerocharge(np.real(MW.glist), MW.xint, MW.zlist,MW.proton)
        massrange=np.nonzero(self.zerofull)
        self.xzero = self.xzero[massrange]
        self.zerofull = self.zerofull[massrange]
        for i in range(len(self.cszero)):
            self.cszero[i] = np.array(self.cszero[i])[massrange]
        self.ax2.plot(self.xzero,self.zerofull)
        for i in range(len(self.cszero)):
            self.ax2.plot(self.xzero,self.cszero[i],color=color[i+1])
        self.score = FT.cszero_score(self.cszero)
        self.labelscore.setText('Score: ' +str(self.score))
        
        self.canvas.draw()
        
    def Lower_Z(self):
        MW.zlist = np.array(MW.zlist)-1
        self.label1.setText('Charge range: ' +str(MW.zlist[-1])+'-'+str(MW.zlist[0]))

        self.ax.clear()
        self.ax2.clear()
        self.ax.set_title('Mass Spectrum')
        self.ax.set_xlabel('m/z')
        self.ax.set_ylabel('Relative Abundance')
        self.ax2.set_title('Zero Charge Spectrum')
        self.ax2.set_xlabel('Mass (Da)')
        self.ax2.set_ylabel('Relative Abundance')

        self.ax.plot(MW.x, MW.y, color='darkgrey')
        for i in range(len(MW.glist)):
            self.ax.plot(MW.xint,np.real(MW.glist[i]),color = color[i+1])
            self.ax.text(MW.xint[np.argmax(MW.glist[i])],max(MW.glist[i]),str(MW.zlist[i]),color=color[i+1],clip_on=True)
        
        self.xzero, self.cszero, self.zerofull = FT.zerocharge(np.real(MW.glist), MW.xint, MW.zlist,MW.proton)
        massrange=np.nonzero(self.zerofull)
        self.xzero = self.xzero[massrange]
        self.zerofull = self.zerofull[massrange]
        for i in range(len(self.cszero)):
            self.cszero[i] = np.array(self.cszero[i])[massrange]
        self.ax2.plot(self.xzero,self.zerofull)
        for i in range(len(self.cszero)):
            self.ax2.plot(self.xzero,self.cszero[i],color=color[i+1])
        self.score = FT.cszero_score(self.cszero)
        self.labelscore.setText('Score: ' +str(self.score))
        
        self.canvas.draw()
        
    def close(self):
        try:
            MW.inclusion
            self.hide()
            self.dialog = iFAMS_STFT(self)
            self.dialog.show()
        except AttributeError:
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(221)
            MW.ax4 = MW.fig.add_subplot(222, sharex=MW.ax)
            MW.ax3 = MW.fig.add_subplot(223, sharex=MW.ax)
            MW.ax2 = MW.fig.add_subplot(224)
        
            MW.ax3.plot(MW.x, MW.y)
            MW.ax3.set_title('Mass Spectrum')
            MW.ax3.set_xlabel('m/z')
            MW.ax3.set_ylabel('Abundance')
        
            MW.ax4.plot(MW.xshort,MW.Xslice,color='orange')
            MW.ax4.set_title('Gabor Spectrogram Slice \nRough Ion mass: '+str(np.round(MW.ruffmass,1))+' Da')
            MW.ax4.set_xlabel('m/z')
            MW.ax4.set_ylabel('Relative Amplitude')
            MW.ax4.set_ylim(0,1.05*np.max(MW.Xslice))
        
            MW.ax.imshow(abs(MW.X.T), origin='lower', aspect='auto',
                         interpolation='nearest', extent=(MW.xshort[0], MW.xshort[-1], MW.yshort[0], MW.yshort[-1]),
                         cmap='jet', vmax=MW.vmax)
            MW.ax.set_title('Gabor Spectrogram')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Frequency')
            MW.ax.set_xlim(MW.xint[0], MW.xint[-1])
            MW.ax.set_ylim(MW.xFT[0]-MW.maxF/MW.winnum/2,MW.xFT[int(np.ceil(len(MW.xFT)/2))]+MW.maxF/MW.winnum/2)
            MW.ax.fill_between(MW.slicex,MW.slicey1,MW.slicey2,color = 'orange',alpha=0.6)
            MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)
            for i in range(len(MW.rlist)):
                rex = Rectangle((MW.rlist[i][0], MW.rlist[i][2]), (MW.rlist[i][1] - MW.rlist[i][0]), (MW.rlist[i][3] - MW.rlist[i][2]),
                                facecolor='none', edgecolor=color[MW.rlist[i][4]], linewidth=3)
                MW.ax.add_artist(rex)
                MW.ax.text(MW.rlist[i][0]+MW.piW,MW.rlist[i][3]+MW.piH,str(MW.zlist[i]),color=color[i],clip_on=True)
                MW.ax4.text(MW.xshort[MW.slicetext[i]],MW.Xslice[MW.slicetext[i]],str(MW.zlist[i]),color=color[i],clip_on=True,ha='center')
        
            MW.ax2.plot(MW.x, MW.y, color='darkgrey')
            if MW.realsel.isChecked() == False:
                for i in range(len(MW.glist)):
                    MW.ax2.plot(MW.xint, abs(MW.glist[i]), color=color[i])
                    MW.ax2.text(MW.xint[np.argmax(abs(MW.glist[i]))],np.max(abs(MW.glist[i])),str(MW.zlist[i]),color=color[i],clip_on=True,ha='center')
            if MW.realsel.isChecked() == True:
                for i in range(len(MW.glist)):
                    MW.ax2.plot(MW.xint, np.real(MW.glist[i]), color=color[i])
                    MW.ax2.text(MW.xint[np.argmax(np.real(MW.glist[i]))],np.max(np.real(MW.glist[i])),str(MW.zlist[i]),color=color[i],clip_on=True,ha='center')
            MW.ax2.set_title('Reconstructed Spectrum')
            MW.ax2.set_xlabel('m/z')
            MW.ax2.set_ylabel('Abundance')
        
            MW.fig.tight_layout()
            MW.canvas.draw()
            self.hide()
            self.dialog = guided_decharger(self)
            self.dialog.show()
       
    def onselect(self, eclick, erelease):
        MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
        MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y2 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y1 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
        
class slice_viewer(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(slice_viewer, self).__init__(parent)
        self.setGeometry(100, 100, 1100, 900)
        self.setWindowTitle("Gabor Slicer")
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))
        
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        ###Layout attempts to match iFAMS STFT layout
        self.ax = self.fig.add_subplot(223) #Raw MS plot bottom left
        self.axGT = self.fig.add_subplot(221, sharex=self.ax) #Gabor spectrogram with slice top left
        self.axS = self.fig.add_subplot(222, sharex=self.ax) #Slice plot top right
        self.axDCN = self.fig.add_subplot(224) #"Working plot" reserved for quick deconvolution bottom right
        self.ax.plot(MW.xint,MW.yint)
        self.ax.set_title('Mass Spectrum',y=0.9)
        self.ax.set_ylabel('Abundance')
        self.ax.set_xlabel('m/z')
        self.ax.set_xlim(MW.xshort[0]-(MW.xshort[2]-MW.xshort[1])/2,MW.xint[-1])
        
        self.axDCN.set_title('Quick Deconvolution',y=0.9)
        
        initx = MW.xshort[int(len(MW.xshort)/2)]
        inity = MW.yshort[int(len(MW.yshort)/8)]
        slicey = 1/MW.xshort * inity*initx
        MW.Xslice = MW.xshort*0
        for i in range(len(MW.Xslice)):
            Fi = int(np.round((slicey[i])/(MW.yshort[2]-MW.yshort[1])))
            try:
                MW.Xslice[i] += np.absolute(MW.X[i][Fi])
            except IndexError:
                pass

        self.dataline = self.axS.plot(MW.xshort,MW.Xslice,color='orange')[0]
        self.axS.set_title('Gabor Spectrogram Slice',y=0.9)
        self.axS.set_xlabel('m/z')
        self.axS.set_ylabel('Relative Magnitude')
        self.axS.set_ylim(0,1.05*np.max(np.absolute(MW.X)))

        self.axGT.imshow(abs(MW.X.T), origin='lower', aspect='auto',
                     interpolation='nearest', extent=(MW.xshort[0], MW.xshort[-1], MW.yshort[0], MW.yshort[-1]),
                     cmap='jet', vmax=MW.vmax)
        self.line = self.axGT.plot(MW.xshort,slicey,color='orange')[0]
        self.axGT.set_title('Gabor Spectrogram',y=0.9,color='w')
        self.axGT.set_ylabel('Frequency')
        self.axGT.set_xlabel('m/z')
        self.axGT.set_xlim(MW.xint[0], MW.xint[-1])
        self.axGT.set_ylim(MW.xFT[0]-MW.maxF/MW.winnum/2,MW.xFT[int(np.ceil(len(MW.xFT)/2))]+MW.maxF/MW.winnum/2)

        for j in range(len(MW.batchstore)):
            for i in range(len(MW.batchstore[j][1])):    
                rex = Rectangle((MW.batchstore[j][1][i][0],MW.batchstore[j][1][i][2]), (MW.batchstore[j][1][i][1] - MW.batchstore[j][1][i][0]), (MW.batchstore[j][1][i][3] - MW.batchstore[j][1][i][2]),
                            facecolor='none', edgecolor=color[-2*(j+1)], linewidth=1)
                self.axGT.add_artist(rex) 

        self.mlabel = QtWidgets.QLabel()
        self.mlabel.setText('Estimated Ion Mass: -')
        self.mlabel.setFixedSize(600,25)
        self.mlabel.setToolTip('Move orange slice to investigate different ion signals. \n'
                               'Release slice to update the ion mass estimates.'
                               'If no mass estimates are showing, try adjusting data \n'
                               'domain and try again.')
        self.addlabel = QtWidgets.QLabel()
        self.addlabel.setText('Estimated (Adduct Mass/Harmonic Number): -')
        self.addlabel.setFixedSize(600,25)

        self.button1 = QtWidgets.QPushButton('Quick Deconvolve')
        self.button1.clicked.connect(self.DCN)
        self.button1.setToolTip("Note: Quick Deconvolve Tool is still in development. \n"
                                "Best practices are to use this tool to verify mass of isotopic signal identified \n"
                                "with the slicer. Results are typically more accurate when using 'High Mass \n"
                                "Resolution' STFT parameters because this tool does not optimize frequencies \n"
                                "(in favor of speed). For accurate deconvolution, use full iFAMS analysis.\n"
                                "Currently, Quick Deconvolve results cannot be saved.")

        self.pressed=False
        self.connect()
        self.layout = QtWidgets.QGridLayout()
        self.mpl_toolbar = NavigationToolbar(self.canvas,parent=None)
        self.layout.addWidget(self.canvas,0,0,12,2)
        self.layout.addWidget(self.mlabel,12,0,1,2)
        self.layout.addWidget(self.addlabel,13,0,1,2)
        self.layout.addWidget(self.button1,14,0,1,2)
        self.layout.addWidget(self.mpl_toolbar)
        self.setLayout(self.layout)
        self.canvas.draw()

        
    def updateLine(self,xin,yin):
        try:
            slicey = 1/MW.xshort * yin*xin
            MW.Xslice = MW.xshort*0
            for i in range(len(MW.Xslice)):
                Fi = int(np.round((slicey[i])/(MW.yshort[2]-MW.yshort[1])))
                try:
                    MW.Xslice[i] += np.absolute(MW.X[i][Fi])
                except IndexError:
                    pass
            self.line.set_data(MW.xshort,slicey)
            self.dataline.set_data(MW.xshort,MW.Xslice)
        except IndexError:
            self.mlabel.setText('Estimated Ion Mass: Undetermined (Slice out of Range)')
            self.addlabel.setText('Estimated (Adduct Mass/Harmonic Number): Undetermined')
            self.canvas.draw()
            self.pressed=False
            return

    def connect(self):
        self.a1=self.fig.canvas.mpl_connect('button_press_event',self.fpressed)
        self.a2=self.fig.canvas.mpl_connect('motion_notify_event',self.fmotion)
        self.a3=self.fig.canvas.mpl_connect('button_release_event',self.fdepressed)
        
    def fpressed(self,event):
        if event.inaxes!=self.axGT: return
        self.line.set_animated(1)
        self.dataline.set_animated(1)
        self.canvas.draw()
        self.bg=self.canvas.copy_from_bbox(self.axGT.bbox)
        self.bg2=self.canvas.copy_from_bbox(self.axS.bbox)
        self.axGT.draw_artist(self.line)
        self.axS.draw_artist(self.dataline)
        self.canvas.blit(self.axGT.bbox)
        self.canvas.blit(self.axS.bbox)
        self.pressed=True
        
    def fmotion(self,event):
        if not self.pressed: return
        if event.inaxes!=self.axGT: return
        self.updateLine(event.xdata,event.ydata)
        self.canvas.restore_region(self.bg)
        self.canvas.restore_region(self.bg2)
        self.axGT.draw_artist(self.line)
        self.axS.draw_artist(self.dataline)
        self.canvas.blit(self.axGT.bbox)
        self.canvas.blit(self.axS.bbox)
    
    def fdepressed(self,event):
        if not self.pressed: return
        self.pressed=False
        self.line.set_animated(0)
        self.dataline.set_animated(0)
        try:
            mass = []
            mscore = []
            slicetemp = MW.Xslice.copy()
            for i in range(4):
                mz1 = MW.xshort[np.argmax(slicetemp)]
                index1 = np.argmax(slicetemp)
                slicetemp[np.argmax(slicetemp)-1:np.argmax(slicetemp)+2] = slicetemp[np.argmax(slicetemp)-1:np.argmax(slicetemp)+2]*0
                mz2 = MW.xshort[np.argmax(slicetemp)]
                index2 = np.argmax(slicetemp)
                if mz1 > mz2:
                    pass
                else:
                    temp = mz1
                    mz1 = mz2
                    mz2 = temp
                z1 = int(np.round(1/(mz1/mz2-1)))
                z2 = z1 + 1
                mass1 = mz1*z1-z1*MW.proton
                mass2 = mz2*z2-z2*MW.proton
                mass.append(np.average((mass1,mass2)))
                sumtemp = 0
                count = 0
                for j in range(z1-2,z2+3):
                    if j <= 0:
                        continue
                    mztemp = (mass[-1]+j*MW.proton)/j
                    mzi = int(np.round((mztemp-MW.xshort[0])/(MW.xshort[2]-MW.xshort[1])))
                    if mzi < 1 or mzi >= len(MW.xshort):
                        continue
                    else:
                        sumtemp += MW.Xslice[mzi]
                        count += 1
                mscore.append(sumtemp/count)
                if i == 0:
                    slicetemp = MW.Xslice.copy()
                    slicetemp[index2-1:index2+2] = slicetemp[index2-1:index2+2]*0
                if i == 1:
                    slicetemp = MW.Xslice.copy()
                    slicetemp[index1-1:index1+2] = slicetemp[index1-1:index1+2]*0
            mass = int(np.round(mass[np.argmax(mscore)]))
            addmass = np.round(mass / (event.xdata*event.ydata),2)
            if addmass >= 1:
                addmass = np.round(mass / (event.xdata*event.ydata),1)
            self.mlabel.setText('Estimated Ion Mass: '+str(mass)+' Da')
            self.addlabel.setText('Estimated (Adduct Mass/Harmonic Number): '+str(addmass)+' Da')
        
            self.axS.set_ylim(0,1.05*np.max(MW.Xslice))
            self.canvas.draw()
            self.mass = mass
        
        except TypeError:
            self.mlabel.setText('Estimated Ion Mass: Undetermined')
            self.addlabel.setText('Estimated (Adduct Mass/Harmonic Number): Undetermined')
        except ValueError:
            self.mlabel.setText('Estimated Ion Mass: Undetermined')
            self.addlabel.setText('Estimated (Adduct Mass/Harmonic Number): Undetermined')
        except IndexError:
            self.mlabel.setText('Estimated Ion Mass: Undetermined')
            self.addlabel.setText('Estimated (Adduct Mass/Harmonic Number): Undetermined')
        except OverflowError:
            self.mlabel.setText('Estimated Ion Mass: Undetermined')
            self.addlabel.setText('Estimated (Adduct Mass/Harmonic Number): Undetermined')
        
    def DCN(self):
        ###Quick Decon
        print('Selecting signal...')
        self.axDCN.clear()
        mass = self.mass
        self.loZ = np.max((int(np.ceil(mass/MW.xint[-1]))+1,2))
        self.hiZ = np.min((int(np.floor(mass/MW.xint[0]))-1,70,int(mass/650)))
        zcheck = np.flip(np.array(range(self.loZ,self.hiZ)))
        mzcheck = (mass+zcheck*MW.proton)/zcheck
        xspace = (MW.xshort[2]-MW.xshort[1])
        mzchecki = np.round((mzcheck-MW.xshort[0])/xspace)
        mzchecki = mzchecki.astype(int)
        mzytemp = MW.Xslice[mzchecki]
        
        value = []
        indices = []
        end = np.min((5,int(np.round(len(mzytemp)/2))))
        try:
            for i in range(0,len(mzytemp)-end,1):
                testlist = mzytemp[i:i+end+1]
                value.append(np.average(testlist)/(np.std(testlist))**(1/2))
                indices.append(np.arange(i,i+end+1))
            indices = indices[np.argmax(value)]
        except ValueError:
            indices = np.arange(0,len(zcheck))
        self.charges = zcheck[indices]
        self.Z = self.charges[np.argmax(mzytemp[indices])]
        self.mzis = mzchecki[indices]
        
        if 10 <= self.Z < 20:
            self.charges = np.arange(self.Z-2,self.Z+2+1)
        elif 20 <= self.Z < 35:
            self.charges = np.arange(self.Z-4,self.Z+4+1)
        elif 35 <= self.Z:
            self.charges = np.arange(self.Z-6,self.Z+6+1)
        zlist = []
        for i in range(len(self.charges)):
            if self.loZ <= self.charges[i] <= self.hiZ:
                zlist.append(self.charges[i])
        zlist = np.array(zlist)
        piW = (MW.xshort[2]-MW.xshort[1])
        piH = (MW.yshort[2]-MW.yshort[1])
        aX = np.absolute(MW.X)
        Xtemp = MW.X.copy()*0
        maxF = np.floor((MW.xFT[-1])/2)
        
        mzcheck = (mass+zlist*MW.proton)/zlist
        picheck = np.round((mzcheck-MW.xshort[0])/piW)
        picheck = picheck.astype(int)
        Rwidths = []
        Lwidths = []
        relmin = np.min(aX[:,0])
        thresh = (np.max(aX[picheck,0])-relmin)*0.3
        for i in range(len(picheck)):
            Ltemp = 2
            Rtemp = 2
            for j in range(1,6):
                if aX[picheck[i]-j][0]-relmin < thresh:
                    Ltemp = j
                    break
            for j in range(1,6):
                if aX[picheck[i]+j][0]-relmin < thresh:
                    Rtemp = j
                    break
            Lwidths.append(Ltemp)
            Rwidths.append(Rtemp)
        Lw = np.average(Lwidths)
        Rw = np.average(Rwidths)+1
        
        Xxi = []
        Xyi = []
        glist = []
        yBL = np.zeros(len(MW.xint))+ min(MW.yint)+ 1
        glistBL = []
        xintBL,yintBL,ypaddBL,po2BL,xpaddBL = load.datapro(MW.xint,yBL)
        XBL = STFT.stft(ypaddBL, MW.winnum, MW.window,10)
        XtempBL = XBL*0
        for i in range(len(picheck)):
            Xxi.append(np.arange(picheck[i]-Lw,picheck[i]+Rw).astype(int))
            tempy = np.arange(0,int(np.ceil(1.8/piH)))
            for j in range(1,5):
                if zlist[i]*j+1 >= maxF:
                    break
                tempy = np.concatenate((tempy,np.arange(int(np.round((zlist[i]*j-0.75)/piH)),int(np.ceil((zlist[i]*j+0.75)/piH)))))
            Xyi.append(tempy.astype(int))
        for i in range(len(picheck)):
            for j in range(len(Xxi[i])):
                Xtemp[Xxi[i][j],Xyi[i]] = MW.X[Xxi[i][j],Xyi[i]]
                XtempBL[Xxi[i][j],Xyi[i]] = XBL[Xxi[i][j],Xyi[i]]
                Xtemp[Xxi[i][j],-Xyi[i]] = MW.X[Xxi[i][j],-Xyi[i]]
                XtempBL[Xxi[i][j],-Xyi[i]] = XBL[Xxi[i][j],-Xyi[i]]
            yrecon = STFT.re_stft(Xtemp,MW.winnum,MW.ypadd,MW.yint,10)*MW.correction_factor
            glist.append(yrecon)
            Xtemp = MW.X.copy()*0
            yreconBL = STFT.re_stft(XtempBL,MW.winnum,ypaddBL,yintBL,10)*MW.correction_factor
            glistBL.append(yreconBL)
            XtempBL = XBL*0
        
        xzero, cszero, zerofull = FT.zerocharge(np.real(glist), MW.xint, zlist,MW.proton)
        xzeroBL, cszeroBL, zerofullBL = FT.zerocharge(np.real(glistBL), MW.xint, zlist,MW.proton)
        massrange=np.nonzero(zerofull)
        xzero = xzero[massrange]
        zerofull = zerofull[massrange]
        xzeroBL = xzeroBL[massrange]
        zerofullBL = zerofullBL[massrange]
        self.BSint = []
        for j in range(len(cszero)):
            cszero[j] = cszero[j][massrange]
            cszeroBL[j] = cszeroBL[j][massrange]
            check = np.max(cszeroBL[j])/4
            if check == 0:
                continue
            mid1 = int(len(xzero)/8)+1
            mid2 = int(len(xzero)*7/8)-1
            for k in range(0,len(cszero[j]),5):
                if cszeroBL[j][k] > check:
                    mid1 = k
                    break
            for k in range(mid1+1,len(cszero[j]),5):
                if cszeroBL[j][k] <= check or cszeroBL[j][k] <= cszeroBL[j][mid1]:
                    mid2 = k
                    break
            int1 = min(sum(cszero[j][0:mid1]),sum(cszero[j][mid2:len(xzero)]))
            int2 = (sum(cszeroBL[j][0:mid1])+sum(cszeroBL[j][mid2:len(xzero)]))/2
            cszero[j] = cszero[j]-cszeroBL[j]*int1/int2
            self.BSint.append(int1/int2)
        BSint = np.average(self.BSint)
        zerofullBL = zerofullBL*BSint
        zerofull = zerofull - zerofullBL
        zero_norm = zerofull/max(zerofull)

        DCNcolor = "mediumseagreen"
        self.axDCN.set_title("Quick Deconvolution",y=0.9)
        self.axDCN.set_ylabel("Normalized Abundance (%)")
        self.axDCN.set_xlabel("Mass (Da)")
        for i in range(len(cszero)):
            self.axDCN.plot(xzero,cszero[i]/max(zerofull)*100)
        self.axDCN.plot(xzero,zero_norm*100,c=DCNcolor)
        self.axDCN.set_ylim(-5,110)
        print(str('Deconvolved spectrum normalized to '+str(np.round(max(zerofull)))))       
        
        self.canvas.draw()        
        
class FileDialog(QtWidgets.QFileDialog):
    def __init__(self, *args):
        QtWidgets.QFileDialog.__init__(self, *args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setFileMode(self.ExistingFiles)
        btns = self.findChildren(QtWidgets.QPushButton)
        self.openBtn = [x for x in btns if 'open' in str(x.text()).lower()][0]
        self.openBtn.clicked.disconnect()
        self.openBtn.clicked.connect(self.openClicked)
        self.tree = self.findChild(QtWidgets.QTreeView)

    def openClicked(self):
        inds = self.tree.selectionModel().selectedIndexes()
        files = []
        for i in inds:
            if i.column() == 0:
                files.append(os.path.join(str(self.directory().absolutePath()),str(i.data())))
        MW.Batchname = files
        self.hide()
        self.close()
        self.dialog = batch_menu(self)
        self.dialog.show()
    
class MW(QtWidgets.QMainWindow):
    """
    Main window
    """
    def __init__(self, parent=None):
        """
        This builds the main window
        it has connects to package that
        runs the script on the backend
        :param parent:
        """
        QtWidgets.QMainWindow.__init__(self, parent)
        self.create_main_frame()
        self.setWindowTitle("iFAMS v6.3 Quant")
        self.setWindowIcon(QtGui.QIcon('icon.png'))
        ##### file menu ####
        self.file_menu = QtWidgets.QMenu('File', self)
        self.file_menu.addAction('Load Ordinary FT',self.FT_load)
        self.file_menu.addAction('Load STFT',self.STFT_load)

        self.file_menu.addAction('Plot STFT From Clipboard',self.STFT_clipboard)
        self.file_menu.addAction('Load Deconvolved or IFFT Spectra',self.deconv_load)
        self.file_menu.addSeparator()
        self.file_menu.addAction('Load Batch Parameters',self.batch_load)
        self.file_menu.addAction('Batch Data Files',self.batch)
        self.file_menu.addAction('Batch Deconvolved Spectra', self.deconbatch)
        self.file_menu.addSeparator()
        self.file_menu.addAction('Export Data',self.export_data)
        self.file_menu.addAction('Export Batch Parameters',self.batch_export_opt)
        self.menuBar().addMenu(self.file_menu)
        ##### file menu #####

        ##### Settings/Data tools #####
        self.data_menu = QtWidgets.QMenu('Settings',self)
        MW.realsel = QtWidgets.QAction('Plot Real Data', self, checkable=True)
        MW.realsel.setChecked(True)
        MW.realsel.setToolTip('If unchecked, the absolute value of the data will be plotted.')
        MW.ionmode = QtWidgets.QAction('Positive Ion Mode', self, checkable=True)
        MW.ionmode.setChecked(True)
        MW.ionmode.setToolTip('If checked, the deconvolution will include a mass correction necessary for \n'
                              'postive ion mode data (subtration of mass of protons)\n'
                              'If unchecked, the deconvolution will include a mass correction necessary for \n'
                              'negative ion mode data (addition of mass of protons)')
        MW.ionmode.triggered.connect(self.ionmode_set)
        MW.MS_interp = QtWidgets.QAction('Show Interpolated MS', self, checkable=True)
        MW.MS_interp.setChecked(False)
        MW.MS_interp.triggered.connect(self.MS_interp_replot)
        MW.normmode = QtWidgets.QAction('Normalize Abundances', self, checkable=True)
        MW.normmode.setChecked(False)
        MW.normmode.setToolTip('If checked, abundances of each spectrum will be displayed as \n'
                               'a percent of its maximum')
        MW.normmode.triggered.connect(self.setting_replot)
        MW.smoothmode = QtWidgets.QAction('Smooth Isotope Resolution', self, checkable=True)
        MW.smoothmode.setChecked(False)
        MW.smoothmode.setToolTip('If checked, all deconvolved spectra will be displayed with isotope \n'
                                 'resolution smoothed for visual purposes only')
        MW.smoothmode.triggered.connect(self.setting_replot)
        self.data_menu.addAction(MW.ionmode)
        self.data_menu.addAction('Adjust Data Domain',self.domainTool)
        self.data_menu.addSeparator()
        self.data_menu.addAction(MW.MS_interp)
        self.data_menu.addAction(MW.normmode)
        self.data_menu.addAction(MW.smoothmode)
        self.menuBar().addMenu(self.data_menu)
        ##### Settings/Data tools #####
        
        ##### Fourier analysis tools #####
        self.fourier_menu = QtWidgets.QMenu('Fourier Analysis Tools',self)
        self.fourier_menu.addAction('Plot Absolute Data',self.plotAbs)
        self.fourier_menu.addAction('Plot Real Data',self.plotReal)
        self.fourier_menu.addSeparator()
        self.fourier_menu.addAction('Open iFAMS Maxima Peak Finder',self.maxima_finder)
        self.fourier_menu.addAction('Manually Enter Charge States and Submass',self.man_input)
        self.fourier_menu.addSeparator()
        self.fourier_menu.addAction('Run iFAMS Analysis',self.harmonic)
        self.fourier_menu.addSeparator()
        self.fourier_menu.addAction('Run Mass Defect Analysis',self.defect)
        self.menuBar().addMenu(self.fourier_menu)
        ##### Fourier analysis tools #####

        ##### Gabor analysis tools #####
        self.Gabor_menu = QtWidgets.QMenu('STFT Analysis Tools',self)
        self.Gabor_menu.addAction('Change STFT Parameters',self.STFT_parm_window)
        self.Gabor_menu.addAction('Open Gabor Slicer',self.STFT_slicer)
        self.Gabor_menu.addSeparator()
        self.Gabor_menu.addAction('Run Guided Search', self.guided_search)
        self.Gabor_menu.addAction('Add Gabor Selection', self.add_gabor,
                                  QtCore.Qt.CTRL + QtCore.Qt.Key_E)
        self.Gabor_menu.addAction('Save Gabor Selections', self.save_gabor_sel,
                                  QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.Gabor_menu.addAction('Delete Previous Gabor Selection', self.del_gabor,
                                  QtCore.Qt.CTRL + QtCore.Qt.Key_Delete)
        self.Gabor_menu.addAction('Open Harmonic Finder', self.harmonic_finder)
        self.Gabor_menu.addSeparator()
        self.Gabor_menu.addAction('Run iFAMS Analysis',self.STFT_iFAMS)
        self.Gabor_menu.addAction('Adjust Charge State Assignments',self.cs_adj_menu)
        self.Gabor_menu.addAction('Open Quantitative Peak Integrator',self.inte)
        self.Gabor_menu.addAction('Run Mass Defect Analysis',self.g_defect)
        self.Gabor_menu.addSeparator()
        self.Gabor_menu.addAction('Show Only Spectrogram',self.show_gabor_only)
        self.Gabor_menu.addAction('Open Noise Calculator',self.noisemenu)
        self.menuBar().addMenu(self.Gabor_menu)
        ##### Gabor analysis tools #####
        
        ##### calibration menu #####
        self.calibration_menu = QtWidgets.QMenu('Calibration Curve Tools')
        self.calibration_menu.addAction('Create New Curve',self.cali)
        self.calibration_menu.addAction('Load Existing Curve',self.load_cali)
        self.calibration_menu.addSeparator()
        self.calibration_menu.addAction('Save Calibration Curve',self.save_cali)
        self.calibration_menu.addAction('Load Files for Concentration Calculation',self.cali2)
        self.menuBar().addMenu(self.calibration_menu)
        ##### calibration menu #####
        
        ##### Isotopic Distribution #####
        self.isotope_menu = QtWidgets.QMenu('Isotopic Distribution',self)
        self.isotope_menu.addAction('Calculate Distribution',self.iso_input)
        self.isotope_menu.addSeparator()
        self.isotope_menu.addAction('Overlay Data Reconstruction',self.data_overlay)
        self.isotope_menu.addAction('Export Calculated Distribution',self.iso_export)
        self.menuBar().addMenu(self.isotope_menu)
        ##### Isotopic Distribution #####

    def create_main_frame(self):
        """
        makes the central plot
        :return: None
        """
        self.main_frame = QtWidgets.QWidget()

        MW.fig = Figure()

        MW.canvas = FigureCanvas(MW.fig)
        MW.canvas.setParent(self.main_frame)
        plt.rcParams.update({'font.size':14})
        MW.file_display = QtWidgets.QLabel()
        MW.file_display.setText('File: - ')   
        MW.file_display.setFixedSize(1600,30)
        MW.ionmode_display = QtWidgets.QLabel()
        MW.ionmode_display.setText('Positive Ion Mode')
        MW.ionmode_display.setFixedSize(1200,30)
        self.mpl_toolbar = NavigationToolbar(MW.canvas, self.main_frame)

        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(MW.canvas)
        vbox.addWidget(MW.file_display)
        vbox.addWidget(MW.ionmode_display)
        vbox.addWidget(self.mpl_toolbar)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
   
        MW.proton = 1.00727647
        #Processing booleans
        MW.ft_pro = False
        MW.stft_pro = False
        MW.stft_selec_man = True
        MW.stft_box_thresh = "NA"
        MW.stft_range_thresh = "NA"
        MW.decon_pro = False
        #Mass spectrum baseline and domain adjustment booleans
        MW.setdatarho = False
        MW.domaintrunc = False
        MW.Nsubtraction = False
        MW.Nlinear_sub = False
        MW.Nseg_sub = False
        #Integration booleans
        MW.NZbased = False
        MW.basecor = False
        MW.based = False
        MW.smoothed = False
        MW.ReInt = False 
        MW.range2 = False
        #Saving booleans
        MW.saveFT = False
        MW.saveIFFT = False
        MW.saveZspec = False
        MW.saveSeries = False
        MW.batchcheck = False
        MW.batchstore = []
        MW.saveMSRT = False
        MW.full_save = False
        MW.saveMS = False
        try:
            del(MW.batchfolder)
        except AttributeError:
            pass
        try:
            del(MW.zerofullstore)
        except AttributeError:
            pass
    
    def deconv_load(self):
        try:
            try:
                del(MW.ax2)
            except AttributeError:
                pass
            #Integration booleans
            MW.NZbased = False
            MW.basecor = False
            MW.based = False
            MW.smoothed = False
            MW.ReInt = False 
            MW.range2 = False
            try:
                del(MW.rlist)
            except AttributeError:
                pass
            try:
                del(MW.zerofullstore)
            except AttributeError:
                pass
            try:
                del(iso.xrot1)
            except AttributeError:
                pass
            try:
                del(MW.x_overlay)
            except AttributeError:
                pass
            try:
                del(MW.xrefcs)
            except AttributeError:
                pass
            try:
                del(MW.inclusion)
            except AttributeError:
                pass

            MW.batchstore = []
            MW.xzeroOL = []   #Lists for spectra to overlay
            MW.zerofullOL = []
            MW.namebaseOL = []
            MW.BoundsOL = []
            MW.csOL = []
            MW.cszeroOL = []
            MW.deconname, extra = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open File(s)')
            MW.deconname[0]
            MW.xzero, MW.zerofull,MW.namebase,MW.Bounds,MW.based,MW.cs,MW.cszero,MW.decon_xtype = load.load_decon(MW.deconname[0])
            for i in range(1,len(MW.deconname)):
                MW.xzeroOL.append(0)
                MW.zerofullOL.append(0)
                MW.namebaseOL.append(0)
                MW.BoundsOL.append(0)
                MW.csOL.append(0)
                MW.cszeroOL.append(0)
                MW.xzeroOL[i-1], MW.zerofullOL[i-1],MW.namebaseOL[i-1],MW.BoundsOL[i-1],based,MW.csOL[i-1],MW.cszeroOL[i-1],domain = load.load_decon(MW.deconname[i])
            
            if MW.Bounds == 0 or len(MW.deconname) > 1 or MW.smoothmode == True:
                deconInt = False
            else:
                deconInt = True
                
            if MW.decon_xtype == 'm/z':
                MW.fig.clf()
                MW.ax2 = MW.fig.add_subplot(111)
                MW.ax2.plot(MW.xzero, MW.zerofull,label='Mass Spectrum',color='darkgray')
                for i in range(len(MW.cs)):
                    if len(MW.csOL) > 0:
                        label = str(MW.cs[i])+' '+os.path.basename(str(MW.namebase))
                    else:
                        label = 'Charge State '+str(MW.cs[i])
                    MW.ax2.plot(MW.xzero,MW.cszero[i],color = color[i],alpha=0.5,label=label)
                    MW.ax2.text(MW.xzero[np.argmax(MW.cszero[i])],np.max(MW.cszero[i]),str(MW.cs[i]),color=color[i],clip_on=True,ha='center')
                for i in range(len(MW.csOL)):
                    for j in range(len(MW.csOL[i])):
                        if i == 0:
                            colori = len(MW.cs)+j
                        else:
                            colori = len(MW.cs)+len(MW.csOL[i-1])+j
                        MW.ax2.plot(MW.xzeroOL[i],MW.cszeroOL[i][j],color = color[colori],alpha=0.5,label=str(MW.csOL[i][j])+' '+os.path.basename(str(MW.namebaseOL[i])))
                        MW.ax2.text(MW.xzeroOL[i][np.argmax(MW.cszeroOL[i][j])],np.max(MW.cszeroOL[i][j]),str(MW.csOL[i][j]),color=color[colori],clip_on=True,ha='center')
                MW.ax2.legend(loc='upper right', title='File Names')
                MW.ax2.set_xlabel('m/z')
                MW.ax2.set_title('Mass Spectrum with Charge-state IFFT Spectra Overlay')
                MW.ax2.set_ylabel('Relative Abundance')
            
            elif deconInt == True:
                MW.xref, MW.yref, MW.sumlist, MW.height, MW.xcent, MW.peakindex = Int.boundsmatch(MW.xzero,MW.zerofull,MW.Bounds)
                if MW.based == True:
                    MW.yref2, MW.sumlist, MW.height = Int.baseline_sub2(MW.xzero, MW.zerofull, MW.xref, MW.yref, MW.peakindex)            
                MW.xpoints = len(np.concatenate(MW.xref))
                MW.minimumX = MW.Bounds[0]-1
                MW.maximumX = MW.Bounds[-1]+1
                integrate.int_show(self)
                normfact = np.max(MW.zerofull)
                for i in range(len(MW.cszero)):
                    MW.ax2.plot(MW.xzero,MW.cszero[i]/normfact,color = color[i],alpha=0.5)
                    MW.ax2.text(MW.xzero[np.argmax(MW.cszero[i])],np.max(MW.cszero[i])/normfact,str(MW.cs[i]),color=color[i],clip_on=True,ha='center')
            else:
                x_plot = []
                y_plot = []
                plot_labels = []
                if MW.smoothmode.isChecked():
                    if MW.normmode.isChecked():
                        ytemp = Int.smooth(MW.xzero,MW.zerofull,1.8)
                        x_plot.append(MW.xzero)
                        y_plot.append(ytemp/max(ytemp)*100)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            ytemp = Int.smooth(MW.xzeroOL[i],MW.zerofullOL[i],1.8)
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(ytemp/max(ytemp)*100)
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Normalized Abundance (%)'
                    else:
                        ytemp = Int.smooth(MW.xzero,MW.zerofull,1.8)
                        x_plot.append(MW.xzero)
                        y_plot.append(ytemp)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            ytemp = Int.smooth(MW.xzeroOL[i],MW.zerofullOL[i],1.8)
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(ytemp)
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Relative Abundance'
                else:
                    if MW.normmode.isChecked():
                        x_plot.append(MW.xzero)
                        y_plot.append(MW.zerofull/max(MW.zerofull)*100)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(MW.zerofullOL[i]/max(MW.zerofullOL[i])*100)
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Normalized Abundance (%)'
                    else:
                        x_plot.append(MW.xzero)
                        y_plot.append(MW.zerofull)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(MW.zerofullOL[i])
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Relative Abundance'
                MW.fig.clf()
                MW.ax2 = MW.fig.add_subplot(111)
                for i in range(len(x_plot)):
                    MW.ax2.plot(x_plot[i], y_plot[i],label=plot_labels[i],alpha=0.8)
                MW.ax2.legend(loc='upper right', title='File Names')
                MW.ax2.set_xlabel('Mass (Da)')
                MW.ax2.set_ylabel(Yaxis_Label)

            print('File: ' + MW.namebase)
            MW.file_display.setText('File: '+MW.namebase) 
            MW.fig.tight_layout()
            MW.canvas.draw()
    
            #Processing booleans
            MW.ft_pro = False
            MW.stft_pro = False
            MW.stft_selec_man = "NA"
            MW.stft_box_thresh = "NA"
            MW.stft_range_thresh = "NA"
            MW.decon_pro = True
            #Saving booleans
            MW.saveFT = False
            MW.saveIFFT = False
            MW.saveZspec = False
            MW.saveSeries = False
            MW.batchcheck = False
            MW.full_save = False
            MW.saveMS = False
            #Noise Determination
            MW.NFrms = 0
            MW.Nrms = 0
            MW.Nh = 0
            MW.N = 0
            MW.GaborFraction = 0
                
        except TypeError:
            MW.ref_num = 0
            print('Either no data loaded or wrong data type. Try .csv or .txt for wrong data type')
        except IndexError:
            print('No file(s) selected')
        except ValueError:
            print('Unable to load file. Check file type and try again')
            MW.message = 'Unable to load file. Check file type and try again'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
   
    def setting_replot(self):
        if MW.ft_pro == False and MW.stft_pro == False and MW.decon_pro == True:
            if MW.decon_xtype == 'm/z':
                pass
            else:
                x_plot = []
                y_plot = []
                plot_labels = []
                if MW.smoothmode.isChecked():
                    if MW.normmode.isChecked():
                        ytemp = Int.smooth(MW.xzero,MW.zerofull,1.8)
                        x_plot.append(MW.xzero)
                        y_plot.append(ytemp/max(ytemp)*100)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            ytemp = Int.smooth(MW.xzeroOL[i],MW.zerofullOL[i],1.8)
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(ytemp/max(ytemp)*100)
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Normalized Abundance (%)'
                    else:
                        ytemp = Int.smooth(MW.xzero,MW.zerofull,1.8)
                        x_plot.append(MW.xzero)
                        y_plot.append(ytemp)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            ytemp = Int.smooth(MW.xzeroOL[i],MW.zerofullOL[i],1.8)
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(ytemp)
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Relative Abundance'
                else:
                    if MW.normmode.isChecked():
                        x_plot.append(MW.xzero)
                        y_plot.append(MW.zerofull/max(MW.zerofull)*100)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(MW.zerofullOL[i]/max(MW.zerofullOL[i])*100)
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Normalized Abundance (%)'
                    else:
                        x_plot.append(MW.xzero)
                        y_plot.append(MW.zerofull)
                        plot_labels.append(os.path.basename(str(MW.namebase)))
                        for i in range(len(MW.xzeroOL)):
                            x_plot.append(MW.xzeroOL[i])
                            y_plot.append(MW.zerofullOL[i])
                            plot_labels.append(os.path.basename(str(MW.namebaseOL[i])))
                        Yaxis_Label = 'Relative Abundance'
                MW.fig.clf()
                MW.ax2 = MW.fig.add_subplot(111)
                for i in range(len(x_plot)):
                    MW.ax2.plot(x_plot[i], y_plot[i],label=plot_labels[i],alpha=0.8)
                MW.ax2.legend(loc='upper right', title='File Names')
                MW.ax2.set_xlabel('Mass (Da)')
                MW.ax2.set_ylabel(Yaxis_Label)
                MW.fig.tight_layout()
                MW.canvas.draw()
        else:
            pass
   
    def ionmode_set(self):
        if MW.ionmode.isChecked():
            MW.proton = 1.00727647
            MW.ionmode_display.setText('Positive Ion Mode')
            MW.canvas.draw()
        else:
            MW.proton = -1.00727647
            MW.ionmode_display.setText('Negative Ion Mode')
            MW.canvas.draw()
    
    def MS_interp_replot(self):
        try:
            if MW.ft_pro == True:
                MW.ax.clear()
                if MW.MS_interp.isChecked():
                    MW.ax.plot(MW.xint,MW.yint)
                    MW.ax.set_title('Interpolated Mass Spectrum')
                else:
                    MW.ax.plot(MW.x,MW.y)
                    MW.ax.set_title('Raw Mass Spectrum')
                MW.ax.set_ylabel('Relative Abundance')
                MW.ax.set_xlabel('m/z')
                MW.fig.tight_layout()
                MW.canvas.draw()
            elif MW.stft_pro == True:
                MW.ax3.clear()
                if MW.MS_interp.isChecked():
                    MW.ax3.plot(MW.xint,MW.yint)
                    MW.ax3.set_title('Interpolated Mass Spectrum')
                else:
                    MW.ax3.plot(MW.x,MW.y)
                    MW.ax3.set_title('Raw Mass Spectrum')
                MW.ax3.set_xlabel('m/z')
                MW.ax3.set_ylabel('Abundance')
                MW.fig.tight_layout()
                MW.canvas.draw()
            else:
                pass
            
        except AttributeError:
            pass
    
    def FT_load(self):
        try:
            global agilentbool
            global agilentpath
            if agilentbool == True:
                MW.x, MW.y, MW.namebase = load.agilent(self, agilentpath)
            if agilentbool == False:
                MW.x, MW.y,MW.namebase = load.load_file_norm(self)
            agilentbool = False
            MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
            MW.datarho = len(MW.xint)
            MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)
            MW.maxF = max(MW.xFT)

            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)#put in upper left hand corner
            MW.ax2 = MW.fig.add_subplot(122)
            if MW.MS_interp.isChecked():
                MW.ax.plot(MW.xint,MW.yint)
                MW.ax.set_title('Interpolated Mass Spectrum')
            else:
                MW.ax.plot(MW.x,MW.y)
                MW.ax.set_title('Raw Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')

            MW.ax2.plot(MW.xFT,abs(MW.yFT))
            MW.ax2.set_title('Fourier Spectrum')
            MW.ax2.set_ylabel('Relative Amplitude')
            MW.ax2.set_xlabel('Frequency')
            MW.ax2.set_xlim(0, MW.xFT[-1] / 2)

            print('File: ' + MW.namebase)
            MW.file_display.setText('File: '+MW.namebase) 
            MW.fig.tight_layout()
            MW.canvas.draw() #shows the graphs
            
            try:
                MW.NFmin = MW.xFT[int(len(MW.xFT)/2 - len(MW.xFT)/10)]
                MW.NFmax = MW.xFT[int(len(MW.xFT)/2)]
                MW.NFrms = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin,MW.NFmax),4)
                MW.Nrms = np.round(MW.NFrms/np.sqrt(len(MW.x)),4)
                print('Frequency Noise RMSD: ' + str(MW.NFrms))
                print('Mass Noise RMSD: ' + str(MW.Nrms))
                
            except TypeError:
                print('Unable to automatically estimate Noise RMSD')
                MW.NFrms = 0
                MW.Nrms = 0

        except TypeError:
            MW.ref_num = 0
            print('Either no data loaded or wrong data type.  Try .csv or .txt for wrong data type')

        try:
            #Processing booleans
            MW.ft_pro = True
            MW.stft_pro = False
            MW.stft_selec_man = "NA"
            MW.stft_box_thresh = "NA"
            MW.stft_range_thresh = "NA"
            MW.decon_pro = False
            #Mass spectrum baseline and domain adjustment booleans
            MW.setdatarho = False
            MW.domaintrunc = False
            MW.Nsubtraction = False
            MW.Nlinear_sub = False
            MW.Nseg_sub = False
            MW.ycopy = MW.y.copy()
            MW.xcopy = MW.x.copy()
            MW.multiple_series = False
            #Integration booleans
            MW.NZbased = False
            MW.basecor = False
            MW.based = False
            MW.smoothed = False
            MW.ReInt = False 
            MW.range2 = False
            #Saving booleans
            MW.saveFT = False
            MW.saveIFFT = False
            MW.saveZspec = False
            MW.saveSeries = False
            MW.batchcheck = False
            MW.full_save = False
            MW.saveMS = False
            try:
                del(MW.batchfolder)
            except AttributeError:
                pass
            try:
                del(MW.zerofullstore)
            except AttributeError:
                pass
            try:
                del(iso.xrot1)
            except AttributeError:
                pass
            try:
                del(MW.x_overlay)
            except AttributeError:
                pass
            try:
                del(MW.xrefcs)
            except AttributeError:
                pass
            try:
                del(MW.inclusion)
            except AttributeError:
                pass
        except AttributeError:
            print('Unable to reset processing parameters')
    
    def data_overlay(self):
        try:
            MW.x_overlay, MW.y_overlay,MW.namebase_overlay = load.load_file_norm(self)
            self.iso_plot()
        except TypeError:
            MW.ref_num = 0
            print('Either no data loaded or wrong data type.  Try .csv or .txt for wrong data type')
    
    def iso_plot(self):
        MW.fig.clf()
        MW.ax = MW.fig.add_subplot(111)
        try:
            MW.ax.plot(iso.xrot1, iso.IFFT_norm1,label='isotope distribution with peak width')
            MW.ax.plot(iso.xrot2, iso.IFFT_norm2,label='isotope-fine distribution')
            text = 'Normalized to the highest peak, with a probability density of ' + str(max(iso.IFFT))
            MW.ax.text(0.5, 0.9,text, horizontalalignment='center', verticalalignment='center',
                        transform=MW.ax.transAxes,color='g',fontsize=12)
        except AttributeError:
            pass
        try:
            MW.ax.plot(MW.x_overlay,(MW.y_overlay/max(MW.y_overlay)),label='imported spectrum')
            MW.file_display.setText('File: '+str(MW.namebase_overlay))
        except AttributeError:
            pass
        MW.ax.set_title('Zero Charge Spectrum')
        MW.ax.set_ylabel('Relative Abundance')
        MW.ax.set_xlabel('Mass (Da)')
        MW.ax.legend(loc='upper right')

        MW.fig.tight_layout()
        MW.canvas.draw()
    
    def iso_export(self):
        try:
            MW.iso_x,MW.iso_y = iso.xrot1, iso.IFFT_norm1
            MW.ReturnFunc = MW.iso_save
            Folder_Select(self).show()

        except AttributeError:
            print('Please calculate isotope distribution first')
            MW.message = 'Please calculate isotope distribution first'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
    
    def iso_save(self):
        directory = os.path.dirname(MW.namebase)
        label = os.path.splitext(MW.namebase)[0] + '_iso_dist.csv'
        temp_path = os.path.join(directory, label)
        
        np.savetxt(temp_path, np.c_[MW.iso_x,MW.iso_y], delimiter=',')
        MW.message = str('Saved isotope distribution to:\n'+temp_path)
        print(MW.message)
        MW.Err = False
        self.dialog = iFAMS_message(self)
        self.dialog.show()
    
    def plotAbs(self):
        """
        plots the absolute value of the Fourier data
        :return:
        """
        try:
            MW.tIFT = np.abs(MW.IFT)
            MW.xzero,MW.cszero,MW.zerofull = FT.zerocharge(MW.tIFT,MW.xint,MW.cs,MW.proton)
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)
            MW.ax.plot(MW.x, MW.y,label = 'Original Spectrum')
            for i in range(len(MW.tIFT)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax.plot(MW.xint,MW.tIFT[i],label = label)
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')
            MW.ax.legend(loc='upper right')

            MW.ax2.plot(MW.xzero, MW.zerofull,label = 'Full Zero',color='mediumseagreen')
            for i in range(len(MW.cszero)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax2.plot(MW.xzero,MW.cszero[i],label = label)
            MW.ax2.set_title('Zero Charge Spectrum')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.ax2.legend(loc = 'upper right')
            MW.fig.tight_layout()
            MW.canvas.draw()

        except AttributeError:
            print('Please run this command after iFAMS analysis')
            MW.message = 'Please run this command after iFAMS analysis'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

    def plotReal(self):
        """
        plots the real value of the Fourier data
        :return:
        """
        try:
            MW.tIFT = np.real(MW.IFT)
            MW.xzero,MW.cszero,MW.zerofull = FT.zerocharge(MW.tIFT,MW.xint,MW.cs,MW.proton)
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(121)
            MW.ax2 = MW.fig.add_subplot(122)
            MW.ax.plot(MW.x, MW.y,label = 'Original Spectrum')
            for i in range(len(MW.tIFT)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax.plot(MW.xint,MW.tIFT[i],label = label)
            MW.ax.set_title('Mass Spectrum')
            MW.ax.set_ylabel('Relative Abundance')
            MW.ax.set_xlabel('m/z')
            MW.ax.legend(loc='upper right')

            MW.ax2.plot(MW.xzero,MW.zerofull,label='Full Zero',color='mediumseagreen')
            for i in range(len(MW.cszero)):
                label = 'Charge State ' + str(int(MW.cs[i]))
                MW.ax2.plot(MW.xzero,MW.cszero[i],label = label)
            MW.ax2.set_title('Zero Charge Spectrum')
            MW.ax2.set_ylabel('Relative Abundance')
            MW.ax2.set_xlabel('Mass (Da)')
            MW.ax2.legend(loc='upper right')
            MW.fig.tight_layout()
            MW.canvas.draw()

        except AttributeError:
            print('Please run this command after iFAMS analysis')
            MW.message = 'Please run this command after iFAMS analysis'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()


    def STFT_load(self):
        try:
            global agilentbool
            global agilentpath
            if agilentbool == True:
                MW.x, MW.y, MW.namebase = load.agilent(self, agilentpath)
            if agilentbool == False:
                MW.x, MW.y, MW.namebase = load.load_file_norm(self)
            agilentbool = False
            self.STFT_Prep()
        except TypeError:
            return
        except ValueError:
            print('Unable to load file. Check file type and try again')
            MW.message = 'Unable to load file. Check file type and try again'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
        
    def STFT_clipboard(self):
        try:
            text = QApplication.clipboard().text()
            dataline = True
            for i in range(len(text)):
                if dataline == False:
                    if text[i] == '\n':
                        dataline = True
                        continue
                    else:
                        continue
                try:
                    int(text[i])
                    start = i
                    break
                except ValueError:
                    dataline = False
                    continue
            columns = 1
            for i in range(start,len(text)):
                if text[i] == '\n':
                    break
                if text[i] == ',' or text[i] == '\t' or text[i] == '':
                    columns += 1
    
            temp = ''
            desig = []            
            x = []
            y = []
            count = 1
            for i in range(start,len(text)):
                if text[i] == ',' or text[i] == '\n' or text[i] == '\t' or text[i] == '':
                    temp = float(temp)
                    if columns == 3:
                        if count % 3 == 0:
                            y.append(temp)
                        elif count % 3 == 2:
                            x.append(temp)
                        else:
                            desig.append(temp)
                    else:
                        if count % 2 == 0:
                            y.append(temp)
                        else:
                            x.append(temp)
                    count+=1
                    temp = ''
                else:
                    temp = temp + str(text[i])
            if len(x) == len(y):
                pass
            else:
                y.append(float(temp))
            MW.x = x
            MW.y = y
            MW.namebase = 'From Clipboard'
            self.STFT_Prep()
        except UnboundLocalError:
            print('Unable to plot clipboard. Re-copy spectrum and try again')
            MW.message = 'Unable to plot clipboard.\nRe-copy spectrum and try again'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
    
    def STFT_Prep(self):
        try:
            if MW.y[0] == 0 and MW.y[-1]==0:
                MW.x = MW.x[1:-1].copy()
                MW.y = MW.y[1:-1].copy()
                     
            MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
            MW.datarho = len(MW.xint)
            MW.xFT,MW.yFT = FT.Fourier(MW.xint,MW.ypadd,MW.po2)
          
            print('File: ' + MW.namebase)
            MW.file_display.setText('File: '+MW.namebase) 
            MW.BL_load = False
            
            try:
                #Estimating Gabor parameters
                MW.maxF = float(MW.xFT[-1])
                MW.maxab = float(max(MW.yint))
                ftW = len(MW.yFT)/MW.maxF
                
                if np.max(np.absolute(MW.yFT[int(8*ftW):int(len(MW.yFT)/2)])) > 10*np.average(np.absolute(MW.yFT[int(8*ftW):int(len(MW.yFT)/2)])):
                    MW.FoI = MW.xFT[int(8*ftW)+np.argmax(np.absolute(MW.yFT[int(8*ftW):int(len(MW.yFT)/2)]))]
                elif np.max(np.absolute(MW.yFT[int(2*ftW):int(10*ftW)])) > 8*np.average(np.absolute(MW.yFT[int(2*ftW):int(10*ftW)])):
                    MW.FoI = MW.xFT[int(2*ftW)+np.argmax(np.absolute(MW.yFT[int(2*ftW):int(10*ftW)]))]
                elif np.max(np.absolute(MW.yFT[int(0.15*ftW):int(5*ftW)])) > 6*np.average(np.absolute(MW.yFT[int(0.15*ftW):int(5*ftW)])):
                    MW.FoI = MW.xFT[int(0.15*ftW)+np.argmax(np.absolute(MW.yFT[int(0.15*ftW):int(5*ftW)]))]
                else:
                    MW.FoI = MW.xFT[int(0.05*ftW)+np.argmax(np.absolute(MW.yFT[int(0.05*ftW):int(len(MW.yFT)/2)]))]
                                
                if MW.FoI >= 8:
                    MW.winnum = int((MW.maxF)*20)
                    MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -5)*len(MW.yFT)/MW.maxF):int((MW.FoI +5)*len(MW.yFT)/MW.maxF)]))/20)
                elif MW.FoI < 8 and MW.FoI > 1:
                    MW.winnum = int((MW.maxF)*30)
                    MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -1)*len(MW.yFT)/MW.maxF):int((MW.FoI +1)*len(MW.yFT)/MW.maxF)]))/10)
                elif MW.FoI <= 1 and MW.FoI >= 0.1:
                    MW.winnum = int((MW.maxF)*80)
                    MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -0.05)*len(MW.yFT)/MW.maxF):int((MW.FoI +0.05)*len(MW.yFT)/MW.maxF)]))/10)
                else:
                    MW.winnum = int((MW.maxF)*1000)
                    MW.vmax = int(np.max(np.absolute(MW.yFT[int((MW.FoI)*len(MW.yFT)/MW.maxF)-3:int((MW.FoI)*len(MW.yFT)/MW.maxF)+3]))/10)
                
                if MW.winnum < 10:
                    MW.winnum = 10
                if len(MW.xint)/MW.winnum < 8:
                    MW.winnum = int(len(MW.xint)/8)
            except IndexError:
                print('Unable to automatically resolve Gabor spectrogram')
                MW.FoI = max(MW.xFT)/3
                MW.winnum = int((MW.maxF)*80)
                MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -0.05)*len(MW.yFT)/MW.maxF):int((MW.FoI +0.05)*len(MW.yFT)/MW.maxF)]))/10)
                MW.message = 'Unable to automatically resolve Gabor spectrogram'
                MW.submessage = 'Manually adjust STFT parameters to improve spectrogram'
                MW.Err = False
                self.dialog = iFAMS_message(self)
                self.dialog.show()
            except ValueError:
                print('Unable to automatically resolve Gabor spectrogram')
                MW.FoI = max(MW.xFT)/3
                MW.winnum = int((MW.maxF)*80)
                MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -0.05)*len(MW.yFT)/MW.maxF):int((MW.FoI +0.05)*len(MW.yFT)/MW.maxF)]))/10)
                MW.message = 'Unable to automatically resolve Gabor spectrogram'
                MW.submessage = 'Manually adjust STFT parameters to improve spectrogram'
                MW.Err = False
                self.dialog = iFAMS_message(self)
                self.dialog.show()
                
            try:
                maxF = int(MW.maxF/2)
                MW.NFmin2 = MW.xFT[int((int(maxF*4/5)+0.4)*len(MW.yFT)/MW.maxF)]
                MW.NFmax2 = MW.xFT[int((int(maxF*4/5)+0.6)*len(MW.yFT)/MW.maxF)]
                MW.NFrms2 = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin2,MW.NFmax2),4)
                MW.NFmin3 = MW.xFT[int((int(maxF*3/5)+0.4)*len(MW.yFT)/MW.maxF)]
                MW.NFmax3 = MW.xFT[int((int(maxF*3/5)+0.6)*len(MW.yFT)/MW.maxF)]
                MW.NFrms3 = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin3,MW.NFmax3),4)
                MW.NFmin4 = MW.xFT[int((int(maxF*2/5)+0.4)*len(MW.yFT)/MW.maxF)]
                MW.NFmax4 = MW.xFT[int((int(maxF*2/5)+0.6)*len(MW.yFT)/MW.maxF)]
                MW.NFrms4 = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin4,MW.NFmax4),4)
                MW.NFrms = np.average([MW.NFrms2,MW.NFrms3,MW.NFrms4])
                MW.Nrms = np.round(MW.NFrms/np.sqrt(len(MW.x)),4)
                MW.NFmin = MW.NFmin3
                MW.NFmax = MW.NFmax3
                print('Frequency Noise RMSD: ' + str(np.round(MW.NFrms)))
                print('Mass Noise RMSD: ' + str(np.round(MW.Nrms)))
            except IndexError:
                print('Unable to automatically estimate Noise RMSD')
                MW.NFrms = 0
                MW.Nrms = 0
            except TypeError:
                print('Unable to automatically estimate Noise RMSD')
                MW.NFrms = 0
                MW.Nrms = 0
            
            MW.OF = 10
                
        except TypeError:
            print('Either no data loaded or wrong data type.  Try .csv or .txt for wrong data type')
            return
        except ValueError:
            print('Unable to load file. Check file type and try again')
            MW.message = 'Unable to load file. Check file type and try again'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return
        
        try:
            #Processing booleans
            MW.ft_pro = False
            MW.stft_pro = True
            MW.stft_selec_man = True
            MW.stft_box_thresh = "NA"
            MW.stft_range_thresh = "NA"
            MW.decon_pro = False
            #Mass spectrum baseline and domain adjustment booleans
            MW.setdatarho = False
            MW.domaintrunc = False
            MW.Nsubtraction = False
            MW.Nlinear_sub = False
            MW.Nseg_sub = False
            MW.ycopy = MW.y.copy()
            MW.xcopy = MW.x.copy()
            MW.multiple_series = False
            #Integration booleans
            MW.NZbased = False
            MW.basecor = False
            MW.based = False
            MW.smoothed = False
            MW.ReInt = False 
            MW.range2 = False
            #Saving booleans
            MW.saveFT = False
            MW.saveIFFT = False
            MW.saveZspec = False
            MW.saveSeries = False
            MW.batchcheck = False
            MW.batchstore = []
            MW.full_save = False
            MW.saveMS = False
            try:
                del(MW.batchfolder)
            except AttributeError:
                pass
            try:
                del(MW.XBL2)
            except AttributeError:
                pass
            try:
                del(MW.zerofullstore)
            except AttributeError:
                pass
            try:
                del(iso.xrot1)
            except AttributeError:
                pass
            try:
                del(MW.x_overlay)
            except AttributeError:
                pass
            try:
                del(MW.xrefcs)
            except AttributeError:
                pass
            try:
                del(MW.inclusion)
            except AttributeError:
                pass
        except AttributeError:
            print('Unable to reset processing parameters')

        self.plot_gabor()


    def plot_gabor(self):
        try:
            MW.fig.clf()
            MW.stdev = int(0.1 * MW.winnum)
            MW.window = str('gaussian')
            MW.X = STFT.stft(MW.ypadd, MW.winnum, MW.window, MW.OF)
            ytemp = STFT.re_stft(MW.X, MW.winnum, MW.ypadd, MW.yint,MW.OF)
            MW.correction_factor = max(MW.yint)/max(ytemp)
            MW.Xtemp = MW.X*0
            MW.Xflat = abs(MW.X.flatten("C")) / 100
            MW.mzspacing = float((MW.xint[-1] - MW.xint[0])/len(MW.xint))
            MW.xshort = np.linspace((MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing + min(MW.xpadd), max(MW.xpadd) - (MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing, len(MW.X))
            MW.yshort = np.linspace(0 - (max(MW.xFT)/len(MW.X[0]))/2, max(MW.xFT) + (max(MW.xFT)/len(MW.X[0]))/2, len(MW.X[0]))
            MW.mzchan = float((MW.xint[-1]-MW.xshort[0])/((MW.xshort[-1]-MW.xshort[0])/len(MW.xshort)))
            
            MW.ax = MW.fig.add_subplot(221)
            MW.ax4 = MW.fig.add_subplot(222, sharey=MW.ax)
            MW.ax3 = MW.fig.add_subplot(223, sharex=MW.ax)
            MW.ax2 = MW.fig.add_subplot(224)

            MW.ax.imshow(abs(MW.X.T), origin='lower', aspect='auto',
                         interpolation='nearest', extent=(MW.xshort[0],MW.xshort[-1],MW.yshort[0],MW.yshort[-1]), cmap='jet',vmax=MW.vmax)
            MW.ax.set_title('Gabor Spectrogram')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Frequency')
            MW.ax.set_xlim(MW.xint[0], MW.xint[-1])
            MW.ax.set_ylim(MW.xFT[0]-MW.maxF/MW.winnum/2,MW.xFT[int(np.ceil(len(MW.xFT)/2))]+MW.maxF/MW.winnum/2)

            MW.ax4.plot(abs(MW.yFT), MW.xFT)
            MW.ax4.set_title('Fourier Spectrum')
            MW.ax4.set_xlabel('Amplitude')
            MW.ax4.set_ylabel('Frequency')

            if MW.MS_interp.isChecked():
                MW.ax3.plot(MW.xint,MW.yint)
                MW.ax3.set_title('Interpolated Mass Spectrum')
            else:
                MW.ax3.plot(MW.x,MW.y)
                MW.ax3.set_title('Raw Mass Spectrum')
            MW.ax3.set_xlabel('m/z')
            MW.ax3.set_ylabel('Abundance')

            MW.ax2.plot(MW.x, MW.y, color='darkgrey')
            MW.ax2.set_title('Reconstructed Spectrum')
            MW.ax2.set_xlabel('m/z')
            MW.ax2.set_ylabel('Abundance')
            
            MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)
            MW.glist = []
            MW.rlist = []
            MW.zlist = []
            
            MW.fig.tight_layout()
            MW.canvas.draw()
            
        except AttributeError:
            print('Unable to plot Gabor Transform')
        except IndexError:
            print('Unable to perform STFT')
            MW.message = 'Insufficient data. Unable to perform STFT.'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            return

    def onselect(self, eclick, erelease):
        if eclick.xdata == erelease.xdata and eclick.ydata == erelease.ydata:
            checki = int(np.round((eclick.xdata-MW.xshort[0])/(MW.xshort[2]-MW.xshort[1])))
            checkiy = int(np.round(eclick.ydata/(MW.yshort[2]-MW.yshort[1])))
            boxleft = checki-2
            boxright = checki+2
            boxbott = checkiy-5
            boxtop = checkiy+5
            AVG = np.min(np.absolute(MW.X[...,checkiy]))
            for j in range(1,20):
                if np.absolute(MW.X[checki-j][checkiy])-AVG < 0.1*(np.absolute(MW.X[checki][checkiy])-AVG):
                    boxleft = checki-j
                    break
                if np.absolute(MW.X[checki-j][checkiy])-AVG > 1.05*(np.absolute(MW.X[checki-j+1][checkiy])-AVG):
                    if np.absolute(MW.X[checki-j-1][checkiy])-AVG < 1.05*(np.absolute(MW.X[checki-j+1][checkiy])-AVG):
                        continue
                    else:
                        boxleft = checki-j+1
                        break
            for j in range(1,20):
                if np.absolute(MW.X[checki+j][checkiy])-AVG < 0.1*(np.absolute(MW.X[checki][checkiy])-AVG):
                    boxright = checki+j
                    break
                if np.absolute(MW.X[checki+j][checkiy])-AVG > 1.05*(np.absolute(MW.X[checki+j-1][checkiy])-AVG):
                    if np.absolute(MW.X[checki+j+1][checkiy])-AVG < 1.05*(np.absolute(MW.X[checki+j-1][checkiy])-AVG):
                        continue
                    else:
                        boxright = checki+j-1
                        break
            for j in range(1,30):
                if np.absolute(MW.X[checki][checkiy-j])-AVG < 0.1*(np.absolute(MW.X[checki][checkiy])-AVG):
                    boxbott = checkiy-j
                    break
            for j in range(1,30):
                if np.absolute(MW.X[checki][checkiy+j])-AVG < 0.1*(np.absolute(MW.X[checki][checkiy+j])-AVG):
                    boxtop = checkiy+j
                    break
            boxbott = checkiy - int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
            boxtop = checkiy + int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
            boxH = boxtop-boxbott
            if eclick.ydata > 4:
                try:
                    beati = np.argmax(np.absolute(MW.X[checki,boxtop:checkiy+boxH]))
                    beattop = np.absolute(MW.X[checki][boxtop+beati])
                    beatbott = np.absolute(MW.X[checki][boxbott-beati])
                    if abs(beattop-beatbott) < 0.01*(np.absolute(MW.X[checki][checkiy])-AVG) and np.average((beattop,beatbott))-AVG > 0.2*(np.absolute(MW.X[checki][checkiy])-AVG):
                        for j in range(1,boxH):
                            if np.absolute(MW.X[checki][boxbott-beati-j])-AVG < 0.2*(np.absolute(MW.X[checki][boxbott-beati])-AVG):
                                boxbott = boxbott-beati-j
                                break
                        for j in range(1,boxH):
                            if np.absolute(MW.X[checki][boxtop+beati+j])-AVG < 0.2*(np.absolute(MW.X[checki][boxtop+beati])-AVG):
                                boxtop = boxtop+beati+j
                                break
                        boxbott = checkiy - int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
                        boxtop = checkiy + int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
                except ValueError:
                    pass
            if boxleft <= 0:
                boxleft = 0
            if boxright >= len(MW.xshort):
                boxright = len(MW.xshort)-1
            MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = MW.xshort[boxleft], MW.yshort[boxbott]
            MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = MW.xshort[boxright], MW.yshort[boxtop]
        else:
            MW.toggle_selector.RS.x1, MW.toggle_selector.RS.y1 = eclick.xdata, eclick.ydata
            MW.toggle_selector.RS.x2, MW.toggle_selector.RS.y2 = erelease.xdata, erelease.ydata
        if MW.toggle_selector.RS.x1 < MW.toggle_selector.RS.x2:
            MW.x1 = MW.toggle_selector.RS.x1
            MW.x2 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y1 < MW.toggle_selector.RS.y2:
            MW.y1 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y2 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)
        if MW.toggle_selector.RS.x2 < MW.toggle_selector.RS.x1:
            MW.x2 = MW.toggle_selector.RS.x1
            MW.x1 = MW.toggle_selector.RS.x2
        if MW.toggle_selector.RS.y2 < MW.toggle_selector.RS.y1:
            MW.y2 = max(MW.toggle_selector.RS.y1,-MW.maxF/MW.winnum/2)
            MW.y1 = max(MW.toggle_selector.RS.y2,-MW.maxF/MW.winnum/2)

    def toggle_selector(event):
        if event.key == 'e':
            print('it works')

    def add_gabor(self):

        try:
            rex = Rectangle((MW.x1, MW.y1), (MW.x2 - MW.x1), (MW.y2 - MW.y1),
                            facecolor='none', edgecolor='r', linewidth=1)
            MW.rlist.append(
                [MW.x1, MW.x2, MW.y1, MW.y2, len(MW.glist)])
            MW.ax.add_artist(rex) 

            MW.ax2.clear()
            piW = MW.xshort[2]-MW.xshort[1]
            piH = MW.yshort[2]-MW.yshort[1]
            i1 = int(np.round((MW.x1-MW.xshort[0])/piW))
            i2 = int(np.round((MW.x2-MW.xshort[0])/piW))+1
            i3 = int(np.round(MW.y1/piH))
            i4 = int(np.round(MW.y2/piH))
            MW.Xtemp[i1:i2,i3:i4+1] = MW.X[i1:i2,i3:i4+1]
            if i3 <= 1:
                MW.Xtemp[i1:i2,-i4:len(MW.X[0])] = MW.X[i1:i2,-i4:len(MW.X[0])]
            else:
                MW.Xtemp[i1:i2,-i4:-i3+1] = MW.X[i1:i2,-i4:-i3+1]
            MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor

            MW.ax2.plot(MW.x, MW.y, color='darkgrey')
            if len(MW.glist) == 0:
                if MW.realsel.isChecked() == False:
                    MW.ax2.plot(MW.xint,abs(MW.yrecon),color = color[len(MW.glist)])
                if MW.realsel.isChecked() == True:
                    MW.ax2.plot(MW.xint,np.real(MW.yrecon),color = color[len(MW.glist)])
            if len(MW.glist) > 0:
                if MW.realsel.isChecked() == False:
                    for i in range(len(MW.glist)):
                        MW.ax2.plot(MW.xint,abs(MW.glist[i]),color = color[i])
                    MW.ax2.plot(MW.xint, abs(MW.yrecon), color=color[len(MW.glist)])
                if MW.realsel.isChecked() == True:
                    for i in range(len(MW.glist)):
                        MW.ax2.plot(MW.xint,np.real(MW.glist[i]),color = color[i])
                    MW.ax2.plot(MW.xint, np.real(MW.yrecon), color=color[len(MW.glist)])
                        
            MW.ax2.set_title('Reconstructed Spectrum')
            MW.ax2.set_xlabel('m/z')
            MW.ax2.set_ylabel('Abundance')
            MW.fig.tight_layout()
            MW.canvas.draw()
        except AttributeError:
            print('No rectangle on screen.  Please draw rectangle first')
            MW.message = 'No rectangle on screen. \n Please draw rectangle first'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            
            
    def save_gabor_sel(self):
        
        try:
            MW.glist.append(MW.yrecon)
            
            for i in range(len(MW.rlist)):
                rex = Rectangle((MW.rlist[i][0], MW.rlist[i][2]), (MW.rlist[i][1] - MW.rlist[i][0]), (MW.rlist[i][3] -
                                                                                                      MW.rlist[i][2]),
                                facecolor='none', edgecolor=color[MW.rlist[i][4]], linewidth=3)
                MW.ax.add_artist(rex)
            MW.Xtemp = MW.X*0
            #print('Rectangle size: ' +str(np.round(MW.rlist[-1][1]-MW.rlist[-1][0],2))+ ' x ' +str(np.round(MW.rlist[-1][3]-MW.rlist[-1][2],4)))
            MW.fig.tight_layout()
            MW.canvas.draw()
            
        except AttributeError:
            print('No rectangle on screen.  Please draw rectangle first')
            MW.message = 'No rectangle on screen. \n Please draw rectangle first'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

    def del_gabor(self):
        try:
            len1 = MW.rlist[-1][4]
            len2 = len(MW.glist)-1
            if len1 == len2:
                MW.glist.pop()
                try:
                    if len(MW.zlist)-len(MW.glist) == 1:
                        MW.zlist.pop()
                    else:
                        MW.zlist = []
                except AttributeError:
                    pass
                try:
                    del(MW.glist_store)
                except AttributeError:
                    pass
                try:
                    MW.glistBL2.pop()
                except AttributeError:
                    pass
                except IndexError:
                    pass
                rlisttemp = []
                for i in range(len(MW.rlist)):
                    if MW.rlist[i][4] <= (len(MW.glist) - 1):
                        rlisttemp.append(MW.rlist[i])
                MW.rlist = rlisttemp
            if len1 > len2:
                rlisttemp = []
                for i in range(len(MW.rlist)):
                    if MW.rlist[i][4] <= (len(MW.glist)-1):
                        rlisttemp.append(MW.rlist[i])
                MW.rlist = rlisttemp
            
            try:
                plottop = MW.ax.get_ylim()[1]
                plotbott = MW.ax.get_ylim()[0]
                plotxL = MW.ax.get_xlim()[0]
                plotxR = MW.ax.get_xlim()[1]
            except AttributeError:
                pass
            MW.ax2.clear()
            MW.ax.clear()
            MW.ax.imshow(abs(MW.X.T), origin='lower', aspect='auto',
                         interpolation='nearest', extent=(MW.xshort[0], MW.xshort[-1], MW.yshort[0], MW.yshort[-1]),
                         cmap='jet', vmax=MW.vmax)
            MW.ax.set_title('Gabor Spectrogram')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Frequency')
            try:
                if plottop > MW.yshort[-1]:
                    MW.ax.set_xlim(MW.xint[0], MW.xint[-1])
                    MW.ax.set_ylim(MW.xFT[0]-MW.maxF/MW.winnum/2,MW.xFT[int(np.ceil(len(MW.xFT)/2))]+MW.maxF/MW.winnum/2)
                else:
                    MW.ax.set_ylim(plotbott,plottop)
                    MW.ax.set_xlim(plotxL,plotxR)
            except UnboundLocalError:
                pass
            MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)

            MW.ax2.plot(MW.x, MW.y, color='darkgrey')
            for i in range(len(MW.rlist)):
                rex = Rectangle((MW.rlist[i][0], MW.rlist[i][2]), (MW.rlist[i][1] - MW.rlist[i][0]), (MW.rlist[i][3] -
                                                                                                      MW.rlist[i][2]),
                                facecolor='none', edgecolor=color[MW.rlist[i][4]], linewidth=3)
                MW.ax.add_artist(rex)
            
            if MW.realsel.isChecked() == False:
                for i in range(len(MW.glist)):
                    MW.ax2.plot(MW.xint, abs(MW.glist[i]), color=color[i])
            if MW.realsel.isChecked() == True:
                for i in range(len(MW.glist)):
                    MW.ax2.plot(MW.xint, np.real(MW.glist[i]), color=color[i])
                       
            MW.ax2.set_title('Reconstructed Spectrum')
            MW.ax2.set_xlabel('m/z')
            MW.ax2.set_ylabel('Abundance')
            MW.fig.tight_layout()
            MW.canvas.draw()

        except AttributeError:
            print('No data loaded')
            MW.message = 'No data loaded'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

        except IndexError:
            print('No more rectangles to delete')
            MW.message = 'No more rectangles to delete'
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
 
        
    def show_gabor_only(self):
        MW.fig.clf()
        MW.ax = MW.fig.add_subplot(111)
        MW.ax.imshow(abs(MW.X.T), origin='lower', aspect='auto',
                     interpolation='nearest', extent=(MW.xshort[0], MW.xshort[-1], MW.yshort[0], MW.yshort[-1]),
                     cmap='jet', vmax=MW.vmax)
        MW.ax.set_title('Gabor Spectrogram')
        MW.ax.set_xlabel('m/z')
        MW.ax.set_ylabel('Frequency')
        MW.ax.set_xlim(MW.xint[0], MW.xint[-1])
        MW.ax.set_ylim(MW.xFT[0]-MW.maxF/MW.winnum/2,MW.xFT[-1]+MW.maxF/MW.winnum/2)
        MW.fig.tight_layout()
        MW.canvas.draw()

    def export_data(self):
        MW.saveFT = False
        MW.saveIFFT = False
        MW.saveZspec = False
        MW.saveSeries = False
        MW.full_save = True
        ### saves parameters file
        try:
            del(MW.batchfolder)
        except AttributeError:
            pass
        try:
            MW.preferences
        except AttributeError:
            try:
                MW.batch_param(self)
            except AttributeError:
                pass
            except ValueError:
                pass
        MW.ReturnFunc = MW.export_Confirm
        Folder_Select(self).show()
        
    def export_Confirm(self):
        MW.save_FT(self)
        try:
            day = str(datetime.datetime.now().day)
            if int(day) < 10:
                day = str(0)+day
            month = str(datetime.datetime.now().month)
            if int(month) < 10:
                month = str(0)+month
            year = str(int(float(datetime.datetime.now().year)-2000))
            label = str(year+month+day)+'_batch_params_from_' + os.path.basename(MW.namebase) + '.csv'
            temp_path = os.path.join(MW.foldername, label)
            np.savetxt(temp_path, np.c_[MW.preferences,MW.preflabel], delimiter=',',fmt='%s')
        except AttributeError:
            pass
        except ValueError:
            print('Error with saving batch parameters')
            MW.submessage = str('Error with saving batch parameters')
        MW.message = str('Data exported to: \n'+str(os.path.basename(MW.foldername)))    
        MW.Err = False
        self.dialog = iFAMS_message(self)
        self.dialog.show()

    def save_FT(self):
        if MW.namebase == 'From Clipboard':
            MW.full_save = True
            MW.saveMS = True
            MW.ReturnFunc = MW.save_FT
            Folder_Select(self)
            return

        if MW.saveZspec == True:
            try:
                MW.xrefcs
            except AttributeError:
                #try:
                MW.deconcsy = []
                MW.deconcsx = []
                if len(MW.batchstore) == 0:
                    for k in range(len(MW.cszero)):
                        cstemp = MW.cszero[k]
                        csBLtemp = MW.cszeroBL2[k]
                        csxtemp = MW.xzero
                        check = np.max(csBLtemp)/4
                        for h in range(0,len(csBLtemp)):
                            if csBLtemp[h] > check:
                                mid1 = h
                                break
                        for h in range(mid1+1,len(csBLtemp)):
                            if csBLtemp[h] <= check or csBLtemp[h] <= csBLtemp[mid1]:
                                mid2 = h
                                break
                        int1 = min(sum(cstemp[0:mid1]),sum(cstemp[mid2:len(cstemp)]))
                        int2 = (sum(csBLtemp[0:mid1])+sum(csBLtemp[mid2:len(csBLtemp)]))/2
                        MW.deconcsy.append(np.asarray(cstemp)-np.asarray(csBLtemp)*int1/int2)
                        MW.deconcsx.append(csxtemp)            
                else:
                    for j in range(len(MW.batchstore)):
                        MW.cszero = MW.batchstore[j][5]
                        MW.cszeroBL2 = MW.batchstore[j][6]
                        csxzero = MW.batchstore[j][3]
                        for k in range(len(MW.cszero)):
                            check = np.max(MW.cszeroBL2[k])/4
                            for h in range(0,len(MW.cszeroBL2[k])):
                                if MW.cszeroBL2[k][h] > check:
                                    mid1 = h
                                    break
                            for h in range(mid1+1,len(MW.cszeroBL2[k])):
                                if MW.cszeroBL2[k][h] <= check or MW.cszeroBL2[k][h] <= MW.cszeroBL2[k][mid1]:
                                    mid2 = h
                                    break
                            int1 = min(sum(MW.cszero[k][0:mid1]),sum(MW.cszero[k][mid2:len(MW.cszero[k])]))
                            int2 = (sum(MW.cszeroBL2[k][0:mid1])+sum(MW.cszeroBL2[k][mid2:len(MW.cszeroBL2[k])]))/2
                            MW.deconcsy.append(np.asarray(MW.cszero[k])-np.asarray(MW.cszeroBL2[k])*int1/int2)
                            MW.deconcsx.append(csxzero)            

                MW.xrefcs = np.zeros(len(MW.cs)).tolist()
                MW.yrefcs = np.zeros(len(MW.cs)).tolist()
                MW.sumlistcs = np.zeros(len(MW.cs)).tolist()
                MW.heightcs = np.zeros(len(MW.cs)).tolist()
                MW.xcentcs = np.zeros(len(MW.cs)).tolist()
                MW.peakindexcs = np.zeros(len(MW.cs)).tolist()
                if MW.paramMatch == True:
                    MW.bounds = []
                    for i in range(len(MW.xref)):
                        MW.bounds.append(MW.xref[i][0])
                        MW.bounds.append(MW.xref[i][-1])
                for i in range(len(MW.cs)):
                    MW.xrefcs[i], MW.yrefcs[i], MW.sumlistcs[i], MW.heightcs[i], MW.xcentcs[i], MW.peakindexcs[i] = Int.boundsmatch(MW.deconcsx[i],MW.deconcsy[i],MW.bounds)
                # except AttributeError:
                #     print('error')

        print('Saving files')
        try:
            MW.batchfolder
        except AttributeError:
            day = str(datetime.datetime.now().day)
            if int(day) < 10:
                day = str(0)+day
            month = str(datetime.datetime.now().month)
            if int(month) < 10:
                month = str(0)+month
            year = str(int(float(datetime.datetime.now().year)-2000))
            count = 1
            folder = str(year+month+day+'_'+os.path.basename(MW.namebase)+'_'+str(count))
            test_path = os.path.exists(os.path.join(os.path.dirname(MW.namebase),folder))
            while test_path == True:
                count += 1
                folder = str(year+month+day+'_'+os.path.basename(MW.namebase)+'_'+str(count))
                test_path = os.path.exists(os.path.join(os.path.dirname(MW.namebase),folder))
            MW.foldername = os.path.join(os.path.dirname(MW.namebase),folder)
            os.mkdir(MW.foldername)
            
        else:
            if MW.saveMSRT == True:
                global RT1
                global RT2
                global RTfull
                if eval(RTfull.title()) == True:
                    MW.foldernameMS = os.path.join(MW.batchfolder,'Extracted_MS_from_fullRT')
                else:
                    MW.foldernameMS = os.path.join(MW.batchfolder,'Extracted_MS_from_RT'+RT1+'-'+RT2)
                if os.path.exists(MW.foldernameMS) == False:
                    os.mkdir(MW.foldernameMS)
            if MW.saveFT == True:
                MW.foldernameFT = os.path.join(MW.batchfolder,'Fourier_spectra')
                if os.path.exists(MW.foldernameFT) == False:
                    os.mkdir(MW.foldernameFT)
            if MW.saveIFFT == True:
                MW.foldernameIFFT = os.path.join(MW.batchfolder,'IFFT_spectra')
                if os.path.exists(MW.foldernameIFFT) == False:
                    os.mkdir(MW.foldernameIFFT)
            if MW.saveZspec == True:
                try:
                    MW.foldernameZ = os.path.join(MW.batchfolder,os.path.basename(MW.namebase)+'_charge_spectra')
                    os.mkdir(MW.foldernameZ)
                except FileExistsError:
                    count = 2
                    while os.path.exists(os.path.join(MW.batchfolder,os.path.basename(MW.namebase)+'_charge_spectra_'+str(count))):
                        count += 1
                    MW.foldernameZ = os.path.join(MW.batchfolder,os.path.basename(MW.namebase)+'_charge_spectra_'+str(count))
                    os.mkdir(MW.foldernameZ)
            MW.foldernamePL = os.path.join(MW.batchfolder,'Peaklists')
            if os.path.exists(MW.foldernamePL) == False:
                os.mkdir(MW.foldernamePL)
            MW.foldernameDecon = os.path.join(MW.batchfolder,'Deconvolved_spectra')
            if os.path.exists(MW.foldernameDecon) == False:
                os.mkdir(MW.foldernameDecon)
            try:
                MW.reconsub
            except AttributeError:
                pass
            else:
                MW.foldernameMMD = os.path.join(MW.batchfolder,'Mass_defect_analysis')
                if os.path.exists(MW.foldernameMMD) == False:
                    os.mkdir(MW.foldernameMMD)
            MW.foldername = MW.batchfolder
        
        if MW.saveMS == True:
            try:
                label = os.path.basename(MW.namebase) + "_MS.csv"
                temp_path = os.path.join(MW.foldername, label)
                np.savetxt(temp_path, np.c_[np.round(MW.x, 4), np.round(MW.y, 4)], delimiter=',')
                print('Saved MS data')
            except AttributeError:
                pass
        if MW.saveMSRT == True:
            try:  
                label = os.path.basename(MW.namebase) + "_MS.csv"
                try:
                    temp_path = os.path.join(MW.foldernameMS, label)
                except AttributeError:
                    temp_path = os.path.join(MW.foldername, label)
                np.savetxt(temp_path, np.c_[np.round(MW.x, 4), np.round(MW.y, 4)], delimiter=',')
                print('Saved extracted MS data')
            except AttributeError:
                pass
        if MW.saveFT == True:
            try:  
                label = os.path.basename(MW.namebase) + "_FT.csv"
                try:
                    temp_path = os.path.join(MW.foldernameFT, label)
                except AttributeError:
                    temp_path = os.path.join(MW.foldername, label)
                np.savetxt(temp_path, np.c_[np.round(MW.xFT, 4), np.round(abs(MW.yFT), 4)], delimiter=',')
                print('Saved FT data')
            except AttributeError:
                print('No FT data to save')
        try:
            if MW.xFT[-1] > 100:
                dec = 4
            else:
                dec = 3
        except AttributeError:
            dec = 3
        
        try:
            self.label = os.path.basename(MW.namebase) + "_peaklist.csv"
            try:    
                self.temp_path = os.path.join(MW.foldernamePL, self.label)
            except AttributeError:
                self.temp_path = os.path.join(MW.foldername, self.label)
            xcent = np.concatenate((['Centroid'],np.round(MW.xcent,dec),['Noise:']))
            sumlist = np.concatenate((['Integration'],np.round(MW.sumlist, 2),[np.round(MW.N,2)]))
            height = np.concatenate((['Height'],np.round(MW.height, 2),[np.round(MW.Nh,2)]))
            peakwidths = np.concatenate((['Width'],np.round(MW.peakwidths,dec),[np.round(MW.integratedx,dec)]))
            np.savetxt(self.temp_path, np.c_[xcent, sumlist, height, peakwidths], delimiter=',',fmt='%s')
            print('Saved peaklist(s)')
        except AttributeError:
            print('No peak list to save')
        except FileNotFoundError:
            pass
            
        try:
            MW.batchlabels.append(self.label)
            MW.batchpeaklists.append(self.temp_path)
        except AttributeError:
            print('No batch started')
            
        if MW.saveIFFT == True:
            try:
                ### to save the inverse FTs of each charge state:
                ifft1 = np.concatenate((['m/z'],np.round(MW.xint, 4)))
                ifft2 = np.concatenate((['Original MS'],np.round(MW.yint,4)))
                ifft3 = np.concatenate((['Charge States:'],np.zeros(len(MW.xint))))
                if len(MW.batchstore) > 0:
                    for i in range(len(MW.batchstore)):
                        label = os.path.basename(MW.namebase) + "_iffts_" + str(np.round(MW.batchstore[i][3][np.argmax(MW.batchstore[i][4])]/1000,1)) + "kDa_series.csv"
                        for j in range(len(MW.batchstore[i][8])):
                            iffttemp = np.concatenate(([MW.batchstore[i][0][j]],MW.batchstore[i][8][j]))
                            if j == 0:
                                ifft4 = iffttemp
                            else:
                                ifft4 = np.c_[ifft4,iffttemp]
                        try:
                            temp_path = os.path.join(MW.foldernameIFFT, label)
                        except AttributeError:
                            temp_path = os.path.join(MW.foldername, label)
                        np.savetxt(temp_path, np.c_[ifft1,ifft2,ifft3,ifft4], delimiter=',',fmt='%s')
                else:
                    label = os.path.basename(MW.namebase) + "_iffts.csv"
                    for i in range(len(MW.cs)):
                        iffttemp = np.concatenate(([np.real(MW.cs[i])],np.real(MW.glist[i])))
                        if i == 0:
                            ifft4 = iffttemp
                        else:
                            ifft4 = np.c_[ifft4,iffttemp]
                    try:
                        temp_path = os.path.join(MW.foldernameIFFT, label)
                    except AttributeError:
                        temp_path = os.path.join(MW.foldername, label)
                    np.savetxt(temp_path, np.c_[ifft1,ifft2,ifft3,ifft4], delimiter=',',fmt='%s')
                print('Saved IFFT spectra')
            except AttributeError:
                print('Error with saving IFFT spectra')
            except IndexError:
                print('Error with saving IFFT spectra')
            
        if MW.saveZspec == True:
            try:
                ### to save individual charge-state zero-charge spectra
                for i in range(len(MW.cs)):
                    label2 = os.path.basename(MW.namebase) +'_decon_selection' +str(i+1)+'_charge'+ str(int(MW.cs[i])) + '.csv'
                    try:
                        temp_path2 = os.path.join(MW.foldernameZ, label2)
                    except AttributeError:
                        temp_path2 = os.path.join(MW.foldername, label2)
                    if MW.based == True:
                        Bounds = [1]
                    if MW.based == False:
                        Bounds = [0]
                    for j in range(0,len(MW.xrefcs[i])):
                        Bounds.append(MW.xrefcs[i][j][0])
                        Bounds.append(MW.xrefcs[i][j][-1])
                    for j in range(int(len(MW.deconcsx[i])-len(Bounds))):
                        Bounds.append(0)
                    np.savetxt(temp_path2, np.c_[np.round(MW.deconcsx[i], 4), np.round(MW.deconcsy[i], 4),np.round(Bounds,4)], delimiter=',')
                    
                    xcentcs = np.concatenate((['Centroid'],np.round(MW.xcentcs[i],dec)))
                    sumlistcs = np.concatenate((['Integration'],np.round(MW.sumlistcs[i], 2)))
                    heightcs = np.concatenate((['Height'],np.round(MW.heightcs[i], 2)))
                    label3 = os.path.basename(MW.namebase) +'_peaklist_selection'+str(i+1)+'_charge' + str(int(MW.cs[i])) + '_mass'+str(np.round((MW.xcentcs[i][np.argmax(MW.sumlistcs[i])])/1000,1))+'kDa.csv'
                    try:
                        temp_path3 = os.path.join(MW.foldernameZ, label3)
                    except AttributeError:
                        temp_path3 = os.path.join(MW.foldername, label3)
                    np.savetxt(temp_path3, np.c_[xcentcs, sumlistcs, heightcs], delimiter=',',fmt='%s')
                print('Saved individual charge state reconstructions')   
            except AttributeError:
                print('Error with saving individual charge state reconstructions')
            except IndexError:
                print('Error with saving individual charge state reconstructions')

        try:
            try:
                if (len(MW.batchstore) > 1 or (len(MW.batchstore) > 0 and MW.smoothed == True)) and MW.saveSeries == True:
                    ### to save each series' full zero-charge spectrum
                    for i in range(len(MW.deconylist)):
                        mass = np.round(MW.deconxlist[i][np.argmax(MW.deconylist[i])]/1000,1)
                        label03 = os.path.basename(MW.namebase) + '_decon_series' +str(i+1)+'_'+str(mass)+'kDa.csv'
                        try:
                            temp_path3 = os.path.join(MW.foldernameDecon, label03)
                        except AttributeError:
                            temp_path3 = os.path.join(MW.foldername, label03)
                        try:
                            if MW.based == True:
                                MW.Bounds = [1]
                            if MW.based == False:
                                MW.Bounds = [0]
                            for j in range(0,len(MW.xref)):
                                MW.Bounds.append(MW.xref[j][0])
                                MW.Bounds.append(MW.xref[j][-1])
                            for j in range(int(len(MW.deconxlist[i])-len(MW.Bounds))):
                                MW.Bounds.append(0)
                                
                            np.savetxt(temp_path3, np.c_[np.round(MW.deconxlist[i], 4), np.round(MW.deconylist[i], 4), np.round(MW.Bounds, 4)], delimiter=',')
                        except AttributeError:
                            np.savetxt(temp_path3, np.c_[np.round(MW.deconxlist[i], 4), np.round(MW.deconylist[i], 4)], delimiter=',')
            except AttributeError:
                pass
            label3 = os.path.basename(MW.namebase) + '_decon_combined.csv'
            try:
                temp_path3 = os.path.join(MW.foldernameDecon, label3)
            except AttributeError:
                temp_path3 = os.path.join(MW.foldername, label3)
            
            try:
                if MW.based == True:
                    MW.Bounds = [1]
                if MW.based == False:
                    MW.Bounds = [0]
                for i in range(0,len(MW.xref)):
                    MW.Bounds.append(MW.xref[i][0])
                    MW.Bounds.append(MW.xref[i][-1])
                for i in range(int(len(MW.xzero)-len(MW.Bounds))):
                    MW.Bounds.append(0)
                decon1 = np.concatenate((['Mass (Da)'],np.round(MW.xzero, 4)))
                decon2 = np.concatenate((['Deconvolved Abundance'],np.round(MW.zerofull, 4)))
                decon3 = np.concatenate((['Integration Bounds'],np.round(MW.Bounds, 4)))
                if len(MW.batchstore) > 0:
                    for i in range(len(MW.batchstore)):
                        for j in range(len(MW.batchstore[i][0])):
                            cszerotemp = np.interp(MW.xzero,MW.batchstore[i][3],MW.batchstore[i][5][j])
                            decontemp = np.concatenate(([MW.batchstore[i][0][j]],np.round(cszerotemp, 4)))
                            if i == 0 and j == 0:
                                deconZs = decontemp
                            else:
                                deconZs = np.c_[deconZs,decontemp]
                else:
                    for i in range(len(MW.cs)):
                        decontemp = np.concatenate(([MW.cs[i]],np.round(MW.cszero[i], 4)))
                        if i == 0:
                            deconZs = decontemp
                        else:
                            deconZs = np.c_[deconZs,decontemp]
                np.savetxt(temp_path3, np.c_[decon1, decon2, decon3, deconZs], delimiter=',',fmt='%s')
            except AttributeError:
                np.savetxt(temp_path3, np.c_[np.round(MW.xzero, 4), np.round(MW.zerofull, 4)], delimiter=',')
            except UnboundLocalError:
                np.savetxt(temp_path3, np.c_[np.round(MW.xzero, 4), np.round(MW.zerofull, 4)], delimiter=',')
            except ValueError:
                np.savetxt(temp_path3, np.c_[np.round(MW.xzero, 4), np.round(MW.zerofull, 4)], delimiter=',')
            
            print('Saved zero charge data')
            
        except AttributeError:
            print('No zero charge data to save')
   
        try:
            try:
                MW.foldernameZMMD = os.path.join(MW.foldernameMMD,os.path.basename(MW.namebase)+'_MMD_charge_spectra')
                os.mkdir(MW.foldernameZMMD)
                for i in range(len(MW.cs)):
                    label = os.path.basename(MW.namebase) + "_MMD" + str(int(MW.cs[i])) + ".csv"
                    temp_path = os.path.join(MW.foldernameZMMD, label)
                    np.savetxt(temp_path, np.c_[np.round(np.real(MW.reconmzsub[i]), 4), np.round(np.real(MW.reconsub[i]), 4)], delimiter=',')
            except AttributeError:
                for i in range(len(MW.cs)):
                    label = os.path.basename(MW.namebase) + "_MMD" + str(int(MW.cs[i])) + ".csv"
                    temp_path = os.path.join(MW.foldername, label)
                    np.savetxt(temp_path, np.c_[np.round(np.real(MW.reconmzsub[i]), 4), np.round(np.real(MW.reconsub[i]), 4)], delimiter=',')
            label2 = os.path.basename(MW.namebase) + '_MMDfullreal.csv' 
            try:
                temp_path2 = os.path.join(MW.foldernameMMD, label2)
            except AttributeError:
                temp_path2 = os.path.join(MW.foldername, label2)
            np.savetxt(temp_path2, np.c_[np.round(np.real(MW.reconmzint), 4), np.round(np.real(MW.recontot), 4)], delimiter=',')
            print('Saved reconstructions and mass defect data')
            print('Finished \n')
        except AttributeError:
            print('No mass defect reconstructions to save')
            print('Finished \n')
    
    def load_cali(self):
        """
        loads the calibration curve and plots it. 
        """
        try:
            MW.Calfiles = []
            MW.conc,MW.concydata,MW.namebase,MW.SN,MW.sigtype,MW.curtype,MW.Calunits,MW.calibrants_used,MW.std_used,self.tol,MW.weight,self.Avg = cal.load_file_cal(self)
            MW.cal_dir = MW.namebase
            self.layout = QtWidgets.QGridLayout() 
            
            try:
                del(Unknown_calculator.unkconc)
            except AttributeError:
                pass
            
            MW.standards = len(MW.std_used)
            MW.tolerance = float(self.tol[2])
            MW.tunit = self.tol[1]
            
            print('Calibration Parameters: ')
            for i in range(len(MW.calibrants_used)):
                  print(' Calibrant Mass '+str(i+1)+' = '+str(MW.calibrants_used[i])+' Da')
            for i in range(len(MW.std_used)):
                  print(' Standard Mass '+str(i+1)+' = '+str(MW.std_used[i])+' Da')
            print(' Tolerance = '+str(MW.tolerance)+' '+MW.tunit)            
            print(' Signal Type = '+MW.sigtype)
            print(' Curve Type = '+MW.curtype)
            print(' Weighting = '+MW.weight)
            
            if self.Avg == True:
                conc = [MW.conc[0]]
                concytemp = [MW.concydata[0]]
                concy = []
                xrep = []
                for i in range(1,len(MW.conc)):
                    if MW.conc[i] == conc[-1]:
                        concytemp.append(MW.concydata[i])
                    else:
                        conc.append(MW.conc[i])
                        concy.append(concytemp)
                        xrep.append(len(concytemp))
                        concytemp = [MW.concydata[i]]
                concy.append(concytemp)
                xrep.append(len(concytemp))
                MW.concAvg = conc
                MW.concydataAvg = np.zeros(len(concy))
                ystddev = np.zeros(len(concy))
                for i in range(len(concy)):
                    MW.concydataAvg[i] = np.mean(concy[i])
                    ystddev[i] = np.std(concy[i])
                MW.concydataAvg = MW.concydataAvg.tolist()
            else:
                MW.concAvg = MW.conc
                MW.concydataAvg = MW.concydata
                ystddev = np.zeros(len(MW.concAvg))
                xrep = 1+np.zeros(len(MW.concAvg))
              
            MW.calx,MW.caly,MW.calparam,MW.calinfo,MW.std_err,MW.confidence,MW.uinfo= cal.calibrate(MW.concAvg,xrep,MW.concydataAvg,MW.weight,MW.curtype,ystddev)
                
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(1,4,(1,3))
            MW.ax.fill_between(MW.confidence[0],MW.confidence[1],MW.confidence[2],color = color[3],alpha=0.3)
            print (MW.calinfo)
            if MW.conc == MW.concAvg:
                if MW.weight == 'no weighting':
                    bars = MW.std_err
                if MW.weight == '1/x':
                    bars = 1/np.array(MW.conc)*MW.std_err
                if MW.weight == '1/x^2':
                    bars = 1/np.array(MW.conc)**2*MW.std_err
                if MW.weight == 'x':
                    bars = MW.std_err*np.array(MW.conc)
                if MW.weight == 'x^2':
                    bars = MW.std_err*np.array(MW.conc)**2
    
                MW.ax.errorbar(MW.conc, MW.concydata, yerr = bars,linestyle='None',marker='o',label='calibrant')
                print(' Standard Deviation About the Regression = ' +str(np.round(MW.std_err,4))+'\n'
                      ' Error bars indicate one standard deviation.')
            else:
                MW.ax.errorbar(MW.conc, MW.concydata, yerr = 0,linestyle='None',marker='o',label='calibrant')
                MW.ax.errorbar(MW.concAvg, MW.concydataAvg, yerr = ystddev,linestyle='None',marker='o',color='r',label='averaged')
                print(' Error bars on averaged data indicate one standard deviation at that concentration.')
                print(' Standard Deviation About the Regression = ' +str(np.round(MW.std_err,4)))
            MW.ax.plot(MW.calx, MW.caly)
            print(' Band indicates 95% confidence interval.')
            MW.ax.set_xlabel('Concentration ('+MW.Calunits+')')
            if MW.sigtype == 'integration':
                siglabel = 'Peak-area Relative Abundance'
            else:
                siglabel = 'Peak-height Relative Abundance'
            MW.ax.set_ylabel(siglabel)
            try:
                err = max(MW.std_err)
            except TypeError:
                err = MW.std_err
            MW.ax.set_ylim(min(MW.concydata)-3*err,max(MW.concydata)+3*err)
            MW.ax.set_xlim(min(MW.conc)-max(MW.conc)/100,max(MW.conc)+max(MW.conc)/100)
            MW.ax.set_title(MW.calinfo)
    
            MW.axN = MW.fig.add_subplot(1,4,4)
            columnsN = ('Calibrant \nConcentration ('+MW.Calunits+')','Response')
            MW.axN.axis('off')
            response = []
            for i in range(len(MW.concydata)):
                response.append(str(np.format_float_scientific(MW.concydata[i],precision=2,unique=False)))
            MW.Ntable = MW.axN.table(cellText=np.c_[np.concatenate((MW.conc,[str('Regression SD')])),np.concatenate((response,[np.format_float_scientific(err,precision=2,unique=False)]))],colLabels=columnsN,loc='center')
            MW.Ntable.auto_set_font_size(False)
            MW.Ntable.set_fontsize(12)
            MW.axN.axis('tight')
            MW.Ntable.scale(1.0,2.2)
    
            MW.file_display.setText('File: '+MW.namebase) 
            MW.canvas.draw()

        except TypeError:
            print('Either no data loaded or wrong data type.  Try .csv for wrong data type')
            MW.message = str('Either no data loaded or wrong data type. \n'
                                  'Please load an iFAMS calibration file.')
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
 
    
    def save_cali(self):
        try:
            print('Saving calibration curve')
            length = max(7,len(MW.calibrants_used),len(MW.std_used))
            caltol = np.zeros(length+1).tolist()
            calmasses = np.zeros(length+1).tolist()
            stdmasses = np.zeros(length+1).tolist()
            
            subheader = []
            subheader.append(MW.Calunits)
            subheader.append(MW.sigtype)
            subheader.append('-')
            
            caltol[0] = 'Tolerance:'
            caltol[1] = MW.tolerance
            caltol[2] = MW.tunit
            caltol[3] = 'Curve Type:'
            caltol[4] = MW.curtype
            caltol[5] = 'Weighting:'
            caltol[6] = MW.weight
            if MW.conc == MW.concAvg:
                caltol[7] = 0
            else:
                caltol[7] = 'Fit averaged data'
            
            calmasses[0] = 'Calibrant Masses:'
            for i in range(0,len(MW.calibrants_used)):
                calmasses[i+1] = MW.calibrants_used[i]
            stdmasses[0] = 'Standard Masses:'
            for i in range(0,len(MW.std_used)):
                stdmasses[i+1] = MW.std_used[i]
            
            calX = np.concatenate((['Concentration',subheader[0]],np.round(MW.conc, 4),calmasses))
            calY = np.concatenate((['Abundance',subheader[1]],np.round(MW.concydata, 7),stdmasses))
            calSN = np.concatenate((['S:N',subheader[2]],MW.SN,caltol))
            
            temp_dir = MW.cal_dir
            
            if MW.weight == 'no weighting':
                callab = 'calcurve_'+str(MW.curtype)+'_with_no_weighting'
            else:
                callab = 'calcurve_'+str(MW.curtype)+'_with_weighting'
            
            if MW.OverwriteCurve == True:
                try:
                    os.remove(os.path.join(temp_dir,MW.prevCurve))
                except AttributeError:
                    pass
            count = 1
            calcurvelabel = str(callab)+'_'+str(count)+'.csv'
            test_path = os.path.exists(os.path.join(temp_dir,calcurvelabel))
            while test_path == True:
                count += 1
                calcurvelabel = str(callab)+'_'+str(count)+'.csv'
                test_path = os.path.exists(os.path.join(temp_dir,calcurvelabel))
            MW.prevCurve = calcurvelabel
                    
            temp_path_calcurve = os.path.join(temp_dir, calcurvelabel)
            np.savetxt(temp_path_calcurve, np.c_[calX,calY,calSN], delimiter=',',fmt='%s')
            print('Saved calibration curve data as ' + str(temp_path_calcurve))
        
            try:
                PeakTable = []
                column1 = ['Peak Type:']
                column2 = ['Peak Masses:']
                for i in range(0,len(MW.calibrants_used)):
                    column1.append('Calibrant '+str(i+1))
                    column2.append(str(MW.calibrants_used[i]))
                for i in range(0,len(MW.std_used)):
                    column1.append('Int STD '+str(i+1))
                    column2.append(str(MW.std_used[i]))
                PeakTable.append(column1)
                PeakTable.append(column2)
                for i in range(0,len(MW.xcent1)):
                    fileColumn1 = [MW.fileconc[MW.concIndex[i]]+MW.Calunits+' Conc. Centroid']
                    fileColumn2 = [MW.fileconc[MW.concIndex[i]]+' Mass Error ('+MW.tunit+')']
                    fileColumn3 = [MW.fileconc[MW.concIndex[i]]+' '+subheader[1]]
                    error1 = []
                    for j in range(0,len(MW.xcent1[i])):
                        if MW.tunit == 'Da':
                            error1.append(np.round(float(MW.xcent1[i][j]-MW.calibrants_used[j]),2))
                        else:
                            error1.append(np.round(float((MW.xcent1[i][j]-MW.calibrants_used[j])/MW.calibrants_used[j]*1000000),2))
                    for j in range(0,len(MW.std_used)):
                        if MW.tunit == 'Da':
                            error1.append(np.round(float(MW.xcent2[i][j]-MW.std_used[j]),2))
                        else:
                            error1.append(np.round(float((MW.xcent2[i][j]-MW.std_used[j])/MW.std_used[j]*1000000),2))
                    if MW.standards > 0 and len(MW.std_used) > 0:
                        PeakTable.append(np.concatenate((fileColumn1,MW.xcent1[i],MW.xcent2[i])).tolist())
                        PeakTable.append(np.concatenate((fileColumn2,error1)).tolist())
                        PeakTable.append(np.concatenate((fileColumn3,MW.ylist1[i],MW.ylist2[i])).tolist())
                    else:
                        PeakTable.append(np.concatenate((fileColumn1,MW.xcent1[i])).tolist())
                        PeakTable.append(np.concatenate((fileColumn2,error1)).tolist())
                        PeakTable.append(np.concatenate((fileColumn3,MW.ylist1[i])).tolist())
                temp_path = os.path.join(temp_dir, 'peak_table_for_'+calcurvelabel)
                np.savetxt(temp_path, np.transpose(PeakTable), delimiter=',',fmt='%s')
                print('Saved Peak Table')
            except AttributeError:
                print('Could not save Peak Table')
    
        except AttributeError:
            print('Could not save')
            MW.message = str('No calibration to save')
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
            
    def deconbatch(self):
        """
        batches deconvolved spectra, returning the peak list for all files selected
        """
        MW.batchpeaklists = []
        MW.batchlabels = []
        
        day = str(datetime.datetime.now().day)
        if int(day) < 10:
            day = str(0)+day
        month = str(datetime.datetime.now().month)
        if int(month) < 10:
            month = str(0)+month
        year = str(int(float(datetime.datetime.now().year)-2000))
        count = 1
        folder = str(year+month+day+'_Batch'+str(count))
        test_path = os.path.exists(os.path.join(os.path.dirname(MW.namebase),folder))
        while test_path == True:
            count += 1
            folder = str(year+month+day+'_Batch'+str(count))
            test_path = os.path.exists(os.path.join(os.path.dirname(MW.namebase),folder))
        MW.batchfolder = os.path.join(os.path.dirname(MW.namebase),folder)
        os.mkdir(MW.batchfolder)
        
        self.save_FT()
        
        try:
            name = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select Files to Batch')
            for n in range(0, len(name[0])):
                MW.xzero, MW.zerofull, MW.namebase,MW.namestr = load.batch_load(name[0][n])
            
                MW.basecor = True
                print('Starting analysis of ' + str(name[0][n]))

                if MW.smoothed == True:
                    x = np.round(np.linspace(min(MW.xzero),max(MW.xzero),len(MW.xzero)),4)
                    FFT = np.fft.fft(MW.zerofull)
                    width = 1.4
                    sig = width/(2*np.sqrt(2*np.log(2)))
                    gaussian = np.exp(-np.power(x-x[0],2) / (2*np.power(sig,2))) + np.exp(-np.power(x-x[-1]-((max(x)-min(x))/len(x)),2) / (2*np.power(sig,2)))
                    FFT = FFT*np.fft.fft(gaussian)
                    MW.zerofull = np.real(np.fft.ifft(FFT))
    
                try:
                    if MW.paramMatch == True:
                        print ('Starting peak integration')
                        MW.ynorm = MW.zerofull / max(MW.zerofull)
                        MW.xlist,MW.ylist,MW.xref,MW.yref,MW.sumlist,MW.peakindex,MW.minL,MW.minR,MW.xcent,MW.height = Int.integrate(MW.ynorm, MW.zerofull,
                                                                                          MW.xzero, MW.deltaY,
                                                                                          MW.minimumY, MW.minimumX, MW.maximumX,MW.ntol)
                        max_value = max(MW.zerofull)
                        for i in range(0, len(MW.sumlist)):
                            MW.sumlist[i] = max_value * MW.sumlist[i]
                            MW.height[i] = max_value*MW.height[i]
                            
                    if MW.paramMatch == False:
                        MW.xref, MW.yref, MW.sumlist, MW.height, MW.xcent, MW.peakindex = Int.boundsmatch(MW.xzero,MW.zerofull,MW.bounds)
                        
                    if MW.based == True:
                        MW.yref2, MW.sumlist, MW.height = Int.baseline_sub2(MW.xzero, MW.zerofull, MW.xref, 
                                                                 MW.yref, MW.peakindex)
                        
                    MW.integratedx = 0 
                    if MW.paramMatch == False:
                        MW.peakwidths = np.zeros(len(MW.sumlist))
                        for i in range(0, len(MW.sumlist)):
                            MW.integratedx += (MW.bounds[2*i+1]-MW.bounds[2*i])*(MW.xzero[1]-MW.xzero[0])
                            MW.peakwidths[i] = float((MW.bounds[2*i+1]-MW.bounds[2*i])*(MW.xzero[1]-MW.xzero[0]))
                    else:
                        MW.peakwidths = np.zeros(len(MW.xref))
                        for i in range(0, len(MW.xref)):
                            MW.integratedx += len(MW.xref[i])
                            MW.peakwidths[i] = len(MW.xref[i])

                except AttributeError:
                    print('No peak integration')
    
                self.save_FT()   
            print('Batch complete')
            MW.message = 'Batch Complete'
            MW.Err = False
            self.dialog = iFAMS_message(self)
            self.dialog.show()

        except TypeError:
            MW.ref_num = 0
            print('Either no data loaded or wrong data type.  Try .csv or .txt for wrong data type')
           
    def batch(self):
        try:
            if MW.namebase == 'From Clipboard':
                MW.full_save = True
                MW.saveMS = True
                MW.ReturnFunc = MW.batch
                Folder_Select(self).show()
                MW.message = 'Please save current data before batching.'
                MW.submessage = 'The batch folder will be generated in the same directory.'
                MW.Err = False
                self.dialog = iFAMS_message(self)
                self.dialog.show()
                return
        except AttributeError:
            pass
        MW.batchpeaklists = []
        MW.batchlabels = []
        
        self.dialog = FileDialog(self)
        self.dialog.show()
        
    def batch_function(self): 
        print('Starting batch')
        try:
            MW.batchfolder
        except AttributeError:
            if MW.batchcheck == False:
                MW.batch_param(self) 
                MW.batch_save(self)
                MW.save_FT(self)
            else:
                MW.namebase = MW.Batchname[0]
                MW.batch_param(self) 
                MW.batch_save(self)
        else:
            MW.save_FT(self)
        
        MW.fig.tight_layout()
        MW.canvas.draw()
        
        MW.batchstore = []
        MW.saveMS = False
        for n in range(0, len(MW.Batchname)):
            try:
                try:
                    del(MW.reconsub)
                    del(MW.foldernameMMD)
                except AttributeError:
                    pass
                print('Starting FT analysis for ' + str(MW.Batchname[n]))
                print(str(datetime.datetime.now()))
                MW.tIFT = []
                MW.namestr = str(MW.Batchname[n])
                MW.namebase = os.path.splitext(MW.namestr)[0]
                MW.x, MW.y, MW.namebase,MW.namestr = load.batch_load(MW.Batchname[n])
                MW.saveMSRT = False
                if MW.y[0] == 0 and MW.y[-1]==0:
                    MW.x = MW.x[1:-1].copy()
                    MW.y = MW.y[1:-1].copy()
                MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
                if MW.setdatarho == True:
                    xtemp = np.linspace(min(MW.x),max(MW.x),MW.datarho)
                    ytemp = np.interp(xtemp,MW.x,MW.y)
                    MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(xtemp,ytemp)
                if MW.domaintrunc == True:
                    xstore = MW.x.copy()
                    ystore = MW.y.copy()
                    MW.x = []
                    MW.y = []
                    for i in range(len(MW.xint)):
                        if (MW.domainminx <= xstore[i]) and (xstore[i] <= MW.domainmaxx):
                            MW.x.append(xstore[i])
                            MW.y.append(ystore[i])
                    MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
                    
                MW.xFT, MW.yFT = FT.Fourier(MW.xint, MW.ypadd, MW.po2)
            except AttributeError:
                print('Unable to perform initial data processing')

            
            try:    #Initial noise determination and MS baseline correction
                MW.maxF = MW.xFT[-1]
                maxF = int(MW.xFT[-1]/2)
                MW.NFmin2 = MW.xFT[int((int(maxF*4/5)+0.4)*len(MW.yFT)/MW.maxF)]
                MW.NFmax2 = MW.xFT[int((int(maxF*4/5)+0.6)*len(MW.yFT)/MW.maxF)]
                MW.NFrms2 = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin2,MW.NFmax2),4)
                MW.NFmin3 = MW.xFT[int((int(maxF*3/5)+0.4)*len(MW.yFT)/MW.maxF)]
                MW.NFmax3 = MW.xFT[int((int(maxF*3/5)+0.6)*len(MW.yFT)/MW.maxF)]
                MW.NFrms3 = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin3,MW.NFmax3),4)
                MW.NFmin4 = MW.xFT[int((int(maxF*2/5)+0.4)*len(MW.yFT)/MW.maxF)]
                MW.NFmax4 = MW.xFT[int((int(maxF*2/5)+0.6)*len(MW.yFT)/MW.maxF)]
                MW.NFrms4 = np.round(cal.noise_calc(MW.xFT,MW.yFT,MW.NFmin4,MW.NFmax4),4)
                MW.NFrms = np.average([MW.NFrms2,MW.NFrms3,MW.NFrms4])
                MW.Nrms = np.round(MW.NFrms/np.sqrt(len(MW.x)),4)
                if MW.Nlinear_sub == True:
                    MW.globalmin1 = min(MW.y[0:int(len(MW.y)/20)])
                    MW.globalmin2 = min(MW.y[-int(len(MW.y)/20):-1])
                    dy = MW.globalmin2-MW.globalmin1
                    dx = MW.x[-1]-MW.x[0]
                    base_m = dy/dx
                    base_b = MW.globalmin1-base_m*MW.x[0]
                    for i in range(0,len(MW.y)):
                        MW.y[i] = MW.y[i] - (base_m*MW.x[i]+base_b)
                    MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
                    MW.xFT, MW.yFT = FT.Fourier(MW.xint, MW.ypadd, MW.po2)
                    
                if MW.Nseg_sub == True:
                    MW.mzminy = []
                    MW.mzminx = []
                    MW.mzslopes = []
                    MW.mzintercepts = []
                    
                    MW.mzminy.append(min(MW.yint[0:int(len(MW.yint)/50)]))
                    MW.mzminx.append(MW.xint[0])
                    for i in range(1,30):
                        MW.mzminy.append(min(MW.yint[int(len(MW.yint)*i/30-len(MW.yint)/100):int(len(MW.yint)*i/30+len(MW.yint)/100)]))
                        MW.mzminx.append(MW.xint[int(len(MW.yint)*i/30)])
                    MW.mzminy.append(min(MW.yint[-int(len(MW.yint)/50):-1]))
                    MW.mzminx.append(MW.xint[-1])
                    for i in range(0,30):
                        dy = MW.mzminy[i+1]-MW.mzminy[i]
                        dx = MW.mzminx[i+1]-MW.mzminx[i]
                        MW.mzslopes.append(float(dy/dx))
                        MW.mzintercepts.append(float(MW.mzminy[i]-(dy/dx)*MW.mzminx[i]))
                    j = 0
                    for i in range(0,len(MW.yint)):
                        if MW.xint[i] > MW.mzminx[j+1]:
                            j += 1
                            if j == 30:
                                j -= 1
                        MW.yint[i] = MW.yint[i] - (MW.mzslopes[j]*MW.xint[i]+MW.mzintercepts[j])
                    
                    MW.y = MW.yint
                    MW.x = MW.xint
                    
                    MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
                    MW.xFT, MW.yFT = FT.Fourier(MW.xint, MW.ypadd, MW.po2)
                    
                if MW.Nsubtraction == True:
                    MW.y = MW.y - MW.RMSDmult*MW.Nrms
                    MW.xint,MW.yint,MW.ypadd,MW.po2,MW.xpadd = load.datapro(MW.x,MW.y)
                    MW.xFT, MW.yFT = FT.Fourier(MW.xint, MW.ypadd, MW.po2)

            except AttributeError:
                print('Unable to determine initial noise')
                MW.NFrms = 0
                MW.Nrms = 0
                MW.Nh = 0
                MW.N = 0
                MW.GaborFraction = 0
                

            MW.X = STFT.stft(MW.ypadd, MW.winnum, MW.window,MW.OF)
            MW.mzspacing = float((MW.xint[-1] - MW.xint[0])/len(MW.xint))
            MW.maxF = float(MW.xFT[-1])
            MW.xshort = np.linspace((MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing + min(MW.xpadd), max(MW.xpadd) - (MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing, len(MW.X))
            MW.yshort = np.linspace(0 - (max(MW.xFT)/len(MW.X[0]))/2, max(MW.xFT) + (max(MW.xFT)/len(MW.X[0]))/2, len(MW.X[0]))
            MW.mzchan = float((MW.xint[-1]-MW.xshort[0])/((MW.xshort[-1]-MW.xshort[0])/len(MW.xshort)))
            ytemp = STFT.re_stft(MW.X, MW.winnum, MW.ypadd, MW.yint,MW.OF)
            MW.correction_factor = max(MW.yint)/max(ytemp)
            print('Starting STFT analysis')
            
            #for baseline model
            MW.yBL2 = np.zeros(len(MW.xint))+ min(MW.yint)+ 1
            MW.glistBL2 = []
            MW.xintBL2,MW.yintBL2,MW.ypaddBL2,MW.po2BL2,MW.xpaddBL2 = load.datapro(MW.xint,MW.yBL2)
            MW.XBL2 = STFT.stft(MW.ypaddBL2, MW.winnum, MW.window, MW.OF)
            MW.XtempBL2 = MW.XBL2*0
            
            MW.Ilist = []
            piW = MW.xshort[2]-MW.xshort[1]
            piH = MW.yshort[2]-MW.yshort[1]
            for j in range(len(MW.rlist)):
                lowmzI = int(np.round((MW.rlist[j][0]-MW.xshort[0])/piW))
                himzI = int(np.round((MW.rlist[j][1]-MW.xshort[0])/piW))
                lowFI = int(np.round(MW.rlist[j][2]/piH))
                hiFI = int(np.round(MW.rlist[j][3]/piH))
                MW.Ilist.append([lowmzI,himzI,lowFI,hiFI,MW.rlist[j][4]])
          
            MW.glist = []
            bstemp = []
            MW.Xtemp = MW.X*0
            Icount = 0
            for j in range(len(MW.Ilist)):
                if MW.Ilist[j][4] > Icount or (MW.Ilist[j][4] < Icount and MW.Ilist[j][4] == 0):
                    MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                    MW.glist.append(MW.yrecon)
                    MW.Xtemp = MW.X*0
                    MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                    MW.glistBL2.append(MW.yreconBL2)
                    MW.XtempBL2 = MW.XBL2*0
                    Icount += 1
                    if MW.Ilist[j][4] == 0:
                        Icount = 0
                        bstemp.append(len(MW.glist))
                if MW.Ilist[j][2] <= 1:
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1]
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:len(MW.X[0])] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:len(MW.X[0])]
                    MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1]
                    MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:len(MW.X[0])] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:len(MW.X[0])]
                else:
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1]
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1]
                    MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1]
                    MW.XtempBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1] = MW.XBL2[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1]
            MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
            MW.glist.append(MW.yrecon)
            MW.Xtemp = MW.X*0
            MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
            MW.glistBL2.append(MW.yreconBL2)
            MW.XtempBL2 = MW.XBL2*0
            bstemp.append(len(MW.glist))
                
           # try:                                
            print('Starting zero-charge analysis')
            MW.deconylist = []
            MW.deconxlist = []
            MW.deconcsy = []
            MW.deconcsx = []
            MW.deconyBL = []
            MW.batchstore = []
            for j in range(len(bstemp)):
                batchstore = []
                if j == 0:
                    MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(np.real(MW.glist[0:bstemp[j]]), MW.xint, MW.cs[0:bstemp[j]],MW.proton)
                    MW.xzeroBL2, MW.cszeroBL2, MW.zerofullBL2 = FT.zerocharge(np.real(MW.glistBL2[0:bstemp[j]]), MW.xint, MW.cs[0:bstemp[j]],MW.proton)
                else:
                    MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(np.real(MW.glist[bstemp[j-1]:bstemp[j]]), MW.xint, MW.cs[bstemp[j-1]:bstemp[j]],MW.proton)
                    MW.xzeroBL2, MW.cszeroBL2, MW.zerofullBL2 = FT.zerocharge(np.real(MW.glistBL2[bstemp[j-1]:bstemp[j]]), MW.xint, MW.cs[bstemp[j-1]:bstemp[j]],MW.proton)
                massrange = np.nonzero(MW.zerofullBL2)
                MW.deconylist.append(MW.zerofull[massrange])
                MW.deconxlist.append(MW.xzero[massrange])
                MW.deconyBL.append(MW.zerofullBL2[massrange])
                MW.BSint = []
                for k in range(len(MW.cszero)):
                    MW.cszero[k] = MW.cszero[k][massrange]
                    MW.cszeroBL2[k] = MW.cszeroBL2[k][massrange]
                    csxtemp = MW.xzero[massrange]
                    check = np.max(MW.cszeroBL2[k])/4
                    if check == 0:
                        continue
                    mid1 = int(len(MW.xzero)/8)+1
                    mid2 = int(len(MW.xzero)*7/8)-1
                    for h in range(0,len(MW.cszeroBL2[k])):
                        if MW.cszeroBL2[k][h] > check:
                            mid1 = h
                            break
                    for h in range(mid1+1,len(MW.cszeroBL2[k])):
                        if MW.cszeroBL2[k][h] <= check or MW.cszeroBL2[k][h] <= MW.cszeroBL2[k][mid1]:
                            mid2 = h
                            break
                    int1 = min(sum(MW.cszero[k][0:mid1]),sum(MW.cszero[k][mid2:len(MW.cszeroBL2[k])]))
                    int2 = (sum(MW.cszeroBL2[k][0:mid1])+sum(MW.cszeroBL2[k][mid2:len(MW.cszeroBL2[k])]))/2
                    MW.BSint.append(int1/int2)
                    MW.deconcsy.append(np.asarray(MW.cszero[k])-np.asarray(MW.cszeroBL2[k])*int1/int2)
                    MW.deconcsx.append(csxtemp)            
                BSint = np.average(MW.BSint)
                MW.deconyBL[j] = np.asarray(MW.deconyBL[j])*BSint 
                MW.deconylist[j] = np.asarray(MW.deconylist[j])-np.asarray(MW.deconyBL[j])  
                if j == 0:
                    batchstore.append(MW.cs[0:bstemp[j]])
                    batchstore.append([0,0,0,0,0])
                    batchstore.append(np.real(MW.glist[0:bstemp[j]]))
                else:
                    batchstore.append(MW.cs[bstemp[j-1]:bstemp[j]])
                    batchstore.append([0,0,0,0,0])
                    batchstore.append(np.real(MW.glist[bstemp[j-1]:bstemp[j]]))
                batchstore.append(MW.deconxlist[j])
                batchstore.append(MW.deconylist[j])
                batchstore.append(MW.cszero)
                batchstore.append(MW.cszeroBL2)
                batchstore.append(bstemp[j])
                batchstore.append(np.real(batchstore[2]))
                MW.batchstore.append(batchstore)
                           
            if MW.NZbased == True:
                piW = MW.xshort[2]-MW.xshort[1]
                piH = MW.yshort[2]-MW.yshort[1]
                MW.NZrex = []
                nzrex = []
                MW.NZcs = []
                nzcs = []
                zcount = -1
                zcount2 = 0
                for i in range(0,len(MW.rlist)):
                    if zcount > MW.rlist[i][4]:
                        zcount2 += zcount+1
                        zcount = -1
                        MW.NZrex.append(nzrex)
                        MW.NZcs.append(nzcs)
                        nzrex = []
                        nzcs = []
                    if zcount == MW.rlist[i][4]:
                        continue
                    else:
                        nzrex.append([int(np.round((MW.rlist[i][0]-MW.xshort[0])/piW)),int(np.round((MW.rlist[i][1]-MW.xshort[0])/piW)),0, int(np.round((MW.rlist[i][3]-((MW.rlist[i][2]+MW.rlist[i][3])/2))/piH)), MW.rlist[i][4]])
                        zcount += 1
                        nzcs.append(MW.cs[zcount+zcount2])
                MW.NZrex.append(nzrex)
                MW.NZcs.append(nzcs)
                MW.Xtemp = MW.X*0
                MW.XtempBL2 = MW.XBL2*0
                MW.NZrecon = []
                MW.NZreconBL2 = []
                for i in range(len(MW.NZrex)):
                    MW.NZrecon.append([])
                    MW.NZreconBL2.append([])
                    for j in range(len(MW.NZrex[i])):
                        MW.Xtemp[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1] = MW.X[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1]
                        MW.Xtemp[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])] = MW.X[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])]
                        MW.XtempBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1] = MW.XBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,0:MW.NZrex[i][j][3]+1]
                        MW.XtempBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])] = MW.XBL2[MW.NZrex[i][j][0]:MW.NZrex[i][j][1]+1,-MW.NZrex[i][j][3]:len(MW.X[0])]
                        MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                        MW.NZrecon[i].append(MW.yrecon/2)
                        MW.Xtemp = MW.X*0
                        MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                        MW.NZreconBL2[i].append(MW.yreconBL2/2)
                        MW.XtempBL2 = MW.XBL2*0
                
                    MW.NZxzero, MW.NZcszero, MW.NZzerofull = FT.zerocharge(np.real(MW.NZrecon[i]), MW.xint, MW.NZcs[i],MW.proton)
                    MW.NZxzeroBL2, MW.NZcszeroBL2, MW.NZzerofullBL2 = FT.zerocharge(np.real(MW.NZreconBL2[i]), MW.xint, MW.NZcs[i],MW.proton)
                    massrange=np.nonzero(MW.NZzerofullBL2)
                    MW.NZxzero = MW.NZxzero[massrange]
                    MW.NZzerofull = MW.NZzerofull[massrange]
                    MW.NZxzeroBL2 = MW.NZxzeroBL2[massrange]
                    MW.NZzerofullBL2 = MW.NZzerofullBL2[massrange]
                    check = np.max(MW.NZzerofullBL2)/4
                    for j in range(len(MW.NZxzero)):
                        if MW.NZzerofullBL2[j] > check:
                            mid1 = j
                            break
                    for j in range(mid1+1,len(MW.NZxzero)):
                        if MW.NZzerofullBL2[j] <= check or MW.NZzerofullBL2[j] <= MW.NZzerofullBL2[mid1]:
                            mid2 = j
                            break
                    int1 = min(sum(MW.NZzerofull[0:mid1]),sum(MW.NZzerofull[mid2:len(MW.NZxzero)]))
                    int2 = (sum(MW.NZzerofullBL2[0:mid1])+sum(MW.NZzerofullBL2[mid2:len(MW.NZxzero)]))/2
                    MW.NZzerofullBL2 = MW.NZzerofullBL2*int1/int2
                    MW.NZzerofull = MW.NZzerofull - MW.NZzerofullBL2
                    MW.deconylist[i] = np.asarray(MW.deconylist[i]) - np.interp(MW.deconxlist[i],MW.NZxzero,MW.NZzerofull)

            min0 = []
            max0 = []
            for j in range(len(MW.deconxlist)):
                min0.append(np.min(MW.deconxlist[j]))
                max0.append(np.max(MW.deconxlist[j]))
            min0 = np.min(min0)
            max0 = np.max(max0)
            space0 = MW.xzero[2]-MW.xzero[1]
            MW.xzero = np.linspace(min0,max0,int(2*np.round((max0-min0)/space0)),endpoint=True)
            y0s = []
            for i in range(len(MW.batchstore)):
                y0s.append(np.interp(MW.xzero,MW.deconxlist[i],MW.deconylist[i]))
            MW.zerofull = np.sum(y0s,axis=0)
            
            if MW.smoothed == True:
                x = np.round(np.linspace(min(MW.xzero),max(MW.xzero),len(MW.xzero)),4)
                FFT = np.fft.fft(MW.zerofull)
                width = 1.4
                sig = width/(2*np.sqrt(2*np.log(2)))
                gaussian = np.exp(-np.power(x-x[0],2) / (2*np.power(sig,2))) + np.exp(-np.power(x-x[-1]-((max(x)-min(x))/len(x)),2) / (2*np.power(sig,2)))
                FFTgauss = np.fft.fft(gaussian)
                FFT = FFT*FFTgauss
                area = np.sum(gaussian)
                MW.zerofull = np.real(np.fft.ifft(FFT))/area
                
            if MW.basecor == True:
                  MW.baseline = STFT.line_baseline(MW.zerofull, MW.xzero, 5)
                  MW.zerofull = MW.zerofull - MW.baseline

            try:
                MW.paramMatch
                print ('Starting integration')
                massrange = np.nonzero(MW.zerofull)
                if len(MW.zerofull[massrange]) == 0:
                    MW.noname
                if MW.paramMatch == True:
                    MW.ynorm = MW.zerofull / max(MW.zerofull)
                    MW.xlist,MW.ylist,MW.xref,MW.yref,MW.sumlist,MW.peakindex,MW.minL,MW.minR,MW.xcent,MW.height = Int.integrate(MW.ynorm[massrange], MW.zerofull[massrange],
                                                                                      MW.xzero[massrange], MW.deltaY,
                                                                                      MW.minimumY, MW.minimumX, MW.maximumX,MW.ntol)
                    max_value = max(MW.zerofull)
                    for i in range(0, len(MW.sumlist)):
                        MW.sumlist[i] = max_value * MW.sumlist[i]
                        MW.height[i] = max_value*MW.height[i]
                    if MW.saveZspec == True:
                        MW.xrefcs = np.zeros(len(MW.cs)).tolist()
                        MW.yrefcs = np.zeros(len(MW.cs)).tolist()
                        MW.sumlistcs = np.zeros(len(MW.cs)).tolist()
                        MW.heightcs = np.zeros(len(MW.cs)).tolist()
                        MW.xcentcs = np.zeros(len(MW.cs)).tolist()
                        MW.peakindexcs = np.zeros(len(MW.cs)).tolist()
                        MW.bounds = []
                        for i in range(len(MW.xref)):
                            MW.bounds.append(MW.xref[i][0])
                            MW.bounds.append(MW.xref[i][-1])
                        for i in range(len(MW.cs)):
                            MW.xrefcs[i], MW.yrefcs[i], MW.sumlistcs[i], MW.heightcs[i], MW.xcentcs[i], MW.peakindexcs[i] = Int.boundsmatch(MW.deconcsx[i],MW.deconcsy[i],MW.bounds)
        
                if MW.paramMatch == False:
                    MW.xref, MW.yref, MW.sumlist, MW.height, MW.xcent, MW.peakindex = Int.boundsmatch(MW.xzero[massrange],MW.zerofull[massrange],MW.bounds)
                    if MW.saveZspec == True:
                        MW.xrefcs = np.zeros(len(MW.cs)).tolist()
                        MW.yrefcs = np.zeros(len(MW.cs)).tolist()
                        MW.sumlistcs = np.zeros(len(MW.cs)).tolist()
                        MW.heightcs = np.zeros(len(MW.cs)).tolist()
                        MW.xcentcs = np.zeros(len(MW.cs)).tolist()
                        MW.peakindexcs = np.zeros(len(MW.cs)).tolist()
                        for i in range(len(MW.cs)):
                            MW.xrefcs[i], MW.yrefcs[i], MW.sumlistcs[i], MW.heightcs[i], MW.xcentcs[i], MW.peakindexcs[i] = Int.boundsmatch(MW.deconcsx[i],MW.deconcsy[i],MW.bounds)
                        
                if MW.based == True:
                    MW.yref, MW.sumlist, MW.height = Int.baseline_sub2(MW.xzero, MW.zerofull, MW.xref, MW.yref, MW.peakindex)            
            
            except AttributeError:
                print('No peak integration or empty deconvolution')

            try:   
                ###Noise Determination
                massrange = np.nonzero(MW.zerofull)
                MW.GaborArea = (MW.xint[-1]-MW.xshort[0])*(MW.xFT[-1]-MW.xFT[0])
                MW.selecArea = 0
                for i in range(0,len(MW.rlist)):
                    MW.selecArea += (MW.rlist[i][1]-MW.rlist[i][0])*(MW.rlist[i][3]-MW.rlist[i][2])*2
                MW.GaborFraction = float(MW.selecArea/MW.GaborArea)
        
                if MW.paramMatch == True:
                    MW.integratedx = 0 
                    MW.peakwidths = np.zeros(len(MW.xref))
                    for i in range(0, len(MW.xref)):
                        MW.integratedx += len(MW.xref[i])
                        MW.peakwidths[i] = len(MW.xref[i])
                    MW.IntFraction = MW.integratedx/len(MW.xzero[massrange])
                else:
                    MW.integratedx = 0 
                    MW.peakwidths = np.zeros(len(MW.sumlist))
                    for i in range(0, len(MW.sumlist)):
                        MW.integratedx += (MW.bounds[2*i+1]-MW.bounds[2*i])*(MW.xzero[1]-MW.xzero[0])
                        MW.peakwidths[i] = float((MW.bounds[2*i+1]-MW.bounds[2*i])*(MW.xzero[1]-MW.xzero[0]))
                    MW.IntFraction = MW.integratedx/len(MW.xzero[massrange])
                                           
                    
                MW.N = np.sqrt(MW.IntFraction)*np.sqrt(MW.GaborFraction)*MW.NFrms*2*np.pi
                MW.Nh = np.sqrt(MW.GaborFraction)*MW.NFrms/np.sqrt(len(MW.xzero[massrange]))*2*np.pi
            except AttributeError:
                print('Unable to determine final noise')
            
            MW.save_FT(self)
            
        try:
            MW.zerofullBL2 = MW.zerofullBL2*sum(MW.zerofullBL)/sum(MW.zerofullBL2)
            temp_path = os.path.join(MW.foldernameDecon, 'Noise Baseline.csv')
            np.savetxt(temp_path, np.c_[np.round(MW.xzeroBL2, 4), np.round(MW.zerofullBL2, 4)], delimiter=',')
        except AttributeError:
            pass
        try:
            temp_path = os.path.join(MW.foldernameDecon, 'Baseline spectrum.csv')
            np.savetxt(temp_path, np.c_[np.round(MW.xzeroBL, 4), np.round(MW.zerofullBL, 4)], delimiter=',')
        except AttributeError:
            pass
        
        print('Batch complete')
        MW.message = 'Batch complete'
        MW.Err = False
        self.dialog = iFAMS_message(self)
        self.dialog.show()
        
        
    def batch_param(self):
        if MW.namebase == 'From Clipboard':
            MW.ReturnFunc = MW.batch_param
            Folder_Select(self).show()
            return
        if MW.batchcheck == False:
            MW.preferences = ["6.3"]
            MW.preflabel = batch.template()
            
            if MW.proton > 0:
                MW.preferences.append('positive')
            else:
                MW.preferences.append('negative')
                
            MW.preferences.append('STFT')
            
            MW.preferences.append(str(int(MW.datarho)))
            global RT1
            global RT2
            global RTfull
            MW.preferences.append(str(RT1))
            MW.preferences.append(str(RT2))
            MW.preferences.append(str(RTfull))
            
            if MW.domaintrunc == True:
                MW.preferences.append(str(True))
                MW.preferences.append(str(np.round(MW.domainminx,7)))
                MW.preferences.append(str(np.round(MW.domainmaxx,7)))
            else:
                MW.preferences.append(str(False))
                MW.preferences.append(str(0))
                MW.preferences.append(str(0))
            
            MW.preferences.append(str(np.round(MW.NFmin,7)))
            MW.preferences.append(str(np.round(MW.NFmax,7)))
            MW.preferences.append(str(MW.Nlinear_sub))
            MW.preferences.append(str(MW.Nseg_sub))
            if MW.Nsubtraction == True:
                MW.preferences.append(str(True))
                MW.preferences.append(str(MW.RMSDmult))
            else:
                MW.preferences.append(str(False))
                MW.preferences.append(str(0))
            
            MW.preferences.append(str(MW.winnum))
            MW.preferences.append(str(MW.OF))
            MW.preferences.append(str(MW.window))
            if len(MW.batchstore) > 0:
                MW.preferences.append('NA')
                MW.preferences.append('NA')
                MW.preferences.append('NA')
            else:
                MW.preferences.append(str(MW.stft_selec_man))
                MW.preferences.append(str(MW.stft_box_thresh))
                MW.preferences.append(str(MW.stft_range_thresh))
            MW.preferences.append(str(True))
            
            MW.preferences.append('NA')
            MW.preferences.append('NA')
            
            MW.preferences.append(str(MW.smoothed))
            MW.preferences.append(str(MW.basecor))
            MW.preferences.append(str(MW.NZbased)) 
            MW.preferences.append(str(MW.based))
            MW.preferences.append(str(MW.paramMatch))
            MW.preferences.append(str(MW.minimumX))
            MW.preferences.append(str(MW.maximumX))
            MW.preferences.append(str(MW.minimumY))
            MW.preferences.append(str(MW.deltaY))
            MW.preferences.append(str(MW.ntol))
            if MW.paramMatch == False:
                MW.preferences.append(str(int(len(MW.bounds)/2)))
            else:
                MW.preferences.append(str(0))
            
            MW.preferences.append(str(len(MW.cs)))
            MW.preferences.append(str(len(MW.rlist)))
            if len(MW.batchstore) > 0:
                MW.preferences.append(str(len(MW.batchstore)))
            else:
                MW.preferences.append(str(1))
            for i in range(len(MW.cs)):
                MW.preferences.append(str(MW.cs[i]))
                MW.preflabel.append('selection '+str(i)+' charge state')
            for i in range(len(MW.rlist)):
                for j in range(5):
                    MW.preferences.append(str(MW.rlist[i][j]))
                    if j == 4:
                        MW.preflabel.append('rectangle '+str(i+1)+' selection number')
                    else:
                        MW.preflabel.append('rectangle '+str(i+1)+' coordinate '+str(j+1))
            
            if MW.paramMatch == False:
                for i in range(0,len(MW.bounds),2):
                    MW.preferences.append(str(MW.bounds[i]))
                    MW.preferences.append(str(MW.bounds[i+1]))
                    MW.preflabel.append('low bound '+str(int(i/2+1)))
                    MW.preflabel.append('high bound '+str(int(i/2+1)))
                  
    def batch_save(self):
        day = str(datetime.datetime.now().day)
        if int(day) < 10:
            day = str(0)+day
        month = str(datetime.datetime.now().month)
        if int(month) < 10:
            month = str(0)+month
        year = str(int(float(datetime.datetime.now().year)-2000))
        
        try:
            MW.batchfolder
            newFolder = False
        except AttributeError:
            newFolder = True
            count = 1
            folder = str(year+month+day+'_Batch'+str(count))
            test_path = os.path.exists(os.path.join(os.path.dirname(MW.namebase),folder))
            while test_path == True:
                count += 1
                folder = str(year+month+day+'_Batch'+str(count))
                test_path = os.path.exists(os.path.join(os.path.dirname(MW.namebase),folder))
            MW.batchfolder = os.path.join(os.path.dirname(MW.namebase),folder)
            os.mkdir(MW.batchfolder)
            if MW.batchcheck == False:
                label = str(year+month+day)+'_batch_params_from_' + os.path.basename(MW.namebase) + '.csv'
            else:
                label = str(year+month+day)+'_batch_params_from_import.csv'
            temp_path = os.path.join(MW.batchfolder, label)
            np.savetxt(temp_path, np.c_[MW.preferences,MW.preflabel], delimiter=',',fmt='%s')
        else:
            if MW.batchcheck == False:
                label = str(year+month+day)+'_batch_params_from_' + os.path.basename(MW.namebase) + '.csv'
            else:
                label = str(year+month+day)+'_batch_params_from_import.csv'
            temp_path = os.path.join(MW.batchfolder, label)
            np.savetxt(temp_path, np.c_[MW.preferences,MW.preflabel], delimiter=',',fmt='%s')
        print('Saved batch parameters as ' + str(temp_path))
        print(' ')
        MW.message = str('Saved batch parameters as:\n' 
                         + str(label))
        if newFolder == True:
            MW.submessage = str('Parameters saved in new folder:\n'+str(os.path.basename(MW.batchfolder)))
        else:
            MW.submessage = str('Parameters saved to batch folder:\n'+str(os.path.basename(MW.batchfolder)))
        MW.Err = False
        self.dialog = iFAMS_message(self)
        self.dialog.show()
        
    def batch_load(self):
        self.Bname = QtWidgets.QFileDialog.getOpenFileName(self,'Open File')
        self.Bnamestr = str(self.Bname[0])
        self.Bnamebase = os.path.splitext(self.Bnamestr)[0]
        self.submess = []
        try:
            del(MW.batchfolder)
        except AttributeError:
            pass
        try:
            MW.preferences,MW.preflabel = np.loadtxt(self.Bnamestr,unpack=True,delimiter=',',dtype=(str,str))
        except UnboundLocalError:
            print('An exception was found')
        except ValueError:
            MW.message('Unable to load batch parameters.\nEither no data loaded or wrong data type.')
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        if MW.preferences[0].lower() == '6.3':
            pass
        else:
            MW.message = 'Warning: batch parameter file might not be in expected format.'
            MW.submessage = 'expected verion: 6.3 | read version: '+str(MW.preferences[0])
            print(MW.message)
            print(MW.submessage)
            MW.Err = False
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
        if MW.preferences[1].lower() == 'positive':
            MW.ionmode.setChecked(True)
            MW.ionmode_set(self)
        elif MW.preferences[1].lower() == 'negative':
            MW.ionmode.setChecked(False)
            MW.ionmode_set(self)
        print('iFAMS process in '+str(MW.preferences[1].lower())+' ion mode')
        self.submess.append('iFAMS process in '+str(MW.preferences[1].lower())+' ion mode')
        
        MW.datarho = int(MW.preferences[3])
        global RT1
        global RT2
        global RTfull
        RT1 = str(MW.preferences[4])
        RT2 = str(MW.preferences[5])
        RTfull = eval(MW.preferences[6].title())
        if RTfull == True:
            print('Extract MS data across all RT')
            self.submess.append('Extract MS data across all RT')
        if float(RT2) > 0:
            print('Extract MS data from RT of '+RT1+'-'+RT2+' min')
            self.submess.append('Extract MS data from RT of '+RT1+'-'+RT2+' min')
        MW.domaintrunc = eval(MW.preferences[7].title())
        MW.domainminx = float(MW.preferences[8])
        MW.domainmaxx = float(MW.preferences[9])
        if MW.domaintrunc == True:
            print('Adjust MS domain to ' +str(MW.domainminx)+'-'+str(MW.domainmaxx)+' m/z')
            self.submess.append(str('Adjust MS domain to ' +str(np.round(MW.domainminx,2))+'-'+str(np.round(MW.domainmaxx,2))+' m/z'))
        MW.NFmin = float(MW.preferences[10])
        MW.NFmax = float(MW.preferences[11])
        print('Calculate noise between frequencies ' +str(np.round(MW.NFmin,2))+' and '+str(np.round(MW.NFmax,2)))
        self.submess.append(str('Calculate noise between frequencies ' +str(np.round(MW.NFmin,2))+' and '+str(np.round(MW.NFmax,2))))
        MW.Nlinear_sub = eval(MW.preferences[12].title())
        if MW.Nlinear_sub == True:
            print('Apply Linear Baseline Subtration to MS')
            self.submess.append(str('Apply Linear Baseline Subtration to MS'))
        MW.Nseg_sub = eval(MW.preferences[13].title())
        if MW.Nseg_sub == True:
            print('Apply Segmented Baseline Subtration to MS')
            self.submess.append(str('Apply Segmented Baseline Subtration to MS'))
        MW.Nsubtraction = eval(MW.preferences[14].title())
        if MW.Nsubtraction == True:
            MW.RMSDmult = float(MW.preferences[15])
            print('Subtract '+str(MW.RMSDmult)+'x noise from MS baseline')
            self.submess.append(str('Subtract '+str(MW.RMSDmult)+'x noise from MS baseline'))
         
        MW.winnum = int(MW.preferences[16])
        MW.OF = int(MW.preferences[17])
        print('STFT made with '+str(MW.winnum)+' frequency channels and '+str(MW.OF)+'x oversampling')
        self.submess.append(str('STFT made with '+str(MW.winnum)+' frequency channels \nand '+str(MW.OF)+'x oversampling'))
        MW.window = str(MW.preferences[18])
            
        self.cslen = int(MW.preferences[36])
        self.rlistlen = int(MW.preferences[37])
        MW.seriesnum = int(MW.preferences[38])
        MW.cs = []
        MW.rlist = []
        base = 39
        for i in range(self.cslen):
            MW.cs.append(int(MW.preferences[base+i]))
        print('Using charge states of '+str(MW.cs))
        for i in range(self.rlistlen):
            MW.rlist.append([0,0,0,0,0])
            for j in range(5):
                MW.rlist[i][j] = float(MW.preferences[base+self.cslen+(5*i)+j])
        print('With ' +str(self.rlistlen)+ ' total rectangles')
        self.submess.append(str('Using ' +str(self.rlistlen)+ ' total rectangles'))
        self.submess.append(str('For '+str(self.cslen)+' total charge states'))
        if MW.seriesnum > 1:
            print('From '+str(MW.seriesnum)+' different series')
            self.submess.append('From '+str(MW.seriesnum)+' different series')
            
        MW.smoothed = eval(MW.preferences[25].title())
        if MW.smoothed == True:
            print('Smooth deconvolved spectrum')
            self.submess.append('Smooth deconvolved spectrum')
        MW.basecor = eval(MW.preferences[26].title())
        if MW.basecor == True:
            print('Perform Minima Baseline Correction before integration')
            self.submess.append(str('Perform Minima Baseline Correction before integration'))
        MW.NZbased = eval(MW.preferences[27].title())
        if MW.NZbased == True:
            print('Perform Fourier Baseline Correction before integration')
            self.submess.append(str('Perform Fourier Baseline Correction before integration'))
        MW.based = eval(MW.preferences[28].title())
        if MW.based == True:
            print('Perform Segmented Baseline Correction to integration')
            self.submess.append(str('Perform Segmented Baseline Correction to integration'))
        MW.paramMatch = eval(MW.preferences[29].title())
        MW.minimumX = float(MW.preferences[30])
        MW.maximumX = float(MW.preferences[31])
        MW.minimumY = float(MW.preferences[32])
        MW.deltaY = float(MW.preferences[33])
        MW.ntol = float(MW.preferences[34])
        if MW.paramMatch == True:
            print('Perform Peak Integration between '+str(MW.minimumX)+' and '+str(MW.maximumX)+'\n'
                  ' with Minimum y of '+str(MW.minimumY)+'\n'
                  '      Minimum Distance Between Peaks of '+str(MW.deltaY)+'\n'
                  '      Noise Tolerance of '+str(MW.ntol))
            self.submess.append(str('Perform Peak Integration between '+str(MW.minimumX)+' and '+str(MW.maximumX)+'\n'
                  ' with Minimum y of '+str(MW.minimumY)+'\n'
                  '      Minimum Distance Between Peaks of '+str(MW.deltaY)+'\n'
                  '      Noise Tolerance of '+str(MW.ntol)))
        if MW.paramMatch == False:
            print('Perform integration with Bounds Match')
            self.submess.append(str('Perform integration with Bounds Match'))
            self.peaknum = int(MW.preferences[35])
            MW.bounds = []
            for i in range(0,self.peaknum*2):
                MW.bounds.append(float(MW.preferences[base+self.cslen+self.rlistlen*5+i]))
    
        MW.batchcheck = True
        print('Batch parameters loaded')
        MW.message = 'Batch parameters loaded:'
        MW.Err = False
        MW.submessage = self.submess[0]
        for i in range(1,len(self.submess)):
            MW.submessage += str('\n'+str(self.submess[i]))
        self.dialog = iFAMS_message(self)
        self.dialog.show()

    def guided_search(self):
        MW.autosearch = False
        MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)
        self.dialog = guided_dialog(self)
        self.dialog.show()
        
    def Go_search(self):
        piW = (MW.xshort[2]-MW.xshort[1])
        piH = (MW.yshort[2]-MW.yshort[1])
        poW = (MW.xint[2]-MW.xint[1])
        try:
            if MW.rangeadjust == True:
                ranthresh1 = 1.25
                ranthresh2 = 0.06
            else:
                ranthresh1 = 10000000
                ranthresh2 = 0.00000001
        except AttributeError:
                ranthresh1 = 1.25
                ranthresh2 = 0.06

        #search by ionmass
        if type(MW.ionmass) == float:
            print('Searching for ions with mass of '+str(MW.ionmass)+' Da')
            loZ = np.max((int(np.ceil(MW.ionmass/MW.xint[-1])),2))
            hiZ = np.min((int(np.floor(MW.ionmass/MW.xint[0])),90,int(MW.ionmass/550)))
            zcheck = np.flip(np.array(range(loZ,hiZ)))
            mzcheck = (MW.ionmass+zcheck*MW.proton)/zcheck
            mzchecki = np.round((mzcheck-MW.xshort[0])/piW)
            mzchecki = mzchecki.astype(int)
            IsoSlice = 1/MW.xshort*MW.ionmass
            IsoSliceX = MW.xshort*0
            if loZ < int((MW.yshort[-1]/2)):
                for i in range(len(MW.xshort)):
                    checki = int(np.round(IsoSlice[i]/piH))
                    if checki >= int(len(MW.yshort)/2):
                        continue
                    IsoSliceX[i] = np.absolute(MW.X[i][checki])
                if np.average(IsoSliceX[mzchecki]) > 1.3*np.average(IsoSliceX):
                    mzpi = mzchecki[np.argmax(IsoSliceX[mzchecki])]
                    Z = zcheck[np.argmax(IsoSliceX[mzchecki])]
                    FoI = IsoSlice[mzpi]
                else:
                    fi1 = int(np.round(0.2/piH))
                    fi2 = int(len(MW.yshort)/4)
                    freqis = fi1 + np.argmax((np.absolute(MW.X[mzchecki,fi1:fi2])),axis=1)
                    Dm = np.average(zcheck/MW.yshort[freqis])
                    IsoSlice = 1/MW.xshort*MW.ionmass/Dm
                    IsoSliceX = MW.xshort*0
                    for i in range(len(MW.xshort)):
                        checki = int(np.round(IsoSlice[i]/piH))
                        IsoSliceX[i] = np.absolute(MW.X[i][checki])
                    mzpi = mzchecki[np.argmax(IsoSliceX[mzchecki])]
                    Z = zcheck[np.argmax(IsoSliceX[mzchecki])]
                    FoI = IsoSlice[mzpi]
            else:
                fi1 = int(np.round(0.2/piH))
                fi2 = int(len(MW.yshort)/4)
                freqis = fi1 + np.argmax((np.absolute(MW.X[mzchecki,fi1:fi2])),axis=1)
                Dm = np.average(zcheck/MW.yshort[freqis])
                IsoSlice = 1/MW.xshort*MW.ionmass/Dm
                IsoSliceX = MW.xshort*0
                for i in range(len(MW.xshort)):
                    checki = int(np.round(IsoSlice[i]/piH))
                    IsoSliceX[i] = np.absolute(MW.X[i][checki])
                mzpi = mzchecki[np.argmax(IsoSliceX[mzchecki])]
                Z = zcheck[np.argmax(IsoSliceX[mzchecki])]
                FoI = IsoSlice[mzpi]
            print('Investigating frequency of '+str(np.round(FoI,2)))
        
            if abs(FoI-MW.FoI) > MW.FoI or abs(FoI-MW.FoI) > FoI:
                MW.FoI = FoI
                if MW.FoI >= 8:
                    MW.winnum = int((MW.maxF)*20)
                    MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -5)*len(MW.yFT)/MW.maxF):int((MW.FoI +5)*len(MW.yFT)/MW.maxF)]))/20)
                elif MW.FoI < 8 and MW.FoI > 1:
                    MW.winnum = int((MW.maxF)*30)
                    MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -1)*len(MW.yFT)/MW.maxF):int((MW.FoI +1)*len(MW.yFT)/MW.maxF)]))/10)
                elif MW.FoI <= 1 and MW.FoI >= 0.1:
                    MW.winnum = int((MW.maxF)*80)
                    MW.vmax = int(max(abs(MW.yFT[int((MW.FoI -0.05)*len(MW.yFT)/MW.maxF):int((MW.FoI +0.05)*len(MW.yFT)/MW.maxF)]))/10)
                else:
                    MW.winnum = int((MW.maxF)*20)
                    MW.vmax = int(np.max(np.absolute(MW.yFT[int((MW.FoI)*len(MW.yFT)/MW.maxF)-3:int((MW.FoI)*len(MW.yFT)/MW.maxF)+3]))/10)
                MW.stdev = int(0.1 * MW.winnum)
                MW.window = str('gaussian')
                MW.X = STFT.stft(MW.ypadd, MW.winnum, MW.window, MW.OF)
                ytemp = STFT.re_stft(MW.X, MW.winnum, MW.ypadd, MW.yint,MW.OF)
                MW.correction_factor = max(MW.yint)/max(ytemp)
                MW.Xtemp = MW.X*0
                MW.Xflat = abs(MW.X.flatten("C")) / 100
                MW.mzspacing = float((MW.xint[-1] - MW.xint[0])/len(MW.xint))
                MW.xshort = np.linspace((MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing + min(MW.xpadd), max(MW.xpadd) - (MW.winnum/2-MW.winnum*0.5/MW.OF) *MW.mzspacing, len(MW.X))
                MW.yshort = np.linspace(0 - (max(MW.xFT)/len(MW.X[0]))/2, max(MW.xFT) + (max(MW.xFT)/len(MW.X[0]))/2, len(MW.X[0]))
                MW.mzchan = float((MW.xint[-1]-MW.xshort[0])/((MW.xshort[-1]-MW.xshort[0])/len(MW.xshort)))
                piW = (MW.xshort[2]-MW.xshort[1])
                piH = (MW.yshort[2]-MW.yshort[1])
                poW = (MW.xint[2]-MW.xint[1])
            mzpi = int(np.round(((MW.ionmass+Z*MW.proton)/Z-MW.xshort[0])/piW))
            mzpi2 = int(np.round(((MW.ionmass+(Z+1)*MW.proton)/(Z+1)-MW.xshort[0])/piW))
            if FoI < 0.2:
                boxbott = 0
                boxtop = 1
            else:
                boxtop = MW.yshort[int(FoI/piH)+8]
                for i in range(1,int(1/piH)):
                    if np.absolute(MW.X[mzpi,int(FoI/piH)+i]) > 0.1*np.absolute(MW.X[mzpi,int(FoI/piH)]):
                        boxtop = np.min((MW.yshort[int(FoI/piH)+i],FoI+3))
                        if boxtop == FoI+3:
                            break
                    else:
                        break
                boxbott = np.max((FoI-(boxtop-FoI),0))
            boxleft = MW.xshort[mzpi-3]
            for i in range(1,int(30/piW)):
                if np.absolute(MW.X[mzpi-i,int(FoI/piH)]) > 0.08*np.absolute(MW.X[mzpi,int(FoI/piH)]):
                    boxleft = MW.xshort[mzpi-i]
                else:
                    break
            boxright = MW.xshort[mzpi+3]
            for i in range(1,int(30/piW)):
                if np.absolute(MW.X[mzpi+i,int(FoI/piH)]) > 0.08*np.absolute(MW.X[mzpi,int(FoI/piH)]):
                    boxright = MW.xshort[mzpi+i]
                else:
                    break
            MW.rlist.append([boxleft,boxright,boxbott,boxtop,0])
            MW.rlist.append([MW.xshort[mzpi2]-(MW.xshort[mzpi]-boxleft),MW.xshort[mzpi2]+(boxright-MW.xshort[mzpi]),np.max((MW.yshort[int(FoI/Z*(Z+1)/piH)]-(boxtop-boxbott)/2,0)),MW.yshort[int(FoI/Z*(Z+1)/piH)]+(boxtop-boxbott)/2,1])
        
            MW.Ilist = []
            for j in range(len(MW.rlist)):
                lowmzI = int(np.round((MW.rlist[j][0]-MW.xshort[0])/piW))
                himzI = int(np.round((MW.rlist[j][1]-MW.xshort[0])/piW))
                lowFI = int(np.round(MW.rlist[j][2]/piH))
                hiFI = int(np.round(MW.rlist[j][3]/piH))
                MW.Ilist.append([lowmzI,himzI,lowFI,hiFI,MW.rlist[j][4]])
            MW.glist = []
            MW.Xtemp = MW.X*0
            for j in range(len(MW.Ilist)):
                if MW.Ilist[j][2] <= 1:
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,0:MW.Ilist[j][3]+1]
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-1]
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-1]
                else:
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,MW.Ilist[j][2]:MW.Ilist[j][3]+1]
                    MW.Xtemp[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1] = MW.X[MW.Ilist[j][0]:MW.Ilist[j][1]+1,-MW.Ilist[j][3]:-MW.Ilist[j][2]+1]
                MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                MW.glist.append(MW.yrecon)
                MW.Xtemp = MW.X*0
        
        # Near-zero frequencies
        if MW.rlist[0][2] <= 0 and MW.rlist[1][2] <= 0:
            try:
                del(MW.Harmrec)
            except AttributeError:
                pass
            if MW.xint[np.argmax(MW.glist[0])] > MW.xint[np.argmax(MW.glist[1])]:
                ifft1 = MW.glist[0]
                ifft2 = MW.glist[1]
                r1 = MW.rlist[0]
                r2 = MW.rlist[1]
            else:
                ifft1 = MW.glist[1]
                ifft2 = MW.glist[0]
                r1 = MW.rlist[1]
                r2 = MW.rlist[0]
            mz1 = MW.xint[np.argmax(ifft1)]
            mz2 = MW.xint[np.argmax(ifft2)]
            z1 = int(np.round(1/(mz1/mz2-1)))
            z2 = z1 + 1
            mass1 = mz1*z1-z1*MW.proton
            mass2 = mz2*z2-z2*MW.proton
            mass = np.average((mass1,mass2))
            print('Searching for ions with rough mass of '+str(np.round(mass,0))+' Da')
            print('Starting from charge states '+str(z1)+' and '+str(z2))
                        
            #exploring distribution
            MW.Xslice = np.sum(np.absolute(MW.X[0:len(MW.xshort),0:2]),axis=1)
            slicex = MW.xshort
            slicey1 = np.zeros(len(MW.xshort))+MW.yshort[1]+piH/2
            slicey2 = np.zeros(len(MW.xshort))
            F1 = 0
            F2 = np.round(MW.yshort[1]+piH/2,3)
            boxW = int(np.round(np.average((MW.rlist[0][1]-MW.rlist[0][0],MW.rlist[1][1]-MW.rlist[1][0]))/piW))
            boxH = int(np.round(np.average((MW.rlist[0][3]-MW.rlist[0][2],MW.rlist[1][3]-MW.rlist[1][2]))/piH))
            pixel2I = int(np.round((mz2-MW.xshort[0])/piW))
            pixel1I = int(np.round((mz1-MW.xshort[0])/piW))
            avgsum = []
            for i in range(0,20):
                avgsum.append(np.min(np.absolute(MW.X[int(i*np.floor(len(MW.xshort)/20)):int((i+1)*np.floor(len(MW.xshort)/20))][0])))
            AVG = np.average(avgsum)
            peak = np.average((np.absolute(MW.X[pixel1I][0]),np.absolute(MW.X[pixel2I][0])))
            peakMS = np.average((np.max(MW.yint[int(np.round((mz1-MW.xint[0])/poW))-2:int(np.round((mz1-MW.xint[0])/poW))+3]),np.max(MW.yint[int(np.round((mz2-MW.xint[0])/poW))-2:int(np.round((mz2-MW.xint[0])/poW))+3])))
            lowz = np.max((int(np.ceil(mass/MW.xint[-1])),int(z1/3)))
            end = z1-lowz+1
            for i in range(1,end):
                checki = int(np.round(((mass+(z1-i)*MW.proton)/(z1-i)-MW.xshort[0])/piW))
                checki = checki-1+np.argmax(np.absolute(MW.X[checki-1:checki+2,0]))
                checkiMS = int(np.round(((mass+(z1-i)*MW.proton)/(z1-i)-MW.xint[0])/poW))
                checkiMS = checki-1+np.argmax(MW.yint[checkiMS-1:checkiMS+2])
                if checki + boxW >= len(MW.xshort):
                    break
                if np.absolute(MW.X[checki][0])-AVG > ranthresh1*3*(peak-AVG):
                    break
                if np.absolute(MW.X[checki][0])-AVG > ranthresh2/2*(peak-AVG) and MW.yint[checkiMS] > 0.15*(peakMS):
                    lowz = z1-i
                else:
                    break
            highz = np.min((int(np.floor(mass/MW.xint[0])),90))
            end = highz-z2+1
            for i in range(1,end):
                checki = int(np.round(((mass+(z2+i)*MW.proton)/(z2+i)-MW.xshort[0])/piW))
                checki = checki-1+np.argmax(np.absolute(MW.X[checki-1:checki+2,0]))
                checkiMS = int(np.round(((mass+(z2+i)*MW.proton)/(z2+i)-MW.xint[0])/poW))
                checkiMS = checki-1+np.argmax(MW.yint[checkiMS-1:checkiMS+2])
                if checki - boxW <= 0:
                    break
                if np.absolute(MW.X[checki][0])-AVG > ranthresh1*3*(peak-AVG):
                    break
                if np.absolute(MW.X[checki][0])-AVG > ranthresh2/2*(peak-AVG) and MW.yint[checkiMS] > 0.15*(peakMS):
                    highz = z2+i
                else:
                    break
            
            #refining distribution
            zcheck = np.flip(np.array(range(lowz,highz+1)))
            peakcheck = []
            for i in range(len(zcheck)):
                checki = int(np.round(((mass+(zcheck[i])*MW.proton)/(zcheck[i])-MW.xshort[0])/piW))
                if checki >= len(MW.xshort):
                    peakcheck.append(0)
                    continue
                try:
                    checki = checki-1+np.argmax(np.absolute(MW.X[checki-1:checki+2,0]))
                except ValueError:
                    pass
                peakcheck.append(np.absolute(MW.X[checki][0]))
            MAX = np.max(peakcheck)
            MAXi = np.argmax(peakcheck)
            MAXz = zcheck[MAXi]
            while MAXz-z2 > 8:
                MAX = np.max(peakcheck[MAXi+2:-1])
                MAXi = np.argmax(peakcheck[MAXi+2:-1])+MAXi+2
                MAXz = zcheck[MAXi]
            while MAXz-z1 < -8:
                MAX = np.max(peakcheck[0:MAXi-2])
                MAXi = np.argmax(peakcheck[0:MAXi-2])
                MAXz = zcheck[MAXi]
            for i in range(1,MAXi):
                if zcheck[MAXi-i] <= z2:
                    highz = z2
                    continue
                if (mass+(zcheck[MAXi-i])*MW.proton)/(zcheck[MAXi-i])-boxW*piW <= MW.xint[0]:
                    highz = zcheck[MAXi-i+1]
                    break
                if peakcheck[MAXi-i]-AVG > (peakcheck[MAXi-i+1]-AVG)*ranthresh1:
                    if zcheck[MAXi-i+1] < MAXz + MAXz/5:
                        pass
                    else:
                        highz = zcheck[MAXi-i+1]
                        break
                if peakcheck[MAXi-i]-AVG < ranthresh2*(MAX-AVG):
                    highz = zcheck[MAXi-i]
                    break
            for i in range(1,len(zcheck)-MAXi):
                if zcheck[MAXi+i] >= z1:
                    lowz = z1
                    continue
                if (mass+(zcheck[MAXi+i])*MW.proton)/(zcheck[MAXi+i])+boxW*piW >= MW.xint[-1]:
                    lowz = zcheck[MAXi+i-1]
                    break
                if peakcheck[MAXi+i]-AVG > (peakcheck[MAXi+i-1]-AVG)*ranthresh1:
                    if zcheck[MAXi+i-1] > MAXz - MAXz/5:
                        pass
                    else:
                        lowz = zcheck[MAXi+i-1]
                        break
                if peakcheck[MAXi+i]-AVG < ranthresh2*(MAX-AVG):
                    lowz = zcheck[MAXi+i]
                    break
            MW.zlist = np.array(range(lowz,highz+1))
            print('Selecting charge states:')
            print(MW.zlist)
            print('Most abundant charge state: '+str(MAXz))
            
            #Determined boxes
            if MW.boxadjust == False:
                lomass1 = r1[0]*z1-z1*MW.proton
                himass1 = r1[1]*z1-z1*MW.proton
                lomass2 = r2[0]*z2-z2*MW.proton
                himass2 = r2[1]*z2-z2*MW.proton
                lomass = np.average((lomass1,lomass2))
                himass = np.average((himass1,himass2))
                boxtop = np.average((r1[3],r2[3]))
                MW.rlist = []
                MW.Ilist = []
                for i in range(len(MW.zlist)):
                    lomz = (lomass+MW.zlist[i]*MW.proton)/MW.zlist[i]
                    lomzi = int(np.round((lomz-MW.xshort[0])/piW))
                    himz = (himass+MW.zlist[i]*MW.proton)/MW.zlist[i]
                    himzi = int(np.round((himz-MW.xshort[0])/piW))
                    MW.rlist.append([lomz,himz,0,boxtop,i])
                    MW.Ilist.append([lomzi,himzi,0,int(np.round(boxtop/piH)),i])
            
            #identifying signal
            else:
                MW.rlist = []
                MW.Ilist = []
                for i in range(len(MW.zlist)):
                    checki = int(np.round(((mass+(MW.zlist[i])*MW.proton)/(MW.zlist[i])-MW.xshort[0])/piW))
                    try:
                        checki = checki-1+np.argmax(np.absolute(MW.X[checki-1:checki+2,0]))
                    except ValueError:
                        pass
                    boxtop = boxH
                    boxleft = checki-int(np.ceil(boxW/3))
                    boxright = checki+int(np.floor(boxW*2/3))
                    for j in range(1,boxW):
                        if np.absolute(MW.X[checki+j][0])-AVG < 0.1*(np.absolute(MW.X[checki][0])-AVG):
                            boxright = checki+j
                            break
                        if np.absolute(MW.X[checki+j][0])-AVG > 1.05*(np.absolute(MW.X[checki+j-1][0])-AVG):
                            if np.absolute(MW.X[checki+j+1][0])-AVG < 1.05*(np.absolute(MW.X[checki+j-1][0])-AVG):
                                continue
                            else:
                                boxright = checki+j-1
                                break
                    for j in range(1,boxW):
                        if np.absolute(MW.X[checki-j][0])-AVG < 0.1*(np.absolute(MW.X[checki][0])-AVG):
                            boxleft = checki-j
                            break
                        if np.absolute(MW.X[checki-j][0])-AVG > 1.05*(np.absolute(MW.X[checki-j+1][0])-AVG):
                            if np.absolute(MW.X[checki-j-1][0])-AVG < 1.05*(np.absolute(MW.X[checki-j+1][0])-AVG):
                                continue
                            else:
                                boxleft = checki-j+1
                                break
                    for j in range(int(boxH/3),3*boxH):
                        if np.absolute(MW.X[checki][j])-AVG < 0.01*(np.absolute(MW.X[checki][1])-AVG):
                            boxtop = j
                            break
                    if boxright >= len(MW.xshort):
                        boxright = len(MW.xshort)-1
                    if boxleft <= 0:
                        boxleft = 0
                    MW.rlist.append([MW.xshort[boxleft],MW.xshort[boxright],MW.yshort[0],MW.yshort[boxtop],i])
                    MW.Ilist.append([boxleft,boxright,0,boxtop,i])

            #extracting signal            
            MW.glist = []
            MW.Xtemp = MW.X*0
            MW.yBL2 = np.zeros(len(MW.xint))+ min(MW.yint)+ 1
            MW.glistBL2 = []
            MW.xintBL2,MW.yintBL2,MW.ypaddBL2,MW.po2BL2,MW.xpaddBL2 = load.datapro(MW.xint,MW.yBL2)
            MW.XBL2 = STFT.stft(MW.ypaddBL2, MW.winnum, MW.window, MW.OF)
            MW.XtempBL2 = MW.XBL2*0
            slicetext = []
            for i in range(len(MW.Ilist)):
                MW.Xtemp[MW.Ilist[i][0]:MW.Ilist[i][1]+1,0:MW.Ilist[i][3]+1] = MW.X[MW.Ilist[i][0]:MW.Ilist[i][1]+1,0:MW.Ilist[i][3]+1]
                MW.Xtemp[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:len(MW.X[0])] = MW.X[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:len(MW.X[0])]
                MW.XtempBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,0:MW.Ilist[i][3]+1] = MW.XBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,0:MW.Ilist[i][3]+1]
                MW.XtempBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:len(MW.X[0])] = MW.XBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:len(MW.X[0])]
                MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                MW.glist.append(MW.yrecon)
                slicetext.append(int(np.argmax(np.absolute(MW.Xtemp))/len(MW.Xtemp[0])))
                MW.Xtemp = MW.X*0
                MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                MW.glistBL2.append(MW.yreconBL2)
                MW.XtempBL2 = MW.XBL2*0
            
        # Fundamentals 
        else:
            if MW.xint[np.argmax(MW.glist[0])] > MW.xint[np.argmax(MW.glist[1])]:
                ifft1 = MW.glist[0]
                ifft2 = MW.glist[1]
                r1 = MW.rlist[0]
                r2 = MW.rlist[1]
            else:
                ifft1 = MW.glist[1]
                ifft2 = MW.glist[0]
                r1 = MW.rlist[1]
                r2 = MW.rlist[0]
            mz1 = MW.xint[np.argmax(ifft1)]
            mz2 = MW.xint[np.argmax(ifft2)]
            z1 = int(np.round(1/(mz1/mz2-1)))
            z2 = z1 + 1
            mass1 = mz1*z1-z1*MW.proton
            mass2 = mz2*z2-z2*MW.proton
            mass = np.average((mass1,mass2))
            if MW.autosearch == True:
                try:
                    MW.foundmass.append(np.round(mass/10,0)*10)
                except AttributeError:
                    pass
            print('Searching for ions with rough mass of '+str(np.round(mass,0))+' Da')
            print('Starting from charge states '+str(z1)+' and '+str(z2))
            Xr1 = np.absolute(MW.X[int(np.round((r1[0]-MW.xshort[0])/piW)):int(np.round((r1[1]-MW.xshort[0])/piW))+1,int(np.round(r1[2]/piH)):int(np.round(r1[3]/piH))+1])            
            Xr2 = np.absolute(MW.X[int(np.round((r2[0]-MW.xshort[0])/piW)):int(np.round((r2[1]-MW.xshort[0])/piW))+1,int(np.round(r2[2]/piH)):int(np.round(r2[3]/piH))+1])            
            freq1 = MW.yshort[np.argmax(Xr1)%len(Xr1[0])+int(np.round(r1[2]/piH))]
            freq2 = MW.yshort[np.argmax(Xr2)%len(Xr2[0])+int(np.round(r2[2]/piH))]
            Deltam = np.average((z1/freq1,z2/freq2))
            sliceI = np.round((1/np.array(MW.xshort)*mass/Deltam)/piH)
            
            #exploring distribution
            MW.Xslice = np.zeros(len(MW.xshort))
            slicex = []
            slicey1 = []
            slicey2 = []
            for i in range(len(MW.xshort)):
                if sliceI[i] > len(MW.yshort)/2 or sliceI[i] < 0:
                    MW.Xslice[i] = 0
                else:
                    MW.Xslice[i] = np.average(np.absolute(MW.X[i,int(sliceI[i]-1):int(sliceI[i]+2)]))
                    slicex.append(MW.xshort[i])
                    slicey1.append(MW.yshort[int(sliceI[i]+1)]+piH/2)
                    slicey2.append(MW.yshort[int(sliceI[i]-1)]-piH/2)
            F1 = np.round(np.max((MW.yshort[int(np.min(sliceI)-1)],0)),2)
            F2 = np.round(MW.yshort[int(np.min((int(np.max(sliceI)+1), len(MW.yshort)/2)))],2)
            boxW = int(np.round(np.average((MW.rlist[0][1]-MW.rlist[0][0],MW.rlist[1][1]-MW.rlist[1][0]))/piW))
            boxH = int(np.round(np.average((MW.rlist[0][3]-MW.rlist[0][2],MW.rlist[1][3]-MW.rlist[1][2]))/piH))
            pixel2I = int(np.round((mz2-MW.xshort[0])/piW))
            pixel1I = int(np.round((mz1-MW.xshort[0])/piW))
            avgsum = []
            for i in range(0,20):
                avgsum.append(np.min(MW.Xslice[int(i*np.floor(len(MW.xshort)/20)):int((i+1)*np.floor(len(MW.xshort)/20))]))
            AVG = np.average(avgsum)
            peak = np.average((MW.Xslice[pixel1I],MW.Xslice[pixel2I]))
            lowz = np.max((int(np.ceil(mass/MW.xshort[-1])),int(z1/3)))
            end = z1-lowz+1
            for i in range(1,end):
                checki = int(np.round(((mass+(z1-i)*MW.proton)/(z1-i)-MW.xshort[0])/piW))
                if checki >= len(MW.xshort):
                    break
                checki = checki-1+np.argmax(MW.Xslice[checki-1:checki+2])
                if checki + boxW >= len(MW.xshort):
                    break
                if MW.Xslice[checki]-AVG > ranthresh1*3*(peak-AVG):
                    break
                if MW.Xslice[checki]-AVG > ranthresh2/2*(peak-AVG):
                    lowz = z1-i
                else:
                    break
            highz = np.min((int(np.floor(mass/MW.xshort[0])),90))
            end = highz-z2+1
            for i in range(1,end):
                checki = int(np.round(((mass+(z2+i)*MW.proton)/(z2+i)-MW.xshort[0])/piW))
                checki = checki-1+np.argmax(MW.Xslice[checki-1:checki+2])
                if checki - boxW <= 0:
                    break
                if MW.Xslice[checki]-AVG > ranthresh1*3*(peak-AVG):
                    break
                if MW.Xslice[checki]-AVG > ranthresh2/2*(peak-AVG):
                    highz = z2+i
                else:
                    break
            
            #refining distribution
            zcheck = np.flip(np.array(range(lowz,highz+1)))
            peakcheck = []
            for i in range(len(zcheck)):
                checki = int(np.round(((mass+zcheck[i]*MW.proton)/(zcheck[i])-MW.xshort[0])/piW))
                checki = checki-1+np.argmax(MW.Xslice[checki-1:checki+2])
                peakcheck.append(MW.Xslice[checki])
            MAX = np.max(peakcheck)
            MAXi = np.argmax(peakcheck)
            MAXz = zcheck[MAXi]
            while MAXz-z2 > 8:
                MAX = np.max(peakcheck[MAXi+2:-1])
                MAXi = np.argmax(peakcheck[MAXi+2:-1])+MAXi+2
                MAXz = zcheck[MAXi]
            while MAXz-z1 < -8:
                MAX = np.max(peakcheck[0:MAXi-2])
                MAXi = np.argmax(peakcheck[0:MAXi-2])
                MAXz = zcheck[MAXi]
            for i in range(1,MAXi):
                if zcheck[MAXi-i] <= z2:
                    highz = z2
                    continue
                if (mass+(zcheck[MAXi-i])*MW.proton)/(zcheck[MAXi-i])-boxW*piW <= MW.xint[0]:
                    highz = zcheck[MAXi-i+1]
                    break
                if peakcheck[MAXi-i]-AVG > (peakcheck[MAXi-i+1]-AVG)*ranthresh1:
                    if zcheck[MAXi-i+1] < MAXz + MAXz/5:
                        pass
                    else:
                        highz = zcheck[MAXi-i+1]
                        break
                if peakcheck[MAXi-i]-AVG < ranthresh2*(MAX-AVG):
                    highz = zcheck[MAXi-i+1]
                    break
                highz = zcheck[MAXi-i]
            for i in range(1,len(zcheck)-MAXi):
                if zcheck[MAXi+i] >= z1:
                    lowz = z1
                    continue
                if (mass+(zcheck[MAXi+i])*MW.proton)/(zcheck[MAXi+i])+boxW*piW >= MW.xint[-1]:
                    lowz = zcheck[MAXi+i-1]
                    break
                if peakcheck[MAXi+i]-AVG > (peakcheck[MAXi+i-1]-AVG)*ranthresh1:
                    if zcheck[MAXi+i-1] > MAXz - MAXz/5:
                        pass
                    else:
                        lowz = zcheck[MAXi+i-1]
                        break
                if peakcheck[MAXi+i]-AVG < ranthresh2*(MAX-AVG):
                    lowz = zcheck[MAXi+i-1]
                    break
                lowz = zcheck[MAXi+i]
            MW.zlist = np.array(range(lowz,highz+1))
            print('Selecting charge states:')
            print(MW.zlist)
            print('Most abundant charge state: '+str(MAXz))
            
            #Determined boxes
            if MW.boxadjust == False:
                lomass1 = r1[0]*z1-z1*MW.proton
                himass1 = r1[1]*z1-z1*MW.proton
                lomass2 = r2[0]*z2-z2*MW.proton
                himass2 = r2[1]*z2-z2*MW.proton
                lomass = np.average((lomass1,lomass2))
                himass = np.average((himass1,himass2))
                loFz = np.average((r1[2]/z1,r2[2]/z2))
                hiFz = np.average((r1[3]/z1,r2[3]/z2))
                MW.rlist = []
                MW.Ilist = []
                harmonic = np.zeros(len(MW.zlist))+2
                topharm = int(np.round(MW.yshort[int(len(MW.yshort)/2)]*Deltam/MW.zlist[0]))
                for i in range(len(MW.zlist)):
                    lomz = (lomass+MW.zlist[i]*MW.proton)/MW.zlist[i]
                    lomzi = int(np.round((lomz-MW.xshort[0])/piW))
                    himz = (himass+MW.zlist[i]*MW.proton)/MW.zlist[i]
                    himzi = int(np.round((himz-MW.xshort[0])/piW))
                    loF = loFz*MW.zlist[i]
                    loFi = int(np.round(loF/piH))
                    hiF = hiFz*MW.zlist[i]
                    hiFi = int(np.round(hiF/piH))
                    checki = int(lomzi+(himzi-lomzi)/2)
                    checkiy = int(loFi+(hiFi-loFi)/2)
                    MW.rlist.append([lomz,himz,loF,hiF,i])
                    MW.Ilist.append([lomzi,himzi,loFi,hiFi,i])
                    for j in range(2,topharm+1):
                        if np.absolute(MW.X[checki][checkiy*j])-AVG < 0.01*(np.absolute(MW.X[checki][checkiy])-AVG):
                            break
                        else:
                            harmonic[i] = j
                MW.Harmrec = int(np.round(np.average(harmonic)))
            
            #identifying signal
            else:
                MW.rlist = []
                MW.Ilist = []
                harmonic = np.zeros(len(MW.zlist))+2
                topharm = int(np.round(MW.yshort[int(len(MW.yshort)/2)]*Deltam/MW.zlist[0]))
                BoxThresh = 0.1
                for i in range(len(MW.zlist)):
                    checki = int(np.round(((mass+(MW.zlist[i])*MW.proton)/(MW.zlist[i])-MW.xshort[0])/piW))
                    checki = checki-1+np.argmax(MW.Xslice[checki-1:checki+2])
                    checkiy = int(sliceI[checki])
                    checkiy = checkiy-2+np.argmax(np.absolute(MW.X[checki,checkiy-2:checkiy+3]))
                    boxleft = checki-int(np.ceil(boxW/3))
                    boxright = checki+int(np.floor(boxW*2/3))
                    boxbott = checkiy-int(np.floor(boxH/2))
                    boxtop = checkiy+int(np.floor(boxH/2))
                    for j in range(1,boxW):
                        if np.absolute(MW.X[checki-j][checkiy])-AVG < BoxThresh*(np.absolute(MW.X[checki][checkiy])-AVG):
                            boxleft = checki-j
                            break
                        if np.absolute(MW.X[checki-j][checkiy])-AVG > 1.05*(np.absolute(MW.X[checki-j+1][checkiy])-AVG):
                            if np.absolute(MW.X[checki-j-1][checkiy])-AVG < 1.05*(np.absolute(MW.X[checki-j+1][checkiy])-AVG):
                                continue
                            else:
                                boxleft = checki-j+1
                                break
                    for j in range(1,boxW):
                        if np.absolute(MW.X[checki+j][checkiy])-AVG < BoxThresh*(np.absolute(MW.X[checki][checkiy])-AVG):
                            boxright = checki+j
                            break
                        if np.absolute(MW.X[checki+j][checkiy])-AVG > 1.05*(np.absolute(MW.X[checki+j-1][checkiy])-AVG):
                            if np.absolute(MW.X[checki+j+1][checkiy])-AVG < 1.05*(np.absolute(MW.X[checki+j-1][checkiy])-AVG):
                                continue
                            else:
                                boxright = checki+j-1
                                break
                    for j in range(1,boxH):
                        if np.absolute(MW.X[checki][checkiy-j])-AVG < BoxThresh*(np.absolute(MW.X[checki][checkiy])-AVG):
                            boxbott = checkiy-j
                            break
                    for j in range(1,boxH):
                        if np.absolute(MW.X[checki][checkiy+j])-AVG < BoxThresh*(np.absolute(MW.X[checki][checkiy+j])-AVG):
                            boxtop = checkiy+j
                            break
                    boxbott = checkiy - int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
                    boxtop = checkiy + int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
                    if Deltam < 10:
                        try:
                            beati = np.argmax(np.absolute(MW.X[checki,boxtop:checkiy+boxH]))
                            beattop = np.absolute(MW.X[checki][boxtop+beati])
                            beatbott = np.absolute(MW.X[checki][boxbott-beati])
                            if abs(beattop-beatbott) < 0.01*(np.absolute(MW.X[checki][checkiy])-AVG) and np.average((beattop,beatbott))-AVG > 0.2*(np.absolute(MW.X[checki][checkiy])-AVG):
                                for j in range(1,boxH):
                                    if np.absolute(MW.X[checki][boxbott-beati-j])-AVG < 0.2*(np.absolute(MW.X[checki][boxbott-beati])-AVG):
                                        boxbott = boxbott-beati-j
                                        break
                                for j in range(1,boxH):
                                    if np.absolute(MW.X[checki][boxtop+beati+j])-AVG < 0.2*(np.absolute(MW.X[checki][boxtop+beati])-AVG):
                                        boxtop = boxtop+beati+j
                                        break
                                boxbott = checkiy - int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
                                boxtop = checkiy + int(np.ceil(np.average((boxtop-checkiy,checkiy-boxbott))))
                        except ValueError:
                            pass
                    for j in range(2,topharm+1):
                        try:
                            if np.absolute(MW.X[checki][checkiy*j])-AVG < 0.01*(np.absolute(MW.X[checki][checkiy])-AVG):
                                break
                            else:
                                harmonic[i] = j
                        except IndexError:
                            break
                    if boxleft <= 0:
                        boxleft = 0
                    if boxright >= len(MW.xshort):
                        boxright = len(MW.xshort)-1
                    MW.rlist.append([MW.xshort[boxleft],MW.xshort[boxright],MW.yshort[boxbott],MW.yshort[boxtop],i])
                    MW.Ilist.append([boxleft,boxright,boxbott,boxtop,i])
                MW.Harmrec = int(np.round(np.average(harmonic)))

            #extracting signal            
            MW.glist = []
            MW.Xtemp = MW.X*0
            MW.yBL2 = np.zeros(len(MW.xint))+ min(MW.yint)+ 1
            MW.glistBL2 = []
            MW.xintBL2,MW.yintBL2,MW.ypaddBL2,MW.po2BL2,MW.xpaddBL2 = load.datapro(MW.xint,MW.yBL2)
            MW.XBL2 = STFT.stft(MW.ypaddBL2, MW.winnum, MW.window, MW.OF)
            MW.XtempBL2 = MW.XBL2*0
            slicetext = []
            for i in range(len(MW.Ilist)):
                MW.Xtemp[MW.Ilist[i][0]:MW.Ilist[i][1]+1,MW.Ilist[i][2]:MW.Ilist[i][3]+1] = MW.X[MW.Ilist[i][0]:MW.Ilist[i][1]+1,MW.Ilist[i][2]:MW.Ilist[i][3]+1]
                MW.Xtemp[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:-MW.Ilist[i][2]+1] = MW.X[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:-MW.Ilist[i][2]+1]
                MW.XtempBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,MW.Ilist[i][2]:MW.Ilist[i][3]+1] = MW.XBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,MW.Ilist[i][2]:MW.Ilist[i][3]+1]
                MW.XtempBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:-MW.Ilist[i][2]+1] = MW.XBL2[MW.Ilist[i][0]:MW.Ilist[i][1]+1,-MW.Ilist[i][3]:-MW.Ilist[i][2]+1]
                MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                MW.glist.append(MW.yrecon)
                slicetext.append(int(np.argmax(np.absolute(MW.Xtemp))/len(MW.Xtemp[0])))
                MW.Xtemp = MW.X*0
                MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                MW.glistBL2.append(MW.yreconBL2)
                MW.XtempBL2 = MW.XBL2*0
        
        if MW.autosearch == False:
            #adding results to plots
            MW.fig.clf()
            MW.ax = MW.fig.add_subplot(221)
            MW.ax4 = MW.fig.add_subplot(222, sharex=MW.ax)
            MW.ax3 = MW.fig.add_subplot(223, sharex=MW.ax)
            MW.ax2 = MW.fig.add_subplot(224)

            MW.ax3.plot(MW.x, MW.y)
            MW.ax3.set_title('Mass Spectrum')
            MW.ax3.set_xlabel('m/z')
            MW.ax3.set_ylabel('Abundance')

            MW.ax4.plot(MW.xshort,MW.Xslice,color='orange')
            MW.ax4.set_title('Gabor Spectrogram Slice \nRough Ion mass: '+str(np.round(mass,1))+' Da')
            MW.ax4.set_xlabel('m/z')
            MW.ax4.set_ylabel('Relative Amplitude')
            MW.ax4.set_ylim(0,1.05*np.max(MW.Xslice))
    
            MW.ax.imshow(abs(MW.X.T), origin='lower', aspect='auto',
                         interpolation='nearest', extent=(MW.xshort[0], MW.xshort[-1], MW.yshort[0], MW.yshort[-1]),
                         cmap='jet', vmax=MW.vmax)
            MW.ax.set_title('Gabor Spectrogram')
            MW.ax.set_xlabel('m/z')
            MW.ax.set_ylabel('Frequency')
            MW.ax.set_xlim(MW.xint[0], MW.xint[-1])
            MW.ax.set_ylim(MW.xFT[0]-MW.maxF/MW.winnum/2,MW.xFT[int(np.ceil(len(MW.xFT)/2))]+MW.maxF/MW.winnum/2)
            MW.ax.fill_between(slicex,slicey1,slicey2,color = 'orange',alpha=0.6)
            MW.toggle_selector.RS = RectangleSelector(MW.ax, self.onselect, drawtype='box', interactive=True)
            for i in range(len(MW.rlist)):
                rex = Rectangle((MW.rlist[i][0], MW.rlist[i][2]), (MW.rlist[i][1] - MW.rlist[i][0]), (MW.rlist[i][3] - MW.rlist[i][2]),
                                facecolor='none', edgecolor=color[MW.rlist[i][4]], linewidth=3)
                MW.ax.add_artist(rex)
                MW.ax.text(MW.rlist[i][0]+piW,MW.rlist[i][3]+piH,str(MW.zlist[i]),color=color[i],clip_on=True)
                MW.ax4.text(MW.xshort[slicetext[i]],MW.Xslice[slicetext[i]],str(MW.zlist[i]),color=color[i],clip_on=True,ha='center')
    
            MW.ax2.plot(MW.x, MW.y, color='darkgrey')
            if MW.realsel.isChecked() == False:
                for i in range(len(MW.glist)):
                    MW.ax2.plot(MW.xint, abs(MW.glist[i]), color=color[i])
                    MW.ax2.text(MW.xint[np.argmax(abs(MW.glist[i]))],np.max(abs(MW.glist[i])),str(MW.zlist[i]),color=color[i],clip_on=True,ha='center')
            if MW.realsel.isChecked() == True:
                for i in range(len(MW.glist)):
                    MW.ax2.plot(MW.xint, np.real(MW.glist[i]), color=color[i])
                    MW.ax2.text(MW.xint[np.argmax(np.real(MW.glist[i]))],np.max(np.real(MW.glist[i])),str(MW.zlist[i]),color=color[i],clip_on=True,ha='center')
            MW.ax2.set_title('Reconstructed Spectrum')
            MW.ax2.set_xlabel('m/z')
            MW.ax2.set_ylabel('Abundance')

            MW.ruffmass = mass
            MW.slicex,MW.slicey1,MW.slicey2 = slicex,slicey1,slicey2
            MW.piW,MW.piH = piW,piH
            MW.slicetext = slicetext
            MW.fig.tight_layout()
            MW.canvas.draw()
            self.dialog = guided_decharger(self)
            self.dialog.show()
        
        elif MW.autosearch == True:
            rtemp = []
            Itemp = []
            gtemp = []
            gBL2temp = []
            MW.Xtemp = MW.X*0
            MW.XtempBL2 = MW.XBL2*0
            for i in range(len(MW.zlist)):
                for j in range(0,MW.Harmrec+1):
                    if j == 0:
                        bottom = 0
                    else:
                        bottom = int(np.round(j*MW.Ilist[i][2]+(j-1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                    top = int(np.round(j*MW.Ilist[i][2]+(j+1)*(MW.Ilist[i][3]-MW.Ilist[i][2])/2))
                    rtemp.append([MW.rlist[i][0],MW.rlist[i][1],MW.yshort[bottom],MW.yshort[top],len(gtemp)])
                    Itemp.append([MW.Ilist[i][0],MW.Ilist[i][1],bottom,top,len(gtemp)])
                    if j == 0:
                        MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1]
                        MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])]
                        MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,0:Itemp[-1][3]+1]
                        MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:len(MW.X[0])]
                    else:
                        MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1]
                        MW.Xtemp[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1] = MW.X[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1]
                        MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,Itemp[-1][2]:Itemp[-1][3]+1]
                        MW.XtempBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1] = MW.XBL2[Itemp[-1][0]:Itemp[-1][1]+1,-Itemp[-1][3]:-Itemp[-1][2]+1]
                MW.yrecon = STFT.re_stft(MW.Xtemp,MW.winnum,MW.ypadd,MW.yint,MW.OF)*MW.correction_factor
                gtemp.append(MW.yrecon)
                MW.yreconBL2 = STFT.re_stft(MW.XtempBL2,MW.winnum,MW.ypaddBL2,MW.yintBL2,MW.OF)*MW.correction_factor
                gBL2temp.append(MW.yreconBL2)
                MW.Xtemp = MW.X*0
                MW.XtempBL2 = MW.XBL2*0
        
            MW.cs = MW.zlist
            MW.rlist = rtemp.copy()
            MW.Ilist = Itemp.copy()
            MW.glist = gtemp.copy()
            MW.glistBL2 = gBL2temp.copy()
            
            if MW.realsel.isChecked() == False:
                MW.tIFT = np.absolute(MW.glist)
            if MW.realsel.isChecked() == True:
                MW.tIFT = np.real(MW.glist)
            
            MW.xzero, MW.cszero, MW.zerofull = FT.zerocharge(MW.tIFT, MW.xint, MW.cs,MW.proton)
            MW.xzeroBL2, MW.cszeroBL2, MW.zerofullBL2 = FT.zerocharge(np.real(MW.glistBL2), MW.xint, MW.cs,MW.proton)
           
            massrange=np.nonzero(MW.zerofullBL2)
            MW.xzero = MW.xzero[massrange]
            MW.zerofull = MW.zerofull[massrange]
            MW.xzeroBL2 = MW.xzeroBL2[massrange]
            MW.zerofullBL2 = MW.zerofullBL2[massrange]
            MW.BSint = []
            for j in range(len(MW.cszero)):
                MW.cszero[j] = MW.cszero[j][massrange]
                MW.cszeroBL2[j] = MW.cszeroBL2[j][massrange]
                check = np.max(MW.cszeroBL2[j])/4
                if check == 0:
                    continue
                mid1 = int(len(MW.xzero)/8)+1
                mid2 = int(len(MW.xzero)*7/8)-1
                for k in range(len(MW.cszero[j])):
                    if MW.cszeroBL2[j][k] > check:
                        mid1 = k
                        break
                for k in range(mid1+1,len(MW.cszero[j])):
                    if MW.cszeroBL2[j][k] <= check or MW.cszeroBL2[j][k] <= MW.cszeroBL2[j][mid1]:
                        mid2 = k
                        break
                int1 = min(sum(MW.cszero[j][0:mid1]),sum(MW.cszero[j][mid2:len(MW.xzero)]))
                int2 = (sum(MW.cszeroBL2[j][0:mid1])+sum(MW.cszeroBL2[j][mid2:len(MW.xzero)]))/2
                MW.BSint.append(int1/int2)
            BSint = np.average(MW.BSint)
            MW.zerofullBL2 = MW.zerofullBL2*BSint
            MW.BL_corrected = (MW.zerofull - MW.zerofullBL2)/max(MW.zerofull)
            MW.zerofull = MW.zerofull - MW.zerofullBL2
    
            MW.glist = np.sum(MW.tIFT,axis=0)
    
            #save selection
            storetemp = []
            storetemp.append(MW.cs)
            storetemp.append(MW.rlist)
            storetemp.append(MW.glist)
            storetemp.append(MW.xzero)
            storetemp.append(MW.zerofull)
            storetemp.append(MW.cszero)
            storetemp.append(MW.cszeroBL2)
            storetemp.append(len(MW.glist))
            storetemp.append(MW.tIFT)
            MW.batchstore.append(storetemp)
        
                
    #### extra windows ######
    def harmonic(self):
        self.dialog = hamu(self)
        self.dialog.show()
        
    def defect(self):
        self.dialog = defect_analysis(self)
        self.dialog.show()
        
    def g_defect(self):
        self.dialog = g_defect_analysis(self)
        self.dialog.show()

    def inte(self):
        self.dialog = integrate(self)
        self.dialog.show()

    def STFT_iFAMS(self):
        MW.stft_selec_man = True
        MW.stft_box_thresh = "NA"
        MW.stft_range_thresh = "NA"
        try:
            del(MW.inclusion)
        except AttributeError:
            pass
        self.dialog = iFAMS_STFT(self)
        self.dialog.show()
     
    def cs_adj_menu(self):
        self.dialog = iFAMS_STFT(self)
        self.dialog.show()

    def STFT_parm_window(self):
        self.dialog = STFT_parm(self)
        self.dialog.show()

    def maxima_finder(self):
        self.dialog = Maxi(self)
        self.dialog.show()

    def man_input(self):
        self.dialog = manu(self)
        self.dialog.show()
        
    def iso_input(self):
        self.dialog = iso(self)
        self.dialog.show()
        
    def cali(self):
        self.dialog = Calibration(self)
        self.dialog.show()
    
    def cali2(self):
        try:
            MW.concydata
        except AttributeError:
            print('Load calibration curve first')
            MW.message = str('Please load calibration curve first')
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

        else:
            self.dialog = Unknown_calculator(self)
            self.dialog.show()
                
    def harmonic_finder(self):
        self.dialog = HarmonicFinder(self)
        self.dialog.show()
        
    def domainTool(self):
        self.dialog = Spec_domain(self)
        self.dialog.show()

    def noisemenu(self):
        self.dialog = Noise(self)
        self.dialog.show()   
        
    def auto_check_open(self):
        self.dialog = autosearch_check(self)
        self.dialog.show()
        
    def STFT_slicer(self):
        try:
            self.dialog = slice_viewer(self)
            self.dialog.show()
        except AttributeError:
            print('No Gabor data to view')
            MW.message = str('No Gabor data to view')
            MW.submessage = str('Please load STFT file and try again')
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()
        
    def batch_export_opt(self):
        try:
            self.batch_param()
            self.batch_save()
        except AttributeError:
            MW.message = str('No batch parameters to save. Please process data first')
            MW.submessage = str('Note that only STFT processing is supported \n'
                                'for batching currently.')
            print(MW.message)
            print(MW.submessage)
            MW.Err = True
            self.dialog = iFAMS_message(self)
            self.dialog.show()

#### extra windows ######
        
def my_excepthook(type, value, tback):
    sys.__excepthook__(type, value, tback)

sys.excepthook = my_excepthook
agilentbool = False
agilentpath = ""
usedFT = False
usedSTFT = False
RT1 = '0'
RT2 = '0'
RTfull = False
def main():
    app = QtWidgets.QApplication(sys.argv) #create iFAMS application
    form = MW()
    form.show() #show GUI on the screen
    app.exec_() #now execute the application, do what I said above

if __name__ == "__main__":
    main()
"""
find whats called main and run it
"""

