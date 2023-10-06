"""
Class Header
"""

import numpy as np
from scipy.special import erf  # we will need this


class Filtered_ADM(object): 
    """
    Filtered ADM testing class before implementation in Fortran. 
    """

    def __init__(self, 
                 xLine, yLine, zLine, 
                 xLoc=0., yLoc=0., zLoc=0., 
                 CT=1.33, diam=1., delta=None, alpha=None, 
                 thickness=1.5):
        
        # initialize a bunch of variables
        self.xLine = xLine
        self.yLine = yLine
        self.zLine = zLine
        self.xLoc = xLoc
        self.yLoc = yLoc
        self.zLoc = zLoc
        self.dx = xLine[1]-xLine[0]
        self.dy = yLine[1]-yLine[0]
        self.dz = zLine[1]-zLine[0]
        self.nx = len(xLine)
        self.ny = len(yLine)
        self.nz = len(zLine)

        self.xG, self.yG, self.zG = np.meshgrid(xLine, yLine, zLine, indexing='ij')
        self.h = np.sqrt(self.dx**2 + self.dy**2 + self.dz**2)
        
        if alpha is None: 
            if delta is None: 
                self.delta = self.h
            else: 
                self.delta = delta
        else: 
            self.delta = alpha * self.h

        self.R = diam/2
        self.s = thickness*self.dx  # thickness

        self.tag_face = np.zeros((self.ny, self.nz)).astype(bool)
        self.tag_face[self.yG[0, :, :]**2 + self.zG[0, :, :]**2 < self.R**2] = 1

        self.npts = np.sum(self.tag_face)

        self.CT = CT
        
    
    def R1(self, ndarray=False): 
        """
        ADM weighting function in cylindrical coordinates. 
        This function is a convolution of a spherical gaussian kernel and 
        the mathematical definition of a finite-thickness actuator disk model. 

        R1 is only the spanwise-weighting component (see R2(), R_xyz())

        Parameters
        ----------
        ndarray (bool) : if True, returns an [nx x ny x nz] 3D array instead. 

        Returns
        -------
        R1(x) (array) : [nx x 1] array of weights; integrates to unity.  
        """

        # duplicate variable pointers for readibility
        delta = self.delta
        s = self.s
        
        if ndarray: 
            x = self.xG
        else: 
            x = self.xLine

        R1 = 1/(2*s) * (erf(np.sqrt(6)/delta*(x+s/2)) - erf(np.sqrt(6)/delta*(x-s/2)))
        return R1

    
    def set_delta(self, delta=None, alpha=None): 
        """
        Resets self.delta parameter (filter smoothing)
        """

        if alpha is None: 
            if delta is None: 
                self.delta = np.sqrt(self.dx**2 + self.dy**2 + self.dz**2)
            else: 
                self.delta = delta
        else: 
            self.delta = alpha * np.sqrt(self.dx**2 + self.dy**2 + self.dz**2)


    def R2(self): 
        """
        ADM weighting function in cylindrical coordinates. 

        R2 is only the radial component, computed via numerical integration of the Green's function. 

        Returns 
        -------
        R2(r) (array) : [ny x nz] array of weights; integrates to unity. 
        """
        
        # duplicate pointers for readability
        R = self.R
        yy = self.yG[0, :, :]
        zz = self.zG[0, :, :]
        delta = self.delta
        B = 6 / (np.pi**2 * R**2 * delta**2)  # normalizing factor

        yFace = yy[self.tag_face]
        zFace = zz[self.tag_face]

        R2 = np.zeros(yy.shape)

        for y, z in zip(yFace, zFace):  # loop through and stamp in the radial gaussian
            R2 += np.exp(-6 * ((y-yy)**2 + (z-zz)**2)/ delta**2) 

        return R2 * B * self.dy * self.dz

    def R_xyz(self): 
        """
        ADM weighting function in three dimensions. 

        Returns
        -------
        R(r, x) (array) : [nx x ny x nz] array of weights, integrates to unity. 
        """

        return self.R1(ndarray=False)[:, np.newaxis, np.newaxis] * self.R2()[np.newaxis, :, :]

    def calc_ud(self): 
        """
        Calculates the disk velocity based on the vortex cylinder model (eq. 18)

        Returns
        -------
        ud (float) : disk velocity as given by the vortex cylinder model
        """

        R = self.R
        CT = self.CT

        integral = np.sum(self.R2()**2) * self.dy * self.dz
        ud = (1 + np.pi * R**2 * CT/4 * integral)

        return 1./ud


    def numerical_M(self): 
        """
        Calculates correction factor M numerically from equation (23). 
        """
        
        R = self.R
        CT = self.CT
        
        integral = np.sum(self.R2()**2) * self.dy * self.dz

        M_inv = 1 + CT/4 * (1 - np.pi * R**2 * integral)
        return 1./M_inv


    def approx_M(self, delta=None, CT=None): 
        """
        Taylor series approximation to M from equation (26)

        Parameters
        ----------
        delta (array) : if provided, uses the given delta instead of self.delta
        CT (array) : if provided, uses the given CT' value instead of self.CT
        """

        if delta is None: 
            delta = self.delta
            
        if CT is None: 
            CT = self.CT
        
        return 1 / (1 + CT/4 * delta/self.R/np.sqrt(3*np.pi))

