# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 16:36:21 2021

@author: hbalsche
"""
import numpy as np                  # for math
import matplotlib.tri as mesher     # for tringulation
from scipy.optimize import leastsq  # for fitting
from stl import mesh                # for stl handling
       

class fitter():
    
    def __init__(self, data = None):
        self.data = data
        self.fitted = False
        self.fittingParams = None
        self.unwrapped = None
            
    
    def getRandomTestzyl(self, D_target, delL, nPhi, nL, fTol):
        ''' creates a point cloud of test cylinder with given dimensions and tolerances which 
            is randomly rotated and got random radial errors'''
            
        pointsX = []
        pointsY = []
        pointsZ = []
        points = []
        actL = 0
        delPhi = 2*np.pi/nPhi
        R = D_target/2.0
        for i in range(nL):
            for j in range(nPhi):
                phi = j *delPhi
                delR = fTol*R*(np.random.rand(1)[0]-0.5)
                x   = (R+delR)*np.cos(phi)
                y   = (R+delR)*np.sin(phi)
                points.append(np.array([x,y,actL]))
            actL += delL
        
        angles = np.random.rand(3)*0.4*np.pi
        
        ca = np.cos(angles[0])
        sa = np.sin(angles[0])
        cb = np.cos(angles[1])
        sb = np.sin(angles[1])
        
        Tx = np.array([[  1,  0,  0],
                       [  0, ca,-sa],
                       [  0, sa, ca]])
        
        Ty = np.array([[ cb,  0, sb],
                       [  0,  1,  0],
                       [-sb,  0, cb]])
        
        for p in points: 
            p = np.dot(p,Tx)
            p = np.dot(p,Ty)
            pointsX.append(p[0])
            pointsY.append(p[1])
            pointsZ.append(p[2])
        
        self.data = [pointsX, pointsY, pointsZ]
        
    def getWavedTestzyl(self, D_target, delL, nPhi, nL, fTol):
        ''' creates a point cloud of test cylinder with given dimensions and tolerances whch 
            is randomly rotated and got random radial errors'''
            
        pointsX = []
        pointsY = []
        pointsZ = []
        points = []
        actL = 0
        delPhi = 2*np.pi/nPhi
        R = D_target/2.0
        for i in range(nL):
            for j in range(nPhi):
                phi = j *delPhi
                delR = fTol*np.sin(phi/2)*np.cos(2*np.pi*i*delL/(nL*delL))
                x   = (R+delR)*np.cos(phi)
                y   = (R+delR)*np.sin(phi)
                points.append(np.array([x,y,actL]))
            actL += delL
        
        angles = np.random.rand(3)*0.4*np.pi
        
        ca = np.cos(angles[0])
        sa = np.sin(angles[0])
        cb = np.cos(angles[1])
        sb = np.sin(angles[1])
        
        Tx = np.array([[  1,  0,  0],
                       [  0, ca,-sa],
                       [  0, sa, ca]])
        
        Ty = np.array([[ cb,  0, sb],
                       [  0,  1,  0],
                       [-sb,  0, cb]])
        
        for p in points: 
            p = np.dot(p,Tx)
            p = np.dot(p,Ty)
            pointsX.append(p[0])
            pointsY.append(p[1])
            pointsZ.append(p[2])
        
        self.data = [pointsX, pointsY, pointsZ]
            
    def computeBestFitZylinder(self, D_guess, p0_guess, d_guess, rot_z_deg = 0, reverse = False):
        ''' computes a best fit cylinder for the given messh based on best guess for 
            Diameter D_guess, and axis defined by point on axis p0_guess and 
            axis direction vector d_guss. 
            an arbitrary rotation around cylinder axis may be aplied via rot_z_deg'''
            
        print('###')
        print('### computing best fit ###')
        print('###')
        
        self.D_guess = D_guess
        # nestet error function for point on cylinder 
        def errorAgainstCyl(params, xp, yp, zp):
            # sarting point
            x0 = params[0]
            y0 = params[1]
            z0 = params[2]
            # direction of axis
            dx = params[3]
            dy = params[4]
            dz = params[5]
            # Rasius
            R = params[6]
     
            cx = xp-x0
            cy = yp-y0
            cz = zp-z0
            
            abs_c_x_d_sq = (cy*dz-cz*dy)**2.0+(cz*dx-cx*dz)**2.0+(cx*dy-cy*dx)**2.0
            abs_d_sq = dx**2.0+dy**2.0+dz**2.0
            
            error = abs_c_x_d_sq/abs_d_sq - R**2
            return error
        
        # pass initial parameter
        pGuess= []
        pGuess.append(p0_guess[0])
        pGuess.append(p0_guess[1])
        pGuess.append(p0_guess[2])
        pGuess.append(d_guess[0])
        pGuess.append(d_guess[1])
        pGuess.append(d_guess[2])
        pGuess.append(D_guess/2.0)
        # measured values 
        xp = self.data[0]
        yp = self.data[1]
        zp = self.data[2]
        # perform fitting
        est_Params, success = leastsq(errorAgainstCyl, pGuess, args=(xp, yp, zp), maxfev=5000)

        self.fittingParams = est_Params
        self.fitted = True 
        self.alignToZ(rot_z_deg, reverse)
    
    def alignToZ(self, rot_z_deg, reverse):
        
        print('###')
        print('### start aligning ###')
        print('###')
        
        
        xd = self.fittingParams[3]
        yd = self.fittingParams[4]
        zd = self.fittingParams[5]
        Ld = (xd**2+yd**2+zd**2)**0.5
        pd = np.array([xd,yd,zd])*(1/Ld)
        
        # first orthogonal vector on d
        x1 = -(yd+zd)/xd
        y1 = 1
        z1 = 1
        L1 = (x1**2+y1**2+z1**2)**0.5
        n1 = np.array([x1,y1,z1])*(1/L1)
        
        # 2nd orthogonal vector on d
        n2 = np.cross(n1,pd)
        x2 = n2[0]
        y2 = n2[1]
        z2 = n2[2]
        L2 = (x2**2+y2**2+z2**2)**0.5
        n2 *= 1/L2
        
        T = np.array([[pd[0], n1[0], n2[0]],
                      [pd[1], n1[1], n2[1]],
                      [pd[2], n1[2], n2[2]]])
        
        revFact = 1.0
        if reverse: 
            revFact = -1.0
        
        cb = np.cos( revFact*0.5*np.pi)
        sb = np.sin( revFact*0.5*np.pi)
        Ty = np.array([[ cb,  0, sb],
                       [  0,  1,  0],
                       [-sb,  0, cb]])
        

        pointsX = []
        pointsY = []
        pointsZ = []
        for i in range(len(self.data[0])):
            
            #create point object from data 
            x = self.data[0][i]
            y = self.data[1][i]
            z = self.data[2][i]
            p = np.array([x,y,z])
            
            # apply rotations 
            p = np.dot(p,T)
            p = np.dot(p,Ty)
            
            pointsX.append(p[0])
            pointsY.append(p[1])
            pointsZ.append(p[2])
        
        self.data = [pointsX, pointsY, pointsZ]

        # apply rotation to direction vector: 
        d = np.array([self.fittingParams[3],
                      self.fittingParams[4],
                      self.fittingParams[5]])       
        d = np.dot(d,T)
        d = np.dot(d,Ty)
        
        self.fittingParams[3] = d[0]
        self.fittingParams[4] = d[1]
        self.fittingParams[5] = d[2]
        
        # apply rotation to axis support point: 
        p0 = np.array([self.fittingParams[0],
                       self.fittingParams[1],
                       self.fittingParams[2]])       
        
        p0 = np.dot(p0,T)
        p0 = np.dot(p0,Ty)
        
        self.fittingParams[0] = p0[0]
        self.fittingParams[1] = p0[1]
        self.fittingParams[2] = p0[2]
        
        
        self.center()
        
        # rotate about z
        cg = np.cos(np.deg2rad(rot_z_deg))
        sg = np.sin(np.deg2rad(rot_z_deg))
        Tz = np.array([[ cg,-sg,  0],
                       [ sg, cg,  0],
                       [  0,  0,  1]])
        
        pointsX = []
        pointsY = []
        pointsZ = []
        for i in range(len(self.data[0])):
            
            #create point object from data 
            x = self.data[0][i]
            y = self.data[1][i]
            z = self.data[2][i]
            p = np.array([x,y,z])
            
            # apply rotations 
            p = np.dot(p,Tz)
            
            pointsX.append(p[0])
            pointsY.append(p[1])
            pointsZ.append(p[2])
            
        self.data = [pointsX, pointsY, pointsZ]
        
        self.unwrap()
    
    
    def center(self):
        
        xp0 = self.fittingParams[0]
        yp0 = self.fittingParams[1]
        zp0 = self.fittingParams[2]       
    
        off_x = xp0
        off_y = yp0
        off_z = min(self.data[2])
        
        for i in range(len(self.data[0])):
            # add some arbitrary offset to axes for testing
            self.data[0][i] -= off_x
            self.data[1][i] -= off_y
            self.data[2][i] -= off_z
        # also apply offsets to axis
        self.fittingParams[0] -= off_x
        self.fittingParams[1] -= off_y
        self.fittingParams[2] -= off_z
            
    
    def unwrap(self):
        # take all points transform them to cylinder coords
        # cast point data to 
        
        print('###')
        print('### start unwrapping surface ###')
        print('###') 
        
        z = []
        phi = []
        h_ref = []
        h_tar = [] 
        R_ref = self.fittingParams[6]
        R_tar = self.D_guess/2.0
        
        for i in range(len(self.data[0])):
            
            x_c = self.data[0][i]
            y_c = self.data[1][i]
            z_c = self.data[2][i]
            
            phi_p = np.arctan2(x_c,y_c)
            r = (x_c**2 + y_c**2)**0.5
            
            z.append(z_c)
            phi.append(phi_p)
            h_ref.append(r-R_ref)
            h_tar.append(r-R_tar)
        
        tria = mesher.Triangulation(phi, z)
        
        self.unwrapped  = {
            'mesh'        : tria,
            'h_reference' : h_ref,
            'h_target'    : h_tar
            }

    
    def readStl(self, pathToSTL):
        
        print('###')
        print('### start reading stl file ###')
        print('###')
        
        givenMesh = mesh.Mesh.from_file(pathToSTL)
              
        # get coordinates of scanned cylinder 
        vertices = givenMesh.vectors.reshape(3*len(givenMesh.vectors), 3)
        vertices = np.unique(vertices, axis=0)
        x, y, z = zip(*vertices)
                
        self.data = [x, y, z]
    
    def getPlottingInputs(self):
        return self.fittingParams, self.D_guess, self.data, self.unwrapped
    
    def removePointsBetween(self, z_min, z_max):
        px, py, pz = [], [], [] 
        for x, y, z  in zip(self.data[0],self.data[1],self.data[2]):
            if z_min <= z <= z_max:
                px.append(x)
                py.append(y)
                pz.append(z)
                
        self.data = [px, py, pz]    
        self.unwrap()
      