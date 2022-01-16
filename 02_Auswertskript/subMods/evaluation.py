# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 11:15:13 2021

@author: hbalsche
"""

import numpy as np
from scipy.spatial import ConvexHull

class pathGenerator:
    def __init__(self, data ):
        self.data  = data 

    
    def getCrossSlice(self, height, cutingWidth):
        z_min = height - 0.5*cutingWidth
        z_max = height + 0.5*cutingWidth
        
        secPoints_x = []
        secPoints_y = []
        secPoints_z = []
        
        self.points2D = []
        self.phi = []
        self.r = []
        for x, y, z  in zip(self.data[0],self.data[1],self.data[2]):
            if z_min <= z <= z_max:
                secPoints_x.append(x)
                secPoints_y.append(y)
                secPoints_z.append(z)
                self.phi.append(np.arctan2(x,y))
                self.r.append((x**2+y**2)**0.5)
                self.points2D.append((x,y))
                
        self.secPoints = [secPoints_x, secPoints_y, secPoints_z]
        
        self.sortToAngles()
        
        return self.secPoints

    def sortToAngles(self):
        sortedList = sorted(zip(self.phi, self.r, self.secPoints[0], self.secPoints[1]))
        phi_sorted, r_sorted, secPoints_x_sorted, secPoints_y_sorted = zip(*sortedList) 
        
        self.phi_sorted = list(phi_sorted)
        self.r_sorted =   list(r_sorted)
        self.secPoints_x_sorted = list(secPoints_x_sorted)
        self.secPoints_y_sorted = list(secPoints_y_sorted)

        
    def getHull(self):
        
        x =  self.secPoints_x_sorted
        x.append( self.secPoints_x_sorted[0])
        y = self.secPoints_y_sorted
        y.append(self.secPoints_y_sorted[0])
        return x,y
    
    def compMinMaxDiameterPoints(self, angularIncrement = np.pi/180):
        ''' sweep a rect along r-phi-pairs and find the verezes with maximum '''
        
        segmentCount = int(np.pi/(angularIncrement))
        
        max_diams = []
        max_ps = []
        min_diams = []
        min_ps = []
        for i_seg in range(segmentCount):

            startAngle = self.phi_sorted[0]+ i_seg*angularIncrement
            endAngle   = startAngle + angularIncrement
            
            r1_samples = []
            r2_samples = []
            p1_samples = []
            p2_samples = []
            for i_r, r in enumerate(self.r_sorted):
                phi = self.phi_sorted[i_r]
                if startAngle <= phi <= endAngle: 
                    r1_samples.append(r)
                    x1 = self.secPoints_x_sorted[i_r]
                    y1 = self.secPoints_y_sorted[i_r]
                    p1_samples.append([x1, y1])
                if (startAngle + np.pi) <= phi <= (endAngle + np.pi): 
                    r2_samples.append(r)
                    x2 = self.secPoints_x_sorted[i_r]
                    y2 = self.secPoints_y_sorted[i_r]
                    p2_samples.append([x2, y2])
            
            # get max/min samples
            r1_max = max(r1_samples)
            r1_min = min(r1_samples)
            r2_max = max(r2_samples)
            r2_min = min(r2_samples)
            
            # get associated points
            p1_max = p1_samples[r1_samples.index(r1_max)]
            p1_min = p1_samples[r1_samples.index(r1_min)]
            
            p2_max = p2_samples[r2_samples.index(r2_max)]
            p2_min = p2_samples[r2_samples.index(r2_min)]
            
            # get asociated "diameters" (dianeter means max chord legth)
            D_max = r1_max + r2_max
            D_min = r1_min + r2_min
            
            max_diams.append(D_max)
            max_ps.append([p1_max, p2_max])
            min_diams.append(D_min)            
            min_ps.append([p1_min, p2_min])
        
        max_diams, max_ps = self.__getSortetLists__(max_diams, max_ps)
        min_diams, min_ps = self.__getSortetLists__(min_diams, min_ps)
        
        pmin = min_ps[0]
        pmax = max_ps[-1]

        return pmin[0], pmin[1], pmax[0], pmax[1], D_min, D_max
    
    def getLongSlice(self, unwrapped,  angle, cuttingAngle, base = 'bestFit'):
        phi_uw = unwrapped['mesh'].x
        z_uw   = unwrapped['mesh'].y
        
        if base == 'bestFit':
            h_uw = unwrapped['h_reference']
        elif base == 'guess':
            h_uw = self.unwrapped['h_target']
        
        angle = np.radians(angle)
        cuttingAngle = np.radians(cuttingAngle)
        
        phi_min = angle - 0.5*cuttingAngle
        phi_max = angle + 0.5*cuttingAngle

        vec_z = []
        vec_phi = []
        vec_h = []
        

        for  z, phi, h  in zip(z_uw, phi_uw, h_uw):
            if phi_min <= phi <= phi_max:
                vec_z.append(z)
                vec_phi.append(phi)
                vec_h.append(h)
        
        
        sorted_pairs  = sorted(zip(vec_z, vec_phi, vec_h))
        tuples  = zip(*sorted_pairs)
        vec_z, vec_phi, vec_h = [ list(tuple) for tuple in  tuples]
        
        
        return np.array(vec_z), np.array(vec_phi), np.array(vec_h)
        
        
          
    def __getSortetLists__(self, list1, list2):
        '''returns sortet lists based on list 1'''     
        
        sorted_pairs  = sorted(zip(list1, list2))
        tuples  = zip(*sorted_pairs)
        list1, list2 = [ list(tuple) for tuple in  tuples]
        
        return list1, list2
    
    def getMaxChord(self, x, y):
        points = [ (xi, yi) for (xi, yi) in zip(x,y) ]
        hull = ConvexHull(points)
    
        hullPoints = [points[i] for i in hull.vertices]
        
        L_max = 0 
        i1, i2 = 0, 0 
        for i, pi in enumerate(hullPoints[:-1]):
            pk = hullPoints[i+1]
            xi, yi = pi[0], pi[1]
            xk, yk = pk[0], pk[1] 
            L_act = 0.5*( (xk-xi)**2 + (yk-yi)**2 )**0.5
            if L_act > L_max:
                L_max = L_act
                i1, i2 = i, i+1
        
        pi = hullPoints[i1]
        pk = hullPoints[i2]
        xi, yi = pi[0], pi[1]
        xk, yk = pk[0], pk[1] 
        
        return [xi, xk], [yi, yk]


    def getImpAmp(self, x, y, xi, xk, yi, yk):
        L = 0 
        i_act = None 
        
        points = [ (z, h) for (z, h) in zip(x,y) ]
        
        xlow  = min(xi, xk)
        xhigh = max(xi, xk)
        points = [p for p in points if xlow <= p[0] <= xhigh]
    
        for i, point in enumerate(points):
            p1 = np.array([xi, yi])
            p2 = np.array([xk, yk])
            p3 = np.array(point)
            L_act = np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1)
            if L_act > L:
                L = L_act
                i_act = i   
                
        x_dent = points[i_act][0]
        y_dent = points[i_act][1]
       
        y_inter = (yk - yi)/(xk-xi)*(x_dent-xi)+yi
        
        Lc  = ((xk - xi)**2 + (yk -  yi)**2)**0.5
            
        return Lc, L, [x_dent, x_dent], [y_dent, y_inter]

            
