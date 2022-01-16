# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 10:36:50 2021

@author: hbalsche
"""

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import gridspec
from scipy.interpolate import griddata
from . import evaluation


class plotter:
    def __init__(self, plottinginputs ):
        self.fittingParams, self.D_guess, self.data, self.unwrapped = plottinginputs
    
    def plotFitting(self, base = 'bestFit', minTol = None, maxTol = None, 
                    colormap = None, nStep = 8, fPath = '.'):
        
        print('###')
        print('### start plot of fitted cylinder ###')
        print('###') 
        
        
        if base == 'bestFit':
            h = self.unwrapped['h_reference']
        elif base == 'guess':
            h = self.unwrapped['h_target']         
        
        if minTol == None:
            minTol = min(h)
        
        if maxTol == None:
            maxTol = max(h)
            
        if colormap == None:
            cmap = plt.get_cmap('Spectral', nStep)
        else:
            cmap = plt.get_cmap(colormap, nStep)
        
        xma = max(self.data[0])
        xmi = min(self.data[0])
        yma = max(self.data[1])
        ymi = min(self.data[1])
        zma = max(self.data[2])
        zmi = min(self.data[2])
        
        skirtX = abs(xma - xmi)*0.1
        skirtY = abs(xma - xmi)*0.1
        skirtZ = abs(xma - xmi)*0.1
        
        fig = plt.figure(figsize=(16,8))
        
        ax1 = fig.add_subplot(131, projection='3d')
        ax1.scatter(self.data[0],
                    self.data[1],
                    self.data[2],
                    marker='.',
                    c = h, alpha = 0.1,
                    vmin = minTol, 
                    vmax = maxTol, 
                    cmap = cmap)
        ax1.set_title('Pointcloud 3D')
        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel('y [mm]')
        ax1.set_zlabel('z [mm]')
        ax1.set_xlim(right = xma+skirtX,left = xmi-skirtX)
        ax1.set_ylim(top = yma++skirtY,bottom = ymi-+skirtY)
        ax1.set_zlim(top = zma,bottom = zmi)
        ax2 = fig.add_subplot(132)
        ax2.set_aspect('equal')
        ax2.scatter(self.data[1],
                    self.data[2],
                    marker='.',  
                    c = h, alpha = 0.1,
                    vmin = minTol, 
                    vmax = maxTol, 
                    cmap = cmap)
        ax2.set_title('Elevation Right')
        ax2.set_xlim(xmax = yma + skirtY,xmin = ymi-  skirtY)
        ax2.set_ylim(ymax = zma + skirtZ,ymin = zmi- skirtZ)
        ax2.set_xlabel('y [mm]')
        ax2.set_ylabel('z [mm]')
        ax2.grid()
        ax3 = fig.add_subplot(133)
        ax3.set_aspect('equal')
        ax3.scatter(self.data[0],self.data[2],
                    marker='.', 
                    c = h, alpha = 0.1,
                    vmin = minTol, 
                    vmax = maxTol, 
                    cmap = cmap)
        ax3.set_xlim(xmax = xma +  skirtX, xmin = xmi -  skirtX)
        ax3.set_ylim(ymax = zma -  skirtZ, ymin = zmi +  skirtZ)
        ax3.set_title('Elevation Front')
        ax3.grid()
        ax3.set_xlabel('x [mm]')
        ax3.set_ylabel('z [mm]')
        fileName = str(fPath) + '/aligned.png'
        plt.savefig(fileName)
        plt.show()
    
    
    def plotUnwrapped(self, base = 'bestFit', minTol = None, maxTol = None, 
                      colormap = None, nStep = 8, reverseSamples = False,
                      fPath = '.'):
        
        print('###')
        print('### start plotting unwrapped shape ###')
        print('###')
                
        if base == 'bestFit':
            h = self.unwrapped['h_reference']
            r = self.fittingParams[6]
        elif base == 'guess':
            h = self.unwrapped['h_target']
            r = self.D_guess/2.0
        
        
        if minTol == None:
            minTol = min(h)
        
        if maxTol == None:
            maxTol = max(h)
            
        if colormap == None:
            cmap = plt.get_cmap('Spectral', nStep)
        else:
            cmap = plt.get_cmap(colormap, nStep)
            
        fig = plt.figure(figsize=(16,8))
        gs = gridspec.GridSpec(1, 3, width_ratios=[10, 6, 1]) 
        ax1 = fig.add_subplot(gs[0], projection='3d')
        ax1.set_xticks([-np.pi*r, -0.5*np.pi*r, 0, 0.5*np.pi*r, np.pi*r] )
        ax1.set_xticklabels(['-180°','-90°','0°','90°','180°'])
        if reverseSamples:
            ax1.scatter((self.unwrapped['mesh'].x*r)[::-1],
                        (self.unwrapped['mesh'].y)[::-1], 
                         h[::-1],
                         marker='.',
                         c = h[::-1], alpha = 0.1,
                         vmin = minTol, 
                         vmax = maxTol, 
                         cmap = cmap)
        else:
            ax1.scatter((self.unwrapped['mesh'].x*r),
                        (self.unwrapped['mesh'].y), 
                         h,
                         marker='.',
                         c = h, alpha = 0.1,
                         vmin = minTol, 
                         vmax = maxTol, 
                         cmap = cmap)
        
        ax1.set_zlabel('$\Delta r$ [mm]')
        ax1.set_ylabel('z [m]')
        ax1.set_xlabel('winding angle')
        ax2 = fig.add_subplot(gs[1])
        ax2.set_aspect('equal')
        ax2.set_xticks([-np.pi*r, -0.5*np.pi*r, 0, 0.5*np.pi*r, np.pi*r] )
        ax2.set_xticklabels(['-180°','-90°','0°','90°','180°'])
        ax2.set_ylabel('z [mm]')
        ax2.set_xlabel('winding angle')
        ax2.set_title('Unwrapped 3D-Scan')
        ax2.grid()
        cnt2 = ax2.scatter(self.unwrapped['mesh'].x*r,
                    self.unwrapped['mesh'].y,
                    marker='.',
                    c = h, vmin = minTol, 
                    vmax = maxTol, 
                    cmap = cmap)
        
        ax3 =  fig.add_subplot(gs[2])     
        ax3.set_title('Legend')
        cb = fig.colorbar(cnt2, ax3)
        cb.set_ticks(np.linspace(minTol, maxTol, nStep+1))
        ax3.set_ylabel('$\Delta r$ [mm]', fontsize = 8)
        fileName = str(fPath) + '/unwrapped.png'
        plt.savefig(fileName, dpi=300)
        plt.show()
    
    def plotCrossSection(self, d_nom = 0, height = 0, cuttingWidth = 10, fPath = '.'):
        pathGen = evaluation.pathGenerator(self.data)
        crosSection = pathGen.getCrossSlice(height, cuttingWidth)
        

        fig = plt.figure(figsize = (12,8))
        gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
        ax1 = fig.add_subplot(gs[0])
        ax1.set_aspect('equal')
        # plot d_nom 
        if d_nom > 0:
            phi = np.linspace(0, 2*np.pi, 50)
            r = d_nom/2
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            ax1.plot(x, y, linestyle = 'dashed', color = 'black')
        # plot x-y data 
        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel('x [mm]')
        # get polygon 
        x,y = pathGen.getHull()
        ax1.plot(x,y, linewidth = 1.5,  color = 'black')
        # get cutting planes for min and max d
        p1_min, p2_min, p1_max, p2_max, D_min, D_max = pathGen.compMinMaxDiameterPoints(np.pi/20)
        x_min = [p1_min[0], p2_min[0]]
        y_min = [p1_min[1], p2_min[1]]
        ax1.plot(x_min,y_min, linewidth = 1,  color = 'black')
        x_max = [p1_max[0], p2_max[0]]
        y_max = [p1_max[1], p2_max[1]]
        ax1.plot(x_max,y_max, linewidth = 1,  color = 'black')
        self.plotLineAnnot(ax1, p1_max, p2_max, '$D_{max}$', 0.3)
        self.plotLineAnnot(ax1, p1_min, p2_min, '$D_{min}$', 0.3)
        ax1.grid()
        # plot annotations for out of roundness
        ax2 = fig.add_subplot(gs[1])
        # hide frame
        ax2.spines['top'].set_color('none')
        ax2.spines['right'].set_color('none')
        ax2.spines['bottom'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.set_xticks([])
        ax2.set_yticks([])
        bbox_props = dict(boxstyle="round", fc="w", ec="black", alpha=1)
        U_r = (D_max-D_min)/d_nom
        text = ('$D_{nom} = $' + str(round(d_nom,1))+ ' mm' +
                '\n$D_{min} = $' + str(round(D_min,1))+ ' mm' +
                '\n$D_{max} = $' + str(round(D_max,))+ ' mm' +
                '\n '+
                '\n$U_{r}$ = ' + str(round(U_r,3)))
        ax2.text(0, 0.8, text, ha="center", va="center",  size=10, bbox=bbox_props)
        ax1.set_title('Cross Section at h = ' + str(height) + 'mm')
        fileName = str(fPath) + '/outOfRoundness_' + str(height) + '.svg'
        plt.savefig(fileName, dpi=300)
        plt.show()
    
    def plotLineAnnot(self, ax, p_start, p_end, text, units):
        x_start = p_start[0]
        y_start = p_start[1]
        x_end   = p_end[0]
        y_end   = p_end[1]
        
        hx = x_end - x_start    
        hy = y_end - y_start
        L = (hx**2 + hy**2)**0.5
        m = hy/hx
        Ln = (m**2+1)**0.5
        
        dx = units*L*1/Ln 
        dy = units*L*m/Ln 
        
        x_pl = x_start + dx
        y_pl = y_start + dy
        bbox_props = dict(boxstyle="round", ec = 'black',  fc="w", alpha=1)
        ax.text(x_pl, y_pl, text, ha="center", va="center",  size=10, bbox=bbox_props)
    

    def plotLongSection(self, angle = 180, amplification = 100,  
                        cuttingAngle = 1, base = 'bestFit', fPath = '.'):
        pathGen = evaluation.pathGenerator(self.data)
        vec_z, vec_phi, vec_h = pathGen.getLongSlice(self.unwrapped, 
                                                     angle, 
                                                     cuttingAngle, 
                                                     base = 'bestFit')
        
        chord_z, chord_h = pathGen.getMaxChord(vec_z, vec_h)
        chord_z, chord_h = np.array(chord_z), np.array(chord_h)
        
        xi, xk = chord_z[0], chord_z[1]
        yi, yk = chord_h[0], chord_h[1]
        L, dw, z_perp, h_perp = pathGen.getImpAmp(vec_z, vec_h, xi, xk, yi, yk)
        z_perp, h_perp = np.array(z_perp), np.array(h_perp)
        
        exp = np.log10(1/amplification)
        
        fig = plt.figure(figsize = (4,10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.9, 0.1]) 
        ax1 = fig.add_subplot(gs[0])
        ax1.plot(h_perp*amplification, z_perp, 
                    color = 'red',
                    linewidth = 0.5)
        ax1.plot(chord_h*amplification, chord_z, 
                    color = 'red',
                    linewidth = 0.5)
        ax1.plot(vec_h*amplification, vec_z, 
                    color = 'black',
                    linewidth = 1.5,
                    label = 'scanned surface')
        ax1.set_xlabel('$\Delta r~[mm \cdot 10^{' + str(int(exp)) + '}]$')
        ax1.set_ylabel('z [mm]') 
        ax1.grid()
        ax1.set_title('at' + str(angle) + ' °')
        
        # plot annotations for out of roundness
        ax2 = fig.add_subplot(gs[1])
        # hide frame
        ax2.spines['top'].set_color('none')
        ax2.spines['right'].set_color('none')
        ax2.spines['bottom'].set_color('none')
        ax2.spines['left'].set_color('none')
        ax2.set_xticks([])
        ax2.set_yticks([])
        bbox_props = dict(boxstyle="round", fc="w", ec="black", alpha=1)
        U_x = dw/L
        text = (r'$L_{chord} = ' + str(round(L,1))+ '~mm $ \n' +
                r'$\Delta w_{max} = ' + str(round(dw,1))+ '~mm $  \n' +
                '\n '+
                r'$U_{x} = ' + str(round(U_x,3)) + '$')
        ax2.text(0.5, 0.3, text, ha="center", va="center",  size=10, bbox=bbox_props)
        
        fileName = str(fPath) + '/longDeviat_' + str(angle) + '.svg'
        plt.savefig(fileName, dpi=300)
        plt.show()
    
    def SurfPlotUnwrapped(self, nPhi = 10, nz = 10, base = 'bestFit'):
        
        print('###')
        print('### start plotting unwrapped shape (surface) ###')
        print('###')
        
                
        if base == 'bestFit':
            h = self.unwrapped['h_reference']
            r = self.fittingParams[6]
        elif base == 'guess':
            h = self.unwrapped['h_target']
            r = self.D_guess/2.0
        
        
        # define regular grid to interpolate scanned data on 
        z_max = max(self.unwrapped['mesh'].y)-5
        z_min = min(self.unwrapped['mesh'].y)+5
        
        x_max = max(self.unwrapped['mesh'].x)*r*0.98
        x_min = min(self.unwrapped['mesh'].x)*r*0.98
        
        px, py = np.linspace(x_min, x_max, nPhi), np.linspace(z_min, z_max, nz)
        px, py = np.meshgrid(px, py)
        
        # interpolate
        x_scan = self.unwrapped['mesh'].x*r
        y_scan = self.unwrapped['mesh'].y
        z_scan = np.array(h)
        
        pz = griddata((x_scan, y_scan), z_scan, (px, py), method = 'cubic')
        
        
        fig = plt.figure(figsize = (16,8))
        ax1 = fig.add_subplot(111, projection = '3d')
        ax1.plot_surface(px,py,pz, cmap = 'Spectral')       # interpolated scan
        plt.show()
        


        