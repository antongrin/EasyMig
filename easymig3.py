# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 13:21:45 2016

@author: GrinevskiyAS
"""

from __future__ import division
import numpy as np
from numpy import sin,cos,tan,pi,sqrt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

#поддержка кариллицы
font = {'family': 'Arial', 'weight': 'normal', 'size':20}
mpl.rc('font', **font)



class Hyperbola:
    def __init__(self, xarray, tarray, x0, v, t0):
        self.x=xarray
        self.x0=x0
        self.v=v
        self.t0=t0
        self.h=t0*v/2
        self.t=sqrt(t0**2 + (2*(xarray-x0)/v)**2)
        xstep=xarray[1]-xarray[0]
        tbegin=tarray[0]
        tend=tarray[-1]
        tstep=tarray[1]-tarray[0]
        
                
        self.x=self.x[ (self.t>=tbegin) & (self.t <= tend) ]        
        self.t=self.t[ (self.t>=tbegin) & (self.t <= tend) ]
        self.amp = 1/(self.t/self.t0)
        self.imgind=((self.x-xarray[0])/xstep).astype(int)       
        self.grid_resample(tarray)
        
        
    def grid_resample(self, tarray):
        tbegin=tarray[0]
        tend=tarray[-1]
        tstep=tarray[1]-tarray[0]
        
        
        self.tind=np.round((self.t-tbegin)/tstep).astype(int)
        self.tind=self.tind[self.tind*tstep<=tarray[-1]]
        self.tgrid=tarray[self.tind]
        self.coord=np.vstack((self.x,tarray[self.tind]))
                
        
    def disp(self):
        fg=plt.figure(figsize=(12,10), facecolor='w')
        ax=fg.add_subplot(111)  
        ax.xaxis.set_ticks_position('top')
        a=ax.plot(self.x,self.t)
        ax.grid(ls=':', alpha=0.5, c=(0.5, 0.5, 0.5))

        ax.set_ylim(np.ceil(max(self.t))+2,0)
        ax.set_xlim(np.floor(min(self.x))-2,np.ceil(max(self.x))+2)
        ax.set_aspect('equal')
        return ax
        
    def disp_grid(self):
        ax=self.disp()
        ax.scatter(self.x,self.tgrid)
    
    def add_to_img(self, img, wavelet):
        maxind=np.size(img,1)        
        wavlen=np.floor(len(wavelet)/2).astype(int)      
        
        self.imgind=self.imgind[self.tind < maxind-1]        
        self.amp=self.amp[self.tind < maxind-1]        
        self.tind=self.tind[self.tind < maxind-1]
        ind_begin=self.tind-wavlen        
        
        
        for i,sample in enumerate(wavelet):
#            if (np.size(img,1) > ind_begin+i):
            img[self.imgind,ind_begin+i]=img[self.imgind,ind_begin+i]+sample*self.amp
        return img
    

class Line:
    def __init__(self, xmin, xmax, tmin, tmax, xarray, tarray):
        self.xmin=xmin
        self.xmax=xmax
        self.tmin=tmin
        self.tmax=tmax
        xstep=xarray[1]-xarray[0]
        tstep=tarray[1]-tarray[0]
        
        xmin=xmin-np.mod(xmin,xstep)
        xmax=xmax-np.mod(xmax,xstep)
        tmin=tmin-np.mod(tmin,tstep)
        tmax=tmax-np.mod(tmax,tstep)
        
        self.x = np.arange(xmin,xmax+xstep,xstep)
        self.t = tmin+(tmax-tmin)*(self.x-xmin)/(xmax-xmin)
                
        self.imgind=((self.x-xarray[0])/xstep).astype(int)
        self.tind=((self.t-tarray[0])/tstep).astype(int)
        
            
    
    def add_to_img(self, img, wavelet):
        
        maxind=np.size(img,1)        
        wavlen=np.floor(len(wavelet)/2).astype(int)      
        
        self.imgind=self.imgind[self.tind < maxind-1]        
        self.tind=self.tind[self.tind < maxind-1]
        ind_begin=self.tind-wavlen                
        

        for i,sample in enumerate(wavelet):
                img[self.imgind,ind_begin+i]=img[self.imgind,ind_begin+i]+sample
        return img
    


def migrate(img,v,aper,xarray,tarray):
    imgmig=np.zeros_like(img)
#    for x0 in xarray[(xarray>xarray[0]+aper) & (xarray<xarray[-1]-aper)]:
    for x0 in xarray:    
        for t0 in tarray[:-1]:
            xmig=xarray[(x0-aper<=xarray) & (xarray<=x0+aper)]
            hi=Hyperbola(xmig,tarray,x0,v,t0)
#            hi.tind=hi.tind[hi.tind*tstep<=tarray[-1]]
#            print "t0 = {0}, x0 = {1}".format(t0,x0)            
#            print hi.tind
            hi.imgind=np.arange(hi.x[0]/xstep,(hi.x[-1]+xstep)/xstep).astype(int)
#            si=np.sum(img[hi.imgind,hi.tind])
            si=np.mean(img[hi.imgind,hi.tind]*hi.amp)
#            si=np.mean(img[hi.imgind,hi.tind])
            imgmig[(x0/xstep).astype(int),(t0/tstep).astype(int)]=si
#            if ( (t0==3 and x0==10) or (t0==7 and x0==17) or (t0==11 and x0==12)  ):
#            if ( (t0==8 and x0==20)):
#                ax_data.plot(hi.x,hi.t,c='m',lw=3,alpha=0.8)
#                ax_data.plot(hi.x0,hi.t0,marker='H', mfc='r', mec='m',ms=5)
#                for xi in xmig:
#                    ax_data.plot([xi,hi.x0],[0,hi.t0],c='#AFFF94',lw=1.5,alpha=1)
#            
    return imgmig
            
        

#создаём сетку, в которой работаем
xstep=1
tstep=1
xarray=np.arange(0, 310, xstep)    
tarray=np.arange(0, 210 , tstep)
img=np.zeros((len(xarray), len(tarray)))

#генерируем синтетику
x1left=1
x1right=1



#строим годограф
#hyp1.disp()
#hyp1.disp_grid()

#генерируем изображение
img=np.zeros((len(xarray), len(tarray)))


#
#hyp1=Hyperbola(xarray, tarray, x0=150, v=2, t0=50)
#hyp2=Hyperbola(xarray, tarray, x0=180, v=2, t0=40)
#img=hyp1.add_to_img(img, [-1,2,-1])
#img=hyp2.add_to_img(img, [-1,2,-1])

hyp1=Hyperbola(xarray, tarray, x0=150, v=2, t0=50)
img=hyp1.add_to_img(img, [-1,2,-1])

hyp2=Hyperbola(xarray, tarray, x0=210, v=2, t0=70)
img=hyp2.add_to_img(img, [-1,2,-1])

hyp3=Hyperbola(xarray, tarray, x0=100, v=2, t0=100)
img=hyp3.add_to_img(img, [-1,2,-1])

#
#line1=Line(50,175,50,150,xarray,tarray)
#img=line1.add_to_img(img, [-1,2,-1])
#line2=Line(100,250,50,150,xarray,tarray)
#img=line2.add_to_img(img, [-1,2,-1])
#line3=Line(150,350,50,150,xarray,tarray)
#img=line3.add_to_img(img, [-1,2,-1])
#line4=Line(200,450,50,150,xarray,tarray)
#img=line4.add_to_img(img, [-1,2,-1])



f2=plt.figure(figsize=(12,10), facecolor='w')
ax_data=f2.add_subplot(111)
q=ax_data.imshow(img.T,interpolation='none',cmap=cm.Greys, vmin=-2,vmax=2, extent=[xarray[0]-xstep/2, xarray[-1]+xstep/2, tarray[-1]+tstep/2, tarray[0]-tstep/2])
#q=ax_data.imshow(img.T,cmap=cm.Greys, vmin=-2,vmax=2, extent=[xarray[0]-xstep/2, xarray[-1]+xstep/2, tarray[-1]+tstep/2, tarray[0]-tstep/2])
ax_data.xaxis.set_ticks_position('top')
ax_data.grid(ls=':', alpha=0.25, lw=1, c='w' ) #(1, 0.72, 0.4)

ax_data.set_ylabel(u'Время, мс',fontsize=16)
ax_data.set_xlabel(u'X, м',fontsize=16)
ax_data.xaxis.set_label_position('top')

#ax_data.plot(hyp1.x,hyp1.t,c='m',lw=3,ls='--')
#ax_data.plot(hyp1.x0,hyp1.t0,marker='H', mfc='r', mec='m',ms=10)

ax_data.set_xlim(xarray[0],xarray[-1])
ax_data.set_ylim(tarray[-1],tarray[0])



##МИГРАЦИЯ
v=2
aper=250
res=migrate(img,v,aper,xarray,tarray)
f3=plt.figure(figsize=(12,10), facecolor='w')
ax_mgr=f3.add_subplot(111)
#q=ax_mgr.imshow(res.T,interpolation='none',cmap=cm.Greys,vmin=-2,vmax=2, extent=[xarray[0]-xstep/2, xarray[-1]+xstep/2, tarray[-1]+tstep/2, tarray[0]-tstep/2])
q=ax_mgr.imshow(res.T,cmap=cm.Greys,vmin=-0.2,vmax=0.2, extent=[xarray[0]-xstep/2, xarray[-1]+xstep/2, tarray[-1]+tstep/2, tarray[0]-tstep/2])
#q=ax_mgr.imshow(res.T,interpolation='none',cmap=cm.Greys, extent=[xarray[0]-xstep/2, xarray[-1]+xstep/2, tarray[-1]+tstep/2, tarray[0]-tstep/2])
ax_mgr.xaxis.set_ticks_position('top')
ax_mgr.grid(ls=':', alpha=0.25, lw=1, c='w' )

ax_mgr.set_ylabel(u'Время, мс',fontsize=16)
ax_mgr.set_xlabel(u'X, м',fontsize=16)
ax_mgr.xaxis.set_label_position('top')
