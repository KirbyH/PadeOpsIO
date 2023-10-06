import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
class PadeOpsViz:
    """
    Class to store 3D data arrays for post processing
    """
    def __init__(self, outputdir_name, runid, tidx, Lx, Ly, Lz, u_normfactor=1., budget=False, numScalars=0):
        self.runid = runid
        self.outputdir_name = outputdir_name
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        info_fname = outputdir_name + '/Run' + '{:02d}'.format(runid) + '_info_t' + '{:06d}'.format(tidx) + '.out'
        self.info = np.genfromtxt(info_fname, dtype=None)
        if budget==True:
            self.outputdir_name = outputdir_name + '/budgets'
        self.time = self.info[0]
        self.nx = int(self.info[1]); self.ny = int(self.info[2]); self.nz = int(self.info[3])
        self.dx = Lx/self.nx
        self.dy = Ly/self.ny
        self.dz = Lz/self.nz
        self.u_normfactor = u_normfactor
        #self.xG,self.yG,self.zG = np.meshgrid(np.linspace(0,Lx-self.dx,self.nx),np.linspace(0,Ly-self.dy,self.ny), np.linspace(self.dz/2,Lz-(self.dz/2),self.nz),indexing='ij')
        self.xLine = np.linspace(0,Lx-self.dx,self.nx)
        self.yLine = np.linspace(0,Ly-self.dy,self.ny)
        self.zLine = np.linspace(self.dz/2,Lz-(self.dz/2),self.nz)
        self.numScalars = numScalars
        if self.numScalars>0:
            self.sc = np.zeros((self.nx,self.ny,self.nz,numScalars))
        print('PadeOpsViz initialized using info files at time:' + '{:.06f}'.format(self.time))
        
    def RealZPlanes_u(self, tidx, zind):
        #Run01_t004200_x00060.plu
        u_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_t' + '{:06d}'.format(tidx) + '_z' +'{:05d}'.format(zind) + '.plu'
        u = np.fromfile(u_fname,dtype=np.dtype(np.float64),count=-1)
        u = self.u_normfactor*u.reshape((self.nx,self.ny),order='F')
        return u
        
    def ReadVelocities_budget(self, tidx, n, budget_number, terms):
        self.budget = np.zeros((self.nx, self.ny, self.nz, np.size(terms)))  # fourth dim. used to be np.size(terms)+1. Changed 12/08/21
        ind = 0
        for i in terms:
            u_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_budget' + '{:01d}'.format(budget_number) + \
            '_term' + '{:02d}'.format(i) + '_t' + '{:06d}'.format(tidx) + '_n' + '{:06d}'.format(n) + '.s3D'
            temp = np.fromfile(u_fname,dtype=np.dtype(np.float64),count=-1)
            self.budget[:,:,:,ind] = self.u_normfactor*temp.reshape((self.nx,self.ny,self.nz),order='F')
            ind+=1
        print('PadeOpsViz loaded the budget fields at time:' + '{:.06f}'.format(tidx))

    def ReadVelocities(self, tidx, readAll=True):
        info_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_info_t' + '{:06d}'.format(tidx) + '.out'
        info = np.genfromtxt(info_fname, dtype=None)
        self.time = info[0]
        if (readAll):
            u_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_uVel_t' + '{:06d}'.format(tidx) + '.out'
            self.u = np.fromfile(u_fname,dtype=np.dtype(np.float64),count=-1)
            self.u = self.u_normfactor*self.u.reshape((self.nx,self.ny,self.nz),order='F')
            v_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_vVel_t' + '{:06d}'.format(tidx) + '.out'
            self.v = np.fromfile(v_fname,dtype=np.dtype(np.float64),count=-1)
            self.v = self.u_normfactor*self.v.reshape((self.nx,self.ny,self.nz),order='F')
            w_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_wVel_t' + '{:06d}'.format(tidx) + '.out'
            self.w = np.fromfile(w_fname,dtype=np.dtype(np.float64),count=-1)
            self.w = self.u_normfactor*self.w.reshape((self.nx,self.ny,self.nz),order='F')
            print('PadeOpsViz loaded the velocity fields at time:' + '{:.06f}'.format(self.time))
        else:
            u_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_uVel_t' + '{:06d}'.format(tidx) + '.out'
            self.u = np.fromfile(u_fname,dtype=np.dtype(np.float64),count=-1)
            self.u = self.u.reshape((self.nx,self.ny,self.nz),order='F')
            print('PadeOpsViz loaded the u velocity field at time:' + '{:.06f}'.format(self.time))
            
    def ReadScalars(self, tidx, readAll=True, scalarNum=1):
        #self.time = self.info[0]; 
        info_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_info_t' + '{:06d}'.format(tidx) + '.out'
        info = np.genfromtxt(info_fname, dtype=None)
        self.time = info[0]
        if self.numScalars==1:
            fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_sc' + '{:02d}'.format(scalarNum) + '_t' + '{:06d}'.format(tidx) + '.out'
            temp = np.fromfile(fname,dtype=np.dtype(np.float64),count=-1)
            self.sc = temp.reshape((self.nx,self.ny,self.nz),order='F')
        else:
            for i in range(self.numScalars):
                fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_sc' + '{:02d}'.format(i+1) + '_t' + '{:06d}'.format(tidx) + '.out'
                temp = np.fromfile(fname,dtype=np.dtype(np.float64),count=-1)
                self.sc[:,:,:,i] = temp.reshape((self.nx,self.ny,self.nz),order='F')
    
    def Read_x_slice(self, tidx, yid, lab):
        fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_t' + '{:06d}'.format(tidx) + '_x'+ '{:05d}'.format(yid) + lab
        f_x_slice = np.fromfile(fname,dtype=np.dtype(np.float64),count=-1)
        f_x_slice = f_x_slice.reshape((self.ny,self.nz),order='F')
        return f_x_slice
    
    def Read_y_slice(self, tidx, yid, lab):
        fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_t' + '{:06d}'.format(tidx) + '_y'+ '{:05d}'.format(yid) + lab
        f_y_slice = np.fromfile(fname,dtype=np.dtype(np.float64),count=-1)
        f_y_slice = f_y_slice.reshape((self.nx,self.nz),order='F')
        return f_y_slice
    
    def Read_z_slice(self, tidx, yid, lab):
        fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_t' + '{:06d}'.format(tidx) + '_z'+ '{:05d}'.format(yid) + lab
        f_z_slice = np.fromfile(fname,dtype=np.dtype(np.float64),count=-1)
        f_z_slice = f_z_slice.reshape((self.nx,self.ny),order='F')
        return f_z_slice
        
    def Read_BPF_Velocities(self, tidx):
        info_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_info_t' + '{:06d}'.format(tidx) + '.out'
        info = np.genfromtxt(info_fname, dtype=None)
        self.time = info[0]
        u_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_uBPF_t' + '{:06d}'.format(tidx) + '.out'
        self.uBPF = np.fromfile(u_fname,dtype=np.dtype(np.float64),count=-1)
        self.uBPF = self.uBPF.reshape((self.nx,self.ny,self.nz),order='F')
        v_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_vBPF_t' + '{:06d}'.format(tidx) + '.out'
        self.vBPF = np.fromfile(v_fname,dtype=np.dtype(np.float64),count=-1)
        self.vBPF = self.vBPF.reshape((self.nx,self.ny,self.nz),order='F')
        w_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_wBPF_t' + '{:06d}'.format(tidx) + '.out'
        self.wBPF = np.fromfile(w_fname,dtype=np.dtype(np.float64),count=-1)
        self.wBPF = self.wBPF.reshape((self.nx,self.ny,self.nz),order='F')
        print('PadeOpsViz loaded the BPF velocity fields at time:' + '{:.06f}'.format(self.time))
        
    def ReadPotTemp(self, tidx):
        info_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_info_t' + '{:06d}'.format(tidx) + '.out'
        info = np.genfromtxt(info_fname, dtype=None)
        self.time = info[0]
        T_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_potT_t' + '{:06d}'.format(tidx) + '.out'
        self.PotT = np.fromfile(T_fname,dtype=np.dtype(np.float64),count=-1)
        self.PotT = self.PotT.reshape((self.nx,self.ny,self.nz),order='F')
        print('PadeOpsViz loaded the Potential Temperature field at time:' + '{:.06f}'.format(self.time))
    
    def ReadPotTemp_restart(self, tidx):
        T_fname = self.outputdir_name + '/RESTART_Run' + '{:02d}'.format(self.runid) + '_T.' + '{:06d}'.format(tidx)
        self.PotT = np.fromfile(T_fname,dtype=np.dtype(np.float64),count=-1)
        self.PotT = self.PotT.reshape((self.nx,self.ny,self.nz),order='F')
        print('PadeOpsViz loaded the Potential Temperature field from the RESTART file')
        
    def ReadU_restart(self, tidx):
        u_fname = self.outputdir_name + '/RESTART_Run' + '{:02d}'.format(self.runid) + '_u.' + '{:06d}'.format(tidx)
        self.u = np.fromfile(u_fname,dtype=np.dtype(np.float64),count=-1)
        self.u = self.u.reshape((self.nx,self.ny,self.nz),order='F')
        print('PadeOpsViz loaded the u field from the RESTART file')
        
    def ReadV_restart(self, tidx):
        v_fname = self.outputdir_name + '/RESTART_Run' + '{:02d}'.format(self.runid) + '_v.' + '{:06d}'.format(tidx)
        self.v = np.fromfile(v_fname,dtype=np.dtype(np.float64),count=-1)
        self.v = self.v.reshape((self.nx,self.ny,self.nz),order='F')
        print('PadeOpsViz loaded the v field from the RESTART file')
    
    def ReadW_restart(self, tidx):
        wE_fname = self.outputdir_name + '/RESTART_Run' + '{:02d}'.format(self.runid) + '_w.' + '{:06d}'.format(tidx)
        self.wE = np.fromfile(wE_fname,dtype=np.dtype(np.float64),count=-1)
        self.wE = self.wE.reshape((self.nx,self.ny,self.nz+1),order='F')
        print('PadeOpsViz loaded the w (edge) field from the RESTART file')
    
    def ReadPressure(self, tidx):
        info_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_info_t' + '{:06d}'.format(tidx) + '.out'
        info = np.genfromtxt(info_fname, dtype=None)
        self.time = info[0]
        P_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_prss_t' + '{:06d}'.format(tidx) + '.out'
        self.Press = np.fromfile(P_fname,dtype=np.dtype(np.float64),count=-1)
        self.Press = self.Press.reshape((self.nx,self.ny,self.nz),order='F')
        print('PadeOpsViz loaded the Pressure field at time:' + '{:.06f}'.format(self.time))
        
    def ReadFringePressure(self, tidx):
        info_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_info_t' + '{:06d}'.format(tidx) + '.out'
        info = np.genfromtxt(info_fname, dtype=None)
        self.time = info[0]
        P_fname = self.outputdir_name + '/Run' + '{:02d}'.format(self.runid) + '_pfrn_t' + '{:06d}'.format(tidx) + '.out'
        self.Pfrn = np.fromfile(P_fname,dtype=np.dtype(np.float64),count=-1)
        self.Pfrn = self.Pfrn.reshape((self.nx,self.ny,self.nz),order='F')
        print('PadeOpsViz loaded the Pressure field at time:' + '{:.06f}'.format(self.time))
        
    def plot_xy(self,field, zid=0, fsize=(4,4), pix=200, colormap='jet', interp='none', 
                printtitle=False, xlim=[0, 40], cbar=True, clim=None, ax=None):
        
        if ax is None:  # modified 01/19/2022
            fig = plt.figure(num=None, figsize=fsize, dpi=pix, facecolor='w', edgecolor='k')
            #plt.pcolormesh(self.xG[:,:,0],self.yG[:,:,0],field[:,:,zid],cmap=colormap)
            ax = plt.gca()
        else: 
            fig = ax.get_figure()
            
        if clim is not None: 
            vmin = clim[0]
            vmax = clim[1]
        else: 
            vmin, vmax = None, None

        im = ax.imshow(np.transpose(field[:,:,zid]), 
                       extent=[self.xLine[0], self.xLine[-1], self.yLine[0], self.yLine[-1]], 
                       cmap=colormap, origin='lower', interpolation=interp, aspect=1., vmin=vmin, vmax=vmax)
#         cbar = fig.colorbar(im)
        ax.set_xlabel('$x/D$')
        ax.set_ylabel('$y/D$')
        ax.set_xlim(xlim)
        
        if colormap is not None: 
            # colorbar: https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph 
            cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
            cbar = fig.colorbar(im, cax=cax)
            if clim is not None: 
                cbar.set_clim(clim)
        
        return fig, ax

    def plot_xy_slice(self,field, vlim, fsize=(4,4), pix=200, colormap='jet', interp='none'):
        
        fig = plt.figure(num=None, figsize=fsize, dpi=pix, facecolor='w', edgecolor='k')
        #plt.pcolormesh(self.xG[:,:,0],self.yG[:,:,0],field[:,:,zid],cmap=colormap)
        ax = plt.gca()
        im = ax.imshow(np.transpose(field), 
                       extent=[self.xLine[0], self.xLine[-1], self.yLine[0], self.yLine[-1]], 
                       cmap=colormap, origin='lower', interpolation=interp, aspect=1.)#, vmin=vlim[0], vmax=vlim[1] )
        cbar = fig.colorbar(im)
        #plt.xlabel('X')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    
    def plot_spacetime_slice(self,field, arrx, xlab, fsize=(4,4), pix=200, colormap='jet', interp='none'):
        
        fig = plt.figure(num=None, figsize=fsize, dpi=pix, facecolor='w', edgecolor='k')
        #plt.pcolormesh(self.xG[:,:,0],self.yG[:,:,0],field[:,:,zid],cmap=colormap)
        ax = plt.gca()
        im = ax.imshow(np.transpose(field), 
                       extent=[arrx[0], arrx[-1], self.times[0], self.times[-1]], 
                       cmap=colormap, origin='lower', interpolation=interp, aspect=1. )
        cbar = fig.colorbar(im)
        ax.set_xlabel(xlab)
        ax.set_ylabel('t')
        
    def plot_xz(self,field, yid=0, fsize=(4,4), pix=200, colormap='jet', 
               printtitle=False, xlim=[0, 40], cbar=True, clim=None, ax=None):
        
        if ax is None:  # modified 02/15/2022
            fig = plt.figure(num=None, figsize=fsize, dpi=pix, facecolor='w', edgecolor='k')
            ax = plt.gca()
        else: 
            fig = ax.get_figure()

#         fig = plt.figure(num=None, figsize=fsize, dpi=pix, facecolor='w', edgecolor='k')
#         plt.pcolormesh(self.xG[:,0,:],self.zG[:,0,:],field[:,yid,:],cmap=colormap)
#         ax = plt.gca()
        if clim is not None: 
            vmin = clim[0]
            vmax = clim[1]
        else: 
            vmin, vmax = None, None

        im = ax.imshow(np.transpose(field[:,yid,:]), 
                       extent=[self.xLine[0], self.xLine[-1], self.zLine[0], self.zLine[-1]], 
                       cmap=colormap, origin='lower', interpolation='none', aspect=1., vmin=vmin, vmax=vmax)
        ax.set_aspect('equal')
        
        ax.set_xlim([0,self.Lx])
        ax.set_ylim([0,self.Lz])
        ax.set_ylabel('$z/D$')
        ax.set_xlabel('$x/D$')
        ax.set_xlim(xlim)
        
        if (colormap is not None) and cbar: 
            # colorbar: https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph 
            cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
            cbar = fig.colorbar(im, cax=cax)
            if clim is not None: 
                cbar.set_clim(clim)
        
        return fig, ax
        
    def plot_yz(self,field, xid=0, fsize=(4,4), pix=200, colormap='jet', clim=None):
        
        fig = plt.figure(num=None, figsize=fsize, dpi=pix, facecolor='w', edgecolor='k')
#         plt.pcolormesh(self.yG[0,:,:],self.zG[0,:,:],field[xid,:,:],cmap=colormap)  # old code
        ax = plt.gca()
    
        if clim is not None: 
            vmin = clim[0]
            vmax = clim[1]
        else: 
            vmin, vmax = None, None
        
        im = ax.imshow(np.transpose(field[xid, :,:]), 
                       extent=[self.yLine[0], self.yLine[-1], self.zLine[0], self.zLine[-1]], 
                       cmap=colormap, origin='lower', interpolation='none', aspect=1., vmin=vmin, vmax=vmax)

        ax.set_aspect('equal')
        ax.set_xlim([0,self.Ly])
        ax.set_ylim([0,self.Lz])
        ax.set_ylabel('$z/D$')
        ax.set_xlabel('$y/D$')
        
        if colormap is not None: 
            # colorbar: https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph 
            cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
            cbar = fig.colorbar(im, cax=cax)
            if clim is not None: 
                cbar.set_clim(clim)
        
        return fig, ax

    def plot_yz_slice(self,field, xid=0, fsize=(4,4), pix=200, colormap='jet'):
        fig = plt.figure(num=None, figsize=fsize, dpi=pix, facecolor='w', edgecolor='k')
        plt.pcolormesh(self.yLine,self.zLine,np.transpose(field),cmap=colormap)
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.xlim([0,self.Ly])
        plt.ylim([0,self.Lz])
        plt.colorbar()
        plt.ylabel('Z')
        plt.xlabel('Y')
        
    def plot_yz_u_and_P(self, xid=0, fsize=(8,4), pix=200, cmap1='jet', cmap2='jet'):
        f, axarr = plt.subplots(1, 2, figsize=(10,5))
        im1 = axarr[0].pcolormesh(self.yG[xid,:,:],self.zG[xid,:,:],self.u[xid,:,:],cmap=cmap1)
        axarr[0].set_xlim([0,self.Ly])
        axarr[0].set_ylim([0,self.Lz])
        axarr[0].set_ylabel('Z')
        axarr[0].set_xlabel('Y')
        axarr[0].set_aspect('equal')
        
        im2 = axarr[1].pcolormesh(self.yG[xid,:,:],self.zG[xid,:,:],self.Press[xid,:,:],cmap=cmap1)
        axarr[1].set_xlim([0,self.Ly])
        axarr[1].set_ylim([0,self.Lz])
        axarr[1].set_ylabel('Z')
        axarr[1].set_xlabel('Y')
        axarr[1].set_aspect('equal')
        
        f.suptitle('x = '+ '{:0.2f}'.format(self.xG[xid,0,0]), fontsize=12)
        
        f.colorbar(im1, ax=axarr[0])
        f.colorbar(im2, ax=axarr[1])
    
    
    def plot_xy_2field(self, f1, f2, zid=0, fsize=(8,4), pix=200, cmap1='jet', cmap2='jet',label1='f1',label2='f2',clim1=[-0.35,0.35],clim2=[-0.35,0.35]):
        f, axarr = plt.subplots(1, 2, figsize=fsize)
        xG,yG = np.meshgrid(self.xLine, self.yLine)
        im1 = axarr[0].pcolormesh(xG,yG,f1[:,:,zid],cmap=cmap1)#,vmin=clim1[0], vmax=clim1[1])
        axarr[0].set_xlim([0,self.Ly])
        axarr[0].set_ylim([0,self.Lz])
        axarr[0].set_ylabel('Z')
        axarr[0].set_xlabel('Y')
        axarr[0].set_aspect('equal')
        axarr[0].set_title(label1)


        im2 = axarr[1].pcolormesh(xG,yG,f2[:,:,zid],cmap=cmap2)#,vmin=clim1[0], vmax=clim1[1])
        axarr[1].set_xlim([0,self.Ly])
        axarr[1].set_ylim([0,self.Lz])
        axarr[1].set_ylabel('Z')
        axarr[1].set_xlabel('Y')
        axarr[1].set_aspect('equal')
        axarr[1].set_title(label2)
                
        f.colorbar(im1, ax=axarr[0])
        f.colorbar(im2, ax=axarr[1])
    
    def plot_yt_2field_slices(self, f1, f2, fsize=(8,4), pix=200, cmap1='jet', cmap2='jet',label1='f1',label2='f2',clim1=[-0.35,0.35],clim2=[-0.35,0.35]):
        f, axarr = plt.subplots(1, 2, figsize=fsize)
        yG,TG = np.meshgrid(self.yLine, self.times,indexing='ij')
        im1 = axarr[0].pcolormesh(yG,TG,f1,cmap=cmap1)#,vmin=clim1[0], vmax=clim1[1])
        axarr[0].set_ylim([0,self.times[-1]])
        axarr[0].set_xlim([0,self.Ly])
        axarr[0].set_ylabel('t')
        axarr[0].set_xlabel('y')
        axarr[0].set_aspect('equal')
        axarr[0].set_title(label1)


        im2 = axarr[1].pcolormesh(yG,TG,f2,cmap=cmap2)#,vmin=clim1[0], vmax=clim1[1])
        axarr[1].set_ylim([0,self.times[-1]])
        axarr[1].set_xlim([0,self.Ly])
        axarr[1].set_ylabel('t')
        axarr[1].set_xlabel('y')
        axarr[1].set_aspect('equal')
        axarr[1].set_title(label2)
                
        f.colorbar(im1, ax=axarr[0])
        f.colorbar(im2, ax=axarr[1])
        
    def plot_xt_2field_slices(self, f1, f2, fsize=(8,4), pix=200, cmap1='jet', cmap2='jet',label1='f1',label2='f2',clim1=[-0.35,0.35],clim2=[-0.35,0.35]):
        f, axarr = plt.subplots(1, 2, figsize=fsize)
        xG,TG = np.meshgrid(self.xLine, self.times,indexing='ij')
        im1 = axarr[0].pcolormesh(xG,TG,f1,cmap=cmap1)#,vmin=clim1[0], vmax=clim1[1])
        axarr[0].set_ylim([0,self.times[-1]])
        axarr[0].set_xlim([0,self.Lx])
        axarr[0].set_ylabel('t')
        axarr[0].set_xlabel('x')
        axarr[0].set_aspect('equal')
        axarr[0].set_title(label1)


        im2 = axarr[1].pcolormesh(xG,TG,f2,cmap=cmap2)#,vmin=clim1[0], vmax=clim1[1])
        axarr[1].set_ylim([0,self.times[-1]])
        axarr[1].set_xlim([0,self.Lx])
        axarr[1].set_ylabel('t')
        axarr[1].set_xlabel('x')
        axarr[1].set_aspect('equal')
        axarr[1].set_title(label2)
                
        f.colorbar(im1, ax=axarr[0])
        f.colorbar(im2, ax=axarr[1])
        
    def mean_fluct_2D(self,u):
        umean = np.zeros((self.nz),dtype=np.float64,order='F')
        ufluct = np.zeros((self.nx,self.ny,self.nz),dtype=np.float64,order='F')
        for i in range(0,self.nz):
            umean[i] = np.mean(u[:,:,i])
            ufluct[:,:,i] = u[:,:,i] - umean[i]
        return umean, ufluct
    
    def mean_fluct_3D(self,u):
        umean = np.mean(u[:,:,:])
        ufluct = u - umean
        return umean, ufluct
    
    def get_stats_global(self):
        umn, fluct = self.mean_fluct_3D(self.u)
        uvar = np.mean(fluct*fluct)
        vmn, fluct = self.mean_fluct_3D(self.v)
        vvar = np.mean(fluct*fluct)
        wmn, fluct = self.mean_fluct_3D(self.w)
        wvar = np.mean(fluct*fluct)
        return uvar, vvar, wvar, umn, vmn, wmn
   
    def get_stats_zprofiles(self):
        umn, fluct = self.mean_fluct_2D(self.u)
        uvar = np.zeros((self.nz),dtype=np.float64,order='F')
        for i in range(0,self.nz):
            uvar[i] = np.mean(fluct[:,:,i]*fluct[:,:,i])
        vmn, fluct = self.mean_fluct_2D(self.u)
        vvar = np.zeros((self.nz),dtype=np.float64,order='F')
        for i in range(0,self.nz):
            vvar[i] = np.mean(fluct[:,:,i]*fluct[:,:,i])
        wmn, fluct = self.mean_fluct_2D(self.u)
        wvar = np.zeros((self.nz),dtype=np.float64,order='F')
        for i in range(0,self.nz):
            wvar[i] = np.mean(fluct[:,:,i]*fluct[:,:,i])
        zplot = self.zG[0,0,:]
        return umn, vmn, wmn, uvar, vvar, wvar, zplot
    
    def get_1pt_cross_corr_2D(self,u,v):
        corr = np.zeros((self.nz),dtype=np.float64,order='F')
        for i in range(0,self.nz):
            umean = np.mean(u[:,:,i])
            ufluct = u[:,:,i] - umean
            vmean = np.mean(v[:,:,i])
            vfluct = v[:,:,i] - vmean
            corr[i] = np.mean(ufluct*vfluct)
        return corr
    
    def get_k_3d(self):
        k1_ind = np.arange(-self.nx/2,self.nx/2)
        k2_ind = np.arange(-self.ny/2,self.ny/2)
        k3_ind = np.arange(-self.nz/2,self.nz/2)
        k1_wav = np.fft.fftshift(k1_ind*(2*np.pi/(self.nx*self.dx)))
        k2_wav = np.fft.fftshift(k2_ind*(2*np.pi/(self.ny*self.dy)))
        k3_wav = np.fft.fftshift(k3_ind*(2*np.pi/(self.nz*self.dz)))
        k1, k2, k3 = np.meshgrid(k1_wav, k2_wav, k3_wav, indexing='ij')
        return k1, k2, k3
    
    def get_xSpectra_from_xy_plane(self,f,zid):
        k1 = np.fft.fftshift(np.arange(-self.nx/2,self.nx/2)*(2*np.pi/(self.nx*self.dx)))
        fhat = np.absolute(sp.fftpack.fft(f[:,:,zid],n=None,axis=0))
        fpow = np.mean(fhat*fhat,axis=1)
        kx = k1[0:np.int32(self.nx/2)]
        dk = 2*np.pi/(self.nx*self.dx)
        fout = fpow[0:np.int32(self.nx/2)]/(dk*(self.nx**2))
        return kx, fout
    
    def TwoPtCorr_xy(self,f,g):
        grot = np.fliplr(np.flipud(g))     
        fhat = np.fft.fft2(f,s=[2*self.nx-1,2*self.ny-1])
        ghat = np.fft.fft2(grot,s=[2*self.nx-1,2*self.ny-1])
        prod = fhat*ghat
        corr2d = np.real(np.fft.ifft2(prod)/(self.nx*self.ny))
        corr1pt = np.mean(f*g)
        corr2d = corr2d/corr1pt
        rx,ry = np.meshgrid(np.linspace(-self.Lx,self.Lx,2*self.nx-1),np.linspace(-self.Ly,self.Ly,2*self.ny-1),indexing='ij')
        return corr2d, rx, ry

    def TwoPtCorr_xt(self,flarge,glarge,steps=1):
        corrsave = np.zeros((2*self.nx-1,2*self.times.size-1),dtype=np.float64,order='F')
        idx = 0
        for yid in np.arange(0,self.ny,steps):
            f = np.squeeze(flarge[:,yid,:])
            g = np.squeeze(glarge[:,yid,:])
            grot = np.fliplr(np.flipud(g))     
            fhat = np.fft.fft2(f,s=[2*self.nx-1,2*self.times.size-1])
            ghat = np.fft.fft2(grot,s=[2*self.nx-1,2*self.times.size-1])
            prod = fhat*ghat
            corr2d = np.real(np.fft.ifft2(prod)/(self.nx*self.times.size))
            corr1pt = np.mean(f*g)
            corr2d = corr2d/corr1pt
            corrsave = corrsave + corr2d
            idx = idx + 1
        corrsave = corrsave/idx  
        rx,rt = np.meshgrid(np.linspace(-self.Lx,self.Lx,2*self.nx-1),np.linspace(-self.times[-1],self.times[-1],2*self.times.size-1),indexing='ij')
        return corrsave, rx, rt
    
    def TwoPtCorr_yt(self,flarge,glarge,steps=1):
        corrsave = np.zeros((2*self.ny-1,2*self.times.size-1),dtype=np.float64,order='F')
        idx = 0
        for xid in np.arange(0,self.nx,steps):
            f = np.squeeze(flarge[xid,:,:])
            g = np.squeeze(glarge[xid,:,:])
            grot = np.fliplr(np.flipud(g))     
            fhat = np.fft.fft2(f,s=[2*self.ny-1,2*self.times.size-1])
            ghat = np.fft.fft2(grot,s=[2*self.ny-1,2*self.times.size-1])
            prod = fhat*ghat
            corr2d = np.real(np.fft.ifft2(prod)/(self.ny*self.times.size))
            corr1pt = np.mean(f*g)
            corr2d = corr2d/corr1pt
            corrsave = corrsave + corr2d
            idx = idx + 1
        corrsave = corrsave/idx  
        ry,rt = np.meshgrid(np.linspace(-self.Ly,self.Ly,2*self.ny-1),np.linspace(-self.times[-1],self.times[-1],2*self.times.size-1),indexing='ij')
        return corrsave, ry, rt
    
    def get_energy_spectrum_3D(self,uMesh, nspect):
        uhat = np.absolute(np.fft.fftn(uMesh)/(self.nx*self.ny*self.nz))
        ek = uhat*uhat
        kline = np.linspace(self.kabs.min(), self.kabs.max(), nspect)
        Ek, k_edges = np.histogram(self.kabs, bins=kline, density=False, weights=ek)
        kline = 0.5*(k_edges[:-1] + k_edges[1:])
        Ek = Ek/(kline[1] - kline[0])
        return kline, Ek
    
    def generate_wavenumbers_3d(self):
        k1_ind = np.arange(-self.nx/2,self.nx/2)
        k2_ind = np.arange(-self.ny/2,self.ny/2)
        k3_ind = np.arange(-self.nz/2,self.nz/2)
        k1_wav = np.fft.fftshift(k1_ind*(2*np.pi/(self.nx*self.dx)))
        k2_wav = np.fft.fftshift(k2_ind*(2*np.pi/(self.ny*self.dy)))
        k3_wav = np.fft.fftshift(k3_ind*(2*np.pi/(self.nz*self.dz)))
        k1, k2, k3 = np.meshgrid(k1_wav, k2_wav, k3_wav, indexing='ij')
        self.kabs = np.sqrt(k1*k1 + k2*k2 + k3*k3)
        print('Wavenumbers generated.')