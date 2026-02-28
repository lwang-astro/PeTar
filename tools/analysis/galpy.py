from sdar.base import *
from sdar.functions import *
from .data import *

HEADER_GALPY_POT_MAP_OFFSET = 16

class GalpyPotentialMapPar(DictNpArrayMix):
    """ Galpy potential measurement parameter file
    keys: 
        time0, initial time
        dt: time step
        nstep: total evolution steps
        dt_out: output interval
        xmin, xmax, nx: x-axis mesh range and number of points
        ymin, ymax, ny: y-axis mesh range and number of points
        zmin, zmax, nz: z-axis mesh range and number of points
    """
    def __init__(self, fpar):
        fp = open(fpar, 'r')
        header = fp.readline()
        fp.close()
        t0, dt, nstep, dt_out, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz = header.split()
        self.t0 = float(t0)
        self.dt = float(dt)
        self.nstep = int(nstep)
        self.dt_out = float(dt_out)
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.nx = int(nx)
        self.ymin = float(ymin)
        self.ymax = float(ymax)
        self.ny = int(ny)
        self.zmin = float(zmin)
        self.zmax = float(zmax)
        self.nz = int(nz)

class GalpyPotentialMapHeader():
    """ Galpy potential map header information
    keys:
        time: snapshot time
        nx: number of x-axis points
        ny: number of y-axis points
    """
    def read(self, filename, snapshot_format='ascii'):
        dtype_header = np.dtype([('time', np.float64), ('nx', np.int32), ('ny', np.int32)])
        if snapshot_format == 'ascii':
            with open(filename, 'r') as f:
                header = f.readline()
                t0, nx, ny = header.split()
                self.time = float(t0)
                self.nx = int(nx)
                self.ny = int(ny)
        else:
            with open(filename, 'rb') as f:
                header_data = np.fromfile(f, dtype=dtype_header, count=1)
            self.time = header_data['time'][0]
            self.nx = header_data['nx'][0]
            self.ny = header_data['ny'][0]

    def __init__(self, filename=None, snapshot_format='ascii'):
        if filename is not None:
            self.read(filename, snapshot_format)
        else:
            self.time = 0.0
            self.nx = 0
            self.ny = 0

class GalpyPotentialMap(DictNpArrayMix):
    """ Galpy measure point
    keys: (class members)
        mass (1D): mass
        pos (2D,3): postion x, y, z
        vel (2D,3): velocity vx, vy, vz
        acc (2D,3): acceleration ax, ay, az
        pot (1D): potential
        den (1D): density
    """
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['mass', np.float64], ['pos', (np.float64, 3)], 
                ['vel', (np.float64, 3)], ['acc',(np.float64, 3)], 
                ['pot',np.float64], ['den',np.float64]]
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
    
    def contour(self, axes, pars, **kwargs):
        """ 
        """

        xyscale=[[pars.xmin, pars.xmax, pars.ymin, pars.ymax],
                 [pars.xmin, pars.xmax, pars.zmin, pars.zmax]]

        count = np.log10(-self.pot)
        vmin = count.min()
        vmax = count.max()

        for key, item in kwargs.items():
            if key == 'vmin':
                vmin = item
            elif key == 'vmax':
                vmax = item

        cset = axes.contour(count, linewidths=2, extent=xyscale[i], 
                            vmin=vmin, vmax=vmax)

        axes.clabel(cset, inline=True, fmt='%1.1f')
        axes.set_xlabel(xylabels[i][0])
        axes.set_ylabel(xylabels[i][1])
    
        im = axes.imshow(count, cmap=pb.cm.RdBu,
                         aspect=(xyscale[i][1]-xyscale[i][0])/(xyscale[i][3]-xyscale[i][2]),
                         interpolation='bilinear', origin='lower', 
                         extent=xyscale[i], vmin=vmin, vmax=vmax)
        #cbaxes = fig.add_axes([0.95, 0.1, 0.03, 0.75]) 
        cbar = plt.colorbar(im, cax = axes) 
        #cbar = plt.colorbar(im)
        cbar.set_label(r'log(Pot)')

