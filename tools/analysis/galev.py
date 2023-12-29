import numpy as np
from .base import *

def toGalevSnap(filename, particles, zmet, time=0, read_binary=False):
    """
    Save particle data for Galev input snapshot file
    Use id, mass, stellar type, luminosity and temperature

    Parameters:
    ----------
    filename (string): filename to save data
    particle (Particle/Binary): petar particle data with star 
    zmet (float): metallicity in [Fe/H]
    time (float): snapshot time (0)
    read_binary (bool): reading snapshot is binary data (False)
    """

    n = particles.size
    def getPars(p):
        star = p.star
        pid = p.id
        mass = star.mass
        temp = 5778*(star.lum/(star.rad*star.rad))**0.25
        
        kw = star.type
        log10L = np.log10(star.lum)
        log10Teff = np.log10(temp)
        
        return np.transpose((pid, kw, mass, log10L, log10Teff))

    data = None
    if (read_binary):
        n *= 2
        data = []
        for p in [particles.p1, particles.p2]:
            data.append(getPars(p))
        data = np.concatenate(data)
    else:
        data = getPars(particles)
    
    np.savetxt(filename, data, 
               fmt='%d %d %.13g %.13g %.13g', 
               header= ('%d %f %f' % (n,time,zmet)),
               comments='   ')

class FilterJohnson(DictNpArrayMix):
    """
    Johnson filter magnitude

    Members (Keys):
    ----------------
    U, B, V, R, I (1D np.ndarray(float64)): magnitude/flux of one filter
    zero_point_mag (dict): zero-point magnitude (vega)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['U', np.float64], ['B', np.float64], ['V', np.float64], 
                ['R', np.float64], ['I', np.float64]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.zero_point_mag = {'U':20.919987, 'B':20.522241, 'V':21.103851, 'R':21.825808, 'I':22.614393}

class FilterSDSS(DictNpArrayMix):
    """
    SDSS filter magnitude

    Members (Keys):
    ----------------
    u, g, r, i, z (1D np.ndarray(float64)): magnitude/flux of one filter
    zero_point_mag (dict): zero-point magnitude (vega)
    zero_point_mode (string): 
       'std': zero_point_mag are the ones with atmosphere at airmass 1.3 (unknown data release) (default)
       'DR6': zero_point_mag are the ones from DR6 without atmosphere
       'Doi2010': zero_point_mag are the ones from Doi+2010
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['u', np.float64], ['g', np.float64], ['r', np.float64], 
                ['i', np.float64], ['z', np.float64]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.zero_point_mode = 'std'
        for key, item in kwargs.items():
            if key == 'zero_point_mode':
                self.zero_point_mode = item
        
        if self.zero_point_mode == 'std':
            self.zero_point_mag = {'u':21.032514, 'g':20.660448, 'r':21.500454, 'i':22.144540, 'z':22.728246}
        elif self.zero_point_mode == 'DR6':
            self.zero_point_mag = {'u':21.042618, 'g':20.647517, 'r':21.496022, 'i':22.147243, 'z':22.732697}
        elif self.zero_point_mode == 'Doi2010':
            self.zero_point_mag = {'u':21.022568, 'g':20.659458, 'r':21.501515, 'i':22.147651, 'z':22.715423}
        else:
            raise ValueError("SDSS zero-point magnitude mode is not supported ", key,", choices: std, DR6, Doi2010")

class FilterHST(DictNpArrayMix):
    """
    HST filter magnitude

    Members (Keys):
    ----------------
    F220W, F250W, F330W, F435W, F555W, F814W (1D np.ndarray(float64)): magnitude/flux of one filter
    zero_point_mag (dict): zero-point magnitude (vega)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['F220W', np.float64], ['F250W', np.float64], ['F330W', np.float64], 
                ['F435W', np.float64], ['F555W', np.float64], ['F814W', np.float64]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.zero_point_mag = {'F220W':20.790323, 'F250W':21.005258, 'F330W':21.081777, 
                               'F435W':20.505857, 'F555W':21.033912, 'F814W':22.393911}


class FilterCSST(DictNpArrayMix):
    """
    CSST filter magnitude

    Members (Keys):
    ----------------
    NUV, u, g, r, i, z, y (1D np.ndarray(float64)): magnitude/flux of one filter
    zero_point_mag (dict): zero-point magnitude (vega)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['NUV', np.float64], ['u', np.float64], ['g', np.float64], 
                ['r', np.float64], ['i', np.float64], ['z', np.float64], ['y', np.float64]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.zero_point_mag = {'NUV':22.97, 'u':24.10, 'g':25.91, 'r':25.46, 'i':25.06, 'z':24.31, 'y':22.91}

class FilterGaia(DictNpArrayMix):
    """
    CSST filter magnitude

    Members (Keys):
    ----------------
    G, Grp, Gbp (1D np.ndarray(float64)): magnitude/flux of one filter
    zero_point_mag (dict): zero-point magnitude (vega)
    """

    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['G', np.float64], ['Grp', np.float64], ['Gbp', np.float64]]

        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.zero_point_mag = {'G':25.6914396869, 'Grp':24.7626744847, 'Gbp':25.3488107670}

class GalevMag(DictNpArrayMix):
    """
    Galev stellar magnitude for different filters

    Members (keys):
    --------
    id   (1D np.ndarray(int64)): particle id
    type (1D np.ndarray(int64)): SSE stellar type
    mass (1D np.ndarray(float64)): particle mass
    lum  (1D np.ndarray(float64)): luminosity
    temp (1D np.ndarray(float64)): temperature
    z    (1D np.ndarray(float64)): metallicity
    filters (Filter class): magnitude of different filters
                            Default filters: ['Johnson', 'HST', 'CSST', 'Gaia']
    mode (string): indicate the current mode of values:
                   abs_mag: absolute magnitudes
                   app_mag: apparent magnitudes
                   abs_flux: flux converted from absolute magnitudes + zero-point magnitudes (vega) of filters
                   abs_flux: flux converted from apparent magnitudes + zero-point magnitudes (vega) of filters

    Keyword arguments:
    ---------
    filter (list): set filter list (['Johnson', 'HST', 'CSST', 'Gaia'])
    mode (string): set initial mode ('abs_mag')
    distance (float): set distance for computing apparent distance [pc] (10)
    """
    
    def __init__(self, _dat=None, _offset=int(0), _append=False, **kwargs):
        """ DictNpArrayMix type initialzation, see help(DictNpArrayMix.__init__)
        """
        keys = [['id', np.int64], ['type', np.int64], ['mass', np.float64], 
                ['lum', np.float64], ['temp', np.float64], ['z', np.float64]]

        mag_keys = [['Johnson', FilterJohnson], ['SDSS', FilterSDSS], ['HST', FilterHST], ['CSST', FilterCSST], ['Gaia', FilterGaia]]
        filters = ['Johnson', 'SDSS', 'HST', 'CSST', 'Gaia']
        distance = 10
        mode = 'abs_mag'

        for key, item in kwargs.items():
            if (key == 'filter'):
                mag_keys = []
                filters = []
                for name in item:
                    if name == 'Johnson':
                        mag_keys.append(['Johnson', FilterJohnson])
                    elif name == 'SDSS':
                        mag_keys.append(['SDSS', FilterSDSS])
                    elif name == 'HST':
                        mag_keys.append(['HST', FilterHST])
                    elif name == 'CSST':
                        mag_keys.append(['CSST', FilterCSST])
                    elif name == 'Gaia':
                        mag_keys.append(['Gaia', FilterGaia])
                    else:
                        raise ValueError('Galev Filter name %s not found!' % name)
                    filters.append(name)
            elif (key == 'mode'):
                mode = item
            elif (key == 'distance'):
                distance = item
            else:
                raise ValueError('Keyword argument %s not found!' % key)   

        keys += mag_keys
        DictNpArrayMix.__init__(self, keys, _dat, _offset, _append, **kwargs)
        self.initargs['filter'] = filters
        self.initargs['mode'] = mode
        self.initargs['distance'] = distance

    def setMode(self, mode):
        """
        Set mode
        
        Parameters:
        -------------
        mode (string): indicate the current mode of values:
                       abs_mag: absolute magnitudes
                       app_mag: apparent magnitudes
                       abs_flux: flux converted from absolute magnitudes + zero-point magnitudes (vega) of filters
                       abs_flux: flux converted from apparent magnitudes + zero-point magnitudes (vega) of filters

        """
        if mode in ['abs_mag', 'app_mag', 'abs_flux', 'app_flux']:
            self.initargs['mode'] = mode
        else:
            raise ValueError('Mode %s is not supported' % mode)

    def getMode(self):
        """
        Return mode
        """
        return self.initargs['mode']

    def setDistance(self, distance):
        """
        set distance in pc
        """
        self.initargs['distance'] = distance

    def getDistance(self):
        """
        Return distance for apparent magnitude
        """
        return self.initargs['distance']

    def convertToApparentMag(self):
        """
        convert absolute magnitudes to apparent magnitudes
        using saved distance value (self.setDistance)
        required mode: 'abs_mag', change to 'app_mag'

        """
        if (self.getMode() == 'abs_mag'):
            distance = self.getDistance()
            for key, item in self.keys:
                if (issubclass(item, DictNpArrayMix)):
                    for ikey, itype in self[key].keys:
                        self[key][ikey] += 5.*np.log10(distance) - 5
            self.setMode('app_mag')
        else:
            raise ValueError("values are not absolute magnitudes, current mode: %s " % self.getMode())

    def convertToAbsoluteMag(self):
        """
        convert absolute magnitudes to apparent magnitudes
        using saved distance value (self.setDistance)
        required mode: 'app_mag', change to 'abs_mag'

        """
        if (self.getMode() == 'app_mag'):
            distance = self.getDistance()
            for key, item in self.keys:
                if (issubclass(item, DictNpArrayMix)):
                    for ikey, itype in self[key].keys:
                        self[key][ikey] -= 5.*np.log10(distance) - 5
            self.setMode('abs_mag')
        else:
            raise ValueError("values are not apparent magnitudes, current mode: %s " % self.getMode())
    

    def convertToFlux(self):
        """
        convert magnitudes to flux using zero-point magnitudes of filters (vega)
        required mode: 'abs_mag' or 'app_mag', change to 'abs_flux' or 'app_flux'

        """
        mode = self.getMode()
        if (mode[-3:] == 'mag'):
            for key, item in self.keys:
                if (issubclass(item, DictNpArrayMix)):
                    for ikey, itype in self[key].keys:
                        mag = self[key][ikey]
                        self[key][ikey] = 10**(-0.4*(mag + self[key].zero_point_mag[ikey]))
            self.setMode(mode[:3]+'_flux')
        else:
            raise ValueError("values are not magnitude, current mode: %s " % mode)

    def convertToMag(self):
        """
        convert flux to magnitudes using zero-point magnitudes of filters (vega)
        required mode: 'abs_flux' or 'app_flux', change to 'abs_mag' or 'app_mag'

        """
        mode = self.getMode()
        if (mode[-4:] == 'flux'):
            for key, item in self.keys:
                if (issubclass(item, DictNpArrayMix)):
                    for ikey, itype in self[key].keys:
                        flux = self[key][ikey]
                        self[key][ikey] = -2.5 * np.log10(flux) - self[key].zero_point_mag[ikey]
            self.setMode(mode[:3]+'_mag')
        else:
            raise ValueError("values are not flux, current mode: %s " % mode)

