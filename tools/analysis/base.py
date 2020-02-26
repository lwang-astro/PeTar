# base class and functions
import numpy as np

class DictNpArrayMix:
    """ the Dictonary with numpy.ndarray function support
    """

    def __getitem__(self, k):
        """ Map getitem to all dictory np.ndarray items
        """
        cls_type = type(self)
        new_dat = cls_type()
        for key, item in self.__dict__.items():
            if (type(item) == np.ndarray):
                new_dat.__dict__[key] = item[k]
                new_dat.size = new_dat.__dict__[key].size
        return new_dat

def InitialDictNpArrayMixMethod(func):
    def initFunc(self, _dat=0):
        """
        _dat: np.ndarray type data reading from np.ndarray data or same class type 
        """
        if (type(_dat) == type(self)):
            self = _dat.copy()
        elif (type(_dat)==np.ndarray):
            keys = func(self)
            icol = 0
            for key, dimension in keys:
                if (dimension==1):
                    self.__dict__[key] = _dat[:,icol]
                else:
                    self.__dict__[key] = _dat[:,icol:icol+dimension]
                icol += dimension
            self.ncols = int(icol)
            self.size  = int(_dat.size/icol)

        elif (_dat==0):
            keys = func(self)
            icol = 0
            for key, dimension in keys:
                self.__dict__[key] = np.empty([0,dimension])
                icol += dimension
            self.ncols = int(icol)
            self.size  = int(0)
        else:
            raise ValueError('Initial fail, date type should be ',type(self),' or np.ndarray, given ',type(_dat))

    return initFunc

def join(*_dat):
    """
    Join multiple data to one
    """
    type0 = type(_dat[0])
    for idat in _dat:
        if (type(idat) != type0):
            raise ValueError('Initial fail, date type not consistent, type [0] is ',type0,' given ',type(idat))
    new_dat = type0()
    for key, item in new_dat.__dict__.items():
        if (type(item) == np.ndarray):
            new_dat.__dict__[key] = np.concatenate(tuple(map(lambda x:x.__dict__[key], _dat)))
    new_dat.size = np.sum(tuple(map(lambda x:x.size, _dat)))
    return new_dat

