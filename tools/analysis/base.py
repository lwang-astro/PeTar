# base class and functions
import numpy as np

class DictNpArrayMix:
    """ the Dictonary with numpy.ndarray function support
    """
    def __init__(self, keys, _dat=None, _offset=int(0), _append=False, **kwargs):
        """
        _dat: np.ndarray type data reading from np.ndarray data or same class type 
        _offset: reading column offset of _dat if it is np.ndarray
        _append: whether append keys and ncols
        """
        self.initargs = kwargs.copy()

        if (_append): self.keys = self.keys + keys
        else: self.keys = keys.copy()
        if (type(_dat) == type(self)):
            self = _dat.copy()
        elif (issubclass(type(_dat), DictNpArrayMix)):
            icol = int(0)
            for key, parameter in self.keys:
                if (type(parameter) == type):
                    if (issubclass(parameter, DictNpArrayMix)):
                        self.__dict__[key] = parameter(_dat.__dict__[key], **kwargs)
                        icol += self.__dict__[key].ncols
                    else:
                        raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given ',parameter)
                elif (type(parameter)==int):
                    self.__dict__[key] = _dat.__dict__[key].copy()
                    icol += parameter
                else:
                    raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
            if (_append): self.ncols += int(icol)
            else: self.ncols = int(icol)
            self.size  = _dat.size
        elif (type(_dat)==np.ndarray):
            icol = _offset
            self.size = int(0)
            for key, parameter in keys:
                if (type(parameter) == type):
                    if (issubclass(parameter, DictNpArrayMix)):
                        self.__dict__[key] = parameter(_dat, icol, False, **kwargs)
                        icol += self.__dict__[key].ncols
                    else:
                        raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given ',parameter)
                    self.size += self.__dict__[key].size*self.__dict__[key].ncols
                elif (type(parameter)==int):
                    if (parameter==1):
                        self.__dict__[key] = _dat[:,icol]
                    else:
                        self.__dict__[key] = _dat[:,icol:icol+parameter]
                    icol += parameter
                    self.size += self.__dict__[key].size
                else:
                    raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
            icol -= _offset
            if (_append): self.ncols += int(icol)
            else: self.ncols = int(icol)
            self.size  = int(self.size/icol)
            if (self.size != _dat.shape[0]):
                raise ValueError('Reading error, final counted size ',self.size,' is not consistent with reading ndarray shape',_dat.shape[0])
        elif (_dat==None):
            icol = int(0)
            for key, parameter in keys:
                if (type(parameter) == type):
                    if (issubclass(parameter, DictNpArrayMix)):
                        self.__dict__[key] = parameter(**kwargs)
                        icol += self.__dict__[key].ncols
                    else:
                        raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given ',parameter)
                elif (type(parameter)==int):
                    if (parameter==1): self.__dict__[key] = np.empty(0)
                    else: self.__dict__[key] = np.empty([0,parameter])
                    icol += parameter
                else:
                    raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
            if (_append): self.ncols += int(icol)
            else: self.ncols = int(icol)
            self.size  = int(0)
        else:
            raise ValueError('Initial fail, date type should be ',type(self),' or np.ndarray, given ',type(_dat))

    def readArray(self, _dat, _offset=int(0),**kwargs):
        icol = _offset
        self.size = int(0)
        for key, parameter in self.keys:
            if (type(parameter) == type):
                if (issubclass(parameter, DictNpArrayMix)):
                    self.__dict__[key].readArray(_dat, icol, **kwargs)
                    icol += self.__dict__[key].ncols
                else:
                    raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given ',parameter)
                self.size += self.__dict__[key].size*self.__dict__[key].ncols
            elif (type(parameter)==int):
                if (parameter==1):
                    self.__dict__[key] = _dat[:,icol]
                else:
                    self.__dict__[key] = _dat[:,icol:icol+parameter]
                icol += parameter
                self.size += self.__dict__[key].size
            else:
                raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
        icol -= _offset
        self.size  = int(self.size/icol)
        if (self.size != _dat.shape[0]):
            raise ValueError('Reading error, final counted size ',self.size,' is not consistent with reading ndarray shape',_dat.shape[0])
        if (self.ncols != icol):
            raise ValueError('Column number inconsistence, self ncols ',self.ncols,' key ncols ', icol)

    def __getitem__(self, k):
        """ Map getitem to all dictory np.ndarray items
        """
        if (type(k)==str):
            return self.__dict__[k]
        else:
            cls_type = type(self)
            new_dat = cls_type(**self.initargs)
            new_dat.ncols = self.ncols
            new_dat.size = int(0)
            new_dat.keys = self.keys.copy()
            icol = int(0)
            for key_type in new_dat.keys:
                key = key_type[0]
                item = self.__dict__[key]
                if (type(item) == np.ndarray):
                    if item.shape[0]!=self.size:
                        raise ValueError('Member ',key,' size/dimension',item.shape, ' is not consistent with the data size',self.size)
                    new_dat.__dict__[key] = item[k]
                    new_dat.size += new_dat.__dict__[key].size
                    if (len(item.shape)>1): icol += item.shape[1]
                    else: icol += 1
                elif (issubclass(type(item), DictNpArrayMix)):
                    new_item = item[k]
                    new_dat.__dict__[key] = new_item
                    new_dat.size += new_item.size*new_item.ncols
                    icol += new_item.ncols
            new_dat.size = int(new_dat.size/new_dat.ncols)
            if (icol != new_dat.ncols):
                raise ValueError('Column number inconsistent, coutned:',icol,' saved ncols:',new_dat.ncols)
            return new_dat

#    def keys(self):
#        return self.__dict__.keys()

    def addNewMember(self, key, member):
        new_key_flag=False
        if (key in self.__dict__.keys()):
            member_old = self.__dict__[key]
            dimension  = int(1)
            if (type(member_old)==np.ndarray):
                if len(member_old.shape)>1:
                    dimension = member_old.shape[1]
            elif (issubclass(type(member_old), DictNpArrayMix)):
                dimension = member_old.ncols
            self.ncols -= dimension
            member_old = member
        else:
            self.__dict__[key] = member
            new_key_flag=True
        dimension = int(1)
        if (type(member)==np.ndarray):
            if len(member.shape)>1:
                dimension = member.shape[1]
            if(new_key_flag): self.keys.append([key,dimension])
        elif (issubclass(type(member), DictNpArrayMix)):
            dimension = member.ncols
            if(new_key_flag): self.keys.append([key,type(member)])
        else:
            raise ValueError('New member type should be np.ndarray or DictNpArrayMix, but given ',type(member))
        self.ncols += dimension
        if (self.size != int(member.size/dimension)):
            raise ValueError('New member has different size: ',member.size/dimension, ' host size: ',self.size)
            
    def getherDataToArray(self):
        """ gether all data to 2D np.ndarray
        """
        dat_out=np.zeros([self.size,self.ncols])
        icol = int(0)
        for key_type in self.keys:
            key = key_type[0]
            member = self.__dict__[key]
            if (type(member)==np.ndarray):
                if len(member.shape)>1:
                    dimension= member.shape[1]
                    if (dat_out.shape[0]!=member.shape[0]):
                        raise ValueError('Member ',key,' size,dimension ',member.shape,' is not consistent with data output size',dat_out.shape[0])
                    for k in range(dimension):
                        dat_out[:,icol] = member[:,k]
                        icol += 1
                else:
                    dat_out[:,icol] = member
                    icol += 1
            elif (issubclass(type(member), DictNpArrayMix)):
                ncols = member.ncols
                dat_out[:,icol:icol+ncols] = member.getherDataToArray()
                icol += ncols
        return dat_out

    def append(self, *_dat):
        for idat in _dat:
            if (type(idat) != type(self)):
                raise ValueError('Initial fail, date type not consistent, type [0] is ',type(self),' given ',type(idat))
        data_with_self = [self]+list(_dat)
        for key, item in self.__dict__.items():
            if (type(item) == np.ndarray):
                self.__dict__[key] = np.concatenate(tuple(map(lambda x:x.__dict__[key], data_with_self)))
            elif(issubclass(type(item), DictNpArrayMix)):
                self.__dict__[key].append(*tuple(map(lambda x:x.__dict__[key], _dat)))
        self.size += np.sum(tuple(map(lambda x:x.size, _dat)))
                
    def savetxt(self, fname, **karg):
        dat_out= self.getherDataToArray()
        np.savetxt(fname, dat_out, **karg)

    def loadtxt(self, fname, **karg):
        dat_int = np.loadtxt(fname, **karg)
        self.readArray(dat_int)
        
def join(*_dat):
    """
    Join multiple data to one
    """
    type0 = type(_dat[0])
    for idat in _dat:
        if (type(idat) != type0):
            raise ValueError('Initial fail, date type not consistent, type [0] is ',type0,' given ',type(idat))
    new_dat = type0(**_dat[0].initargs)
    for key, item in _dat[0].__dict__.items():
        if (type(item) == np.ndarray):
            new_dat.__dict__[key] = np.concatenate(tuple(map(lambda x:x.__dict__[key], _dat)))
        elif(issubclass(type(item), DictNpArrayMix)):
            new_dat.__dict__[key] = join(*tuple(map(lambda x:x.__dict__[key], _dat)))
        else:
            new_dat.__dict__[key] = _dat[0].__dict__[key]
    new_dat.size = np.sum(tuple(map(lambda x:x.size, _dat)))
    return new_dat

# vector dot of x, y 
vec_dot = lambda x,y: np.sum(x*y,axis=1)
