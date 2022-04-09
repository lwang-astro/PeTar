# base class 
import numpy as np
import warnings

class DictNpArrayMix:
    """ The basic class of data structure
        The member functions are initialized by provided keys in initial function
        Member functions can be accessed by using the stype of either Dictonary or numpy.ndarray
    """
    def __init__(self, keys, _dat=None, _offset=int(0), _append=False, **kwargs):
        """
        Parameters
        ----------
        keys: list of class member name and the corresponding types or numpy.ndarray shape
            Class members list description. Defined by inherited types
            For exmaple: keys=[['mass',numpy.float64],['pos',(numpy.float64,3)],['sub1',typename],['sub2',(typename,kwargs)]], will provide class members: mass (1D numpy.ndarray), pos ( 2D numpy.ndarray with a shape of (*,3)), sub1 (a class instance with the type of typename) and sub2 (a type based on DictNpArrayMix with additional keyword arguments, kwargs)
        _dat: numpy.ndarray | type(self) | None
            If it is 2D numpy.ndarray type data, read data as readArray function
            If it is the same class type, copy the data 
            If it is None, initial class with empty data
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        _append: bool (False)
            If true, append keys and ncols to the current class instead of create new class members
        kwargs: dict ()
            keyword arguments, defined by inherited types
        """
        self.initargs = kwargs.copy()

        if (_append): self.keys = self.keys + keys
        else: self.keys = keys.copy()
        if (issubclass(type(_dat), DictNpArrayMix)):
            icol = int(0)
            for key, parameter in keys:
                if (type(parameter) == type):
                    if (issubclass(parameter, DictNpArrayMix)):
                        self.__dict__[key] = parameter(_dat.__dict__[key], **kwargs)
                        icol += self.__dict__[key].ncols
                    else:
                        self.__dict__[key] = _dat.__dict__[key].copy()
                        icol += 1
                        #raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given ',parameter)
                elif (type(parameter)==tuple):
                    if (type(parameter[0]) == type) & (type(parameter[1]) == int):
                        self.__dict__[key] = _dat.__dict__[key].copy()
                        icol += parameter[1]
                    elif (type(parameter[0]) == type) & (type(parameter[1])==dict):
                        if(issubclass(parameter[0], DictNpArrayMix)):
                            self.__dict__[key] = parameter[0](_dat.__dict__[key], **{**kwargs, **parameter[1]})
                            icol += self.__dict__[key].ncols
                        else:
                            self.__dict__[key] = _dat.__dict__[key].copy()
                            icol += 1
                    else:
                        raise ValueError('Initial fail, unknown key type ',parameter[0],' and column count ', parameter[1] )
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
                        self.size += self.__dict__[key].size*self.__dict__[key].ncols
                    else:
                        self.__dict__[key] = _dat[:,icol].astype(parameter)
                        icol += 1
                        self.size += self.__dict__[key].size
                        #raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given ',parameter)
                elif (type(parameter)==tuple):
                    if (type(parameter[0]) == type) & (type(parameter[1]) == int):
                        self.__dict__[key] = _dat[:,icol:icol+parameter[1]].astype(parameter[0])
                        icol += parameter[1]
                        self.size += self.__dict__[key].size
                    elif (type(parameter[0]) == type) & (type(parameter[1])==dict):
                        if(issubclass(parameter[0], DictNpArrayMix)):
                            self.__dict__[key] = parameter[0](_dat, icol, False, **{**kwargs, **parameter[1]})
                            icol += self.__dict__[key].ncols
                            self.size += self.__dict__[key].size*self.__dict__[key].ncols
                        else:
                            self.__dict__[key] = _dat[:,icol].astype(parameter[0])
                            icol += 1
                            self.size += self.__dict__[key].size
                    else:
                        raise ValueError('Initial fail, unknown key type ',parameter[0],' and column count ', parameter[1] )
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
                        self.__dict__[key] = np.empty(0).astype(parameter)
                        icol += 1
                        #raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given b',parameter)
                elif (type(parameter)==tuple):
                    if (type(parameter[0]) == type) & (type(parameter[1]) == int):
                        self.__dict__[key] = np.empty([0,parameter[1]]).astype(parameter[0])
                        icol += parameter[1]
                    elif (type(parameter[0]) == type) & (type(parameter[1])==dict):
                        if(issubclass(parameter[0], DictNpArrayMix)):
                            self.__dict__[key] = parameter[0](**{**kwargs, **parameter[1]})
                            icol += self.__dict__[key].ncols
                        else:
                            self.__dict__[key] = np.empty(0).astype(parameter[0])
                            icol += 1
                    else:
                        raise ValueError('Initial fail, unknown key type ',parameter[0],' and column count ', parameter[1] )
                else:
                    raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
            if (_append): self.ncols += int(icol)
            else: self.ncols = int(icol)
            self.size  = int(0)
        else:
            raise ValueError('Initial fail, date type should be ',type(self),' or np.ndarray, given ',type(_dat))

    def readArray(self, _dat, _offset=int(0), ncol_check=True, **kwargs):
        """ Read class member data from a 2D numpy.ndarray
        Parameters
        ----------
        _dat: numpy.ndarray 
            Read 2D array, rows are the event, columns are members. The class members are filled in the order of items in keys provided in the initial function.
            For exmaple: if keys are [['mass',1],['pos',3]], the member mass = _dat[:,_offset] and pos = _dat[:,_offset+1:_offset+3]
        _offset: int (0)
            Reading column offset of _dat if it is 2D np.ndarray
        ncol_check: bool (True)
            If True, check whether self.ncols (the number of columns in class) - _offset == _dat.shape[1] (the number of columns). If not equal, output warning
        kwaygs: dict ()
            keyword arguments
        """
        if (self.ncols + _offset != _dat.shape[1]) & (ncol_check):
            warnings.warn('The reading data shape[1] or the number of columns (%d) mismatches the number of columns defined in the class instance (%d)! Make sure whether this is intended and whether you properly choose the correct keyword arguments for the class instance initialziation' % (_dat.shape[1], self.ncols+_offset))

        icol = _offset
        self.size = int(0)
        for key, parameter in self.keys:
            if (type(parameter) == type):
                if (issubclass(parameter, DictNpArrayMix)):
                    self.__dict__[key].readArray(_dat, icol, False, **kwargs)
                    icol += self.__dict__[key].ncols
                    self.size += self.__dict__[key].size*self.__dict__[key].ncols
                else:
                    self.__dict__[key] = _dat[:,icol].astype(parameter)
                    icol += 1
                    self.size += self.__dict__[key].size
                    #raise ValueError('Initial fail, unknown key type, should be inherience of  DictNpArrayMix, given ',parameter)
            elif (type(parameter)==tuple):
                if (type(parameter[0]) == type) & (type(parameter[1]) == int):
                    self.__dict__[key] = _dat[:,icol:icol+parameter[1]].astype(parameter[0])
                    icol += parameter[1]
                    self.size += self.__dict__[key].size
                elif (type(parameter[0]) == type) & (type(parameter[1])==dict):
                    if(issubclass(parameter[0], DictNpArrayMix)):
                        self.__dict__[key].readArray(_dat, icol, False, **{**kwargs,**parameter[1]})
                        icol += self.__dict__[key].ncols
                        self.size += self.__dict__[key].size*self.__dict__[key].ncols
                    else:
                        self.__dict__[key] = _dat[:,icol].astype(parameter[0])
                        icol += 1
                        self.size += self.__dict__[key].size
                else:
                    raise ValueError('Initial fail, unknown key type ',parameter[0],' and column count ', parameter[1] )
            else:
                raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
        icol -= _offset
        self.size  = int(self.size/icol)
        if (self.size != _dat.shape[0]):
            raise ValueError('Reading error, final counted size ',self.size,' is not consistent with reading ndarray shape',_dat.shape[0])
        if (self.ncols != icol):
            raise ValueError('Column number inconsistence, self ncols ',self.ncols,' key ncols ', icol)

    def __getitem__(self, k):
        """ Map getitem to all members generated from the keys in the initial function, and return a new data filtered by k
        If the member is an inherited type of DictNpArrayMix, also map all sub-members if it.

        Parameters
        ----------
        k: filter
            The same type of arguments for numpy.ndarray.__getitem__
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
                raise ValueError('Column number inconsistent, counted:',icol,' saved ncols:',new_dat.ncols,'keys:',new_dat.keys,'fileter: ',k,' original size:',self.size,' original ncols:',self.ncols)
            return new_dat

    def __setitem__(self, k, data):
        """ Map setitem to all members generated from the keys in the initial function
        If the member is an inherited type of DictNpArrayMix, also map all sub-members if it.

        Parameters
        ----------
        k: filter
            The same type of arguments for numpy.ndarray.__getitem__
        data: numpy.ndarray | DictNpArrayNix
            The new data to set
        """
        if (type(k)==str):
            self.__dict__[k] = data
        else:
            for key_type in self.keys:
                key = key_type[0]
                self.__dict__[key][k] = data[key]

#    def keys(self):
#        return self.__dict__.keys()

    def addNewMember(self, key, member):
        """ Add a new class member
            The ncols is updated also.
            Be careful if the target for adding members is a sub member, the ncols of its parent is not updated.
            This can cause issue for the parent when the size of data is needed to calculate.
            Thus after calling of this function, please also increase the ncols of parents for consistence.

        Parameters
        ----------
        key: string
            new member name
        member: numpy.ndarray | DictNpArrayNix
            data binding to the member, should be the same size as existing members in the class
        """
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
                if(new_key_flag): self.keys.append([key,(type(member[:,0]),dimension)])
            else:
                if(new_key_flag): self.keys.append([key,type(member)])
        elif (issubclass(type(member), DictNpArrayMix)):
            dimension = member.ncols
            if(new_key_flag): self.keys.append([key,type(member)])
        else:
            raise ValueError('New member type should be np.ndarray or DictNpArrayMix, but given ',type(member))
        self.ncols += dimension
        if (self.size != int(member.size/dimension)):
            raise ValueError('New member has different size: ',member.size/dimension, ' host size: ',self.size)
            
    def getherDataToArray(self):
        """ gether all data to a 2D numpy.ndarray and return it
        An inverse function to readArray
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

    def printTable(self, column_format):
        """
        print Table with defined column list and formats

        Parameters
        column_format: a list of column label (class member name) and format, enclosed by tuple, for sub-member, use . to access
                       For exmaple: [(key1,'%s'), (key2,'%12.7f'), (key3.subkey1,'%d'), (key3.subkey2,'%e')]
        """
        if self.size==0: return
        table=[]
        fmt_list=''
        for key, fmt in column_format:
            keylst = key.split('.')
            dat_key=self
            for ikey in keylst:
                if '[' in ikey:
                    key_name, index = ikey.split('[')
                    index = int(index.split(']')[0])
                    dat_key = dat_key[key_name][:,index]
                else:
                    dat_key = dat_key[ikey]
            table.append(dat_key)
            fmt_list += fmt
        table=np.transpose(np.array(table))
        for line in table:
            print(fmt_list % tuple(line))

    def append(self, *_dat):
        """ Map the numpy.append function to each member
        Append the numpy.ndarray of each member in a group of input data to the corresponding member in self

        Parameters
        *_dat: inherited DictNpArrayMix
            The data should contain all members existing in the self 
        """
        #for idat in _dat:
        #    if (type(idat) != type(self)):
        #        raise ValueError('Initial fail, date type not consistent, type [0] is ',type(self),' given ',type(idat))
        data_with_self = [self]+list(_dat)
        for key, item in self.__dict__.items():
            if (type(item) == np.ndarray):
                if (len(item.shape)!=len(_dat[0][key].shape)):
                    raise ValueError('Appending data member ',key,' has shape',_dat[0][key].shape,' but self data has shape',item.shape)
                self.__dict__[key] = np.concatenate(tuple(map(lambda x:x.__dict__[key], data_with_self)))
            elif(issubclass(type(item), DictNpArrayMix)):
                self.__dict__[key].append(*tuple(map(lambda x:x.__dict__[key], _dat)))
        self.size += np.sum(tuple(map(lambda x:x.size, _dat)))
                
    def savetxt(self, fname, **kwargs):
        """ Save class member data to a file
        Use the getherDataToArray and then numpy.savetxt

        Parameters
        ----------
        fname: string
            name of the output file
        kwargs: dict
            keyword arguments for numpy.savetxt
        """
        dat_out= self.getherDataToArray()
        np.savetxt(fname, dat_out, **kwargs)

    def loadtxt(self, fname, **kwargs):
        """ Load class member data from a file
        Use numpy.loadtxt to read data and then use readArray

        Parameters
        ----------
        fname: string
            name of the input file
        kwargs: dict
            keyword arguments for numpy.loadtxt
        """
        dat_int = np.loadtxt(fname, ndmin=2, **kwargs)
        if (dat_int.size>0):
            self.readArray(dat_int, **kwargs)

    def collectDtype(self):
        """ Collect dtype from keys iteratively for reading BINARY format
        For member with type of DictNpArrayMix, use column name with the prefix of member name + '.'
        """
        dt=[]
        for key, parameter in self.keys:
            if (type(parameter) == type):
                if (issubclass(parameter, DictNpArrayMix)):
                    dt_sub = self.__dict__[key].collectDtype()
                    for item in dt_sub:
                        dt.append((key+'.'+item[0], item[1]))
                else:
                    dt.append((key, parameter))
            elif (type(parameter) == tuple):
                if (type(parameter[0]) == type) & (type(parameter[1]) == int):
                    dt.append((key, parameter))
                elif (type(parameter[0]) == type) & (type(parameter[1])==dict):
                    if(issubclass(parameter[0], DictNpArrayMix)):
                        dt_sub = self.__dict__[key].collectDtype()
                        for item in dt_sub:
                            dt.append((key+'.'+item[0], item[1]))
                    else:
                        dt.append((key, parameter[0]))
                else:
                    raise ValueError('Initial fail, unknown key type ',parameter[0],' and column count ', parameter[1] )
            else:
                raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
        return dt

    def readArrayWithName(self, _dat, _prefix='', **kwargs):
        """ Read class member data from a numpy.ndarray with names
        Parameters
        ----------
        _dat: numpy.ndarray 
            Read array with names of columns (key/member name of class)
        _prefix: string ('')
            The prefix add in front of the name of columns to read. This is used when the current class instance is a sub member (member name is consistent with prefix)
        kwaygs: dict ()
            keyword arguments
        """
        icol = int(0)
        self.size = int(0)
        for key, parameter in self.keys:
            if (type(parameter) == type):
                if (issubclass(parameter, DictNpArrayMix)):
                    self.__dict__[key].readArrayWithName(_dat, _prefix+key+'.', **kwargs)
                    icol += self.__dict__[key].ncols
                    self.size += self.__dict__[key].size*self.__dict__[key].ncols
                else:
                    self.__dict__[key] = _dat[_prefix+key]
                    icol += 1
                    self.size += self.__dict__[key].size
            elif (type(parameter) == tuple):
                if (type(parameter[0]) == type) & (type(parameter[1]) == int):
                    self.__dict__[key] = _dat[_prefix+key]
                    icol += parameter[1]
                    self.size += self.__dict__[key].size
                elif (type(parameter[0]) == type) & (type(parameter[1])==dict):
                    if(issubclass(parameter[0], DictNpArrayMix)):
                        self.__dict__[key].readArrayWithName(_dat, _prefix+key+'.', **kwargs,**parameter[1])
                        icol += self.__dict__[key].ncols
                        self.size += self.__dict__[key].size*self.__dict__[key].ncols
                    else:
                        self.__dict__[key] = _dat[_prefix+key]
                        icol += 1
                        self.size += self.__dict__[key].size
                else:
                    raise ValueError('Initial fail, unknown key type ',parameter[0],' and column count ', parameter[1] )
            else:
                raise ValueError('Initial fail, unknown key parameter, should be DictNpArrayMix type name or value of int, given ',parameter)
        self.size = int(self.size/icol)
        if (self.size != _dat.size):
            raise ValueError('Reading error, final counted size ',self.size,' is not consistent with reading ndarray shape',_dat.size)
        if (self.ncols != icol):
            raise ValueError('Column number inconsistence, self ncols ',self.ncols,' key ncols ', icol)

    def fromfile(self, fname, **kwargs):
        """ Load clas member data from a file using BINARY format
        Use numpy.fromfile to read data, the dtype is defined by keys (members)
        Notice if the first line is header, offset counts in byte should be used

        Parameters
        ----------
        fname: string
            name of the input file
        kwargs: dict
            keyword arguments for numpy.fromfile, notice dtype is already defined, do not provide that
        """
        dt = self.collectDtype()
        dat_int = np.fromfile(fname, dtype=dt, **kwargs)
        self.readArrayWithName(dat_int, '', **kwargs)

    def printSize(self):
        """ print size of each member
        Print the shape of each members, used for testing whether the members have consistent size
        """
        for key_type in self.keys:
            key = key_type[0]
            member = self.__dict__[key]
            if (type(member)==np.ndarray):
                print(key,member.shape)
            elif (issubclass(type(member), DictNpArrayMix)):
                member.printSize()
                
                
        
def join(*_dat):
    """ Join multiple data to one
    For a list of data with the same type of inherited DictNpArrayNix, this function join all data to one 
    
    Parameters
    ----------
    *_dat: inherited DictNpArrayNix
        a group of data to join

    return: new joined data
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

