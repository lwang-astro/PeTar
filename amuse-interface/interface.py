from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravityFieldInterface
from amuse.community.interface.gd import GravityFieldCode

class petarInterface(CodeInterface,
                     LiteratureReferencesMixIn,
                     GravitationalDynamicsInterface,
                     GravityFieldInterface):

    """
    Parallel, Particle-Particle & Particle-Tree & Few-body integration module
    
    .. [#] Iwasawa M., Tanikawa A., Hosono N., Nitadori K., Muranushi T., Makino J., 2016, PASJ, 68, 54
    .. [#] Iwasawa M., Portegies Zwart S., Makino J., 2015, ComAC, 2, 6
    .. [#] Wang, L., Nitadori, K., Makino, J. (2019, pre.)
    .. [#] Wang, L., Nitadori, K., Iwasawa, M., Makino, J. (2020, pre.)
    """
    
    include_headers = ['interface.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="petar_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        
    
class petar(GravitationalDynamics, GravityFieldCode):

    def __init__(self, convert_nbody = None, **keyword_arguments):
        GravitationalDynamics.__init__(self, 
                                       petarInterface(**keyword_arguments), 
                                       convert_nbody, 
                                       **keyword_arguments)