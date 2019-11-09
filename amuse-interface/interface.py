from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamics
from amuse.community.interface.gd import GravitationalDynamicsInterface

class petarInterface(CodeInterface,
                     GravitationalDynamicsInterface):
    
    include_headers = ['interface.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="petar_worker", **keyword_arguments)
        
    
class petar(GravitationalDynamics):

    def __init__(self, convert_nbody = None, **keyword_arguments):
        GravitationalDynamics.__init__(self, 
                                       petarInterface(**keyword_arguments), 
                                       convert_nbody, 
                                       **keyword_arguments)
