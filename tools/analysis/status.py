# analysis status data
from .base import *

class Status(DictNpArrayMix):
    """ Status data
    """
    def __init__(self, _dat=0):
        keys = [["time",1],["n_real_loc",1],["n_real_glb",1],["n_all_loc",1],["n_all_glb",1],["dE",1],["Ekin",1],["Epot",1],["Etot",1],["dE_SD",1],["Ekin_SD",1],["Epot_SD",1],["Etot_SD",1],["d|L|",1],["Lx",1],["Ly",1],["Lz",1],["|L|",1],["dE_hard",1],["dE_SD_hard",1],["CM.mass",1],["CM.pos.x",1],["CM.pos.y",1],["CM.pos.z",1],["CM.vel.x",1],["CM.vel.y",1],["CM.vel.z",1]]
        DictNpArrayMix.__init__(self, keys, _dat)
