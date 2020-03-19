import numpy as np
from .base import *
from .data import *

class Escaper:
    """ escaper information
    """

    def __init__(self):
        self.single = Particle()
        self.single.addNewMember('r2',np.empty(0))
        self.single.addNewMember('ekin',np.empty(0))
        self.single.addNewMember('etot',np.empty(0))
        self.single.addNewMember('time',np.empty(0))
        self.binary = Binary(member_particle_type=Particle)
        self.binary.addNewMember('r2',np.empty(0))
        self.binary.addNewMember('pot',np.empty(0))
        self.binary.addNewMember('ekin',np.empty(0))
        self.binary.addNewMember('etot',np.empty(0))
        self.binary.addNewMember('time',np.empty(0))
        self.rcut = 0

    def calcRCutIsolate(self,rh):
        self.rcut = 20*rh

    def findEscaper(self, time, single, binary, G=1):
        single.calcR2()
        single.calcEkin()
        single.calcEtot()

        binary.calcR2()
        binary.calcPot(G)
        binary.calcEkin()
        binary.calcEtot()

        rcut2 = self.rcut*self.rcut
        ssel = ((single.r2>rcut2) & (single.etot>0.0))
        single_esc = single[ssel]
        nssel = ssel.sum()
        single_esc.addNewMember('time',np.ones(nssel)*time)
        self.single.append(single_esc)
        unid, index= np.unique(self.single.id, return_index=True)
        self.single = self.single[index]

        bsel = (binary.r2>rcut2) & (binary.etot>0.0)
        binary_esc = binary[bsel]
        nbsel = bsel.sum()
        binary_esc.addNewMember('time',np.ones(nbsel)*time)
        self.binary.append(binary_esc)
        unid, index= np.unique(self.binary.p1.id, return_index=True)
        self.binary =self.binary[index]

def joinEscaper(*esc_list):
    esc_merge = Escaper()
    for ei in esc_list:
        esc_merge.rcut = ei.rcut
        esc_merge.single.append(ei.single)
        esc_merge.binary.append(ei.binary)
    unid, index= np.unique(esc_merge.single.id, return_index=True)
    esc_merge.single = esc_merge.single[index]

    unid, index= np.unique(esc_merge.binary.p1.id, return_index=True)
    esc_merge.binary = esc_merge.binary[index]

    return esc_merge
