import collections
import numpy as np
from .base import vec_dot

def Lagrangian(single, binary, radius_fraction, calc_mode):
    """ Calculate Langrangian radii and related properties
    """    
    n_frac = int(0)
    if (type(radius_fraction)==list):
        n_frac = len(radius_fraction)
    elif (type(radius_fraction)==np.ndarray):
        n_frac = radius_fraction.size

    # lagrangian radii
    lagr = collections.OrderedDict()

    keys      = ['r','m','n','vel','sigma'] # radius, mass, number, velocity, dispersion
    vel_keys  = ['3d','x','y','z','rad','tan','rot'] # all, x, y, z, radial, tangential, rotational
    comp_keys = ['s','b','all'] # single, binary, all

    for key in keys:
        lagr[key] = collections.OrderedDict()
        for ckey in comp_keys:
            lagr[key][comp_keys] = np.zeros(n_frac)

    for key in ['vel','sigma']:
        for vkey in vel_keys:
            lagr[key][vkey] = collections.OrderedDict()
            for ckey in comp_keys:
                lagr[key][vkey][ckey] = np.zeros(n_frac)

    total_single_mass = single.mass.sum()
    total_binary_mass = binary.m1.sum()+binary.m2.sum()
    
    r2_single = vec_dot(single.pos,single.pos)
    r2_binary = vec_dot(binary.pos,binary.pos)
        rb2 = 0
        r2 = rs2

        if (bflag): 
            tbmass = float(bm1.sum() + bm2.sum())
            rb2 = bxc1*bxc1 + bxc2*bxc2 + bxc3*bxc3
            r2 = np.append(r2,rb2)

        rm2 = 0
        if (mflag):
            tmmass = float(mm1.sum() + mm2.sum() + mm3.sum())
            rm2 = mxc1*mxc1 + mxc2*mxc2 + mxc3*mxc3
            r2 = np.append(r2,rm2)
#        rb2 = np.append(rb2,rm2)

# Total mass
        tmass = tsmass + tbmass + tmmass
# Treat triple as binary
        tbmass = tbmass + tmmass

#  mass limit for different R_lagr
        rmass = tmass*rfrac
        rsmass = tsmass*rfrac
        rbmass = tbmass*rfrac
        
#  Get distance sorting index
        idx = r2.argsort()

# counter for different R_lagr bins
        kk = 0
        kks = 0
        kkb = 0

#   Number counter
        nc = 0
        ncs = 0
        ncb = 0
        #   Previous counter
        ncprev = 0
        ncsprev = 0
        ncbprev = 0

#   Mass counter
        cmass = 0
        bmass = 0
        smass = 0
        #   Previous counter
        mcprev = 0
        mcsprev = 0
        mcbprev = 0

#  Single and binary number
        N_SB = N_SINGLE + N_BINARY
#  binary and merger
        N_BM = N_BINARY + N_MERGER
        # Resolved case
        if (fbres): N_BM = N_BINARY*2 + N_MERGER*3
#  Total number 
        N_TOT = N_SB + N_MERGER
        # For total counter, unresolved and resolved case 
        N_TOTR = N_BM + N_SINGLE
#  Initialize the velocity vectors
        # velocity, x,y,z,r,t,rot
        vx = np.zeros(N_TOT)
        vy = np.zeros(N_TOT)
        vz = np.zeros(N_TOT)
        # The radial velocity uses the spherical coordinate system, thus not a vector
        vr = np.zeros(N_TOT)
        # Notice the tangential velocity is also vector in (x,y,z) coordinate system
        vt = np.zeros((N_TOT,3))
        # Average tangetial velocity
        vtave = np.zeros((rfrac.size,3))
        vtbave = np.zeros((rfrac.size,3))
        vtsave = np.zeros((rfrac.size,3))
        # rotational velocity
        vrot = np.zeros(N_TOT)
#   Initialize mass array
        mmb = np.zeros(N_TOT)

        for j in idx:
#  Initialization
            # position
            rx = 0.0
            ry = 0.0
            rz = 0.0
            ri = math.sqrt(r2[j])

#   Binary/merger case
            if (j>=N_SINGLE):
                # increase binary counter by one
                ncb += 1
                if (fbres):
                    # increase total counter by one in binary case
                    nc += 1
                    # increase binary counter by one for resolved case
                    ncb += 1
#   Merger case
                
                if (j>=N_SB):
                    if (fbres):
                        # increase total counter by one in merger case
                        nc += 1
                        # increase binary counter by one for resolved case
                        ncb += 1
                    inx = j - N_SB
                    mmb[j] = mm1[inx]+mm2[inx]+mm3[inx]
                    vx[j] = mvc1[inx]
                    vy[j] = mvc2[inx]
                    vz[j] = mvc3[inx]
                    rx    = mxc1[inx]
                    ry    = mxc2[inx]
                    rz    = mxc3[inx]
#   Binary case
                else:
                    inx = j - N_SINGLE
                    mmb[j] = bm1[inx]+bm2[inx]
                    vx[j] = bvc1[inx]
                    vy[j] = bvc2[inx]
                    vz[j] = bvc3[inx]
                    rx    = bxc1[inx]
                    ry    = bxc2[inx]
                    rz    = bxc3[inx]

#########--debug-------------------------
###            btt=bxc1[inx]*bxc1[inx] + bxc2[inx]*bxc2[inx]+ bxc3[inx]*bxc3[inx]
###            if (btt == r2[j]): 
###                print "Inconsistence: inx=%d, j=%d, r2=%e, r2s=%e, " % (inx,j,btt,r2[j])
#########--debug-------------------------

#   Single case
            else:
                # increase number counter
                ncs += 1
                mmb[j] = mass[j]
                rx     = x1[j]
                ry     = x2[j]
                rz     = x3[j]
                vx[j]  = v1[j]
                vy[j]  = v2[j]
                vz[j]  = v3[j]

            #  increase total counter
            nc += 1

#   Get velocity information
            # x-y plane radial velocity * rxy
            rvxy = rx*vx[j] + ry*vy[j]
            # radial velocity value
            vr[j] = (rvxy + rz*vz[j])/ri
            # tangential velocity vector
            vt[j,0] = vx[j] - vr[j]*rx/ri
            vt[j,1] = vy[j] - vr[j]*ry/ri
            vt[j,2] = vz[j] - vr[j]*rz/ri
            # x-y plane radial position square
            rxy2 = rx*rx + ry*ry
            # rotational velocity
            vrot1 = vx[j] - rvxy*rx/rxy2
            vrot2 = vy[j] - rvxy*ry/rxy2
            vrot[j] = math.sqrt(vrot1*vrot1+vrot2*vrot2)
            # rotational direction sign
            vrotd = vrot1*ry - vrot2*rx
            if (vrotd<0.0): vrot[j] = -vrot[j]

#   Check whether need to reach next bin
            if (j>=N_SINGLE):
                #  Add mass
                bmass += mmb[j]

                # increase mass/number counter for binaries in R_lagr
                if (kk < rfrac.size):
                    msblagr[kk] += mmb[j]
                    nsblagr[kk] += 1
                    if (fbres): nsblagr[kk] += 1

                    # primordial binareis
                    if (j<N_SB):
                        if (abs(bn1[j-N_SINGLE]-bn2[j-N_SINGLE])==1):
                            mspblagr[kk] += mmb[j]
                            nspblagr[kk] += 1
                            if (fbres): nspblagr[kk] += 1

                if (kkb < rfrac.size):
                    #  average velocity
                    vxblagr[kkb] += mmb[j]*vx[j] 
                    vyblagr[kkb] += mmb[j]*vy[j] 
                    vzblagr[kkb] += mmb[j]*vz[j] 
                    vrblagr[kkb] += mmb[j]*vr[j] 
                    vtbave[kkb] += mmb[j]*vt[j] 
                    vrotblagr[kkb] += mmb[j]*vrot[j]

                    # Go to next bin if mass reach the R_Lagr limit
                    if ((bmass >= rbmass[kkb]) | ((kkb == rfrac.size-1) & (ncb == N_BM))):
                        # update mass
                        rbmass[kkb] = bmass
                        # Get R_lagr for binary
                        rblagr[kkb] = ri
                        # Get number for binary
                        nblagr[kkb] = ncb
                        # For shell cases:
                        if ((fshell) & (kkb>0)):
                            rbmass[kkb] -= mcbprev
                            nblagr[kkb] -= ncbprev
                        mcbprev = bmass
                        ncbprev = ncb
                        # Increase bin index
                        kkb += 1
                        # initial next bins
                        if (kkb < rfrac.size):
                            if ((not fshell) & (nblagr[kkb] == 0)):
                                vxblagr[kkb] = vxblagr[kkb-1]
                                vyblagr[kkb] = vyblagr[kkb-1]
                                vzblagr[kkb] = vzblagr[kkb-1]
                                vrblagr[kkb] = vrblagr[kkb-1]
                                vtbave[kkb] = vtbave[kkb-1]
                                vrotblagr[kkb] = vrotblagr[kkb-1]
###########---debug--------------
#                if(kkb > rfrac.size - 1): print j , ncb, N_BINARY, bmass, tbmass, rbmass[rfrac.size-1],time
###########---debug--------------
                        # Avoid overflow of bin index
                        # kkb = min(kkb,rfrac.size-1)
            else:
                # Add mass
                smass += mmb[j]

                if (kks < rfrac.size):
                    #  average velocity
                    vxslagr[kks] += mmb[j]*vx[j] 
                    vyslagr[kks] += mmb[j]*vy[j] 
                    vzslagr[kks] += mmb[j]*vz[j] 
                    vrslagr[kks] += mmb[j]*vr[j] 
                    vtsave[kks] += mmb[j]*vt[j] 
                    vrotslagr[kks] += mmb[j]*vrot[j]

                    # Go to next bin if mass reach the R_lagr limit
                    if ((smass >= rsmass[kks]) | ((kks == rfrac.size-1) & (ncs == N_SINGLE))):
                        # update mass
                        rsmass[kks] = smass
                        # Get R_lagr for single
                        rslagr[kks] = ri
                        # Get number for single
                        nslagr[kks] = ncs
                        # For shell cases:
                        if ((fshell) & (kks>0)):
                            rsmass[kks] -= mcsprev 
                            nslagr[kks] -= ncsprev
                        mcsprev = smass
                        ncsprev = ncs
                        # increase bin index
                        kks += 1
                        # initial next bins
                        if (kks < rfrac.size):
                            if((not fshell) & (nslagr[kks] == 0)):
                                vxslagr[kks] = vxslagr[kks-1]
                                vyslagr[kks] = vyslagr[kks-1]
                                vzslagr[kks] = vzslagr[kks-1]
                                vrslagr[kks] = vrslagr[kks-1]
                                vtsave[kks] = vtsave[kks-1]
                                vrotslagr[kks] = vrotslagr[kks-1]
                        # Avoid overflow of bin index
                        # kks = min(kks,rfrac.size-1)
            
#   Go to next R_lagr if mass reach limit
            cmass += mmb[j]

            if (kk < rfrac.size):
                #  average velocity
                vxlagr[kk] += mmb[j]*vx[j] 
                vylagr[kk] += mmb[j]*vy[j] 
                vzlagr[kk] += mmb[j]*vz[j] 
                vrlagr[kk] += mmb[j]*vr[j] 
                vtave[kk] += mmb[j]*vt[j] 
                vrotlagr[kk] += mmb[j]*vrot[j]

                if ((cmass >= rmass[kk]) | ((kk == rfrac.size-1) & (nc == N_TOTR))):
                    # update mass
                    rmass[kk] = cmass
                    # Get R_lagr 
                    rlagr[kk] = ri
                    # Get number
                    nlagr[kk] = nc
                    # For shell cases:
                    if ((fshell) & (kk>0)):
                        rmass[kk] -= mcprev
                        nlagr[kk] -= ncprev
                    mcprev = cmass
                    ncprev = nc
                    # increase bin index
                    kk += 1
                    # Get initial value for next bin 
                    if (kk < rfrac.size ):
                        if (not fshell):
                            # binary counter
                            if (nsblagr[kk] == 0):
                                msblagr[kk] = msblagr[kk-1]
                                nsblagr[kk] = nsblagr[kk-1]
                            if (nspblagr[kk] == 0):
                                mspblagr[kk] = mspblagr[kk-1]
                                nspblagr[kk] = nspblagr[kk-1]
                            # total counter
                            if (nlagr[kk] == 0):
                                vxlagr[kk] = vxlagr[kk-1]
                                vylagr[kk] = vylagr[kk-1]
                                vzlagr[kk] = vzlagr[kk-1]
                                vrlagr[kk] = vrlagr[kk-1]
                                vtave[kk] = vtave[kk-1]
                                vrotlagr[kk] = vrotlagr[kk-1]
                    # Avoid overflow of bin index
                    # kk = min(kk,rfrac.size-1)

#   Fill empty bins with neighbor bin values
        if (not fshell):
            #   Total
            kn = kk - 1
            while (kk < rfrac.size):
                rlagr[kk] = rlagr[kn]
                nlagr[kk] = nlagr[kn]
                nsblagr[kk] = nsblagr[kn]
                msblagr[kk] = msblagr[kn]
                nspblagr[kk] = nspblagr[kn]
                mspblagr[kk] = mspblagr[kn]
                vxlagr[kk] = vxlagr[kn]
                vylagr[kk] = vylagr[kn]
                vzlagr[kk] = vzlagr[kn]
                vrlagr[kk] = vrlagr[kn]
                vtave[kk] = vtave[kn]
                vrotlagr[kk] = vrotlagr[kn]
                kk += 1
            #   Single
            ksn = kks - 1
            while (kks < rfrac.size):
                rslagr[kks] = rslagr[ksn]
                nslagr[kks] = nslagr[ksn]
                vxslagr[kks] = vxslagr[ksn]
                vyslagr[kks] = vyslagr[ksn]
                vzslagr[kks] = vzslagr[ksn]
                vrslagr[kks] = vrslagr[ksn]
                vtsave[kks] = vtsave[ksn]
                vrotslagr[kks] = vrotslagr[ksn]
                kks += 1
            #   Binary + Merger
            kbn = kkb - 1
            while (kkb < rfrac.size):
                rblagr[kkb] = rblagr[kbn]
                nblagr[kkb] = nblagr[kbn]
                vxblagr[kkb] = vxblagr[kbn]
                vyblagr[kkb] = vyblagr[kbn]
                vzblagr[kkb] = vzblagr[kbn]
                vrblagr[kkb] = vrblagr[kbn]
                vtbave[kkb] = vtbave[kbn]
                vrotblagr[kkb] = vrotblagr[kbn]
                kkb += 1

#   Average mass
        mlagr  = np.array(map(fxovery,rmass,nlagr)  ) 
        mslagr = np.array(map(fxovery,rsmass,nslagr)) 
        mblagr = np.array(map(fxovery,rbmass,nblagr)) 
#   Average velocity
        # total
        vxlagr  = np.array(map(fxovery,vxlagr  ,rmass)) 
        vylagr  = np.array(map(fxovery,vylagr  ,rmass)) 
        vzlagr  = np.array(map(fxovery,vzlagr  ,rmass)) 
        vrlagr  = np.array(map(fxovery,vrlagr  ,rmass)) 
        vtave   = np.array(map(fxovery,vtave   ,rmass)) 
        vrotlagr= np.array(map(fxovery,vrotlagr,rmass)) 
        vlagr   = np.sqrt(vxlagr*vxlagr + vylagr*vylagr + vzlagr*vzlagr)
        vtlagr  = np.array(map(lambda x: math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]), vtave))
        #single
        vxslagr  = np.array(map(fxovery,vxslagr  ,rsmass)) 
        vyslagr  = np.array(map(fxovery,vyslagr  ,rsmass)) 
        vzslagr  = np.array(map(fxovery,vzslagr  ,rsmass)) 
        vrslagr  = np.array(map(fxovery,vrslagr  ,rsmass)) 
        vtsave   = np.array(map(fxovery,vtsave   ,rsmass)) 
        vrotslagr= np.array(map(fxovery,vrotslagr,rsmass)) 
        vslagr   = np.sqrt(vxslagr*vxslagr + vyslagr*vyslagr + vzslagr*vzslagr)
        vtslagr  = np.array(map(lambda x: math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]), vtsave))
        #binary/merger
        if(bflag):
            vxblagr  = np.array(map(fxovery,vxblagr  ,rbmass)) 
            vyblagr  = np.array(map(fxovery,vyblagr  ,rbmass)) 
            vzblagr  = np.array(map(fxovery,vzblagr  ,rbmass)) 
            vrblagr  = np.array(map(fxovery,vrblagr  ,rbmass)) 
            vtbave   = np.array(map(fxovery,vtbave   ,rbmass)) 
            vrotblagr= np.array(map(fxovery,vrotblagr,rbmass)) 
            vblagr   = np.sqrt(vxblagr*vxblagr + vyblagr*vyblagr + vzblagr*vzblagr)
            vtblagr  = np.array(map(lambda x: math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]), vtbave))
    
#   Loop again to get velocity dispersion
#   counter for different R_lagr bins
        kk = 0
        kks = 0
        kkb = 0
#   Number counter
        nc = 0
        ncs = 0
        ncb = 0
        #   Previous counter
        ncprev = 0
        ncsprev = 0
        ncbprev = 0
#
        for j in idx:
            #  increase total counter
            nc += 1
#   Binary/merger case
            if (j>=N_SINGLE):
                # increase binary counter by two 
                ncb += 1
                if (fbres):
                    # increase total counter by one in binary case
                    nc += 1
                    # increase binary counter by one for resolved case
                    ncb += 1
                    if (j>=N_SB):
                        # For merger case
                        nc  += 1
                        ncb += 1
                # x,y,z
                dx = vx[j] - vxblagr[kkb]
                dy = vy[j] - vyblagr[kkb]
                dz = vz[j] - vzblagr[kkb]
                dr = vr[j] - vrblagr[kkb]
                dt = vt[j] - vtbave[kkb]
                drot = vrot[j] - vrotblagr[kkb]
                # mass weighted square
                dx2 = mmb[j]*dx*dx
                dy2 = mmb[j]*dy*dy
                dz2 = mmb[j]*dz*dz
                dr2 = mmb[j]*dr*dr
                dt2 = mmb[j]*(dt[0]*dt[0] + dt[1]*dt[1] + dt[2]*dt[2])
                drot2= mmb[j]*drot*drot
                # velocity value
#                dv2 = dx2 + dy2 + dz2
                # dispersion
                sigxblagr[kkb] += dx2
                sigyblagr[kkb] += dy2
                sigzblagr[kkb] += dz2
#                sigblagr[kkb] += dv2
                sigrblagr[kkb] += dr2
                sigtblagr[kkb] += dt2
                sigrotblagr[kkb] += drot2

                # check whether need to increase kkb
                if (ncb-ncbprev==nblagr[kkb]):
                    if (fshell): ncbprev += nblagr[kkb]
                    kkb += 1
                    if ((not fshell) & (kkb < rfrac.size)):
                        sigxblagr[kkb] = sigxblagr[kkb-1]
                        sigyblagr[kkb] = sigyblagr[kkb-1]
                        sigzblagr[kkb] = sigzblagr[kkb-1]
                        sigrblagr[kkb] = sigrblagr[kkb-1]
                        sigtblagr[kkb] = sigtblagr[kkb-1]
                        sigrotblagr[kkb] = sigrotblagr[kkb-1]
                
#   Single case
            else:
                # increase number counter
                ncs += 1
                # x,y,z
                dx = vx[j] - vxslagr[kks]
                dy = vy[j] - vyslagr[kks]
                dz = vz[j] - vzslagr[kks]
                dr = vr[j] - vrslagr[kks]
                dt = vt[j] - vtsave[kks]
                drot = vrot[j] - vrotslagr[kks]
                # mass weighted square
                dx2 = mmb[j]*dx*dx
                dy2 = mmb[j]*dy*dy
                dz2 = mmb[j]*dz*dz
                dr2 = mmb[j]*dr*dr
                dt2 = mmb[j]*(dt[0]*dt[0] + dt[1]*dt[1] + dt[2]*dt[2])
                drot2= mmb[j]*drot*drot
                # velocity value
#                dv2 = dx2 + dy2 + dz2
                # dispersion
                sigxslagr[kks] += dx2
                sigyslagr[kks] += dy2
                sigzslagr[kks] += dz2
#                sigslagr[kks] += dv2
                sigrslagr[kks] += dr2
                sigtslagr[kks] += dt2
                sigrotslagr[kks] += drot2

                # check whether need to increase kks
                if (ncs-ncsprev==nslagr[kks]):
                    if (fshell): ncsprev += nslagr[kks]
                    kks += 1
                    if ((not fshell) & (kks < rfrac.size)):
                        sigxslagr[kks] = sigxslagr[kks-1]
                        sigyslagr[kks] = sigyslagr[kks-1]
                        sigzslagr[kks] = sigzslagr[kks-1]
                        sigrslagr[kks] = sigrslagr[kks-1]
                        sigtslagr[kks] = sigtslagr[kks-1]
                        sigrotslagr[kks] = sigrotslagr[kks-1]

#   Total                 
            # x,y,z
            dx = vx[j] - vxlagr[kk]
            dy = vy[j] - vylagr[kk]
            dz = vz[j] - vzlagr[kk]
            dr = vr[j] - vrlagr[kk]
            dt = vt[j] - vtave[kk]
            drot = vrot[j] - vrotlagr[kk]
            # mass weighted square
            dx2 = mmb[j]*dx*dx
            dy2 = mmb[j]*dy*dy
            dz2 = mmb[j]*dz*dz
            dr2 = mmb[j]*dr*dr
            dt2 = mmb[j]*(dt[0]*dt[0] + dt[1]*dt[1] + dt[2]*dt[2])
            drot2= mmb[j]*drot*drot
            # velocity value
 #           dv2 = dx2 + dy2 + dz2
            # dispersion
            sigxlagr[kk] += dx2
            sigylagr[kk] += dy2
            sigzlagr[kk] += dz2
#            siglagr[kk] += dv2
            sigrlagr[kk] += dr2
            sigtlagr[kk] += dt2
            sigrotlagr[kk] += drot2

            # check whether need to increase kk
            if (nc-ncprev==nlagr[kk]):
                if (fshell): ncprev += nlagr[kk]
                kk += 1
                if ((not fshell) & (kk < rfrac.size)):
                    sigxlagr[kk] = sigxlagr[kk-1]
                    sigylagr[kk] = sigylagr[kk-1]
                    sigzlagr[kk] = sigzlagr[kk-1]
                    sigrlagr[kk] = sigrlagr[kk-1]
                    sigtlagr[kk] = sigtlagr[kk-1]
                    sigrotlagr[kk] = sigrotlagr[kk-1]

        if (not fshell):
            kn = kk - 1
            while (kk < rfrac.size):
                sigxlagr[kk] = sigxlagr[kn]
                sigylagr[kk] = sigylagr[kn]
                sigzlagr[kk] = sigzlagr[kn]
                sigrlagr[kk] = sigrlagr[kn]
                sigtlagr[kk] = sigtlagr[kn]
                sigrotlagr[kk] = sigrotlagr[kn]
                kk += 1
            ksn = kks - 1
            while (kks < rfrac.size):
                sigxslagr[kks] = sigxslagr[ksn]
                sigyslagr[kks] = sigyslagr[ksn]
                sigzslagr[kks] = sigzslagr[ksn]
                sigrslagr[kks] = sigrslagr[ksn]
                sigtslagr[kks] = sigtslagr[ksn]
                sigrotslagr[kks] = sigrotslagr[ksn]
                kks += 1
            kbn = kkb - 1
            while (kkb < rfrac.size):
                sigxblagr[kkb] = sigxblagr[kbn]
                sigyblagr[kkb] = sigyblagr[kbn]
                sigzblagr[kkb] = sigzblagr[kbn]
                sigrblagr[kkb] = sigrblagr[kbn]
                sigtblagr[kkb] = sigtblagr[kbn]
                sigrotblagr[kkb] = sigrotblagr[kbn]
                kkb += 1

# Divide mass
        # total
#        siglagr   = siglagr   /(rmass*3.0)
        sigxlagr  = np.array(map(fxovery,sigxlagr  ,rmass)) 
        sigylagr  = np.array(map(fxovery,sigylagr  ,rmass)) 
        sigzlagr  = np.array(map(fxovery,sigzlagr  ,rmass)) 
        sigrlagr  = np.array(map(fxovery,sigrlagr  ,rmass)) 
        sigtlagr  = np.array(map(fxovery,sigtlagr  ,rmass)) 
        sigrotlagr= np.array(map(fxovery,sigrotlagr,rmass))
        siglagr   = sigxlagr + sigylagr + sigzlagr
        #single
#        sigslagr   = sigslagr   /(rsmass*3.0)
        sigxslagr  = np.array(map(fxovery,sigxslagr  ,rsmass)) 
        sigyslagr  = np.array(map(fxovery,sigyslagr  ,rsmass)) 
        sigzslagr  = np.array(map(fxovery,sigzslagr  ,rsmass)) 
        sigrslagr  = np.array(map(fxovery,sigrslagr  ,rsmass)) 
        sigtslagr  = np.array(map(fxovery,sigtslagr  ,rsmass)) 
        sigrotslagr= np.array(map(fxovery,sigrotslagr,rsmass)) 
        sigslagr   = sigxslagr + sigyslagr + sigzslagr
        #binary/merger
#        sigblagr   = sigblagr   /(rbmass*3.0)
        if(bflag):
            sigxblagr  = np.array(map(fxovery,sigxblagr  ,rbmass)) 
            sigyblagr  = np.array(map(fxovery,sigyblagr  ,rbmass)) 
            sigzblagr  = np.array(map(fxovery,sigzblagr  ,rbmass)) 
            sigrblagr  = np.array(map(fxovery,sigrblagr  ,rbmass)) 
            sigtblagr  = np.array(map(fxovery,sigtblagr  ,rbmass)) 
            sigrotblagr= np.array(map(fxovery,sigrotblagr,rbmass)) 
            sigblagr   = sigxblagr + sigyblagr + sigzblagr

    
