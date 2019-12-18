import sys
sys.path.append('/home/lwang/code/amuse-git/src')
import numpy as np
import getopt
import amuse
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation, rc
from matplotlib import cm
from IPython.display import HTML

from amuse.units import nbody_system
from amuse.units import units
from amuse.community.petar.interface import petar
from amuse.community.sse.interface import SSE
from amuse.test.amusetest import get_path_to_results

from amuse.io import write_set_to_file

from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

plt.style.use('dark_background')

def generate_cluster(number_of_stars, mass_of_cluster, radius_of_cluster):
    # numpy.random.seed(1)

    if (mass_of_cluster>0|units.MSun):
        number_of_stars = int(mass_of_cluster.value_in(units.MSun))*5
    
    salpeter_masses = new_salpeter_mass_distribution(number_of_stars)
    imf = salpeter_masses

    if (mass_of_cluster>0|units.MSun):
        print("sorted sampling salpeter IMF, total mass:",mass_of_cluster)
        imf_cumsum = imf.cumsum()
        pick_up_number = (imf_cumsum<=mass_of_cluster).sum()
        if (pick_up_number<number_of_stars):
            imf_pick_up = imf[:pick_up_number+1]
            sorted_imf_pick_up = np.sort(imf_pick_up)
            if (sorted_imf_pick_up[pick_up_number]-mass_of_cluster < mass_of_cluster-sorted_imf_pick_up[pick_up_number-1]):
                pick_up_number+=1
            imf = sorted_imf_pick_up[:pick_up_number]
            number_of_stars = pick_up_number
    total_mass = imf.sum()
        
    convert_nbody = nbody_system.nbody_to_si(total_mass, radius_of_cluster | units.parsec)

    print("generate plummer model, N=",number_of_stars)
    particles = new_plummer_model(number_of_stars, convert_nbody)

    print("setting masses of the stars")
    particles.radius = 0.0 | units.RSun
    particles.mass = imf

    #particles.mass[0] = 40.0 |units.MSun
    #particles.mass[1] = 50.0 |units.MSun
    
    return particles, convert_nbody


def evolve_cluster(particles, convert_nbody, end_time=40 | units.Myr, dt=0.25 | units.Myr, R=1.0, plot_HRdiagram=True, framescale=20.0):

    #gravity = petar(convert_nbody,redirection='none')
    gravity = petar(convert_nbody)
    
    stellar_evolution = SSE()

    print("initializing the particles")
    stellar_evolution.particles.add_particles(particles)
    from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(particles)
    from_stellar_evolution_to_model.copy_attributes(["mass"])

    print("centering the particles")
    particles.move_to_center()
    print("scaling particles to viridial equilibrium")
    particles.scale_to_standard(convert_nbody)

    gravity.particles.add_particles(particles)
    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
       
    total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy

    time = 0.0 | units.Myr
    particles.savepoint(time)
 
    print("evolving the model until t = " + str(end_time))


    ptcls=[]
    ncol= 2
    xsize=16
    fig, axe = plt.subplots(1,ncol,figsize=(xsize, 8))
    boxsize=2.0*R
    axe[0].set_xlim(-boxsize,boxsize)
    axe[0].set_ylim(-boxsize,boxsize)
    axe[0].set_aspect(1.0)
    axe[0].set_xlabel('x[pc]')
    axe[0].set_ylabel('y[pc]')

    if plot_HRdiagram:
        axe[1].set_xlim(30e3,1e3);
        axe[1].set_ylim(1e-5,1e5);
        axe[1].set_yscale('log');
        axe[1].set_xscale('log');
        axe[1].set_xlabel(r'$T_{eff}[K]$')
        axe[1].set_ylabel(r'$L[L_\odot]$')
    else :
        boxsize /= framescale
        axe[1].set_xlim(-boxsize,boxsize)
        axe[1].set_ylim(-boxsize,boxsize)
        axe[1].set_aspect(1.0)
        axe[1].set_xlabel('x[pc]')
        axe[1].set_ylabel('y[pc]')

    nlayer_cross=5
    nlayer_point=10
    nlayer = nlayer_cross + nlayer_point
    alpha_amplifier = 2.5
    alphascale=np.linspace(1,nlayer,nlayer)*2.0/(nlayer*(nlayer+1))*alpha_amplifier
    print('Alpha layer sequence:',alphascale,' sum:',alphascale.sum())
    sizescale=np.logspace(0,3,nlayer)[::-1]
    print('Size layer sequence:',sizescale)
    #print(sizescale)
    for i in range(nlayer_cross):
        pt =axe[0].scatter([],[],marker='+',alpha=alphascale[i],edgecolors='none')
        ptcls.append(pt)
    for i in range(nlayer_point):
        pt =axe[0].scatter([],[],alpha=alphascale[i],edgecolors='none')
        ptcls.append(pt)

    for i in range(nlayer_cross):
        pt =axe[1].scatter([],[],marker='+',alpha=alphascale[i],edgecolors='none')
        ptcls.append(pt)
    for i in range(nlayer_point):
        pt =axe[1].scatter([],[],alpha=alphascale[i],edgecolors='none')
        ptcls.append(pt)


    itext=2*nlayer
    ptcls.append(axe[0].text(.05, 0.95, '', transform = axe[0].transAxes, color='white'))
    ptcls.append(axe[0].text(.35, 0.95, '', transform = axe[0].transAxes, color='white'))
    ptcls.append(axe[0].text(.55, 0.95, '', transform = axe[0].transAxes, color='white'))

    #fig.patch.set_facecolor('black')
    #for i in range(2):
    #    axe[i].patch.set_facecolor('black')
    #    axe[i].spines['bottom'].set_color('white')
    #    axe[i].spines['top'].set_color('white') 
    #    axe[i].spines['right'].set_color('white')
    #    axe[i].spines['left'].set_color('white')
    #    axe[i].get

    luminosity = stellar_evolution.particles.luminosity.value_in(units.LSun)
    lum_min = np.min(luminosity)
    lum_max = np.max(luminosity)
    radius = stellar_evolution.particles.radius.value_in(units.RSun)
    temperature_eff = 5778*(luminosity/(radius*radius))**0.25
    temp_max = np.max(temperature_eff)
    temp_min = np.min(temperature_eff)

    def getCM(x,y,m,boxsize):
        sel=(x>-boxsize) & (x<boxsize) & (y>-boxsize) & (y<boxsize)
        mtot =m[sel].sum()
        nsel=sel.sum()
        xcm = (x[sel]).sum()/float(nsel)
        ycm = (y[sel]).sum()/float(nsel)
        return xcm, ycm

    def init():

        x_values = particles.position.x.value_in(units.parsec)
        y_values = particles.position.y.value_in(units.parsec)
        mass_values = particles.mass.value_in(units.MSun)
        #xcm,ycm = getCM(x_values,y_values,mass_values,R);
        #x_values = x_values - xcm
        #y_values = y_values - ycm
        axe[0].set_title('T=%f Myr' % 0)
        ptcls[itext  ].set_text(r"$M_{tot}: %f M_\odot$" % (particles.mass.sum().value_in(units.MSun)))
        ptcls[itext+1].set_text(r"$m_{max}: %f M_\odot$" % (particles.mass.max().value_in(units.MSun)))       
        axe[1].set_title('dE = %f ' % 0)

        luminosity = stellar_evolution.particles.luminosity.value_in(units.LSun)
        radius = stellar_evolution.particles.radius.value_in(units.RSun)
        temperature_eff = 5778*(luminosity/(radius*radius))**0.25
        LogTeff=(np.log10(temperature_eff)-np.log10(temp_min))/(np.log10(temp_max)-np.log10(temp_min))
        colors=cm.rainbow(1.0-LogTeff)

        for i in range(nlayer):
            #sizes = mass_values * (100*(1-float(i)/nlayer))
            sizes = mass_values*sizescale[i]
            #sizes = (np.log10(luminosity)-np.log10(lum_min)+1)
            #sizes = luminosity
            ptcls[i].set_offsets(np.array([x_values,y_values]).transpose())
            ptcls[i].set_sizes(sizes)
            ptcls[i].set_color(colors)

        if plot_HRdiagram:
            #types = stellar_evolution.particles.stellar_type.value_in(units.NO_UNIT)

            ptcls[2*nlayer-1].set_offsets(np.array([temperature_eff, luminosity]).transpose())
            #ptcls[2].set_color(types)
        else:
            for i in range(nlayer):
                sizes = mass_values*framescale*sizescale[i]
                ptcls[nlayer+i].set_offsets(np.array([x_values,y_values]).transpose())
                ptcls[nlayer+i].set_sizes(sizes)
                ptcls[nlayer+i].set_color(colors)

        return ptcls

    def animate(k):

        time = (dt*k)
        print("Evolve to time: ",time.as_quantity_in(units.Myr))
        gravity.evolve_model(time)
        stellar_evolution.evolve_model(time)

        from_gravity_to_model.copy()
        from_stellar_evolution_to_model.copy_attributes(["mass", "radius"])
        
        total_energy_at_this_time = gravity.kinetic_energy + gravity.potential_energy

        from_model_to_gravity.copy_attributes(["mass"])

        particles.savepoint(time)

        x_values = particles.position.x.value_in(units.parsec)
        y_values = particles.position.y.value_in(units.parsec)
        mass_values = particles.mass.value_in(units.MSun)
        #xcm,ycm = getCM(x_values,y_values,mass_values,R);
        #x_values = x_values - xcm
        #y_values = y_values - ycm
        axe[0].set_title('T = %f Myr' % time.value_in(units.Myr))
        de = (total_energy_at_this_time - total_energy_at_t0) /total_energy_at_t0
        axe[1].set_title('dE = %f ' % de)
        ptcls[itext  ].set_text(r"$M_{tot}: %f M_\odot$" % (particles.mass.sum().value_in(units.MSun)))
        ptcls[itext+1].set_text(r"$m_{max}: %f M_\odot$" % (particles.mass.max().value_in(units.MSun)))       

        luminosity = stellar_evolution.particles.luminosity.value_in(units.LSun)
        radius = stellar_evolution.particles.radius.value_in(units.RSun)
        temperature_eff = 5778*(luminosity/(radius*radius))**0.25
        LogTeff=(np.log10(temperature_eff)-np.log10(temp_min))/(np.log10(temp_max)-np.log10(temp_min))
        colors=cm.rainbow(1.0-LogTeff)
        #sizes = (np.log10(luminosity)-np.log10(lum_min)+1)*100
        #sizes = luminosity

        for i in range(nlayer):
            sizes = mass_values*sizescale[i]
            #sizes = (np.log10(luminosity)-np.log10(lum_min)+1)
            #sizes = luminosity
            ptcls[i].set_offsets(np.array([x_values,y_values]).transpose())
            ptcls[i].set_sizes(sizes)
            ptcls[i].set_color(colors)

        if plot_HRdiagram:
            #types = stellar_evolution.particles.stellar_type.value_in(units.NO_UNIT)

            ptcls[2*nlayer-1].set_offsets(np.array([temperature_eff, luminosity]).transpose())
            #ptcls[2].set_color(types)
        else:
            for i in range(nlayer):
                sizes = mass_values*framescale*sizescale[i]
                ptcls[nlayer+i].set_offsets(np.array([x_values,y_values]).transpose())
                ptcls[nlayer+i].set_sizes(sizes)
                ptcls[nlayer+i].set_color(colors)

        return ptcls

    n_frame = int(end_time/dt)

    anime = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=n_frame, interval=80, blit=True)
    
    #gravity.stop()
    #stellar_evolution.stop()

    return anime


if __name__ == '__main__':

    N = 50
    M = 0
    R = 1.0
    t_end = 10
    dt = 0.25
    output_file = 'test'
    plot_HRdiagram = False
    framescale = 20.0

    def usage():
        print("option:")
        print("  -h(--help): help")
        print("  -N: number of stars")
        print("  -M: mass of cluster in Msun (sorted sampling)")
        print("  -R: radius of cluster")
        print("  -t: finishing time in Myr")
        print("  -o(--output-filename): output filename")
        print("  -s(--dt-output): output time interval in Myr")
        print("  --frame-scale: second frame x,y enlarged ratio")

    try:
        shortargs = 'N:M:R:t:o:s:h'
        longargs = ['output-filename=','plot-HRdiagram', 'frame-scale=','help','frame-scale=']
        opts,args= getopt.getopt( sys.argv[1:], shortargs, longargs)

        for opt,arg in opts:
            if opt in ('-h','--help'):
                usage()
                sys.exit(1)
            elif opt in ('-N'):
                N = int(arg)
            elif opt in ('-M'):
                M = float(arg)
            elif opt in ('-R'):
                R = float(arg)
            elif opt in ('-t'):
                t_end = float(arg)
            elif opt in ('-o','--output-filename'):
                output_file = arg
            elif opt in ('--plot-HRdiagram'):
                plot_HRdiagram=True
            elif opt in ('-s','--dt-output'):
                dt = float(arg)
            elif opt in ('--frame-scale'):
                framescale = float(arg)
            else:
                assert False, "unhandeld option"

    except getopt.GetoptError:
        print('getopt error!')
        usage()
        sys.exit(1)

    particles, convert_nbody = generate_cluster(N, M |units.MSun, R)

    anime = evolve_cluster(particles, convert_nbody, t_end| units.Myr, dt |units.Myr, R, plot_HRdiagram, framescale)

    write_set_to_file(particles, output_file+".hdf5", "amuse", append_to_file=False)    

    anime.save(output_file+'.mp4', fps=30, extra_args=['-vcodec', 'libx264'])    
