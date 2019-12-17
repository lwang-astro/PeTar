import sys
sys.path.append('/home/lwang/code/amuse-git/src')
import numpy as np
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

def generate_cluster(number_of_stars, radius_of_cluster):
    # numpy.random.seed(1)

    salpeter_masses = new_salpeter_mass_distribution(number_of_stars)
    total_mass = salpeter_masses.sum()
        
    convert_nbody = nbody_system.nbody_to_si(total_mass, radius_of_cluster | units.parsec)

    print("generate plummer model, N=",number_of_stars)
    particles = new_plummer_model(number_of_stars, convert_nbody)

    print("setting masses of the stars")
    particles.radius = 0.0 | units.RSun
    particles.mass = salpeter_masses

    #particles.mass[0] = 40.0 |units.MSun
    #particles.mass[1] = 50.0 |units.MSun
    
    return particles, convert_nbody


def evolve_cluster(particles, convert_nbody, end_time=40 | units.Myr, dt=0.25 | units.Myr, plot_HRdiagram=True):
    

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
    boxsize=2
    axe[0].set_xlim(-boxsize,boxsize)
    axe[0].set_ylim(-boxsize,boxsize)
    axe[0].set_aspect(1.0)

    framescale=20.0

    if plot_HRdiagram:
        axe[1].set_xlim(30e3,1e3);
        axe[1].set_ylim(1e-5,1e5);
        axe[1].set_yscale('log');
        axe[1].set_xscale('log');
    else :
        boxsize /= framescale
        axe[1].set_xlim(-boxsize,boxsize)
        axe[1].set_ylim(-boxsize,boxsize)
        axe[1].set_aspect(1.0)

    nlayer=20
    alphascale=np.linspace(1,nlayer,nlayer)*2.0/(nlayer*(nlayer+1))*2.0
    print(alphascale,alphascale.sum())
    sizescale=np.logspace(-1,2,nlayer)[::-1]
    #print(sizescale)
    for i in range(nlayer):
        pt =axe[0].scatter([],[],cmap='rainbow',alpha=alphascale[i],edgecolors='none')
        ptcls.append(pt)

    for i in range(nlayer):
        pt =axe[1].scatter([],[],cmap='rainbow',alpha=alphascale[i],edgecolors='none')
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

    def init():

        x_values = particles.position.x.value_in(units.parsec)
        y_values = particles.position.y.value_in(units.parsec)
        mass_values = particles.mass.value_in(units.MSun)
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
    R = 1.0
    t_end = 10
    dt = 0.25
    output_file = 'test'
    plot_HRdiagram=True

    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    if len(sys.argv) > 2:
        R = float(sys.argv[2])
    if len(sys.argv) > 3:
        t_end = float(sys.argv[3])
    if len(sys.argv) > 4:
        dt = float(sys.argv[4])
    if len(sys.argv) > 5:
        output_file = sys.argv[5]
    if len(sys.argv) > 6:
        plot_HRdiagram = False

    particles, convert_nbody = generate_cluster(N, R)

    anime = evolve_cluster(particles, convert_nbody, t_end| units.Myr, dt |units.Myr, plot_HRdiagram)

    write_set_to_file(particles, output_file+".hdf5", "amuse", append_to_file=False)    

    anime.save(output_file+'.mp4', fps=30, extra_args=['-vcodec', 'libx264'])    
