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

def generate_cluster(number_of_stars):
    # numpy.random.seed(1)

    salpeter_masses = new_salpeter_mass_distribution(number_of_stars)
    total_mass = salpeter_masses.sum()
        
    convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)

    print "generate plummer model, N=",number_of_stars
    particles = new_plummer_model(number_of_stars, convert_nbody)

    print "setting masses of the stars"
    particles.radius = 0.0 | units.RSun
    particles.mass = salpeter_masses

    #particles.mass[0] = 40.0 |units.MSun
    #particles.mass[1] = 50.0 |units.MSun
    
    return particles, convert_nbody


def evolve_cluster(particles, convert_nbody, end_time=40 | units.Myr, dt=0.25 | units.Myr):
    

    #gravity = petar(convert_nbody,redirection='none')
    gravity = petar(convert_nbody)
    
    stellar_evolution = SSE()

    print "initializing the particles"
    stellar_evolution.particles.add_particles(particles)
    from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(particles)
    from_stellar_evolution_to_model.copy_attributes(["mass"])
    #print stellar_evolution.particles

    print "centering the particles"
    particles.move_to_center()
    print "scaling particles to viridial equilibrium"
    particles.scale_to_standard(convert_nbody)

    gravity.particles.add_particles(particles)
    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
       
    total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy

    time = 0.0 | units.Myr
    #particles.savepoint(time)
 
    print "evolving the model until t = " + str(end_time)


    ptcls=[]
    fig, axe = plt.subplots(1,2,figsize=(16, 8))
    pt =axe[0].scatter([],[],marker='o')
    axe[0].set_xlim(-2,2)
    axe[0].set_ylim(-2,2)
    ptcls.append(pt)
    pt =axe[1].scatter([],[],marker='o')
    axe[1].set_xlim(30e3,1e3);
    axe[1].set_ylim(1e-5,1e5);
    axe[1].set_yscale('log');
    axe[1].set_xscale('log');
    ptcls.append(pt)    
    ptcls.append(axe[0].text(.05, 0.95, '', transform = axe[0].transAxes))
    ptcls.append(axe[0].text(.15, 0.95, '', transform = axe[0].transAxes))
    ptcls.append(axe[0].text(.25, 0.95, '', transform = axe[0].transAxes))
    
    def init():
        x_values = particles.position.x.value_in(units.parsec)
        y_values = particles.position.y.value_in(units.parsec)
        mass_values = particles.mass.value_in(units.MSun)
        sizes = mass_values * 100.0
        ptcls[0].set_offsets(np.array([x_values,y_values]).transpose())
        ptcls[0].set_sizes(sizes)
        axe[0].set_title('T=%f Myr' % 0)
        axe[1].set_title('dE = %f ' % 0)
        #print  gravity.total_mass.as_quantity_in(units.MSun)
        ptcls[2].set_text("M_tot: %f Msun" % (particles.mass.sum().value_in(units.MSun)))
        ptcls[2].set_text("m_max: %f Msun" % (particles.mass.max().value_in(units.MSun)))       

        luminosity = stellar_evolution.particles.luminosity.value_in(units.LSun)
        radius = stellar_evolution.particles.radius.value_in(units.RSun)
        temperature_eff = 5778*(luminosity/(radius*radius))**0.25
        #types = stellar_evolution.particles.stellar_type.value_in(units.NO_UNIT)

        ptcls[1].set_offsets(np.array([temperature_eff, luminosity]).transpose())
        #ptcls[1].set_color(types)

        return ptcls

    def animate(k):
        time = (dt*k)
        print "Evolve to time: ",time.as_quantity_in(units.Myr)
        #print "gravity evolve step starting"
        gravity.evolve_model(time)
        #print "gravity evolve step done"

        #print "stellar evolution step starting"
        stellar_evolution.evolve_model(time)
        #print "stellar evolution step done"

        from_gravity_to_model.copy()
        from_stellar_evolution_to_model.copy_attributes(["mass", "radius"])
        
        total_energy_at_this_time = gravity.kinetic_energy + gravity.potential_energy

        from_model_to_gravity.copy_attributes(["mass"])

        #particles.savepoint(time)

        x_values = particles.position.x.value_in(units.parsec)
        y_values = particles.position.y.value_in(units.parsec)
        mass_values = particles.mass.value_in(units.MSun)
        sizes = mass_values * 100.0
        ptcls[0].set_offsets(np.array([x_values,y_values]).transpose())
        ptcls[0].set_sizes(sizes)
        axe[0].set_title('T = %f Myr' % time.value_in(units.Myr))
        de = (total_energy_at_this_time - total_energy_at_t0) /total_energy_at_t0
        axe[1].set_title('dE = %f ' % de)
        ptcls[2].set_text("M_tot: %f Msun" % (particles.mass.sum().value_in(units.MSun)))
        ptcls[2].set_text("m_max: %f Msun" % (particles.mass.max().value_in(units.MSun)))       

        luminosity = stellar_evolution.particles.luminosity.value_in(units.LSun)
        radius = stellar_evolution.particles.radius.value_in(units.RSun)
        temperature_eff = 5778*(luminosity/(radius*radius))**0.25
        #types = stellar_evolution.particles.stellar_type.value_in(units.NO_UNIT)

        ptcls[1].set_offsets(np.array([temperature_eff, luminosity]).transpose())
        #ptcls[1].set_color(types)

        return ptcls

    n_frame = int(end_time/dt)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=n_frame, interval=80, blit=True)
    
    #gravity.stop()
    #stellar_evolution.stop()

    return anim    


if __name__ == '__main__':

    N = 50
    t_end = 10
    dt = 0.25
    output_file = 'test'

    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    if len(sys.argv) > 2:
        t_end = float(sys.argv[2])
    if len(sys.argv) > 3:
        dt = float(sys.argv[3])
    if len(sys.argv) > 4:
        output_file = sys.argv[4]

    particles_init, convert_nbody = generate_cluster(N)

    ani = evolve_cluster(particles_init, convert_nbody, t_end| units.Myr, dt |units.Myr)

    #ani=plot_particles(particles)

    ani.save(output_file+'.mp4', fps=30, extra_args=['-vcodec', 'libx264'])    
