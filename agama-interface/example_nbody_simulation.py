#!/usr/bin/python
"""
This example demonstrates the use of Agama as a source of external gravitational potential
in several N-body simulation codes: GyrfalcON (from NEMO), Gadget4, Arepo, and PeTar.
The first of these codes uses the agama.so shared library directly as a plug-in, whereas
the others need some patches and recompilation; this script attempts to do this automatically.
In additional, a "restricted N-body simulation" method that uses the built-in tools from Agama
can be used as a faster (but more approximate) alternative: it evolves the collection of orbits
in the moving potential of the cluster on a prescribed trajectory in the host galaxy, updating
this cluster potential periodically to account for the mass loss due to tidal stripping.
This method is also illustrated by the script "example_tidal_stream.py".
The host galaxy potential mimics the Milky Way with a rotating bar, and is created separately
using the script "example_mw_potential_hunter24.py". The cluster is placed on an orbit within
the bar region, which leads to the stripping of ~1/3 of its mass by the end of the simulation
(which takes a few minutes).
"""
import sys, agama, numpy, os, subprocess, shutil, multiprocessing, zipfile, platform
if sys.version_info.major == 2:
    from urllib import urlretrieve
else:
    from urllib.request import urlretrieve

agama.setUnits(length=1, velocity=1, mass=1)  # length = 1 kpc, velocity = 1 km/s, mass = 1 Msun

# create an isolated star cluster
c_pot = agama.Potential(type='plummer', scaleRadius=0.01, mass=1000)
c_df  = agama.DistributionFunction(type='quasispherical', potential=c_pot)
c_gm  = agama.GalaxyModel(c_pot, c_df)
N_DEFAULT = 1000
xv, m = c_gm.sample(N_DEFAULT)

# shift it to some initial point in the Galaxy
c_center = numpy.hstack([2.0, 0, 0, 0, -100, 50])
xv   += c_center
sim_time = 0.25  # total simulation time in units of kpc/(km/s) = 0.98/5 Gyr

potfile = 'MWPotentialHunter24_rotating.ini'
try:
    g_pot = agama.Potential(potfile)
except RuntimeError:
    print('You need to create the barred Milky Way potential by running the script "example_mw_potential_hunter24.py"')
    exit(1)
print('''Initial conditions for a cluster in the galaxy are created, what to do next?
  1 - show the orbit of the cluster center as a test particle
  2 - run the simulation of a disrupting cluster using the Agama-only "restricted N-body" approach
  3 - run the simulation using the GyrfalcON code (from existing NEMO installation)
  4 - run the simulation using the Gadget4 code (downloading and compiling if necessary)
  5 - run the simulation using the Arepo code (downloading and compiling if necessary)
  6 - run the simulation using the PeTar code (Smart detection, N=1000 optimized IC, no auto-download)
  any other key - exit''')

def runSingleOrbit():
    import matplotlib.pyplot as plt
    o = agama.orbit(potential=g_pot, ic=c_center, time=sim_time, trajsize=500)[1]
    plt.figure(figsize=(10,10), dpi=75)
    plt.plot(o[:,0], o[:,1], 'k')
    plt.xlim(-2.5, 2.5)
    plt.ylim(-2.5, 2.5)
    plt.show()

def runRestrictedNbody():
    import matplotlib.pyplot as plt
    plt.ion()
    plt.figure(figsize=(10,10), dpi=75)
    ax = plt.axes([0.08, 0.08, 0.9, 0.9])
    num_intervals = 32
    num_subint = 16
    interval = sim_time / num_intervals
    time_center, orbit_center = agama.orbit(potential=g_pot, ic=c_center, time=sim_time,
        trajsize=num_intervals * num_subint + 1)
    cpot = c_pot
    snap = xv.copy()
    print('Playing animation as the simulation progresses...')
    for i in range(num_intervals+1):
        c_bound = c_pot.potential(snap[:,0:3] - orbit_center[i*num_subint, 0:3]) + 0.5 * numpy.sum((snap[:,3:6] - orbit_center[i*num_subint, 3:6])**2, axis=1) < 0
        time = i * interval
        ax.cla()
        ax.scatter(snap[:,0], snap[:,1], c=c_bound, cmap='bwr_r', vmin=0, vmax=1, s=2, linewidths=0)
        ax.plot(orbit_center[0:i*num_subint, 0], orbit_center[0:i*num_subint, 1], color='k')
        ax.text(0.01, 0.99, 'time=%.4f, bound fraction=%.3f' % (time, numpy.sum(c_bound)*1./N_DEFAULT), ha='left', va='top', transform=ax.transAxes)
        grid = numpy.linspace(-2.5, 2.5, 101)
        xyz = numpy.column_stack([numpy.repeat(grid,len(grid)), numpy.tile(grid,len(grid)), numpy.zeros(len(grid)**2)])
        den = g_pot.density(xyz, t=time).reshape(len(grid), len(grid)).T
        ax.contour(grid, grid, numpy.log10(den), levels=numpy.linspace(8.1, 12.1, 17), cmap='earth_r', zorder=2)
        ax.set_xlim(min(grid), max(grid))
        ax.set_ylim(min(grid), max(grid))
        plt.draw()
        plt.pause(0.01)
        if i == num_intervals: continue
        t_pot = agama.Potential(g_pot,
                agama.Potential(potential=cpot, center=numpy.column_stack((time_center, orbit_center))))
        snap = numpy.vstack(agama.orbit(ic=snap, potential=t_pot, time=interval, timestart=time, trajsize=1, accuracy=1e-5, verbose=False)[:,1])
        cpot = agama.Potential(type='multipole', particles=(snap[:,0:3] - orbit_center[(i+1)*num_subint, 0:3], m), symmetry='s')
    filename = 'example_nbody_simulation_last.nemo'
    agama.writeSnapshot(filename, (snap,m), 'nemo')
    print('Saved the final snapshot to %s' % filename)
    plt.ioff()
    plt.show()

def runGyrfalcon():
    import os, subprocess
    path = 'gyrfalcon'
    if not os.path.isdir(path):
        os.mkdir(path)
    os.chdir(path)
    infile = 'IC.nemo'
    print('*** Writing the initial conditions into %s' % infile)
    agama.writeSnapshot(infile, (xv, m), 'nemo')
    outfile = 'output.nemo'
    endfile = 'output_last.nemo'
    if os.path.isfile(outfile):
        os.remove(outfile)
    if os.path.isfile(endfile):
        os.remove(endfile)
    runstr = ('gyrfalcON %s %s eps=0.001 kmax=14 Nlev=3 fea=0.5 step=0.001953125 tstop=0.25 logstep=32 '
        'Grav=%g accname=agama accpars=0,%g accfile=%s') % (infile, outfile, agama.G, agama.G, '../'+potfile)
    print('*** Launching the simulation (end time: %g):\n%s' % (sim_time, runstr))
    sub = subprocess.check_call(runstr, shell=True)
    if os.path.isfile(outfile):
        print('*** Done! You may play the simulation movie using "glnemo %s/%s"' % (path,outfile))
        subprocess.check_call('s2s %s %s times=%g' % (outfile, endfile, sim_time), shell=True)
        showSnapshot(endfile)
    else:
        print('No output files produced, try running the above command manually to see the error messages')

def runArepoOrGadget(code):
    import os, platform, subprocess, zipfile, shutil, multiprocessing
    if sys.version_info.major == 2:
        from urllib import urlretrieve
    else:
        from urllib.request import urlretrieve
    path = code
    if os.path.isdir(path):
        print('*** Using the existing installation of %s' % code)
    else:
        print('*** Downloading %s' % code)
        filename = code + '.zip'
        if code == 'arepo':
            url = 'https://gitlab.mpcdf.mpg.de/vrs/arepo/-/archive/73d6bc4821daece021f02c0ae5da834c1d05b3c7/arepo-73d6bc4821daece021f02c0ae5da834c1d05b3c7.zip'
        elif code == 'gadget4':
            url = 'https://gitlab.mpcdf.mpg.de/vrs/gadget4/-/archive/1e171a4a679d30ac1e6accabe8a76a037ccbacac/gadget4-1e171a4a679d30ac1e6accabe8a76a037ccbacac.zip'
        urlretrieve(url, filename)
        if not os.path.isfile(filename):
            raise RuntimeError('Cannot find downloaded file %s' % filename)
        with zipfile.ZipFile(filename, 'r') as f:
            for info in f.infolist():
                extr = f.extract(info)
                attr = info.external_attr >> 16
                if attr:
                    os.chmod(extr, attr)
        os.remove(filename)
        folder = url[url.rfind('/')+1:-4]
        os.rename(folder, path)

    with open(path + '/Template-Config.sh', 'r') as f:
        lines = ''.join(f.readlines())
        if 'AGAMA' in lines:
            print('*** Code is already patched')
        else:
            print('*** Patching %s' % code)
            subprocess.check_call('patch -r -u -N -d %s -p 1 < example_nbody_simulation_%s.patch' % (path, code), shell=True)
            agama_root = os.path.abspath('../')
            if not (os.path.isfile(agama_root+'/Makefile.local') and os.path.isfile(agama_root+'/agama.a')):
                raise RuntimeError('Cannot find Makefile.local and agama.a in %s, exiting...' % agama_root)
            allflags = []
            with open(agama_root+'/Makefile.local', 'r') as mf:
                for line in mf.readlines():
                    if (line.startswith('LINK_FLAGS_ALL') or line.startswith('LINK_FLAGS_LIB_AND_EXE_STATIC') or
                            code=='arepo' and line.startswith('LINK_FLAGS_EXE_STATIC_NONCPP')):
                        flags = line[line.find('=')+1:].strip().split(' ')
                        for i in range(len(flags)):
                            if flags[i].startswith('extras/'):
                                flags[i] = agama_root + '/' + flags[i]
                            if flags[i].startswith('-Lextras/'):
                                flags[i] = '-L' + agama_root + '/' + flags[i][2:]
                        allflags += flags
            with open(path + '/Makefile', 'r') as mf:
                lines = mf.readlines()
                for i in range(len(lines)):
                    if lines[i].startswith('ifeq (EXTERNALGRAVITY_AGAMA'):
                        lines = lines[:i+1] + ['LIBS += %s/agama.a ' % agama_root + ' '.join(allflags) + '\n'] + lines[i+1:]
                        break
            with open(path + '/Makefile', 'w') as mf:
                mf.write(''.join(lines))

    os.chdir(path)
    executable = code[0].upper() + code[1:]
    if not os.path.isfile(executable):
        if not os.path.isfile('Makefile.systype'):
            with open('Makefile.systype', 'w') as f:
                if platform.uname()[0]=='Darwin':
                    systype = 'Darwin'
                else:
                    systype = 'Generic-gcc' if code == 'gadget4' else 'Ubuntu'
                f.write('SYSTYPE="%s"\n' % systype)
        macros = ['SELFGRAVITY', 'EXTERNALGRAVITY', 'EXTERNALGRAVITY_AGAMA']
        if code == 'gadget4':
            macros += ['GADGET2_HEADER']
        elif code == 'arepo':
            macros += ['GRAVITY_NOT_PERIODIC']
        with open('Config.sh', 'w') as f:
            f.write('\n'.join(macros)+'\n')
        print('*** I am going to compile the executable %s/%s, ' % (path, executable) +
            'please adjust the files Config.sh, Makefile, Makefile.systype in the %s directory as needed.\n' % path +
            'In particular, make sure that GSL and HDF5 libraries are available, if necessary, '
            'adding -I/path_to_gsl_or_hdf5_include to CFLAGS and -L/path_to_gsl_or_hdf5_lib to LIBS.\n'
            'Press any key when ready...')
        input()
        subprocess.check_call('make', shell=True)

    infile = 'IC.snap'
    print('*** Writing the initial conditions into %s' % infile)
    agama.writeSnapshot(infile, (xv, m), 'gadget')

    with open('agama_potential.ini', 'w') as f:
        f.write('[Potential]\nfile=../%s\n' % potfile)

    shutil.copyfile('../example_nbody_simulation_%s.param' % code, 'param.txt')

    numproc = max(1, min(8, multiprocessing.cpu_count()))

    runstr = 'mpirun -np %i ./%s param.txt' % (numproc, executable)
    print('*** Launching the simulation (end time: %g):\n%s' % (sim_time, runstr))
    subprocess.check_call(runstr + ' | grep Time:', shell=True)

    if os.path.isdir('output'):
        outputfiles = sorted(['output/'+name for name in os.listdir('output') if name.startswith('snapshot')])
    else:
        outputfiles = []
    if outputfiles:
        with open('output_files.txt', 'w') as f:
            f.write('\n'.join(outputfiles) + '\n')
        print('*** Done! You may play the simulation movie using "glnemo %s/output_files.txt"' % path)
        showSnapshot(outputfiles[-1])
    else:
        print('No output files produced, try running the above command manually to see the error messages')

def runPeTar():
    """ 
    Automated PeTar Runner with Smart Executable Detection. 
    Optimized for N=1000 IC with no stellar evolution.
    """
    import os, subprocess, shutil
    
    # Smart Executable Detection Logic

    def scan_dir_for_petar_agama(dir_path):
        matches = []
        try:
            for name in os.listdir(dir_path):
                if ('petar' in name) and ('agama' in name) and (not 'bse' in name) and (not 'format' in name) and (not 'test' in name) and (not '.g' in name) and (not 'debug' in name):
                    full_path = os.path.join(dir_path, name)
                    if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
                        matches.append(full_path)
        except OSError:
            pass
        return matches

    search_dirs = []
    for p in os.environ.get('PATH', '').split(os.pathsep):
        if p:
            search_dirs.append(p)
    search_dirs.append('../petar/bin')

    scanned_paths = []
    for d in search_dirs:
        scanned_paths += scan_dir_for_petar_agama(d)

    all_candidates = []
    for path in scanned_paths:
        if path not in all_candidates:
            all_candidates.append(path)

    def feature_score(exe_path):
        name = os.path.basename(exe_path).lower()
        score = 0
        score += 1000 * int('avx512' in name)
        score += 500  * int('avx2' in name) 
        score += 300  * int('omp' in name)
        score += 200  * int('mpi' in name)
        return score

    all_candidates.sort(key=feature_score, reverse=True)
    petar_exe = all_candidates[0] if all_candidates else None
            
    if not petar_exe:
        print("\n" + "!"*70)
        print("*** ERROR: No suitable PeTar executable with Agama support found.")
        print("Please compile PeTar manually following the README of PeTar.")
        print("!"*70 + "\n")
        return

    print("\n" + "="*70)
    print(f"*** Auto-detected optimal PeTar executable: {petar_exe}")
    print("="*70 + "\n")

    original_dir = os.getcwd()
    path = 'petar_run'
    if not os.path.isdir(path): 
        os.mkdir(path)
    os.chdir(path)
    
    petar_time_myr = sim_time * 977.7922
    output_interval = petar_time_myr / 10.0

    c_xv = xv - c_center
    numpy.savetxt('raw_ic.txt', numpy.hstack((m.reshape(-1, 1), c_xv)), fmt='%.8e')

    center_str = ",".join(map(str, c_center))
    init_cmd = f"petar.init -f petar.input -r 1000.0 -v kms2pcmyr -t -c {center_str} raw_ic.txt"
    
    run_cmd_prefix = "mpirun -np 1 " if "mpi" in petar_exe else ""
    runstr = (f"{run_cmd_prefix}{petar_exe} -u 1 -t {petar_time_myr} -o {output_interval} "
              f" -a 0 --agama-conf-file ../{potfile} petar.input")

    log_filename = 'petar_simulation.log'
    print(f'*** Initializing and Launching PeTar...')
    print(f'*** All PeTar console outputs are redirected to: petar_run/{log_filename}')
    print(f'*** Please wait. Simulation is running in the background...')
    
    with open(log_filename, 'w') as log_file:
        subprocess.check_call(init_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        subprocess.check_call(runstr, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
    
    print('*** PeTar simulation completed successfully!')

    output_files = sorted([f for f in os.listdir('.') if f.startswith('data.') and f[5:].isdigit()])
    if output_files:
        print(f"*** Output data saved in 'petar_run/'. Last snapshot: {output_files[-1]}")
        pos = numpy.loadtxt(output_files[-1], usecols=(1,2,3), skiprows=1)
        vel = numpy.loadtxt(output_files[-1], usecols=(4,5,6), skiprows=1)
        snap = numpy.hstack((pos/1000, vel/1.022712165045695))
        agama.writeSnapshot('final_snapshot.nemo', (snap, m), 'nemo')
        print('Saved the final snapshot to final_snapshot.nemo')
        showSnapshot('final_snapshot.nemo')
    else:
        print('No output files produced.')

    os.chdir(original_dir)

def showSnapshot(filename):
    import matplotlib.pyplot as plt
    pos = agama.readSnapshot(filename)[0]
    plt.figure(figsize=(10,10), dpi=75)
    plt.axes([0.08, 0.08, 0.9, 0.9])
    plt.scatter(pos[:,0], pos[:,1], c='b', s=2, linewidths=0)
    grid = numpy.linspace(-2.5, 2.5, 101)
    xyz = numpy.column_stack([numpy.repeat(grid,len(grid)), numpy.tile(grid,len(grid)), numpy.zeros(len(grid)**2)])
    den = g_pot.density(xyz, t=sim_time).reshape(len(grid), len(grid)).T
    plt.contour(grid, grid, numpy.log10(den), levels=numpy.linspace(8.1, 12.1, 17), cmap='earth_r', zorder=2)
    plt.xlim(min(grid), max(grid))
    plt.ylim(min(grid), max(grid))
    plt.show()



if sys.version_info.major == 2:
    input = raw_input
if len(sys.argv)>1:
    choice = sys.argv[1]
    print(choice+'\n')
else:
    choice = input()
if choice == '1': runSingleOrbit()
if choice == '2': runRestrictedNbody()
if choice == '3': runGyrfalcon()
if choice == '4': runArepoOrGadget('gadget4')
if choice == '5': runArepoOrGadget('arepo')
if choice == '6': runPeTar()