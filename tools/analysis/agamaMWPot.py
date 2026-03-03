import os
import numpy as np
import scipy.integrate
from scipy.integrate import quad
import scipy.special
import agama

AGAMA_TIME_CONVERSION = 977.8131076864  # Myr to kpc/(km/s)

def createMWPotentialHunter24ConstRotationConfig(
    file_name="MW_potential.ini",
    mode="forward",
    evolve_duration=1000.0,     # Evolution duration (Myr)
    initial_time=0.0,           # Initial time (Myr) for PeTar simulation starting snapshot
    current_angle_bar=-0.44,    # Current angle of the bar at the end of forward evolution (rad)
    current_angle_spiral=-0.44, # Current angle of the spiral at the end of forward evolution (rad)
    v_bar=-37.5,                # Unit: km/s/kpc
    v_spiral=-22.5,             # Unit: km/s/kpc
    main_pot_file="MWPotentialHunter24_full.ini",
    spiral_pot_file="MWPotentialHunter24_spiral.ini"
):
    """
    Generate Agama potential configuration file for Milky Way potential of (Hunter et al. 2024) with rotating bar and spiral arm using constant speed 
    
    Parameters:
    -----------
    file_name: str (default: "MW_potential.ini")
        Name of the output configure file for Agama.
    mode: str (default: "forward")
        "forward" for simulating forward evolution, "backward" for simulating backward evolution.
    evolve_duration: float (default: 1000.0)
        Total duration of the evolution in Myr for bar/spiral rotation.
    initial_time: float (default: 0.0)
        Starting time of the evolution in Myr for PeTar input snapshot.
    current_angle_bar: float (default: -0.44)
        The angle of the bar at the end of the evolution (present day) in radians.
    current_angle_spiral: float (default: -0.44)
        The angle of the spiral arms at the end of the evolution (present day) in radians.
    v_bar: float (default: -37.5)
        The pattern speed of the bar in km/s/kpc (negative for clockwise rotation).
    v_spiral: float (default: -22.5)
        The pattern speed of the spiral arms in km/s/kpc (negative for clockwise rotation).
    main_pot_file: str (default: "MWPotentialHunter24_full.ini")
        The base potential file for the main component (e.g., bar + axisymmetric).
    spiral_pot_file: str (default: "MWPotentialHunter24_spiral.ini")
        The base potential file for the spiral arms. If None, switch off spiral arms.

    Example usage
    -------------
    # Forward: Simulate evolution from 8.8 Gyr (PeTar input snapshot time) to 11.8 Gyr (present day)
    createMWPotentialHunter24ConstRotationConfig(mode="forward", evolve_duration=3000, initial_time=8800, file_name="MW_forward.ini")

    # Backward: Same 3 Gyr time span, but the starting and ending angles are swapped
    createMWPotentialHunter24ConstRotationConfig(mode="backward", evolve_duration=3000, initial_time=8800, file_name="MW_backward.ini")

    --------------------------------------------------------------------------
    [Detailed Explanation of Units and Physical Logic]
    
    1. Time Conversion:
       Agama uses a unit system where G = 1 by default. If the unit of length 
       is kpc and the unit of velocity is km/s, the internal time unit (UnitTime) 
       is defined as 1 kpc / (1 km/s).
       Conversion factor = (3.08567758e16 / 1.0) / 3.15576e13 ~ 977.813106
       Therefore, AGAMA_TIME_CONVERSION = 977.8131076864 converts Myr to kpc/(km/s).
       
    2. Pattern Speed:
       The units for v_bar and v_spiral are km/s/kpc. Numerically, this is 
       equivalent to the angular velocity Omega (rad/UnitTime).
       Because: (km/s/kpc) * (kpc/(km/s)) = dimensionless (radians).
       
    3. Angle Calculation:
       Change in angle Delta_Phi = Omega * Delta_T_internal.
       Since current_angle is set as the end of the evolution (t_now), 
       the starting angle is derived as:
       past_angle = current_angle - (v * dt_agama).
       
    4. Forward vs. Backward Logic:
       - Forward: Time axis increases, angle evolves from past to current.
       - Backward: Time axis increases, but the starting angle is set to 
         current and the ending angle to past. This simulates a "reverse 
         rotation" effect within the same absolute time frame.
    --------------------------------------------------------------------------
    """
    
    # Conversion to Agama time units
    t1 = initial_time / AGAMA_TIME_CONVERSION
    t2 = (initial_time + evolve_duration) / AGAMA_TIME_CONVERSION
    dt_agama = evolve_duration / AGAMA_TIME_CONVERSION
    
    # Calculate the change in angle Delta_Theta = v (km/s/kpc) * dt (kpc/(km/s))
    delta_angle_bar = v_bar * dt_agama
    delta_angle_spiral = v_spiral * dt_agama
    
    # Extrapolate the angle at the starting time (past angle = current angle - change)
    past_angle_bar = current_angle_bar - delta_angle_bar
    past_angle_spiral = current_angle_spiral - delta_angle_spiral
    
    # Set rotation parameters based on forward/backward mode
    if mode.lower() == "forward":
        # Forward: From [past time, past angle] to [current time, current angle]
        rot_main = f"[[{t1:.8f}, {past_angle_bar:.8f}], [{t2:.8f}, {current_angle_bar:.8f}]]"
        rot_spiral = f"[[{t1:.8f}, {past_angle_spiral:.8f}], [{t2:.8f}, {current_angle_spiral:.8f}]]"
        header = f"# Forward Evolution: From {initial_time} Myr to {initial_time + evolve_duration} Myr"
    else:
        # Backward: Time points remain the same, but initial and final angles are swapped
        rot_main = f"[[{t1:.8f}, {current_angle_bar:.8f}], [{t2:.8f}, {past_angle_bar:.8f}]]"
        rot_spiral = f"[[{t1:.8f}, {current_angle_spiral:.8f}], [{t2:.8f}, {past_angle_spiral:.8f}]]"
        header = f"# Backward Simulation: Angles swapped to reverse the rotation effect"

    # Format the configure file content
    content = [
        header,
        "# Bar + Axisymmetric Potential (Main Component)",
        "[Potential main]",
        f"file={main_pot_file}",
        f"rotation={rot_main}\n"]
    if spiral_pot_file is not None:
        content += [
            "# Spiral Arms Potential",
            "[Potential spiral]",
            f"file={spiral_pot_file}",
            f"rotation={rot_spiral}\n"
        ]

    # Output file
    with open(file_name, "w") as f:
        f.write("\n".join(content))
    
    print(f"File saved: {file_name} (Mode: {mode})")


def createMWPotentialHunter24VariableRotationConfig(
    file_name="MW_decelerating_bar.ini",
    mode="forward",            
    evolve_duration=3000.0,     # Evolution duration (Myr) for bar/spiral rotation
    initial_time=0.0,           # Initial time (Myr) for PeTar simulation starting snapshot
    nstep = 1000,               # Number of time steps for sampling the rotation sequence (affects interpolation accuracy and file size)
    current_bar_angle=-0.44,    # Current angle of the bar at the end of forward evolution (rad)
    current_angle_spiral=-0.44, # Current angle of the spiral at the end of forward evolution (rad)
    v_start_bar=-50.0,          # Initial pattern speed of the bar in km/s/kpc
    v_end_bar=-35.0,            # Final pattern speed of the bar in km/s/kpc at the end of evolution
    v_start_spiral=-30.0,       # Initial pattern speed of the spiral arms in km/s/kpc
    v_end_spiral=-20.0,         # Final pattern speed of the spiral arms in km/s/kpc at the end of evolution
    main_pot_file="MWPotentialHunter24_full.ini", 
    spiral_pot_file="MWPotentialHunter24_spiral.ini"
):
    """
    Generate Agama potential configuration file for Milky Way potential of (Hunter et al. 2024) with a decelerating bar using a linearly evolving rotation sequence.
    
    Parameters:
    -----------
    file_name: str (default: "MW_decelerating_bar.ini")
        Name of the output configure file for Agama.
    mode: str (default: "forward")
        "forward" for simulating forward evolution, "backward" for simulating backward evolution.
    evolve_duration: float (default: 3000.0)
        Total duration of the evolution in Myr for bar/spiral rotation.
    initial_time: float (default: 0.0)
        The starting time of the evolution in Myr for PeTar input snapshot.
    nstep: int (default: 1000)
        Number of time steps for sampling the rotation sequence. Affects interpolation accuracy and file size.
    current_bar_angle: float (default: -0.44)
        The angle of the bar at the end of the evolution (present day) in radians.
    current_angle_spiral: float (default: -0.44)
        The angle of the spiral at the end of the evolution (present day) in radians.
    v_start_bar: float (default: -50.0)
        The initial pattern speed of the bar in km/s/kpc at the start of evolution.
    v_end_bar: float (default: -35.0)
        The final pattern speed of the bar in km/s/kpc at the end of evolution.
    v_start_spiral: float (default: -30.0)
        The initial pattern speed of the spiral arms in km/s/kpc at the start of evolution.
    v_end_spiral: float (default: -20.0)
        The final pattern speed of the spiral arms in km/s/kpc at the end of evolution
    main_pot_file: str (default: "MWPotentialHunter24_full.ini")
        The base potential file for the main component (e.g., bar + axisymmetric).
    spiral_pot_file: str (default: "MWPotentialHunter24_spiral.ini")
        The base potential file for the spiral arms. If None, switch off spiral arms.
    
    Example usage
    -------------
    createMWPotentialHunter24VariableRotationConfig(evolve_duration=3000, current_bar_angle=-0.44)
    --------------------------------------------------------------------------

    [Detailed Explanation of Variable Rotation Configuration]

    For potentials with a rotation speed that changes over time (variable rotation), 
    Agama cannot define the rotation using just the initial and final time points 
    like it does for uniform rotation. In this case, we must provide a "serialized table" 
    directly in the configuration file.

    For example, a typical .ini file content looks like this:
    [Potential main]
    file=MWPotentialHunter24_full_NFW_84.ini
    rotation=[[t1, angle1], [t2, angle2], ...]

    Detailed parameter description:
    1. file: The base potential file to be rotated.
    2. rotation: Represents the rotation of this potential around the z-axis. 
    In the table that follows, the first value is the time point (in Agama's 
    internal time units), and the second value is the corresponding angle (in radians).
    3. Interpolation mechanism: For times not explicitly given in the nodes, Agama uses 
    a regularized cubic spline for automatic interpolation. This ensures not only 
    smooth angle transitions but also continuous angular velocity.
    ======================================================================
    """

    t1 = initial_time / AGAMA_TIME_CONVERSION
    dt_agama = evolve_duration / AGAMA_TIME_CONVERSION

    # The Angular Velocity Function over Time, assuming the bar's rotation speed linearly decreases
    def fw(t, v_start, v_end):
        return v_start + (v_end - v_start) * (t / dt_agama)

    # Create Time Sequence and Obtain the Integrated Angle
    tlist = np.linspace(0, dt_agama, nstep)
    bar_integral_values = np.zeros_like(tlist)
    spiral_integral_values = np.zeros_like(tlist)

    # Calculate the cumulative angle rotated at each time point (integral of w over t)
    for i, t in enumerate(tlist):
        integral, _ = quad(fw, 0, t, args=(v_start_bar, v_end_bar))
        bar_integral_values[i] = integral
        integral, _ = quad(fw, 0, t, args=(v_start_spiral, v_end_spiral))
        spiral_integral_values[i] = integral

    tlist += t1  # Shift time points to start from initial_time in Agama's internal units

    bar_theta = np.array(bar_integral_values)
    spiral_theta = np.array(spiral_integral_values)
    # Adjust the real angle sequence: Current Integral - Final Integral + Actual Present Angle
    # This ensures that at the final time point (t_max_agama), the angle is exactly -0.44
    bar_theta_real = bar_theta - bar_theta[-1] + current_bar_angle
    spiral_theta_real = spiral_theta - spiral_theta[-1] + current_angle_spiral

    if mode.lower() == "forward":
        header = f"# Forward Evolution: From 0 Myr to {evolve_duration} Myr with decelerating rotation"
        # Format into the table structure required by Agama: [[t1, angle1], [t2, angle2], ...]
        bar_rotation_points = [f"[{t:.8f},{ang:.8f}]" for t, ang in zip(tlist, bar_theta_real)]
        bar_rotation_str = "[" + ",".join(bar_rotation_points) + "]"
        spiral_rotation_points = [f"[{t:.8f},{ang:.8f}]" for t, ang in zip(tlist, spiral_theta_real)]
        spiral_rotation_str = "[" + ",".join(spiral_rotation_points) + "]"
    else:
        header = f"# Backward Simulation: From {evolve_duration} Myr to 0 Myr with decelerating rotation (angles reversed)"
        # For backward, we reverse the order of the angles to simulate backward evolution
        bar_rotation_points = [f"[{t:.8f},{ang:.8f}]" for t, ang in zip(tlist, bar_theta_real[::-1])]
        bar_rotation_str = "[" + ",".join(bar_rotation_points) + "]"
        spiral_rotation_points = [f"[{t:.8f},{ang:.8f}]" for t, ang in zip(tlist, spiral_theta_real[::-1])]
        spiral_rotation_str = "[" + ",".join(spiral_rotation_points) + "]"

    # --------------------------------------------------------------------
    # 5. Output to Agama .ini Configuration File
    # --------------------------------------------------------------------
    content = [
        header,
        "# Agama potential configuration with a decelerating bar (main component)",
        "[Potential main]",
        f"file={main_pot_file}",
        f"rotation={bar_rotation_str}\n",
    ]
    if spiral_pot_file is not None:
        content += [
            "# Spiral Arms Potential",
            "[Potential spiral]",
            f"file={spiral_pot_file}",
            f"rotation={spiral_rotation_str}\n"
        ]

    # Save to file
    with open(file_name, "w") as f:
        f.write("\n".join(content))

    print(f"File saved: {file_name} (Mode: {mode})")


def appendLMCPotentialConfig(
    target_ini_file="MW_potential.ini",
    mode="forward",
    evolve_duration=3000.0,
    dt_myr=1.0,
    initial_time=0.0
):
    """
    Compute the past trajectory of the Large Magellanic Cloud (LMC) and append 
    its time-dependent potential configuration (along with the Milky Way reflex 
    acceleration) to an existing Agama .ini configuration file.
    
    Parameters:
    -----------
    target_ini_file: str (default: "MW_potential.ini")
        The path of the target .ini file. The LMC configurations will be appended 
        to the end of this file without overwriting existing content.
    mode: str (default: "forward")
        "forward" for simulating forward evolution (LMC trajectory moves from past to present),
        "backward" for simulating backward evolution (LMC trajectory moves from present to past).
    evolve_duration: float (default: 3000.0)
        Total duration of the evolution in Myr.
    dt_myr: float (default: 1.0)
        Time step in Myr for sampling the trajectory. Smaller values yield larger .txt files
        but higher precision for Agama's spline interpolation.
    initial_time: float (default: 0.0)
        Starting time of the evolution in Myr in the simulation timeframe.

    Example usage
    -------------
    # Step 1: Generate the main Milky Way potential
    createMWPotentialHunter24VariableRotationConfig(file_name="MW_final.ini", mode="forward")
    
    # Step 2: Append the LMC potential and reflex acceleration to the same file
    appendLMCPotentialConfig(target_ini_file="MW_final.ini", mode="forward")

    --------------------------------------------------------------------------
    [Detailed Explanation of LMC Integration and Time Mapping Logic]
    
    1. Physical Integration (Backward in time):
       We only have observational constraints for the LMC's position and velocity 
       at the *present day* (t=0). Therefore, to compute its history, we must 
       integrate the orbit backwards in time from 0 to -evolve_duration using 
       the Milky Way potential and dynamical friction (Chandrasekhar formula).
       
    2. Time Sequence Generation:
       The ODE solver computes the states corresponding to t = [0, -dt, ..., -duration].
       We then generate a simulation time array mapped from `initial_time` to 
       `initial_time + evolve_duration`.

    3. Forward vs. Backward Mapping:
       - Forward: The simulation moves from past to present. However, our ODE 
         solution starts at the present and goes to the past. Therefore, we 
         reverse (`[::-1]`) the ODE solution arrays so that the trajectory 
         evolves naturally from past to present alongside the simulation time.
       - Backward: The simulation moves from present to past (simulating reverse 
         time flow). This exactly matches the direction of our ODE solution, 
         so we keep the array as is.
    --------------------------------------------------------------------------
    """
    print(f"Generating LMC trajectory for {evolve_duration} Myr (Mode: {mode})...")
    
    # 1. Set Agama units
    agama.setUnits(length=1, velocity=1, mass=1)
    
    # 2. Initialize LMC current observation coordinates and velocity (t=0)
    ra, dec, dist, pmra, pmdec, vlos = 81.28, -69.78, 49.6, 1.858, 0.385, 262.2
    l, b, pml, pmb = agama.transformCelestialCoords(
        agama.fromICRStoGalactic, ra * np.pi/180, dec * np.pi/180, pmra, pmdec)
    posvelLMC = agama.getGalactocentricFromGalactic(l, b, dist, pml*4.74, pmb*4.74, vlos)

    # 3. Create a background MW model strictly for the ODE orbit integration
    # (Computes gravitational pull and dynamical friction acting on the LMC)
    paramBulge = dict(type='Spheroid', mass=1.2e10, scaleRadius=0.2, outerCutoffRadius=1.8, gamma=0.0, beta=1.8)
    paramDisk  = dict(type='Disk', mass=5.0e10, scaleRadius=3.0, scaleHeight=-0.4)
    paramHalo  = dict(type='Spheroid', densityNorm=1.35e7, scaleRadius=14, outerCutoffRadius=300, cutoffStrength=4, gamma=1, beta=3)
    densMWhalo = agama.Density(paramHalo)
    potMW      = agama.Potential(paramBulge, paramDisk, paramHalo)
    potMWsph   = agama.Potential(type='Multipole', potential=potMW, lmax=0, rmin=0.01, rmax=1000)
    gmHalo     = agama.GalaxyModel(potMWsph, agama.DistributionFunction(type='quasispherical', density=densMWhalo, potential=potMWsph))
    
    rgrid    = np.logspace(1, 3, 16)
    xyzgrid  = np.column_stack([rgrid, rgrid*0, rgrid*0])
    sigmafnc = agama.Spline(rgrid, gmHalo.moments(xyzgrid, dens=False, vel=False, vel2=True)[:,0]**0.5)

    # 4. Create LMC self-potential model
    massLMC    = 1.5e11
    radiusLMC  = (massLMC/1e11)**0.6 * 8.5
    bminCouLog = radiusLMC * 2.0
    potLMC     = agama.Potential(type='spheroid', mass=massLMC, scaleradius=radiusLMC, outercutoffradius=radiusLMC*10, gamma=1, beta=3)

    # 5. Define ODE system
    def difeq(vars, t):
        x0, v0 = vars[0:3], vars[3:6]
        x1, v1 = vars[6:9], vars[9:12]
        dx, dv = x1-x0, v1-v0
        dist = sum(dx**2)**0.5
        vmag = sum(dv**2)**0.5
        f0 = potLMC.force(-dx)
        f1 = potMW.force(dx)
        rho = potMW.density(dx)
        sigma = sigmafnc(dist)
        couLog= max(0, np.log(dist / bminCouLog)**0.5)
        X = vmag / (sigma * 2**.5)
        drag = -(4*np.pi * rho * dv / vmag *
            (scipy.special.erf(X) - 2/np.pi**.5 * X * np.exp(-X*X)) *
            massLMC * agama.G**2 / vmag**2 * couLog)
        return np.hstack((v0, f0, v1, f1 + drag))

    # 6. Compute physical trajectory (integrating backwards from t=0 to t=-evolve_duration)
    t_phys_start = -evolve_duration / AGAMA_TIME_CONVERSION
    steps = int(round(evolve_duration / dt_myr))
    t_ode = np.linspace(0, t_phys_start, steps + 1)
    
    ic = np.hstack((np.zeros(6), posvelLMC))
    sol = scipy.integrate.odeint(difeq, ic, t_ode)

    # 7. Map to simulation time based on forward/backward mode
    t1 = initial_time / AGAMA_TIME_CONVERSION
    t2 = (initial_time + evolve_duration) / AGAMA_TIME_CONVERSION
    t_sim = np.linspace(t1, t2, steps + 1)
    
    if mode.lower() == "forward":
        # Reverse ODE solution so it evolves from past -> present alongside simulation time
        mapped_sol = sol[::-1]
    else:
        # Keep ODE solution as is (present -> past) for backward simulation
        mapped_sol = sol

    # 8. Extract relative trajectory and inertial acceleration
    trajLMC = np.column_stack([t_sim, mapped_sol[:, 6:12] - mapped_sol[:, 0:6]])
    
    trajMWx = agama.Spline(t_sim, mapped_sol[:, 0], der=mapped_sol[:, 3])
    trajMWy = agama.Spline(t_sim, mapped_sol[:, 1], der=mapped_sol[:, 4])
    trajMWz = agama.Spline(t_sim, mapped_sol[:, 2], der=mapped_sol[:, 5])
    accMW   = np.column_stack([t_sim, -trajMWx(t_sim, 2), -trajMWy(t_sim, 2), -trajMWz(t_sim, 2)])

    # 9. Determine output paths and export data
    target_dir = os.path.dirname(os.path.abspath(target_ini_file))
    if not target_dir: 
        target_dir = os.getcwd()
        target_ini_file = os.path.join(target_dir, target_ini_file)

    file_static = os.path.join(target_dir, 'LMC_static.ini')
    file_traj   = os.path.join(target_dir, 'trajLMC.txt')
    file_acc    = os.path.join(target_dir, 'accMW.txt')

    potLMC.export(file_static)
    np.savetxt(file_traj, trajLMC, fmt='%.8e')
    np.savetxt(file_acc, accMW, fmt='%.8e')

    # 10. Append configurations to the target .ini file
    with open(target_ini_file, "a") as f:
        f.write("\n\n# =====================================================================\n")
        f.write(f"# LMC and Reflex Acceleration Perturbations (Mode: {mode})\n")
        f.write(f"# Evolution: {evolve_duration} Myr, dt: {dt_myr} Myr\n")
        f.write("# =====================================================================\n")
        
        f.write("[Potential lmc_moving]\n")
        f.write(f"file = {file_static}\n")
        f.write(f"center = {file_traj}\n\n")
        
        f.write("[Potential reflex_acceleration]\n")
        f.write("type = UniformAcceleration\n")
        f.write(f"file = {file_acc}\n")

    print(f"LMC configuration successfully appended to: {target_ini_file}")