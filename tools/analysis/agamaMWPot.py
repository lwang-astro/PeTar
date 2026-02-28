import numpy as np
from scipy.integrate import quad

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
    generate_agama_config(mode="forward", evolve_duration=3000, initial_time=8800, file_name="MW_forward.ini")

    # Backward: Same 3 Gyr time span, but the starting and ending angles are swapped
    generate_agama_config(mode="backward", evolve_duration=3000, initial_time=8800, file_name="MW_backward.ini")

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
    
    # Conversion factor: 1 InternalTime = 977.8131076864 Myr
    
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

    # Format the .ini content
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
            f"rotation={rot_spiral}"
        ]

    # Output file handling
    with open(file_name, "w") as f:
        f.write("\n".join(content))
    
    print(f"File saved: {file_name} (Mode: {mode})")


def createMWPotentialHunter24VariableRotationConfig(
    file_name="MW_decelerating_bar.ini",
    mode="forward",            
    evolve_duration=3000.0,     # Evolution duration (Myr) for bar/spiral rotation
    initial_time=0.0,           # Initial time (Myr) for PeTar simulation starting snapshot
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
    createMWPotentialHunter24VariableRotationConfig(time_myr=3000, current_bar_angle=-0.44)
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

    Below is a demonstration of how to generate this sequence using Python, 
    taking a simple monotonically decelerating Bar model as an example.
    ======================================================================
    """


    # --------------------------------------------------------------------
    # 1. Base Parameters and Unit Conversion
    # --------------------------------------------------------------------
    # Agama time conversion factor: from Myr to internal units (kpc / (km/s))
    
    t1 = initial_time / AGAMA_TIME_CONVERSION
    t2 = (initial_time + evolve_duration) / AGAMA_TIME_CONVERSION
    dt_agama = evolve_duration / AGAMA_TIME_CONVERSION

    # --------------------------------------------------------------------
    # 2. Define the Angular Velocity Function over Time
    # --------------------------------------------------------------------
    def fw(t, v_start, v_end):
        """
        Example of a single decelerating angular velocity (Unit: km/s/kpc):
        Assuming the bar's rotation speed linearly decreases from -50 to -35 
        during the evolution time.
        You can replace this with your own complex piecewise deceleration function.
        """
        return v_start + (v_end - v_start) * (t / dt_agama)

    # --------------------------------------------------------------------
    # 3. Create Time Sequence and Obtain the Integrated Angle
    # --------------------------------------------------------------------
    # Generate enough points for the sequence so Agama can perform an accurate 
    # regularized cubic spline interpolation.
    # For a 3000 Myr span, 1000 points provide very high precision.
    tlist = np.linspace(0, dt_agama, 1000)
    bar_integral_values = np.zeros_like(tlist)
    spiral_integral_values = np.zeros_like(tlist)

    # Calculate the cumulative angle rotated at each time point (integral of w over t)
    for i, t in enumerate(tlist):
        integral, _ = quad(fw, 0, t, args=(v_start_bar, v_end_bar))
        bar_integral_values[i] = integral
        integral, _ = quad(fw, 0, t, args=(v_start_spiral, v_end_spiral))
        spiral_integral_values[i] = integral

    tlist += t1  # Shift time points to start from initial_time in Agama's internal units

    # --------------------------------------------------------------------
    # 4. Align Angles and Combine into Agama Sequence
    # --------------------------------------------------------------------
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
            f"rotation={spiral_rotation_str}"
        ]

    # Save to file
    with open(file_name, "w") as f:
        f.write("\n".join(content))

    print(f"File saved: {file_name} (Mode: {mode})")
