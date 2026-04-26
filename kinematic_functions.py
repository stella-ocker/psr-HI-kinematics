import numpy as np

distances = np.linspace(0, 20, 100000) # initialize distances

# Reid et al. parameters
# Constants from Table 3 column A5
R0 = 8.15
U = 10.6
V = 10.7
W = 7.6
U_s = 6.0
V_s = -4.3
A2 = 0.96
A3 = 1.62
THETA_0 = 236

def calculate_galactocentric_radius(longitude):
    """ Calculates the galactocentric radius at a given longitude and uses discretized
        distances.

    Args:
        longitude (float): The longitude at a specific line of sight.

    Returns:
        Numpy Array: An array that represents a range of galactocentric radii.
    """
    return np.sqrt(R_0**2 + distances**2 - 2* R_0 *distances * np.cos(longitude))

def calculate_best_fit(radius):
    """ Calculates the best fit for the rotation curve which returns us a range of 
        circular velocities at the line of sight that was calculated for the radius.
        This function uses the power law best-fit.

    Args:
        radius (Numpy Array): An array of galactocentric radii along a line of sight.

    Returns:
        Numpy Array: An array of radial velocities from the best fit line
    """
    return 109.190 + 108.201 * (radius**0.0042069)
    # return (221.461 - 0.44286 * radius)

def angular_velocities_translation(theta, radius):
    """ This translates our calculated radial velocities and radii to angular velocity.

    Args:
        theta (Numpy Array): Radial velocities along a line of sight.
        radius (Numpy Array): Galactocentric radii along a line of sight.

    Returns:
        Numpy Array: The angular velocities along the line of sight of the pulsar.
    """
    return theta / radius

def calc_radial_velocities(omega, longitude, latitude):
    """ Solves for the radial velocities along a line of sight.

    Args:
        omega (Numpy Array): The angular velocities along a line of sight.
        longitude (float): The longitude of the pulsar.
        latitude (float): The latitude of the pulsar.

    Returns:
        Numpy Array: An array of radial velocities along a line of sight.
    """
    return (omega - W_0) * (R_0 * np.sin(longitude) * np.cos(latitude))

def calculate_intersection(observed_velocity, radial_velocities):
    """ Calculatees the intersection between an observed velocities and the rotation curve
        along the line of sight of a specific pulsar.

    Args:
        observed_velocity (float): The observed velocity of a bound of a pulsar.
        radial_velocities (Numpy Array): The radial velocities calculated along the line 
        of sight of a pulsar.

    Returns:
        tuple: The distance and radial velocity at which it intersects on the rotation 
        curve along the line of sight.
    """
    observed_velocity_list = np.full((1, len(radial_velocities)), observed_velocity)
    differences = np.abs(observed_velocity_list - radial_velocities)
    max_threshold_index = np.argmax(radial_velocities)
    min_threshold_index = np.argmin(radial_velocities)
    index = 0
    if radial_velocities[0] < observed_velocity and not max_threshold_index == 0:
        for difference in differences[0, :max_threshold_index]:
            if difference <= np.min(differences[0, :max_threshold_index]):
                break
            index+=1
        return (distances[index], radial_velocities[index])
    elif radial_velocities[0] < observed_velocity and not max_threshold_index == 0 and min_threshold_index < max_threshold_index:
        for difference in differences[0, min_threshold_index:max_threshold_index]:
            if difference <= np.min(differences[0, min_threshold_index:max_threshold_index]):
                break
            index+=1
        return (distances[index], radial_velocities[index])
    elif radial_velocities[0] > observed_velocity and not min_threshold_index == 0:
        for difference in differences[0, :min_threshold_index]:
            if difference <= np.min(differences[0, :min_threshold_index]):
                break
            index+=1
        return (distances[index], radial_velocities[index])
    elif radial_velocities[0] > observed_velocity and not min_threshold_index == 0 and max_threshold_index < min_threshold_index:
        for difference in differences[0, max_threshold_index:min_threshold_index]:
            if difference <= np.min(differences[0, max_threshold_index:min_threshold_index]):
                break
            index+=1
        return (distances[index], radial_velocities[index])
        
def calculate_universal_rotation_curve(radius):
    """Calculating the new universal rotation curve presented in the appendix of the Reid et al. 2019 paper

    Args:
        radius (list): Galactocentric Radii

    Returns:
        list: Radial velocities at the given Galactocentric Radii
    """
    radius = np.array(radius)
    lam = np.power((A3 / 1.5), 5)
    ropt = A2 * R0
    rho = radius / ropt
    log_lam = np.log10(lam)
    term1 = 200 * np.power(lam, 0.41)
    top = 0.75 * np.exp(-0.4 * lam)
    bot = 0.47 + 2.25 * np.power(lam, 0.4)
    term2 = np.sqrt(0.80 + 0.49 * log_lam + (top/bot))
    top = 1.97 * np.power(rho, 1.22)
    bot = np.power(rho**2 + 0.61, 1.43)
    term3 = (0.72 + 0.44 * log_lam) * (top/bot)
    top = np.power(rho, 2)
    bot = np.power(rho, 2) + 2.25 * np.power(lam, 0.4)
    term4 = 1.6 * np.exp(-0.4 * lam) * (top/bot)
    return (term1/term2) * np.sqrt(term3 + term4)

def get_radial_velocities(omega, l, b):
    """Calculates the radial velocities along a particular line of sight

    Args:
        omega (Numpy Array): The angular velocities along a line of sight.
        longitude (float): The longitude of the pulsar.
        latitude (float): The latitude of the pulsar.

    Returns:
        list: Radial velocities along the line of sight
    """
    return (omega - (THETA_0 / R0)) * (R0 * np.sin(l) * np.cos(b)) 

def reid_galactocentric_radius(longitude):
    """ Calculates the galactocentric radius at a given longitude and uses discretized
        distances.

    Args:
        longitude (float): The longitude at a specific line of sight.

    Returns:
        Numpy Array: An array that represents a range of galactocentric radii.
    """
    return np.sqrt(R0**2 + distances**2 - 2* R0 *distances * np.cos(longitude))