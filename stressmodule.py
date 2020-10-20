import numpy as np

def StrainRosette(gauge_directions, gauge_measurements):
    """
        Calculates the canonical axial stresses on the element based 
        on the strain gauges measurements.

        Parameters:
        -----------
        gauge_directions: numpy array-like
            1D array that contains the strain gauges directions.
        
        gauge_measurements: numpy array-like
            1D array that contains the strain gauges measurements.

        -------
        Return: numpy array-like
            1D array that contains the stresses in the canonical
            directions.  
    """
    gauge_directions = (np.pi/180)*gauge_directions
    coefList = [[(np.cos(angle))**2, (np.sin(angle))**2, np.sin(2*angle)/2] for angle in gauge_directions]
    coefMatrix = np.array(coefList)
    return np.linalg.solve(coefMatrix, gauge_measurements)

def StrainToStress(canonical_strains, young_module, poisson_module):
    """
        Calculates the normal stresses on the element based on it's 
        the canonical strains.

        Parameters:
        -----------
        canonical_strains: numpy array-like
            1D array that contains the normal strains.
        
        young_module: float
            The material Young's Module.

        poisson_module: float
            The material Poisson's Module.

        -------
        Return: numpy array-like
            1D array that contains the normal stresses.   
    """
    coefMatrix = ( -1*poisson_module*np.ones((3,3)) + (poisson_module+1)*np.identity(3) )/young_module
    return np.linalg.solve(coefMatrix, canonical_strains)

def PrincipalStresses(stress_tensor):
    """
        Calculates the normal stresses on the element based on it's 
        the canonical strains.

        Parameters:
        -----------
        stress_tensor: numpy array-like
            3x3 array that contains the Cauchy's stress tensor for the element.

        -------
        Return: numpy array-like
            1D array that contains the principal stresses for the element.   
    """
    return np.linalg.eigvals(stress_tensor)

def VonMisesStress(stress_tensor):
    """
        Calculates the Von Mises stress on the element based on it's 
        Cauchy's stress tensor.

        Parameters:
        -----------
        canonical_strains: numpy array-like
            3x3 array that contains the Cauchy's stress tensor.

        -------
        Return: float
            The Von Mises Stress Tensor.   
    """
    deviatoric_tensor = (stress_tensor - (np.trace(stress_tensor)/3)*np.identity(3))
    von_mises = np.multiply(deviatoric_tensor, deviatoric_tensor) 
    return np.sqrt(1.5*np.sum(von_mises))

