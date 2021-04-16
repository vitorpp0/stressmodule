import numpy as np
import pandas as pd

def strain_rosette(gauge_directions, gauge_measurements):
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

def strain_to_stress(canonical_strains, young_module, poisson_module):
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
    dimension = canonical_strains.shape[0]

    coefMatrix = ( -1*poisson_module*np.ones((dimension,dimension)) + (poisson_module+1)*np.identity(dimension) )/young_module
    return np.linalg.solve(coefMatrix, canonical_strains)

def principal_stresses(stress_tensor):
    """
        Calculates the normal stresses on the element based on it's 
        the canonical strains.

        Parameters:
        -----------
        stress_tensor: numpy array-like
            2x2 or 3x3 array that contains the Cauchy's stress tensor for the element.

        -------
        Return: numpy array-like
            1D array that contains the principal stresses for the element.   
    """
    return np.linalg.eigvals(stress_tensor)

def stress_state(stress_tensor, rotation_angle, decimal_places=6):
    """
        Calculates the normal stresses on the element based on it's 
        the canonical strains.

        Parameters:
        -----------
        stress_tensor: numpy array-like
            2x2 or 3x3 array that contains the Cauchy's stress tensor for the element.

        rotation_angle: float or numpy array-like [degrees]
            If the stress_tensor is 2X2 then the rotation_angle must be angle of
            a counter-clockwise rotation.
            If the stress_tensor is 3x3 the the rotation_angle must be a 1x3 array 
            with the yaw, pitch, roll, in this order, angles of a intrinsic rotation
            about the principal axis.     
        
        decimal_places: Int, default=6
            The answers number of decimal places. 
        -------
        Return: numpy array-like
            Pandas DataFrame with the element stress state.   
    """
    try:
        if not isinstance(decimal_places,int):
            raise TypeError
    except TypeError:
        print("Error: decimal_places must be a integer.")
        raise
    try:
        if not isinstance(stress_tensor,np.matrix):
            raise TypeError
    except TypeError:
        print("Error: stress_tensor must be numpy-matrix.")
        raise
    try:
        if stress_tensor.shape[0] != stress_tensor.shape[1]:
            raise AttributeError
        elif (stress_tensor.shape[0] <=1 or stress_tensor.shape[0] >=4):
            raise AttributeError
        else:
            pass
    except AttributeError:
        print("Error: stress_tensor must be a 2x2 or 3x3 square matrix.")
        raise
    try:
        if (stress_tensor.shape[0]==3 and 
            not isinstance(rotation_angle, np.ndarray) ):
            raise TypeError
    except TypeError:
        print("Error: The rotation_angle must be 3x1 numpy-ndarray.")
        raise
    try:
        if (stress_tensor.shape[0]==2 and 
            not ( isinstance(rotation_angle, int) or isinstance(rotation_angle, float) )
            ):
            raise TypeError
    except TypeError:
        print("Error: The rotation_angle must be integer or float.")
        raise

    rotation_angle = np.pi*rotation_angle/180

    space_dimension = stress_tensor.shape[0]
    if space_dimension==2:
        rotation_matrix = np.matrix([
            [np.cos(rotation_angle), -np.sin(rotation_angle)],
            [np.sin(rotation_angle), np.cos(rotation_angle)]
        ])
    else:
        rotation_matrix = np.matrix([
            [np.cos(rotation_angle[0]),-np.sin(rotation_angle[0]),0],
            [np.sin(rotation_angle[0]),np.cos(rotation_angle[0]),0],
            [0,0,1]
        ]).dot(
            np.matrix([
                [np.cos(rotation_angle[1]),0,np.sin(rotation_angle[1])],
                [0,1,0],
                [-np.sin(rotation_angle[1]),0,np.cos(rotation_angle[1])]
            ])
        ).dot(
            np.matrix([
                [1,0,0],
                [0,np.cos(rotation_angle[2]),-np.sin(rotation_angle[2])],
                [0,np.sin(rotation_angle[2]),np.cos(rotation_angle[2])]
            ])
        )
    
    stress_vector = stress_tensor.dot(rotation_matrix)
    normal_stresses = np.sum(np.multiply(stress_vector,rotation_matrix), axis=0)

    if space_dimension==2:
        shear_stresses = np.sum(np.multiply(stress_vector,rotation_matrix[:,[1,0]]), axis=0)
        index, columns = ["normal(dir:x',y')", "shear(dir:y',x')"], ["x'","y'"]
    else:
        shear_stresses =  np.sum(np.multiply(stress_vector,rotation_matrix[:,[2,0,1]]), axis=0)
        index, columns = ["normal(dir:x',y',z')]", "shear(dir:z',x',y')"], ["x'","y'","z'"]
    
    stress_state = np.vstack([normal_stresses,shear_stresses])
    stress_state = pd.DataFrame(stress_state, index=index, columns=columns)

    return stress_state.round(decimal_places)

def von_mises_stress(stress_tensor):
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

