
def mobility_low_field(N, muMin, muMax, Ng, l):
    m = (muMin + (muMax - muMin) / (1 + ((N / 1e6) / Ng) **l))/10000

    return m

#def calculate_mobility(material, holes, N, x=0.0, y=0.0, T=300):
    """ Calculates the mobility using the model by Sotoodeh et al. If the material is not in the database, then the function returns the mobility for GaAs at that temperature, T, and impurity concentration, N.

    :param material: A string with the material name
    :param holes: If calculation should be done for electrons (holes=0) or holes (holes=1)
    :param N: Impurity concentration
    :param x: The fractional composition in the case of ternaries
    :param y: The other fractional composition in the case of quaternaries
    :param T: Temperature
    :return: The calculated mobility
    """
# To convert it from cm2 to m2

#m = mobility_low_field(N / 1e6, muMin, muMax, Ng, l) / 10000  

