def deltaPhi(phi1,phi2):
    phi = phi1-phi2
    while phi >= np.pi: phi -= np.pi*2.
    while phi < -np.pi: phi += np.pi*2.
    return phi

def delta_R(eta1, phi1, eta2, phi2):
    return np.sqrt(deltaPhi(phi1,phi2)**2+(eta1-eta2)**2)

def gaussian_fit(target):
    mean = np.average(target)
    _sigma = 0
    for i in range(len(target)):
        _sigma += (target[i] - mean)**2
    sigma = np.sqrt(_sigma/len(target))

    return mean, sigma 

def to_matrix(target_1, target_2):
    len_parton = len(target_1)
    len_jet = len(target_2)
    _matrix = np.zeros(len_parton*len_jet)
    matrix =  _matrix.reshape((len_parton, len_jet))
    for x in range(len_parton):
        _value = int(target_1[x])
        matrix[x][_value] = 1
            
    return matrix


