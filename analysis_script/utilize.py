import numpy as np
import itertools
import pandas as pd

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


#Define Chi square minimizer
def chi_square_minimizer( jet_pt_chi2, jet_eta_chi2, jet_phi_chi2, jet_btag_chi2, jet_mass_chi2):
    
    num_of_btag = np.sum(np.array(jet_btag_chi2) ==1)

    class jet_cand_properties():
        def __init__(self, idx):
            self.idx = idx
            self.pt = jet_pt_chi2[self.idx]
            self.eta = jet_eta_chi2[self.idx]
            self.phi = jet_phi_chi2[self.idx]
            self.mass = jet_mass_chi2[self.idx]
            self.px = self.pt*np.cos(self.phi)
            self.py = self.pt*np.sin(self.phi)
            self.pz = self.pt*np.sinh(self.eta)
            self.e = np.sqrt( (self.px**2 + self.py**2 + self.pz**2) + self.mass**2 )

    def cal_W_inv(jet1, jet2):
        part_1 = (jet1.e + jet2.e)**2
        part_2 = (jet1.px + jet2.px)**2
        part_3 = (jet1.py + jet2.py)**2
        part_4 = (jet1.pz + jet2.pz)**2
        return np.sqrt( part_1 - part_2 - part_3 - part_4 )

    def cal_top_inv(jet1, jet2, jet3):
        part_1 = (jet1.e + jet2.e + jet3.e)**2
        part_2 = (jet1.px + jet2.px + jet3.px)**2
        part_3 = (jet1.py + jet2.py + jet3.py)**2
        part_4 = (jet1.pz + jet2.pz + jet3.pz)**2
        return np.sqrt( part_1 - part_2 - part_3 - part_4 )

    _bjet_list = []
    
    
    bjet = []
    _jet_list = []

    min_chi2 = -1
    m_W = 80.9
    sigma_W = 17.77
    sigma_t = 27.49

    jet_idx_list = np.array(['Nan', 'Nan', 'Nan', 'Nan', 'Nan', 'Nan'])
    
    _jet_index = []
    for i in range(len(jet_pt_chi2)):
        _jet_index.append(i)

    for i in range(len(jet_btag_chi2)):
        if jet_btag_chi2[i] == 1:
            _bjet_list.append(i)
        else :
            _jet_list.append(i)
            
    _bjet = itertools.combinations(_bjet_list, 2)
    
    for a in _bjet:
        bjet.append(a)

    bjet = np.array(bjet, dtype='object')

    #print(_bjet_list)
    #print(bjet)
    jet_index_candidate = []

    for i in range(len(bjet)):
        #print(bjet[i])

        jet = []
        
        tmp_jet_index = _bjet_list.copy()
        for c in range(len(bjet[i])): 
            
            _tmp = bjet[i][c]
            tmp_jet_index.remove(_tmp)
        
        tmp_jet_list = _jet_list.copy()

        for d in tmp_jet_index:
            tmp_jet_list.append(d)
        
        _jet = itertools.permutations(tmp_jet_list, 4)

        for b in _jet:
            jet.append(b)

        jet = np.array(jet, dtype='object')
        #print( jet)
        
        for j in range(len(jet)):
        
            _jet_index_candidate = []
            _jet_index_candidate.append(bjet[i][0])
            _jet_index_candidate.append(bjet[i][1])
            _jet_index_candidate.append(jet[j][0])
            _jet_index_candidate.append(jet[j][1])
            _jet_index_candidate.append(jet[j][2])
            _jet_index_candidate.append(jet[j][3])
            jet_index_candidate.append(_jet_index_candidate)
    #print(jet_index_candidate)   

    for i in range(len(jet_index_candidate)):

        b_1_idx = jet_index_candidate[i][0]
        b_2_idx = jet_index_candidate[i][1]

        j_1_idx = jet_index_candidate[i][2]
        j_2_idx = jet_index_candidate[i][3]
        j_3_idx = jet_index_candidate[i][4]
        j_4_idx = jet_index_candidate[i][5]

        bjet_1 = jet_cand_properties(b_1_idx)
        bjet_2 = jet_cand_properties(b_2_idx)

        jet_1 = jet_cand_properties(j_1_idx)
        jet_2 = jet_cand_properties(j_2_idx)
        jet_3 = jet_cand_properties(j_3_idx)
        jet_4 = jet_cand_properties(j_4_idx)
        
        W_1_inv = cal_W_inv(jet_1, jet_2)
        W_2_inv = cal_W_inv(jet_3, jet_4)
        top_1_inv = cal_top_inv(bjet_1, jet_1, jet_2)
        top_2_inv = cal_top_inv(bjet_2, jet_3, jet_4)

        chi2_part_1 = (top_1_inv - top_2_inv)**2
        chi2_part_2 = (W_1_inv - m_W)**2
        chi2_part_3 = (W_2_inv - m_W)**2
        
        chi2_tmp = chi2_part_1/(2*(sigma_t**2)) + chi2_part_2/sigma_W**2 + chi2_part_3/sigma_W**2
        #print(chi2_tmp, b_1_idx, j_1_idx, j_2_idx, b_2_idx, j_3_idx, j_4_idx)
        if (min_chi2 < 0 or chi2_tmp < min_chi2 ):
            min_chi2 = chi2_tmp
            jet_1_best_idx = j_1_idx
            jet_2_best_idx = j_2_idx
            jet_3_best_idx = j_3_idx
            jet_4_best_idx = j_4_idx
            b_1_best_idx = b_1_idx
            b_2_best_idx = b_2_idx
            jet_idx_list = np.array([b_1_best_idx, jet_1_best_idx, jet_2_best_idx, b_2_best_idx, jet_3_best_idx, jet_4_best_idx])
        else: 
            pass
        
    return min_chi2, jet_idx_list


def chi_square_minimizer_old( jet_pt_chi2, jet_eta_chi2, jet_phi_chi2, jet_btag_chi2, jet_mass_chi2):
            
    
    if (np.sum(np.array(jet_btag_chi2, dtype='object') == 1) > 2):

        length_of_btag = np.sum( np.array(jet_btag_chi2, dtype='object') == 1 )
        a = 0
        while a < length_of_btag-1:
            b = a + 1
            while b < length_of_btag:

                m_W = 80.9
                sigma_W = 9999
                sigma_t = 1288.28
                _btag_idx_tmp = []
                _jet_idx_tmp = []
                W_inv_cand = []
                W_minus_inv_cand = [] 
                top_inv_cand = []
                top_bar_inv_cand = []

                chi_square_value = []
                min_chi2 = -1

                jet_idx_list = np.array(['Nan', 'Nan', 'Nan', 'Nan', 'Nan', 'Nan'])
                for i in range(len(jet_pt_chi2)):
                    if jet_btag_chi2[i] == 1:
                        _btag_idx_tmp.append(i)
                    else :
                        _jet_idx_tmp.append(i)

                _btag_idx_tmp_re = _btag_idx_tmp
                


                _btag_idx_tmp_prime = []
                _btag_idx_tmp_prime.append(_btag_idx_tmp_re[a])
                _btag_idx_tmp_prime.append(_btag_idx_tmp_re[b])

                _a = _btag_idx_tmp_re[a]
                _b = _btag_idx_tmp_re[b]

                _btag_idx_tmp_re.remove(_a)
                _btag_idx_tmp_re.remove(_b)


                for x in range(len(_btag_idx_tmp)):
                    _jet_idx_tmp.append(_btag_idx_tmp_re[x])

                _length = len(_jet_idx_tmp) 

                for i in range(0,_length):

                    _jet_1_idx = _jet_idx_tmp[i]
                    _jet_1_pt = jet_pt_chi2[_jet_1_idx]
                    _jet_1_eta = jet_eta_chi2[_jet_1_idx]
                    _jet_1_phi = jet_phi_chi2[_jet_1_idx]
                    _jet_1_mass = jet_mass_chi2[_jet_1_idx]
                    _jet_1_px = _jet_1_pt*np.cos( _jet_1_phi )
                    _jet_1_py = _jet_1_pt*np.sin( _jet_1_phi )
                    _jet_1_pz = _jet_1_pt*np.sinh( _jet_1_eta )
                    _jet_1_e = np.sqrt( ( _jet_1_px**2 + _jet_1_py**2 + _jet_1_pz**2) + _jet_1_mass**2 )
                    for j in range(i+1, _length):

                        _jet_2_idx = _jet_idx_tmp[j]
                        _jet_2_pt = jet_pt_chi2[_jet_2_idx]
                        _jet_2_eta = jet_eta_chi2[_jet_2_idx]
                        _jet_2_phi = jet_phi_chi2[_jet_2_idx]
                        _jet_2_mass = jet_mass_chi2[_jet_2_idx]
                        _jet_2_px = _jet_2_pt*np.cos( _jet_2_phi )
                        _jet_2_py = _jet_2_pt*np.sin( _jet_2_phi )
                        _jet_2_pz = _jet_2_pt*np.sinh( _jet_2_eta )
                        _jet_2_e = np.sqrt( ( _jet_2_px**2 + _jet_2_py**2 + _jet_2_pz**2) + _jet_2_mass**2 )
                        for k in range(0, _length):

                            _jet_3_idx = _jet_idx_tmp[k]
                            if (_jet_3_idx == _jet_1_idx or _jet_3_idx ==_jet_2_idx):
                                _jet_3_found = 0
                            else:
                                _jet_3_found = 1
                                _jet_3_pt = jet_pt_chi2[_jet_3_idx]
                                _jet_3_eta = jet_eta_chi2[_jet_3_idx]
                                _jet_3_phi = jet_phi_chi2[_jet_3_idx]
                                _jet_3_mass = jet_mass_chi2[_jet_3_idx]
                                _jet_3_px = _jet_3_pt*np.cos( _jet_3_phi )
                                _jet_3_py = _jet_3_pt*np.sin( _jet_3_phi )
                                _jet_3_pz = _jet_3_pt*np.sinh( _jet_3_eta )
                                _jet_3_e = np.sqrt( ( _jet_3_px**2 + _jet_3_py**2 + _jet_3_pz**2) + _jet_3_mass**2 )
                            for l in range(k+1, _length):

                                _jet_4_idx = _jet_idx_tmp[l]
                                if (_jet_4_idx == _jet_1_idx or _jet_4_idx ==_jet_2_idx):
                                    _jet_4_found = 0
                                else :
                                    _jet_4_found = 1
                                    _jet_4_pt = jet_pt_chi2[_jet_4_idx]
                                    _jet_4_eta = jet_eta_chi2[_jet_4_idx]
                                    _jet_4_phi = jet_phi_chi2[_jet_4_idx]
                                    _jet_4_mass = jet_mass_chi2[_jet_3_idx]
                                    _jet_4_px = _jet_4_pt*np.cos( _jet_4_phi )
                                    _jet_4_py = _jet_4_pt*np.sin( _jet_4_phi )
                                    _jet_4_pz = _jet_4_pt*np.sinh( _jet_4_eta )
                                    _jet_4_e = np.sqrt( ( _jet_4_px**2 + _jet_4_py**2 + _jet_4_pz**2) + _jet_4_mass**2 )

                                for m in range(len(_btag_idx_tmp_prime)):

                                    _b_cand_idx_1 = _btag_idx_tmp_prime[m]
                                    _b_cand_pt_1 = jet_pt_chi2[_b_cand_idx_1]
                                    _b_cand_eta_1 = jet_eta_chi2[_b_cand_idx_1]
                                    _b_cand_phi_1 = jet_phi_chi2[_b_cand_idx_1]
                                    _b_cand_mass_1 = jet_mass_chi2[_b_cand_idx_1]
                                    _b_cand_px_1 = _b_cand_pt_1 * np.cos(_b_cand_phi_1)
                                    _b_cand_py_1 = _b_cand_pt_1 * np.sin(_b_cand_phi_1)
                                    _b_cand_pz_1 = _b_cand_pt_1 * np.sinh(_b_cand_eta_1)
                                    _b_cand_e_1 = np.sqrt( ( _b_cand_px_1**2 + _b_cand_py_1**2 + _b_cand_pz_1**2) + _b_cand_mass_1**2 )

                                    for n in range(m+1, len(_btag_idx_tmp_prime)):

                                        if _btag_idx_tmp[m] == _btag_idx_tmp_prime[n] :
                                            _btag_jet_found = 0
                                        else : 
                                            _btag_jet_found = 1
                                            _b_cand_idx_2 = _btag_idx_tmp_prime[n]
                                            _b_cand_pt_2 = jet_pt_chi2[_b_cand_idx_2]
                                            _b_cand_eta_2 = jet_eta_chi2[_b_cand_idx_2]
                                            _b_cand_phi_2 = jet_phi_chi2[_b_cand_idx_2]
                                            _b_cand_mass_2 = jet_mass_chi2[_b_cand_idx_2]
                                            _b_cand_px_2 = _b_cand_pt_2 * np.cos(_b_cand_phi_2)
                                            _b_cand_py_2 = _b_cand_pt_2 * np.sin(_b_cand_phi_2)
                                            _b_cand_pz_2 = _b_cand_pt_2 * np.sinh(_b_cand_eta_2)
                                            _b_cand_e_2 = np.sqrt( ( _b_cand_px_2**2 + _b_cand_py_2**2 + _b_cand_pz_2**2) + _b_cand_mass_2**2 )

                                        if (_jet_3_found == 1 and _jet_4_found == 1 and _btag_jet_found == 1):
                                            _W_inv_mass = np.sqrt( (_jet_1_e + _jet_2_e)**2 
                                                                - (_jet_1_px + _jet_2_px)**2 
                                                                - (_jet_1_py + _jet_2_py)**2 
                                                                - (_jet_1_pz + _jet_2_pz)**2 )
                                            _W_inv_mass_1 = np.sqrt( (_jet_3_e + _jet_4_e)**2 
                                                                - (_jet_3_px + _jet_4_px)**2  
                                                                - (_jet_3_py + _jet_4_py)**2 
                                                                - (_jet_3_pz + _jet_4_pz)**2 )
                                            _top_inv_mass = np.sqrt( (_jet_1_e + _jet_2_e + _b_cand_e_1)**2 
                                                                - (_jet_1_px + _jet_2_px + _b_cand_px_1)**2 
                                                                - (_jet_1_py + _jet_2_py + _b_cand_py_1)**2 
                                                                - (_jet_1_pz + _jet_2_pz + _b_cand_pz_1)**2 )
                                            _top_inv_mass_1 = np.sqrt( (_jet_3_e + _jet_4_e + _b_cand_e_2)**2 
                                                                - (_jet_3_px + _jet_4_px + _b_cand_px_2)**2 
                                                                - (_jet_3_py + _jet_4_py + _b_cand_py_2)**2 
                                                                - (_jet_3_pz + _jet_4_pz + _b_cand_pz_2)**2 )
                                            W_inv_cand.append(_W_inv_mass)
                                            W_minus_inv_cand.append(_W_inv_mass_1)
                                            top_inv_cand.append(_top_inv_mass)
                                            top_bar_inv_cand.append(_top_inv_mass_1)

                                            chi_square_tmp = 0        
                                            chi_square_tmp_1 = (_top_inv_mass - _top_inv_mass_1)**2 
                                            chi_square_tmp_2 = (_W_inv_mass - m_W)**2
                                            chi_square_tmp_3 = (_W_inv_mass_1 - m_W)**2
                                            #print(type(chi_square_tmp_1), type(chi_square_tmp_2), type(chi_square_tmp_3))
                                            chi_square_tmp = chi_square_tmp_1/(2*sigma_t**2) + chi_square_tmp_2/(sigma_W**2) + chi_square_tmp_3/(sigma_W**2)
                                            #print(W_inv_cand, W_minus_inv_cand, top_inv_cand, top_bar_inv_cand, chi_square_tmp)  
                                            if (min_chi2 < 0 or chi_square_tmp < min_chi2 ):
                                                min_chi2 = chi_square_tmp
                                                jet_1_best_idx = _jet_1_idx
                                                jet_2_best_idx = _jet_2_idx
                                                jet_3_best_idx = _jet_3_idx
                                                jet_4_best_idx = _jet_4_idx
                                                b_1_best_idx = _b_cand_idx_1
                                                b_2_best_idx = _b_cand_idx_2
                                                jet_idx_list = np.array([b_1_best_idx, jet_1_best_idx, jet_2_best_idx, b_2_best_idx, jet_3_best_idx, jet_4_best_idx])
                                            else: 
                                                pass
                    b += 1
                a += 1

    else: 
        length_of_btag = np.sum( np.array(jet_btag_chi2, dtype='object') == 1 )
        #print("length_of_btag: {0}".format(length_of_btag))
        m_W = 80.9
        sigma_W = 197.4
        sigma_t = 1288.28
        _btag_idx_tmp = []
        _jet_idx_tmp = []
        W_inv_cand = []
        W_minus_inv_cand = [] 
        top_inv_cand = []
        top_bar_inv_cand = []

        chi_square_value = []
        set_info = []

        min_chi2 = -1

        jet_idx_list = np.array(['Nan', 'Nan', 'Nan', 'Nan', 'Nan', 'Nan'])
        for i in range(len(jet_pt_chi2)):
            if jet_btag_chi2[i] == 1:
                _btag_idx_tmp.append(i)
            else :
                _jet_idx_tmp.append(i)
        #print(_btag_idx_tmp, _jet_idx_tmp)
        _length = len(_jet_idx_tmp)
        for i in range(0,_length):

            _jet_1_idx = _jet_idx_tmp[i]
            _jet_1_pt = jet_pt_chi2[_jet_1_idx]
            _jet_1_eta = jet_eta_chi2[_jet_1_idx]
            _jet_1_phi = jet_phi_chi2[_jet_1_idx]
            _jet_1_mass = jet_mass_chi2[_jet_1_idx]
            _jet_1_px = _jet_1_pt*np.cos( _jet_1_phi )
            _jet_1_py = _jet_1_pt*np.sin( _jet_1_phi )
            _jet_1_pz = _jet_1_pt*np.sinh( _jet_1_eta )
            _jet_1_e = np.sqrt( ( _jet_1_px**2 + _jet_1_py**2 + _jet_1_pz**2) + _jet_1_mass**2 )
            #print("q1 has been found.")
            for j in range(i+1, _length):

                _jet_2_idx = _jet_idx_tmp[j]
                _jet_2_pt = jet_pt_chi2[_jet_2_idx]
                _jet_2_eta = jet_eta_chi2[_jet_2_idx]
                _jet_2_phi = jet_phi_chi2[_jet_2_idx]
                _jet_2_mass = jet_mass_chi2[_jet_2_idx]
                _jet_2_px = _jet_2_pt*np.cos( _jet_2_phi )
                _jet_2_py = _jet_2_pt*np.sin( _jet_2_phi )
                _jet_2_pz = _jet_2_pt*np.sinh( _jet_2_eta )
                _jet_2_e = np.sqrt( ( _jet_2_px**2 + _jet_2_py**2 + _jet_2_pz**2) + _jet_2_mass**2 )
                #print("q2 has benn found.")
                for k in range(0, _length):

                    _jet_3_idx = _jet_idx_tmp[k]
                    if (_jet_3_idx == _jet_1_idx or _jet_3_idx ==_jet_2_idx):
                        _jet_3_found = 0
                    else:
                        _jet_3_found = 1
                        _jet_3_pt = jet_pt_chi2[_jet_3_idx]
                        _jet_3_eta = jet_eta_chi2[_jet_3_idx]
                        _jet_3_phi = jet_phi_chi2[_jet_3_idx]
                        _jet_3_mass = jet_mass_chi2[_jet_3_idx]
                        _jet_3_px = _jet_3_pt*np.cos( _jet_3_phi )
                        _jet_3_py = _jet_3_pt*np.sin( _jet_3_phi )
                        _jet_3_pz = _jet_3_pt*np.sinh( _jet_3_eta )
                        _jet_3_e = np.sqrt( ( _jet_3_px**2 + _jet_3_py**2 + _jet_3_pz**2) + _jet_3_mass**2 )
                        #print("q3 has been found.")
                    for l in range(k+1, _length):

                        _jet_4_idx = _jet_idx_tmp[l]
                        if (_jet_4_idx == _jet_1_idx or _jet_4_idx ==_jet_2_idx):
                            _jet_4_found = 0
                        else :
                            _jet_4_found = 1
                            _jet_4_pt = jet_pt_chi2[_jet_4_idx]
                            _jet_4_eta = jet_eta_chi2[_jet_4_idx]
                            _jet_4_phi = jet_phi_chi2[_jet_4_idx]
                            _jet_4_mass = jet_mass_chi2[_jet_3_idx]
                            _jet_4_px = _jet_4_pt*np.cos( _jet_4_phi )
                            _jet_4_py = _jet_4_pt*np.sin( _jet_4_phi )
                            _jet_4_pz = _jet_4_pt*np.sinh( _jet_4_eta )
                            _jet_4_e = np.sqrt( ( _jet_4_px**2 + _jet_4_py**2 + _jet_4_pz**2) + _jet_4_mass**2 )
                            #print("q4 has been found.")

                        for m in range(len(_btag_idx_tmp)):

                            _b_cand_idx_1 = _btag_idx_tmp[m]
                            _b_cand_pt_1 = jet_pt_chi2[_b_cand_idx_1]
                            _b_cand_eta_1 = jet_eta_chi2[_b_cand_idx_1]
                            _b_cand_phi_1 = jet_phi_chi2[_b_cand_idx_1]
                            _b_cand_mass_1 = jet_mass_chi2[_b_cand_idx_1]
                            _b_cand_px_1 = _b_cand_pt_1 * np.cos(_b_cand_phi_1)
                            _b_cand_py_1 = _b_cand_pt_1 * np.sin(_b_cand_phi_1)
                            _b_cand_pz_1 = _b_cand_pt_1 * np.sinh(_b_cand_eta_1)
                            _b_cand_e_1 = np.sqrt( ( _b_cand_px_1**2 + _b_cand_py_1**2 + _b_cand_pz_1**2) + _b_cand_mass_1**2 )
                            #print("b1 has been found")

                            for n in range(m+1, len(_btag_idx_tmp)):

                                if _btag_idx_tmp[m] == _btag_idx_tmp[n] :
                                    _btag_jet_found = 0
                                else : 
                                    _btag_jet_found = 1
                                    _b_cand_idx_2 = _btag_idx_tmp[n]
                                    _b_cand_pt_2 = jet_pt_chi2[_b_cand_idx_2]
                                    _b_cand_eta_2 = jet_eta_chi2[_b_cand_idx_2]
                                    _b_cand_phi_2 = jet_phi_chi2[_b_cand_idx_2]
                                    _b_cand_mass_2 = jet_mass_chi2[_b_cand_idx_2]
                                    _b_cand_px_2 = _b_cand_pt_2 * np.cos(_b_cand_phi_2)
                                    _b_cand_py_2 = _b_cand_pt_2 * np.sin(_b_cand_phi_2)
                                    _b_cand_pz_2 = _b_cand_pt_2 * np.sinh(_b_cand_eta_2)
                                    _b_cand_e_2 = np.sqrt( ( _b_cand_px_2**2 + _b_cand_py_2**2 + _b_cand_pz_2**2) + _b_cand_mass_2**2 )
                                    #print("b2 has been found.")

                                if (_jet_3_found == 1 and _jet_4_found == 1 and _btag_jet_found == 1):
                                    _W_inv_mass = np.sqrt( (_jet_1_e + _jet_2_e)**2 
                                                        - (_jet_1_px + _jet_2_px)**2 
                                                        - (_jet_1_py + _jet_2_py)**2 
                                                        - (_jet_1_pz + _jet_2_pz)**2 )
                                    _W_inv_mass_1 = np.sqrt( (_jet_3_e + _jet_4_e)**2 
                                                        - (_jet_3_px + _jet_4_px)**2  
                                                        - (_jet_3_py + _jet_4_py)**2 
                                                        - (_jet_3_pz + _jet_4_pz)**2 )
                                    _top_inv_mass = np.sqrt( (_jet_1_e + _jet_2_e + _b_cand_e_1)**2 
                                                        - (_jet_1_px + _jet_2_px + _b_cand_px_1)**2 
                                                        - (_jet_1_py + _jet_2_py + _b_cand_py_1)**2 
                                                        - (_jet_1_pz + _jet_2_pz + _b_cand_pz_1)**2 )
                                    _top_inv_mass_1 = np.sqrt( (_jet_3_e + _jet_4_e + _b_cand_e_2)**2 
                                                        - (_jet_3_px + _jet_4_px + _b_cand_px_2)**2 
                                                        - (_jet_3_py + _jet_4_py + _b_cand_py_2)**2 
                                                        - (_jet_3_pz + _jet_4_pz + _b_cand_pz_2)**2 )
                                    W_inv_cand.append(_W_inv_mass)
                                    W_minus_inv_cand.append(_W_inv_mass_1)
                                    top_inv_cand.append(_top_inv_mass)
                                    top_bar_inv_cand.append(_top_inv_mass_1)
                                            
                                    chi_square_tmp_1 = (_top_inv_mass - _top_inv_mass_1)**2 
                                    chi_square_tmp_2 = (_W_inv_mass - m_W)**2
                                    chi_square_tmp_3 = (_W_inv_mass_1 - m_W)**2
                                    chi_square_tmp = chi_square_tmp_1/(2*sigma_t**2) + chi_square_tmp_2/(sigma_W**2) + chi_square_tmp_3/(sigma_W**2)
                                    #print(W_inv_cand, W_minus_inv_cand, top_inv_cand, top_bar_inv_cand, chi_square_tmp) 
                                        
                                    if (min_chi2 < 0 or chi_square_tmp < min_chi2 ):
                                        min_chi2 = chi_square_tmp
                                        jet_1_best_idx = _jet_1_idx
                                        jet_2_best_idx = _jet_2_idx
                                        jet_3_best_idx = _jet_3_idx
                                        jet_4_best_idx = _jet_4_idx
                                        b_1_best_idx = _b_cand_idx_1
                                        b_2_best_idx = _b_cand_idx_2
                                        jet_idx_list = np.array([b_1_best_idx, jet_1_best_idx, jet_2_best_idx, b_2_best_idx, jet_3_best_idx, jet_4_best_idx])
                                    else: 
                                        pass

        
    return min_chi2, jet_idx_list