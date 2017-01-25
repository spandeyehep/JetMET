from math import pi

def deltaPhi(phi1, phi2):
    dphi = phi2-phi1
    if  dphi > pi:
        dphi -= 2.0*pi
    if dphi <= -pi:
        dphi += 2.0*pi
    return abs(dphi)

def deltaR2(l1, l2):
    return deltaPhi(l1['phi'], l2['phi'])**2 + (l1['eta'] - l2['eta'])**2

def jetID(j):
    if abs(j.eta())<3.0:
        jetId = True if abs(j.eta())>2.4 else j.chargedHadronEnergyFraction()>0 and j.chargedMultiplicity()>0 and j.chargedEmEnergyFraction()<0.99
        return jetId and j.neutralHadronEnergyFraction()<0.99 and j.neutralEmEnergyFraction()<0.99 and j.chargedMultiplicity()+j.neutralMultiplicity()>1
    else:
        return j.neutralEmEnergyFraction()<0.9 and j.neutralMultiplicity()>10

