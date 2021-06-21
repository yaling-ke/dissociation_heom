from scipy.constants import physical_constants as phy_const

hartree = phy_const["Hartree energy in eV"][0]                            # 1 hartree=27.21138602 eV
hartree_to_K = phy_const["hartree-kelvin relationship"][0]                # 1 hartree=315775.13 K
hartree_to_cm_1 = phy_const["hartree-inverse meter relationship"][0]/100  # 1 hartree=219474.6313702 cm^{-1}
k_B = phy_const["Boltzmann constant in eV/K"][0]                          # k_B in units of eV/K
au_to_fs = phy_const["atomic unit of time"][0]*1e15                       # 1 atomic unit of time = 2.41888432e-2 fs
h_planck=phy_const["Planck constant in eV s"][0]

hbar_eVs=phy_const["Planck constant over 2 pi in eV s"][0]
hbar_Js=phy_const["Planck constant over 2 pi"][0]
eV_J=phy_const["electron volt-joule relationship"][0]
au2mass=phy_const["atomic unit of mass"][0]
au2momentum=phy_const["atomic unit of mom.um"][0]
au2eV=phy_const["atomic unit of electric potential"][0]
au2m=phy_const["atomic unit of length"][0]
amu=phy_const["atomic mass constant"][0]
Bohr=phy_const["Bohr radius"][0]

