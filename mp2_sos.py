#!/usr/bin/env python
import numpy as np
from pyscf.pbc import gto, scf, df
from pyscf.pbc import lib as pbc_lib
# from pyscf.pbc.mp import kmp2_time_idan as kmp2
# from pyscf.pbc.mp import kmp2_new as kmp2
from numpy import array
from pyscf import lib
# from pyscf.pbc.mp import kmp2_lap_sos_n4_tami_time as kmp2_lap_sos
# from pyscf.pbc.mp import kmp2_lap_sos
from pyscf.pbc.tools.pbc import super_cell
from pyscf.pbc.tools import pyscf_ase
from pyscf.pbc.tools import lattice
from pyscf.gto.basis import parse_nwchem
# from pyscf.pbc.mp import kmp2_sep_reg_k
# from pyscf.pbc.mp import kmp2_SCS_reg
# import pyscf.pbc.tools.lattice as tlattice
from pbcflow import sc40 as tlattice
import re
import warnings

warnings.filterwarnings("ignore")

import sys
import time
from importlib import import_module
# from pyscf import gto as mol_gto
from pyscf import scf as mol_scf
from parse_vasp_file_aug_ghost import parse_vasp_file

# Import relevant modules from pyscf library
kmp2_lap_sos = import_module("pyscf.pbc.mp.kmp2_lap_sos_" + str(sys.argv[8]))
if str(sys.argv[10]) == '2.0.1':
    kmp2 = import_module("pyscf.pbc.mp.kmp2_new_time")
else:
    # kmp2= import_module("pyscf.pbc.mp.dfkmp2")
    kmp2 = import_module("pyscf.pbc.mp.kmp2_time_idan")

##############################
# Create a "Basis" from Hong-zhouye - in relation to HZ from mul_script, lc is an acronym for large core
##############################
if str(sys.argv[12]) == 'yes':
    basname = str(sys.argv[5])

    if sys.argv[16] and sys.argv[16] != '':
        print("arg 16 (path_to_basis_data) is ", sys.argv[16])
        path_to_basis_data = sys.argv[16]
        fbas = f'{path_to_basis_data}/{basname}.dat'
        fbas_aug = f'{path_to_basis_data}/aug-{basname}.dat'
else:
    raise "No valid path to basis sets was provided"

compounds = str(sys.argv[1])

# Filter atoms, augmented atoms using RegEx
arrays = re.findall('[A-Z][^A-Z]*', compounds)
arrays_aug = [atm + '*' for atm in arrays]

atms = arrays
atms_aug = arrays_aug

# Create the basis sets for regular atoms and augmented atoms
basis_reg = {atm: parse_nwchem.load(fbas, atm) for atm in atms}
basis_aug = {atm: parse_nwchem.load(fbas_aug, atm) for atm in atms_aug}

# Merge dictionaries (Python 3.5 and above)
unified_basis = basis_reg | basis_aug
print("Unified basis is:")
print(unified_basis)

##############################
# Create a "Cell"               - we need to change here, according to vasp's unit(bohr\atomic\angstrom)
##############################

cell = gto.Cell()
cell.unit = 'Angstrom'
if sys.argv[14] != 'default':
    # Custom vasp file requested by user - display variables for sanity checks
    print("arg 14 (relevant_vasp_file) is ", sys.argv[14])
    relevant_vasp_file = sys.argv[14]
    print("arg 15 (relevant_vasp_path) is ", sys.argv[15])
    relevant_vasp_path = sys.argv[15]
    lattice_vectors, cell_structure, ghost_dict = parse_vasp_file(f'{relevant_vasp_path}/{relevant_vasp_file}.vasp')
    print(ghost_dict)

else:
    # default structure - this line to be removed or edited in the future #TODO #FIXME
    lattice_vectors, cell_structure, ghost_dict = parse_vasp_file(
        '/private/projects/MgO_CO/supporting_data/2023/arXiv:2309.14651/geom/surf+mol/2l1a_331.vasp')

# Display lattice vectors and cell structure for sanity check
print("Lattice vectors are:")
print(lattice_vectors)
print("Cell structure is:")
print(cell_structure)
cell.a = lattice_vectors
cell.atom = cell_structure

# Take care of creating the unique basis sets for regular and augmented ghost atoms
if ghost_dict:
    # Create the ghost key:value pairs to be inserted to a dictionary later
    atm_types_ghost = list(ghost_dict.keys())
    X_map_ghost = list(ghost_dict.values())
    ghost_aug_dic = {key: value + '*' for key, value in ghost_dict.items()}
    atm_types_ghost_aug = list(ghost_aug_dic.keys())
    X_aug_map_ghost = list(ghost_aug_dic.values())

    # NOTE: Ghost atoms are marked with the format of "X" prefix, some number suffix

    # Use dict comprehension and zip attributes to allow pyscf library to parse relevant basis sets
    basis_ghost = {ghost_atm: gto.basis.load(fbas, atm_type_ghost) for ghost_atm, atm_type_ghost in
                   zip(X_map_ghost, atm_types_ghost)}
    basis_ghost_aug = {ghost_aug_atm: gto.basis.load(fbas_aug, atm_type_ghost_aug) for ghost_aug_atm, atm_type_ghost_aug
                       in zip(X_aug_map_ghost, atm_types_ghost_aug)}
    unified_basis = unified_basis | basis_ghost | basis_ghost_aug
    print("Unified basis with ghost is:")
    print(unified_basis)

cell.basis = unified_basis
cell.pseudo = "gth-hf-rev"

cell.precision = 1e-16
cell.verbose = 7
# cell.max_memory = 1000000
# cell.max_memory = 40000
cell.incore_anyway = True
cell.build()

# Change lattice constants and atomic positions proportionally
# rescale = float(sys.argv[4])
# cell.a = cell.lattice_vectors() * rescale
# cell.atom = cell.format_atom(cell.atom, unit=cell.unit)
# cell.atom = [(atom[0], list(np.array(atom[1]) * rescale)) for atom in cell.atom]


# cell.build()

##############################
#  KRHF with Gaussian density fitting (GDF)
##############################

nk = float(sys.argv[2])
c_nk = float(sys.argv[3])
# number of kpts for each axis
# nmp = [7, 7, 7]
nmp = [nk, nk, nk]
# nmp = [nk, nk, 1]
# nmp= [nk, 1, 1]
# set the origin of kpt mesh
scaled_center = [c_nk, c_nk, c_nk]
# scaled_center=[c_nk, c_nk, c_nk]
print("nmp", nmp, "scaled_center", scaled_center)
kpts = cell.make_kpts(nmp, scaled_center=scaled_center)

# set up k-point scf (khf)
kmf = scf.KRHF(cell, kpts=kpts, exxdiv="ewald")
kmf.max_cycle = 100
# mymf = scf.KRHF(cell, kpts=kpts, exxdiv="None")
# kmf = kmf.density_fit()
# auxbasis_i = {atm: mol_gto.basis.load('/private/ccgto/aux/cc-pvdz-fit.dat', atm) for atm in atms}

if str(sys.argv[6]) == "ri":
    kmf = kmf.density_fit(auxbasis='cc-pv' + str(sys.argv[9]) + 'z-ri')
    # kmf=kmf.rs_density_fit(auxbasis=auxbasis_i)
    # kmf.with_df =df.df.DF(cell,kpts)
    # kmf.with_df.auxbasis = 'cc-pv'+str(sys.argv[9])+'z-ri'
    # kmf.with_df =df.rsdf.RSGDF(cell,kpts)
    # kmf.with_df.auxbasis = 'cc-pv'+str(sys.argv[9])+'z-ri'
elif str(sys.argv[6]) == "custom":
    basname = str(sys.argv[5])
    # change here the path to relevant "supporting_data"
    # if 'lc' in basname:
    #    fbas = '/private/ccgto/basis/gth-hf-rev/%s.dat' % basname
    # else:
    #    fbas = '/private/ccgto/basis/ccecp/%s.dat' % basname

    # fbas = '/private/projects/MgO_CO/supporting_data/2023/arXiv:2309.14651/basis/%s-jkfit.dat' % basname
    # aug_fbas = '/private/projects/MgO_CO/supporting_data/2023/arXiv:2309.14651/basis/aug-%s-jkfit.dat' % basname

    if sys.argv[16] and sys.argv[16] != '':
        path_to_basis_data = sys.argv[16]
        # fbas = f'/private/projects/MgO_CO/supporting_data/2023/arXiv:2309.14651/basis/%s.dat' % basname
        # fbas_aug = f'/private/projects/MgO_CO/supporting_data/2023/arXiv:2309.14651/basis/aug-%s.dat' % basname
        fbas = f'{path_to_basis_data}/{basname}-jkfit.dat'
        aug_fbas = f'{path_to_basis_data}/aug-{basname}-jkfit.dat'
    else:
        raise "No valid path to basis sets was provided"

    compounds = str(sys.argv[1])
    arrays = re.findall('[A-Z][^A-Z]*', compounds)
    # arrays = ['Mg', 'O']
    arrays_aug = [atm + '*' for atm in arrays]
    atms = arrays
    atms_aug = arrays_aug
    # atms = [str(sys.argv[1])]
    aux_basis = {atm: parse_nwchem.load(fbas, atm) for atm in atms}
    aug_aux_basis = {atm: parse_nwchem.load(aug_fbas, atm) for atm in atms_aug}
    unified_aux_basis = aux_basis | aug_aux_basis
    print("Auxiliary basis set is:")
    print(unified_aux_basis)
    if ghost_dict:
        aux_basis_ghost = {ghost_atm: gto.basis.load(fbas, atm_type_ghost) for ghost_atm, atm_type_ghost in
                           zip(X_map_ghost, atm_types_ghost)}
        aux_basis_ghost_aug = {ghost_aug_atm: gto.basis.load(aug_fbas, atm_type_ghost_aug) for
                               ghost_aug_atm, atm_type_ghost_aug in zip(X_aug_map_ghost, atm_types_ghost_aug)}
        unified_aux_basis = unified_aux_basis | aux_basis_ghost | aux_basis_ghost_aug
        print("Ghost auxiliary basis set is:")
        print(unified_aux_basis)

    kmf = kmf.density_fit(auxbasis=unified_aux_basis)


else:
    kmf = kmf.density_fit()

if float(sys.argv[7]) > 0.0:
    kmf.with_df.exp_to_discard = float(sys.argv[7])

# with_df = df.GDF(cell, kpts)
# with_df._prefer_ccdf = True
# kmf = scf.KRHF(cell, kpts, exxdiv="ewald").density_fit(with_df=with_df)

if str(sys.argv[11]) == 'yes':
    kmf = mol_scf.addons.remove_linear_dep_(kmf)

ekrhf = kmf.kernel()

##############################
# Save SCF
##############################

# chkfile = '333_qz_no_pad_results.chk'
# kmf.chkfile = chkfile
# scf.chkfile.dump_scf(cell, chkfile, kmf.e_tot, kmf.mo_energy, kmf.mo_coeff, kmf.mo_occ)

##############################
# KMP2
##############################

# mymp = mp.KMP2(mymf, frozen=0)
# e_mp, t2 = mymp.kernel(with_t2=False)

t0_mp2 = time.time()

# reference KMP2
if str(sys.argv[10]) == '2.0.1':
    kmp = kmp2.KMP2_SOS(kmf)
    emp2, emp2_ref, t2 = kmp.kernel(with_t2=False)

else:
    t0_mp2 = time.time()
    kmp = kmp2.KMP2(kmf)
    t1_mp2 = time.time()
    emp2_ref = kmp.kernel(with_t2=False)[0]
    t1_kernel = time.time()
    emp2_ref = kmp.e_corr_os
    t1_emp2_ref = time.time()
    print("t_mp2", t1_mp2 - t0_mp2)
    print("t_kernel", t1_kernel - t1_mp2)
    print("t_mp2", t1_emp2_ref - t1_kernel)

# mymp_sep_reg=kmp2_sep_reg_k.KMP2_reg_sep(mymf)
# mymp_sep_reg.kappa=[1., 1.45, 2.]
# mymp_sep_reg.kappa=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
# emp2, kmp2_reg, emp2j, emp2_ex, t2 = mymp_sep_reg.kernel(with_t2=False)
# print('E(MP2scs) = %.9g' % emp2scs)
t1_mp2 = time.time()
dt_mp2 = t1_mp2 - t0_mp2

t0_sos = time.time()
# laplace KMP2
kmp_laplace = kmp2_lap_sos.LaplaceSOSKMP2(kmf)

energies = []

# taus = range(6,21)
taus = range(6, 8)

print(taus)
for tau in taus:
    emp2_new = kmp_laplace.kernel(
        with_t2=False, tau=tau
    )[0]
    emp2_new = kmp_laplace.e_corr_os
    diff = emp2_new - emp2_ref
    energies.append([tau, emp2_ref, emp2_new, diff])

    print(f"\n{'Tau':>8} {'EMP2_ref':>20} {'EMP2_new':>20} {'Diff':>20}")

for tau in taus:
    tau, emp2_ref, emp2_new, diff = energies.pop(0)
    print(
        f"{tau:>8} {emp2_ref:>20.12f} {emp2_new:>20.12f} {diff:>20.12f}"
    )
t1_sos = time.time()
dt_sos = t1_sos - t0_sos

print("time_mp2", dt_mp2, "time_sos", dt_sos, "dt", dt_mp2 - dt_sos)
# print( "kmp2:", emp2, "kmp2_reg", kmp2_reg, "emp2j", emp2j, "emp2_ex", emp2_ex )

