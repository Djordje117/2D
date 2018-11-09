from pymatgen.io.vasp.inputs import Incar
from mpinterfaces.utils import get_magmom_string
from pymatgen.core.structure import Structure
import os
from functions import update_sub_dict,submit_builder, update_INCAR_dict
import shutil

def GW_step1(name,spin_polarized=True, submit = True,incar_dict={}, force_overwrite = False):
    #name is the name of the job that will be used in submit file
    incar = Incar.from_string("""ISMEAR = 1 ; SIGMA = 0.2
ALGO = Exact  ; NELM = 1               # exact diagonalization one step suffices
EDIFF = 1E-8                           # high precision for groundstate calculation
NBANDS = 96                            # need for a lot of bands in GW
LOPTICS = .TRUE.                       # we need d phi/ d k  for GW calculations
     
KPAR = 2""")

    if spin_polarized:
        if os.path.exists('POSCAR'):
            incar.update({'LORBIT':11,'ISPIN':2,'MAGMOM':get_magmom_string(Structure.from_file('POSCAR'))})
        else:
            raise IOError('No POSCAR found at '+os.getcwd())
    
    if os.path.isdir('step1'):
        os.chdir('step1')
    else:
        os.mkdir('step1')
        os.chdir('step1')
    
    if force_overwrite:
        os.system('rm -rf *')
        
    incar.update(incar_dict)
    incar.write_file('INCAR')
    os.system('cp ../{KPOINTS,POSCAR,CONTCAR,POTCAR} .')
    
    if os.path.exists('CONTCAR') and os.stat('CONTCAR').st_size > 10:
        shutil.move('CONTCAR','POSCAR')
    
    submit_builder(update_dict={'--mem-per-cpu':'3gb','--job-name':name,'--nodes':2,'--ntasks':32})
    
    
    if submit:
        os.system('sbatch submit')
    
    os.chdir('../')
        
def GW_step2(name,submit=True,spin_polarized=True,incar_dict={},force_overwrite = False):
    #Calculation is supposed to be run in parent directory (step1 should be there as well)
    #incar_dict should not have NBANDS tag, since it will be extracted from step1/INCAR
    incar = Incar.from_string("""
    ISMEAR = -5
NBANDS = 96                            # need for a lot of bands in GW
ALGO = GW0                             # 
NELM = 1                               
PRECFOCK = Fast                        
ENCUTGW = 100                         
NOMEGA = 200                           
KPAR = 2""")
    
    if not os.path.isdir('step1'):
        raise ValueError('It is not possible to run step2 before step1. Check if you are in right directory.')
    
    os.chdir('step1')
    
    if os.path.exists('WAVECAR') and os.path.exists('WAVEDER') and os.stat('WAVECAR').st_size > 100 and os.stat('WAVEDER').st_size > 100:
        os.chdir('../')
        try:
            os.mkdir('step2')
            os.chdir('step2')
        except:
            os.chdir('step2')
    else:
        raise ValueError('step1 calculation did not finish properly. WAVECAR and/or WAVEDER files are not written properly.')
    
    if force_overwrite:
         os.system('rm -rf *')
    os.system('cp ../step1/{WAVE*,KPOINTS,POSCAR,POTCAR} .')
    NBANDS = Incar.from_file('../step1/INCAR')['NBANDS'] 
    incar.update({'NBANDS':NBANDS})
    if spin_polarized:
        if os.path.exists('POSCAR'):
            incar.update({'LORBIT':11,'ISPIN':2,'MAGMOM':get_magmom_string(Structure.from_file('POSCAR'))})
        else:
            raise IOError('No POSCAR found at '+os.getcwd())
   
    incar.write_file('INCAR')     
    submit_builder(update_dict={'--mem-per-cpu':'3gb','--job-name':name,'--nodes':2,'--ntasks':32})
    
    if submit:
        os.system('sbatch submit')

    os.chdir('../')
