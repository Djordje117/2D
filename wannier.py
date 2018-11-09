from pymatgen.symmetry.bandstructure import HighSymmKpath
import os
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints, Incar
from mpinterfaces.utils import remove_z_kpoints
import subprocess
from collections import defaultdict
import re
from functions import *


def find_kpath(f = 'KPOINTS_temp'):
    try:
        A = open(f,'r').read()
    except:
        raise IOError('No file "'+f+'" in '+os.getcwd())
        
    L = A.split('\n')
    temp = []
    temp2 = []
    
    for l in L:
        if '!' in l:
            temp.append(l)
        
    for t in temp:
        t=re.sub('!','',t)
        t=re.sub('\\\\','',t)
        temp2.append(t)

    temp3 = []
    
    for t in temp2:
        a = t.split('  ')
        b = a[1]+' '+a[0]
        temp3.append(b)
    
    final = list(set(temp3))
    string = ''
    for i in range(len(final)-1):
        string +=final[i]+' '+final[i+1]+'\n'
    
    string+= final[i+1]+' '+final[0]
    
    return string

def get_kpoints_wannier(structure, n_kpts=1, dim=2):
    kpath = HighSymmKpath(structure)
    os.system('mv KPOINTS K_temp')
    Kpoints.automatic_linemode(n_kpts, kpath).write_file('KPOINTS')
    if dim == 2:
        remove_z_kpoints()
    path = find_kpath(f='KPOINTS')
    os.system('rm KPOINTS')
    os.system('mv K_temp KPOINTS')
    return path
    
            
def write_wannier_input_file(S,spin_polarized=True,write_incar = True):        
        string = """Begin Projections
random
End Projections
guiding_centres=true
bands_plot = true

begin kpoint_path\n"""

        string +=   get_kpoints_wannier(S)+'\n'    
        string += 'end kpoint_path\n'
        
        if spin_polarized:
            string+='spinors=true\n'
            
        w = open('wannier90.win','w')
        w.write(string)
        w.close()
        
        #write INCAR
        if write_incar:
            incar=Incar.from_string("""ISMEAR = -5
        # usefull energy range for density of states    
ALGO = None ; NELM = 1                 # exact diagonalization one step suffices
NBANDS = 96                            # need for a lot of bands in GW    
LORBIT = 11
LWRITE_MMN_AMN = TRUE
LWANNIER90_RUN = .TRUE.""")
            try:
                NBANDS = Incar.from_file('../../step1/INCAR')['NBANDS'] 
            except:
                NBANDS = Incar.from_file('../INCAR')['NBANDS'] 
            
            if spin_polarized:
                if os.path.exists('POSCAR'):
                    incar.update({'LORBIT':11,'ISPIN':2,'MAGMOM':get_magmom_string(Structure.from_file('POSCAR'))})
                else:
                    raise IOError('No POSCAR found at '+os.getcwd())
            incar.update({'NBANDS':NBANDS})
            incar.write_file('INCAR')

def start_wannier(name,spin_polarized=True, submit = True,incar_dict={}, force_overwrite = False):
    if force_overwrite:
        if os.path.exists('wannier'):
            os.system('rm -rf wannier')
            start_wannier(name,spin_polarized=spin_polarized, submit = submit,incar_dict=incar_dict, force_overwrite = False)
            return
            
    if os.path.exists('wannier'):
        print 'Wannier already done'
        return
    else:
        os.mkdir('wannier')
        os.chdir('wannier')
        os.system('cp ../{PO*,KPOINTS,INCAR,W*} .')
        write_wannier_input_file(Structure.from_file('POSCAR'))
        submit_builder(update_dict={'--mem-per-cpu':'3gb','--job-name':name,'--nodes':1,'--ntasks':4,'--time':'00:30:00'})
        if submit:
            os.system('sbatch submit')
    
