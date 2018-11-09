import time
import os
import os.path
import re
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import operator
from operator import itemgetter
from mpinterfaces.utils import add_vacuum_padding
from pymatgen.io.vasp.inputs import Incar
from decimal import Decimal
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.matproj.rest import MPRester
from pymatgen.io.vasp.inputs import Kpoints
from openpyxl import Workbook
from openpyxl.chart import (ScatterChart,Reference,Series)
import numpy as np

pseudos = {'H':'H','He':'He','Li':'Li_sv','Be':'Be','B':'B','C':'C','N':'N',
'O':'O','F':'F','Ne':'Ne','Na':'Na_pv','Mg':'Mg','Al':'Al','Si':'Si','P':'P','S':'S'
,'Cl':'Cl','Ar':'Ar','K':'K_sv','Ca':'Ca_sv','Sc':'Sc_sv','Ti':'Ti_sv','V':'V_sv'
,'Cr':'Cr_pv','Mn':'Mn_pv','Fe':'Fe','Co':'Co','Ni':'Ni','Cu':'Cu','Zn':'Zn','Ga':'Ga_d','Ge':'Ge_d',
'As':'As','Se':'Se','Br':'Br','Kr':'Kr','Rb':'Rb_sv','Sr':'Sr_sv','Y':'Y_sv','Zr':'Zr_sv','Nb':'Nb_sv','Mo':'Mo_sv',
'Tc':'Tc_pv','Ru':'Ru_pv','Rh':'Rh_pv','Pd':'Pd','Ag':'Ag','Cd':'Cd','In':'In_d','Sn':'Sn_d','Sb':'Sb','Te':'Te','I':'I',
'Xe':'Xe','Cs':'Cs_sv','Ba':'Ba_sv', 'La':'La','Ce':'Ce','Pr':'Pr_3','Nd':'Nd_3','Pm':'Pm_3','Sm':'Sm_3','Eu':'Eu_2','Gd':'Gd_3','Tb':'Tb_3',
'Dy':'Dy_3','Ho':'Ho_3','Er':'Er_3','Tm':'Tm_3','Yb':'Yb_2','Lu':'Lu_3','Hf':'Hf_pv','Ta':'Ta_pv','W':'W_pv','Re':'Re'
,'Os':'Os','Ir':'Ir','Pt':'Pt','Au':'Au','Hg':'Hg','Tl':'Tl_d','Pb':'Pb_d','Bi':'Bi_d','Po':'Po_d','At':'At_d','Rn':'Rn','Fr':'Fr_sv'
,'Ra':'Ra_sv','Ac':'Ac','Th':'Th','Pa':'Pa','U':'U','Np':'Np','Pu':'Pu','Am':'Am','Cm':'Cm'}

class Log:
    
    def __init__(self, log_dir = '/ufrc/hennig/djordje.gluhovic/log',file_name='log'):
        cur_dir = os.getcwd()
        if os.path.isdir(log_dir):
            os.chdir(log_dir)
            self.log_dir = log_dir
            if os.path.exists(file_name):
                self.file_name = file_name
                os.chdir(cur_dir)
            else:
                os.chdir(cur_dir)
		raise IOError('File '+file_name+' does not exist') 

    def current_time(self):
        Time = time.strftime(("%H:%M:%S"))
        Date = (time.strftime("%d:%m:%Y")).replace(':','/')
        return ' '.join([Date,Time])
        
    def append_data(self,data):
       Current_dir = os.getcwd()
       os.chdir(self.log_dir)
       File = open(self.file_name,'a')
      # current_time1 = (self.current_time().split())[1] 
       File.write(data)
       os.chdir(Current_dir)
       File.close()
       
    def check_if_a_new_day(self):
	cur_dir = os.getcwd()
	os.chdir(self.log_dir)
        os.system('grep "Date:  " '+self.file_name+' > temp.txt') # grep to substitute FINDSTR on hipergator
        temp = open('temp.txt','r').read()
        current_date = (self.current_time()).split()[0]
        if current_date in temp:
            os.system('rm temp.txt') # rm on hipergator
	    os.chdir(cur_dir)
            return False
        os.system('rm temp.txt') # rm on hipergator
	os.chdir(cur_dir)
        return True
        
        
      
    def write_event(self,comment):
       check = self.check_if_a_new_day()
       if check:
            date = (self.current_time().split())
            A = ('Date:  '+date[0])
            self.append_data('\n\n')
            self.append_data(A)
            self.append_data(comment+'  '+date[1])
       else:
            self.append_data(comment+'  '+date[1])
    
            
    def check_if_done(self):
        if os.path.exists('OUTCAR'):
            OUT = open('OUTCAR','r').readlines()
            l =len(OUT)
            for i in xrange(l-1,int(90*l/100),-1):
                    if 'reached required accuracy' in OUT[i]:
                                return True
            return 0
        return 0
    
    def proper(self):
        if os.path.exists('OSZICAR'):
            OSZ = open('OSZICAR','r').readlines()
            if len(OSZ) < 3:
                print 'OSZICAR empty'
                return 0
            numb = re.findall('^ *[0-9]+', OSZ[len(OSZ)-1])
            number_of_ionic_steps = int(re.sub(' *','',numb[0]))
            if number_of_ionic_steps > 5:
                return 0
            else: return 1
        else:
            return 0                
                                                                                        
    def get_name(self):
        if os.path.exists('POSCAR'):
            POS = open('POSCAR','r').readlines()
            name = POS[5].split()
            num = POS[6].split()
            if int(num[1]) > int(num[0]):
                Name = name[0]+num[0]+name[1]+num[1]
            else:
                Name = name[1]+num[1]+name[0]+num[0]
            return Name
  
    def extract_energy(self,DIR = os.getcwd()):
        OSZ = open('OSZICAR','r').readlines()
        sake = re.sub('.*E0= ','',OSZ[len(OSZ)-1])
        sake = re.sub(' *d .*\n','',sake)
        sake = float(sake)
        return sake
        
        
    def formation_energy_anti_fero_magnetism(self,twoD_DIR,threeD_DIR):
        twoList = []
        threeList = []

        os.chdir(twoD_DIR)
        Titel = 'MPs 2D_energy 3D_energy formation_energy\n'
        form = open('formation_energies_anti_fero_magnetism','w')
        form.write(Titel)
        two_info = open('energy_info_anti_fero_magnetism','r').readlines()
        for i in range(len(two_info)):
            twoList.append(two_info[i].split())
    
        os.chdir(threeD_DIR)
        three_info = open('energy_info_anti_fero_magnetism','r').readlines()
        for i in range(len(three_info)):
            threeList.append(three_info[i].split())

        for i in range(len(twoList)):
            for j in range(len(threeList)):
                if twoList[i][0] == threeList[j][0]:
                    form_energy = (float(twoList[i][1])-float(threeList[j][1]))/num(os.path.join(twoD_DIR,twoList[i][0]))*1000
                    form.write(twoList[i][0]+'  '+twoList[i][1]+'   '+threeList[j][1]+'  '+str(form_energy)+'\n')
        form.close()
    
    def create_energy_dic_magnetism(self,DIR):
        DIC = {}
        for child in os.listdir(DIR):
            path = os.path.join(DIR,child,'magnetism')
            if os.path.isdir(path):
                os.chdir(path)
                if self.check_if_done() and self.proper():
                    DIC.update({child: self.extract_energy()})
        return DIC
    
    def formation_energy_magnetism(self,twoD_DIR,threeD_DIR):
        twoList = []
        threeList = []

        os.chdir(twoD_DIR)
        Titel = 'MPs 2D_energy 3D_energy formation_energy\n'
        form = open('formation_energies_magnetism','w')
        form.write(Titel)
        two_info = open('energy_info_magnetism','r').readlines()
        for i in range(len(two_info)):
            twoList.append(two_info[i].split())

        os.chdir(threeD_DIR)
        three_info = open('energy_info_magnetism','r').readlines()
        for i in range(len(three_info)):
            threeList.append(three_info[i].split())

        for i in range(len(twoList)):
            for j in range(len(threeList)):
                if twoList[i][0] == threeList[j][0]:
                    form_energy = (float(twoList[i][1])-float(threeList[j][1]))/num(os.path.join(twoD_DIR,twoList[i][0]))*1000
                    form.write(twoList[i][0]+'  '+twoList[i][1]+'   '+threeList[j][1]+'  '+str(form_energy)+'\n')
        form.close()

    def write_energy_from_dic(self,choice,DIC = {}):
        if choice == 0:
            raise SystemExit

        if choice ==1:
            energy = open('energy_info','w')
        if choice ==2:
            energy = open('energy_info_magnetism','w')
        if choice ==3:
            energy = open('energy_info_anti_fero_magnetism','w')
        if choice ==4:
            energy = open('energy_info_spin_x','w')
        
        if choice ==5:
            energy = open('energy_info_spin_z','w')
        
        for key,value in DIC.iteritems():
            energy.write(key+'  '+str(value)+'\n')
        energy.close()
    
    def extract_magnetisation(self):
        if self.check_if_done():
            OUT = open('OSZICAR','r').readlines()
            l = len(OUT)
            OSZ = OUT[l-1]
	    mag = re.findall('mag=.*\n',OSZ)[0].strip()
	    try:
	        mag= mag.split('     ')[1]
	    except:
		mag= mag.split('    ')[1]
	    return float(mag)
	
        
    def get_KPOINTS_mesh(self,KPO='KPOINTS'):
        try:
            KP = open(KPO,'r').read()
            if KP:
                return KP
            else:
                print 'Empty KPOINTS file found'
                return 'Empty KPOINTS file found'
        except:
            raise ValueError('KPOINTS file does not exists in '+os.getcwd())
    
    
    def get_coordinates(self,CON = 'CONTCAR'):
        fpath = os.path.join(os.getcwd(),CON)
        fpath_POS = os.path.join(os.getcwd(),'POSCAR')
        
        if os.path.isfile(fpath) and os.path.getsize(fpath) > 0:
            CON = open('CONTCAR','r').read()
            return CON
       
        elif os.path.isfile(fpath_POS) and os.path.getsize(fpath_POS) > 0:
            POS = open('POSCAR','r').read()
            return POS
        else:
            raise ValueError('No POSCAR or CONTCAR found in '+os.getcwd())
        
    def get_INCAR(self,input_file = 'INCAR'):
        fpath = os.path.join(os.getcwd(),input_file)
        if os.path.isfile(fpath) and os.path.getsize(fpath) > 0:
            IN = open(input_file,'r').read()
            return IN
        else:
            raise ValueError('No '+input_file+' found in '+os.getcwd())
            
    def get_submit(self,submit='submit'):
        fpath = os.path.join(os.getcwd(),submit)
        if os.path.isfile(fpath) and os.path.getsize(fpath) > 0:
            SUB = open(submit,'r').read()
            return SUB
            
        elif os.path.exists('runjob'):
            SUB = open('runjob','r').read()
            return SUB
        else:
            raise ValueError('No submit script found in '+os.getcwd())
    
    def get_POTCAR_order(self):
        if not os.path.exists('POTCAR'):
            raise ValueError('No POTCAR found in '+os.getcwd())
        os.system('grep TITEL POTCAR >> temp')
        TEMP = open('temp','r').readlines()
        for i in range(len(TEMP)):
            TEMP[i] = TEMP[i].strip(' ')
        os.system('rm temp')
        return ''.join(TEMP)
        
    def get_mp_id(self): 
        DIR = os.getcwd()
        mp = re.search('mp.+?/',DIR).group(0)
        mp = mp.strip('/')
        return mp
        
    def write_comment(self,comment):
        L = len(comment)
        space = ' '*4
        self.append_data('*'*L+'**********\n')
        self.append_data('*'+space+comment+space+'*\n')
        self.append_data('*'*L+'**********\n')
        
    def store_all_info_of_single_structure(self,comment =''): #Needs to be in the directory where the files are located
        END_DEL = '##END##\n'
        s = ' '*4
        INCAR = self.get_INCAR().split('\n')
        Coordinates = self.get_coordinates().split('\n')
        KPOINTS = self.get_KPOINTS_mesh().split('\n')
        submit = self.get_submit().split('\n')
        if '/magnetism' in os.getcwd():
            magnetism = self.extract_magnetisation()
        POTCAR = self.get_POTCAR_order().split('\n')

        self.append_data('Date:   '+(self.current_time()).split()[0]+'  '+(self.current_time()).split()[1]+'\n')
        self.append_data('Structure: '+os.getcwd()+'\n')
        self.append_data(s+'INCAR: \n')
        self.append_data(''.join(('\t'+IN+'\n' for IN in INCAR)))
        self.append_data('\n')
        self.append_data(s+'KPOINTS: \n')
        self.append_data(''.join(('\t'+KP+'\n' for KP in KPOINTS)))
        self.append_data('\n')
        self.append_data(s+'POSCAR: \n')
        self.append_data(''.join(('\t'+Coor+'\n' for Coor in Coordinates)))
        self.append_data('\n')
        self.append_data(s+'POTENTIALS: \n')
        self.append_data(''.join(('\t'+POT+'\n' for POT in POTCAR)))
        self.append_data('\n')
        self.append_data(s+'Submit script: \n')
        self.append_data(''.join(('\t'+sub+'\n' for sub in submit)))
        self.append_data('\n')
        if ('hse' in os.getcwd()) or ('pbe' in os.getcwd()):
            if self.check_if_HSE_done():
                energy = self.extract_energy();
                time_of_computation = self.extract_time_of_computation()
                self.append_data(s+'Energy computed:  ')
                self.append_data(str(energy)+'eV'+'\n')
                self.append_data(s+'Time of computation:  ')
                self.append_data(time_of_computation+' s'+'\n')
                if '/magnetism' in os.getcwd():
                    self.append_data(s+'Magnetism computed:  ')
                    magnetism = self.extract_magnetisation()
                    self.append_data(str(magnetism)+'\n')
        else:
            if self.check_if_done():
                energy = self.extract_energy();
                time_of_computation = self.extract_time_of_computation()
                self.append_data(s+'Energy computed:  ')
                self.append_data(str(energy)+'\n')
                self.append_data(s+'Time of computation:  ')
                self.append_data(time_of_computation+' s'+'\n')
                if '/magnetism' in os.getcwd():
                    self.append_data(s+'Magnetism computed:  ')
                    magnetism = self.extract_magnetisation()
                    self.append_data(str(magnetism)+'\n')
 
        if comment != '':
            self.write_comment(comment)
        self.append_data(END_DEL)
        
    def search_structure(self,structure):
            Dates = []
            start_index = 0
            cur_dir = os.getcwd()
            os.chdir(self.log_dir)
            Log = open(self.file_name,'r').read()
            Str_info = re.findall(structure+'(.*?)##END##',Log,re.DOTALL)
            if len(Str_info) == 0:
                print 'No information about that structure exists'
                return 

            X=[x.start() for x in re.finditer(structure, Log)] #Find indices of all matches for a given structure
            end_index = int(X[0])
            for i in range(len(X)):
                date_index = Log.find('Date:',start_index,end_index) #Find where date starts
                if date_index >-1:
                    new_line_index = Log.find('\n',date_index,date_index+40) # within next 40 characters a newline must exists, due to previous formating
                    date = Log[date_index:new_line_index] # String  = 'Data   mm/dd/yyyy'
                    date = date.split('   ')[1]
                    Dates.append(date)
                    print date
                else:
                    print Dates[len(Dates)-1]
        	print structure        
                print Str_info[i]
                start_index = int(X[i])
		if i < len(X)-1:
                    end_index = int(X[i+1])
               
            os.chdir(cur_dir)
    def extract_data_to_file(self,structure,name_of_file=''):
    	    Dates = []
            cur_dir = os.getcwd()
            os.chdir(self.log_dir)
            start_index = 0
            if name_of_file == '':
                mp = re.search('mp.+',structure).group(0)
                mp = mp.strip('/')
		mp = re.sub('/','_',mp) 
                name_of_file = mp
            if os.path.exists(name_of_file):
                struct = open(name_of_file,'a')
            else:    
                struct = open(name_of_file,'w')
            Log = open(self.file_name,'r').read()
            Str_info = re.findall(structure+'(.*?)##END##',Log,re.DOTALL)
            if len(Str_info) == 0:
                print 'No information about that structure exists'
                struct.close()
                return

            X=[x.start() for x in re.finditer(structure, Log)] #Find indices of all matches for a given structure
            end_index = int(X[0])
            for i in range(len(X)):
                date_index = Log.find('Date:',start_index,end_index) #Find where date starts
                if date_index >-1:
                    new_line_index = Log.find('\n',date_index,date_index+40) # within next 40 characters a newline must exists, due to previous formating
                    date = Log[date_index:new_line_index] # String  = 'Data   mm/dd/yyyy'
                    date = date.split('   ')[1]
                    Dates.append(date)
                    struct.write(date+'\n')
                else:
                    struct.write(Dates[len(Dates)-1]+'\n')
                struct.write(structure)
                struct.write(Str_info[i])
		struct.write('\n')
                start_index = int(X[i])
                if i < len(X)-1:
                    end_index = int(X[i+1])
            struct.close()
            os.chdir(cur_dir)

    def extract_time_of_computation(self,OUTCAR = 'OUTCAR'):
        if os.path.exists(OUTCAR):
            OUT = open(OUTCAR,'r').readlines()
            L = len(OUT)
            for i in xrange(L-1,L-200,-1):
                if 'Total CPU time used (sec):' in OUT[i]:
                    string = OUT[i]
                    time = re.findall('\d+\.\d+',string)
                    if len(time) > 0:
                        return time[0]
                    else:
                        print 'Something went wrong while extractinf time of computation in '+ os.getcwd()
                        return 'No data'
        else:
            print 'No '+OUTCAR +' present'
            return 'No data'

    def check_if_HSE_done(self, OUTCAR = 'OUTCAR'):
        if os.path.exists(OUTCAR):
            OUT = open(OUTCAR,'r').readlines()
            l =len(OUT)
            for i in xrange(l-1,int(98*l/100),-1):
                    if 'General timing and accounting informations for this job:' in OUT[i]:
                                return True
            return 0
        return 0
