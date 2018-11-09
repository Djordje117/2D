from collections import Counter
from pymatgen.core.composition import Composition
from functions import pseudos
from pymatgen.core.structure import *
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

Elem_radius = {i:Element(i).atomic_radius for i in pseudos.keys()}

def cluster(S, tol = 1.1,seed_index=0, write_poscar_from_cluster=False):
# S = Structure object
# tol = tolerance bonding

        S = SpacegroupAnalyzer(S,0.1).get_conventional_standard_structure()
   
  
        if len (S) < 20: S.make_supercell(2)
        Distance_matrix = S.distance_matrix 
        #np.fill_diagonal(Distance_matrix,100) 
		
        radii = [Elem_radius[site.species_string] for site in S.sites]
        radiiT = np.array(radii)[np.newaxis].T #transpose of radii list
        radii_matrix = radii + radiiT*tol
        temp = Distance_matrix -radii_matrix
        binary_matrix = (temp < 0).astype(int)
        
        original=list(range(len(S)))
        final=[]
        difference=[]
        counter = 0
        while(len(original)!=len(final)):
            seed = set((np.where( binary_matrix[seed_index]==1 ))[0])
            cluster  = seed
            NEW = seed
            check = True
            while check:
                Ss = set()
                for n in NEW:
                    Ss.update(set(np.where( binary_matrix[n]==1 )[0]))
                    #print n,set(np.where( binary_matrix[n]==1 )[0])
                
                if Ss.issubset(cluster):
                    check = False
                else:
                    NEW = Ss - cluster
                    cluster.update(Ss)  
            L = []
            cl = list(cluster)
            sites = [S[i] for i in cl]
            S.from_sites(sites).to('POSCAR', 'POSCAR{}'.format(str(counter)))
            counter +=1
            final = set(final)
            final.update(cluster)
            final = list(final)
            difference = list(set(original)-set(final))
            if(len(difference)>0):
                seed_index = difference[0]
        
