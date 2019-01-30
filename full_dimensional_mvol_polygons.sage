load("polytope_tuples.sage")
load("volume_classification.sage")

import os.path
from collections import defaultdict

FILE_NAME_MV='data/dim_%d_mv_%d.txt'

# ================================================
# Front-end functions
# ================================================

def update_dim_2_mixed_volume_database(mvBound):
    # the files storing polytopes are created in the data subfolder
    if not os.path.exists('data'):
        os.mkdir('data')

    # let's see what files are missing
    missingMixedVolumes = [mv for mv in range(1,mvBound+1) if not os.path.isfile(FILE_NAME_MV % (2,mv))]

    if len(missingMixedVolumes)>0:
        # we should run mixed volume classification.
        
        # but first, we adjust mvBound
        mvBound = max(missingMixedVolumes)
        
        result=maximal_pairs_classification(mvBound)
        # result[vl] contains all lattice polytopes with Volume equal to vl
        
        for mv in missingMixedVolumes: 
            f=open(FILE_NAME_MV % (2,mv),'w')
            print >>f, [[map(tuple,P.vertices()) for P in t] for t in result[mv]]
            f.close()


def maximal_pairs_of_polygons_with_given_mixed_volume(mv):
    """
        That's the main function for users of this 'module'. It returns the list of all maximal pairs of lattice polygons with mixed volume mv
    """
    update_dim_2_mixed_volume_database(mv)
    # now, we can read the list polytopes from the corresponding file and return them
    f=open(FILE_NAME_MV % (2,mv),'r')
    L=eval(f.read().replace('\n',' '))
    f.close()
    return [[Polyhedron(P),Polyhedron(Q)] for P,Q in L]

def pairs_of_polygons_with_given_mixed_volume(mv):
    # pairs is a dictionary of dictionaries: 
    # it provides access to pairs as follows:
    # volume of B in pair -> normal form of pair (A,B) -> the pair (A,B)
    
    # initialization of the pairs dictionary
    pairs=defaultdict(dict)
    for A,B in maximal_pairs_of_polygons_with_given_mixed_volume(mv):
        NF=cayley_normal_form([A,B])
        totalVolume=Volume(A)+Volume(B)
        pairs[totalVolume][NF]=[A,B]
    
    result=[]
    while len(pairs)>0:
        print sorted(pairs.keys())
        totalVolume=max(pairs.keys())
        for A,B in pairs[totalVolume].values():
            for v in A.vertices():
                smallerA=Polyhedron([z for z in A.integral_points() if tuple(z)!=tuple(v)])
                newTotalVolume=Volume(smallerA)+Volume(B)
                newPair=[smallerA,B]
                if smallerA.dim()==2 and mixed_volume(newPair)==mv:
                    NF=cayley_normal_form(newPair)
                    pairs[newTotalVolume][NF]=newPair
            for v in B.vertices():
                smallerB=Polyhedron([z for z in B.integral_points() if tuple(z)!=tuple(v)])
                newTotalVolume=Volume(A)+Volume(smallerB)
                newPair=[A,smallerB]
                if smallerB.dim()==2 and mixed_volume(newPair)==mv:
                    NF=cayley_normal_form(newPair)
                    pairs[newTotalVolume][NF]=newPair
        result+=pairs[totalVolume].values()
        del pairs[totalVolume]
            
    return result


# ===================================================
# Back-end (=auxiliary) functions
# ===================================================

def maximal_second_polygons_up_to_mv(A,max_mv,d_min):
    
    maximal_polygons = set([])

    for mv in range(1,max_mv+1):
        maximal_polygons = maximal_polygons.union(maximal_second_polygons(A,mv,d_min))

    return maximal_polygons


def maximal_pairs_classification(max_mv):
    """
        returns a dictionary of all maximal pairs A,B of polygons 
        with V(A,B)<=max_mv
        the keys are values, 1,..,max_mv
        For each mv, the value is a list of pairs A,B with V(A,B)=mv
    """
    
    # we need to be sure that we have all polygons with normalized volume at most max_mv
    update_volume_classification_database(2,max_mv)
    
    result={}
    
    for vl in range(1,max_mv+1):
        for A in lattice_polytopes_with_given_dimension_and_volume(2,vl):
            A = translative_normal_form(A)
            maximal_polygons = maximal_second_polygons_up_to_mv(A,max_mv,2)
            for B in maximal_polygons:
                if is_maximal_pair([A,B]):
                    mv=mixed_volume([A,B])
                    if mv not in result.keys():
                        result[mv]={}
                    NF=cayley_normal_form([A,B])
                    if NF not in result[mv].keys():
                        result[mv][NF]=[A,B]

    for mv in result.keys():
        result[mv]=result[mv].values()

    return result
    
                


