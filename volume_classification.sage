load("polytopes.sage")
#load("volume_classification.sage.py")

import logging
import os.path
import sys

# Using the logging package one can conveniently turn off and on the auxiliary messages  

logging.basicConfig(format='%(message)s',stream=sys.stdout,level=logging.INFO)
# After modifing the level from, say, logging.INFO to logging.WARNING , the change will come into force only after _restarting the sage session_ and reloading

# Sandwich is a pair of lattice polytopes A,B with A being a subset of B. 
# The volume gap of a sandwich A,B is the difference Volume(B)-Volume(A).

# that's the template for names of files, in which we store polytopes
FILE_NAME_VOLUME='data/dim_%d_volume_%d.txt'


def prepare_sandwiches(d,volumeBound):
    for A in empty_simplices(d,volumeBound):
        # first, we let B to be an homothetic copy of A as follows
        scalingFactor=(d+1)*(volumeBound/Volume(A)-1)+1
        B=scalingFactor*A+(1-scalingFactor)*A.center() # note: A.center() is the barycenter of the vertex set of A
        # the above polyhedron B is not necessarily integral, so we replace it by its integer hull:
        B=Polyhedron(B.integral_points()) 
        # still, B can contain some redundant integral points, and so we do:
        B=iterative_pealing_off(A,B,volumeBound)
        yield A,B

def pealed_off_container(A,B,volumeBound):
    """
        For a given sandwich (A,B) and a volumeBound
        the function
        returns a polytope
        obtained by taking the lattice points of B and removing all vertices v of B 
        with the property that if v is added to A, then after adding the volume of A exceeds the volumeBound.
        
        NOTE: we could have pealed off the lattice points of B in one step. That might be more efficient.
    """
    to_peal_off=[]
    for v in B.vertices():
        blow_up_of_A=Polyhedron(list(A.vertices())+[v])
        if Volume(blow_up_of_A)>volumeBound:
            to_peal_off.append(vector(v))
    Z=[vector(z) for z in B.integral_points()]
    return Polyhedron([z for z in Z if z not in to_peal_off])

def iterative_pealing_off(A,B,volumeBound):
    pealed_off_B=pealed_off_container(A,B,volumeBound)
    while pealed_off_B != B:
        B=pealed_off_B
        pealed_off_B=pealed_off_container(A,B,volumeBound)
    return B

def layered_polytope_from_sandwich(A,B):
    """ 3*B is embeded into height 0, two copies of 3*A are embedded into heights 1 and -1.
        Then, one generates a polytope based on these three layers at heights -1,0 and 1
    """ 
    middleLayer=[tuple(3*vector(v))+(0,) for v in B.vertices()]
    upperLayer=[tuple(3*vector(v))+(1,) for v in A.vertices()]
    lowerLayer=[tuple(3*vector(v))+(-1,) for v in A.vertices()]
    return Polyhedron(middleLayer+upperLayer+lowerLayer)

def sandwich_normal_form(A,B):
    """
        returns data that allows to distinguish two sandwiches (A,B) 
        (A',B') up to affine unimodular transformations.
    """
    return affine_normal_form(layered_polytope_from_sandwich(A,B))

# Sandwich factory is used to store sandwiches up to affine unimodular transformations.
# A sandwich factory is a dictionary of dictionaries. For each possible volume gap, a storage
# for sandwiches with this volume gap is created. The latter storage
# is a dictionary with key,value pairs such that the value is a sandwich and 
# the respective key is the sandwich normal form of this sandwich.

def append_sandwich(sf,A,B):
    """
        If no affine unimodular image of the sandwich (A,B) is in the sandwich factory sf,
        the sanwich (A,B) is appended to sf.
    """
    vlGap=Volume(B)-Volume(A)
    SNF=sandwich_normal_form(A,B)
    if vlGap not in sf.keys():
        sf[vlGap]={}
    if SNF not in sf[vlGap].keys():
        sf[vlGap][SNF]=[A,B]

            
def new_sandwich_factory(d,volumeBound):
    sandwich_factory={}
    for A,B in prepare_sandwiches(d,volumeBound):
        B=iterative_pealing_off(A,B,volumeBound)
        append_sandwich(sandwich_factory,A,B)
    return sandwich_factory


def sandwich_factory_statistics(sf):
    logging.info("Maximum volume gap in sandwiches: %d",max(sf.keys()))
    logging.info("Number of sandwiches: %d",sum([len(sf[vlGap]) for vlGap in sf.keys() if vlGap!=0]))
    if 0 in sf.keys():
        logging.info("Number of polytopes found: %d", len(sf[0]))
    logging.info(50*"-")

def volume_classification(d,m):
    sf=new_sandwich_factory(d,m)
    maxVolGap=max(sf.keys())
    
    while maxVolGap>0:
        
        sandwich_factory_statistics(sf)
        
        for SNF in sf[maxVolGap].keys():
            A,B=sf[maxVolGap][SNF]
            
            for v in B.vertices(): # pick a vertex of B which is not in A
                if v not in A:
                    break
            
            blow_up_of_A=Polyhedron(list(A.vertices())+[v])
            reduction_of_B=Polyhedron([z for z in B.integral_points() if vector(z)!=vector(v)])
            
            append_sandwich(sf,blow_up_of_A,iterative_pealing_off(blow_up_of_A,B,m))
            append_sandwich(sf,A,reduction_of_B)
            
        del sf[maxVolGap]
        maxVolGap=max(sf.keys())

    sandwich_factory_statistics(sf)
        
    result={}
    for vl in xrange(1,m+1):
        result[vl]=[]
    
    for A,B in sf[0].values():
        result[Volume(A)].append(A)

    return result


def update_volume_classification_database(d,volumeBound):
    # the files storing polytopes are created in the data subfolder
    if not os.path.exists('data'):
        os.mkdir('data')

    # let's see what files are missing
    missingVolumes = [vl for vl in range(1,volumeBound+1) if not os.path.isfile(FILE_NAME_VOLUME % (d,vl))]

    if len(missingVolumes)>0:
        # we should run volume classification.
        
        # but first, we adjust volumeBound
        VolumeBound = max(missingVolumes)
        
        result=volume_classification(d,volumeBound)
        # result[vl] contains all lattice polytopes with Volume equal to vl
        
        for vl in missingVolumes: 
            f=open(FILE_NAME_VOLUME % (d,vl),'w')
            print >>f, [[tuple(p) for p in P.vertices()] for P in result[vl]]
            f.close()
        
def lattice_polytopes_with_given_dimension_and_volume(d,vl):
    """
        That's the main function for users of this module. It returns the list of all d-dimensional lattice polytopes with Volume vl.
        
        For efficiency reasons, it may be advantageous to make an explicit call of update_volume_classification(d,volumeBound).
        For example, if you know you'd work with 3-dimensional lattice polytopes of Volume at most 5,
        first call update_volume_classification_database(3,5) and then go on using lattice_polytopes_with_given_dimension_and_volume(d,vl).
    """
    # first, we update the database of lattice polytopes with a given volume
    update_volume_classification_database(d,vl)
    # now, we can read the list polytopes from the corresponding file and return them
    f=open(FILE_NAME_VOLUME % (d,vl),'r')
    L=eval(f.read().replace('\n',' '))
    f.close()
    return [Polyhedron(P) for P in L]

    
