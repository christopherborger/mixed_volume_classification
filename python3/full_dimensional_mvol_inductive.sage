load("polytope_tuples.sage")
load("volume_classification.sage")
from itertools import *



def subtriples_up_to_equivalence(container_triples):
    """
        container_triples - list of triples
        
        return, up to equivalence (!), all full-dimensional triples (P_1,P_2,P_3) that are contained in triples (Q_1,Q_2,Q_3) from the container_triples,
        where under containment we understand the inclusions P_i \subseteq Q_i for i=1,2,3.
    """
    print("subtriples_up_to_equivalence")
    if len(container_triples)==0:
        return [] 

    b = max([total_number_of_points_in_triple(t) for t in container_triples])
    
    # here we'll save the triples. subtriples[i] will contain the list of triples with the total number of points equal to i.
    subtriples={}
    
    # we build the normal forms of cayley(2*A,2*B,2*C) of triples (A,B,C) to test if it's already found, up to equivalence.
    cayleys={}

    # The minimum total number of points in a triple is 4+4+4=12
    # In the beginning, the dictionaries contain empty lists and empty sets, respectively.
    for i in range(12,b+1):
        # Looking for subtriples with exactly i total number of points
        subtriples[i]=[]
        cayleys[i] = set([])
        
    for t in container_triples:
        myCayley=cayley_normal_form(t)
        myNumberOfPoints=total_number_of_points_in_triple(t)
        # just in case, container_triples have some repetitions, up to equivalence.
        if myCayley not in cayleys[myNumberOfPoints]:
            subtriples[myNumberOfPoints].append(t)
            cayleys[myNumberOfPoints].add(myCayley)
    
    # Iteratively pealing off the vertices, we go through all subtriples.
    # We go descedingly starting from the largest possible total number of points.
    
    for p in range(b,12,-1):
        print("The case of the total number of points equal to",p)
        for t in subtriples[p]:
            # go through the polytopes in the triple
            for i in range(0,3):
                P = t[i]
                # go through the vertices of the current polytope from the triple
                for v in P.vertices():
                    # peal off this vertex like that:
                    L= list(P.integral_points())
                    w = vector(v)
                    L.remove(w)
                    # this is the polytope with a vertex v pealed off
                    Q = Polyhedron([list(l) for l in L])
                    # it might happen that Q.dim() <= 2. We do not deal with such cases here.
                    if Q.dim()>=3:
                        # smallerTriple is derived from t by replacing P[i] with Q
                        smallerTriple=t[:]
                        smallerTriple[i]=Q
                        cayley = cayley_normal_form(smallerTriple)
                        # this smallerTriple has p-1 points in total
                        if not cayley in cayleys[p-1]:
                            subtriples[p-1].append(smallerTriple)
                            cayleys[p-1].add(cayley)

    # Converting from the dictionary of subtriples to a list
    subtriples_list = []
    for k in subtriples.keys():
        for t in subtriples[k]:
            subtriples_list.append(t)
    
    return subtriples_list

def empty_subsimplices_up_to_symmetry(container):
    """
        container is a lattice polytope. 
        
        returns the list of empty simplices S contained in container up to the following equivalence:
        
        
        Two empty simplices S and T are equivalent if there exists an affine unimodular transformation 
        that preserves container and sends S to T.
    """
    print("empty_subsimplices_up_to_symmetry(container)")

    # the normal forms that we introduce will enable to check the above equivalence.
    empty_subsimplices = []
    normal_forms = set([])

    cntr=1
    
    # trying to choose a vertex set for an empty simplex...
    for V in combinations(container.integral_points(),4):
        P = Polyhedron(vertices=V)
        # if P is an empty simplex...
        if P.dim()==3 and len(P.integral_points())==4:
            # we can see if it's new...
            snf = sandwich_normal_form(P,container)       
            if not snf in normal_forms:
                # if P is a new empty simplex, we append it...
                empty_subsimplices.append(P)
                normal_forms.add(snf)
                print("Found a new empty simplex in container",cntr)
                cntr+=1

    return empty_subsimplices

def empty_subsimplices(container):
    """
        container - lattice polytope to determine the empty subsimplices of

        return a list of all empty subsimplices of a container
    """
    empty_subsimplices = []

    for V in combinations(container.integral_points(),4):
        P = Polyhedron(vertices=V)
        if P.dim()==3 and len(P.integral_points())==4:
            empty_subsimplices.append(P)

    return empty_subsimplices

def append_sandwich_up_to_translation(sf,A,B):
    """
        (A,B) is a sandwich, that is, A, B are polytpoes with A \subseteq B
        
        sf is 'sandwich factory' used to keep sandwiches up to translations. That means,
        (A,B) and (C,D) are equivalent if A+u=C and B+u=D holds for some translation vector u.
        
        sf is a dictionary
        where the keys represent the volume gaps between B and A for the
        pairs they contain and whose values are sets of tuples of A,B)
        
        (A,B) is appended to sf if it is new (up to translations)
    
        NOTE: in contrast to the sandwiches in the volume classification we consider
        the tuples A,B up to common translations but not up to general common
        affine unimodular transformations
    """
    vlGap=Volume(B)-Volume(A)
    if vlGap not in sf.keys():
        sf[vlGap]=set([])
    sf[vlGap].add(common_translative_normal_form(tuple([A,B])))

def create_translative_sandwich_factory(initial_polytopes,container,m):
    """
        initial_polytopes - list of initial values for A
        container - prescribed container for B
        m - upper bound on the volume of B

        return a sandwich factory containing tuples A,B where A is in the
        list given by initial_polytopes and B is a subset of container such that the union of A with every point in B has volume at most m.
    """
    sf = {}
    for A in initial_polytopes:
        B=iterative_pealing_off(A,container,volumeBound=m)
        append_sandwich_up_to_translation(sf,A,B)
    return sf                

def subpolytopes_using_sandwiches(container,fixed_polytope,m):
    """
        container - bounding box in which to search for subpolytopes
        m - bound both for the volume of the subpolytopes and for V(subpolytope,subpolytope,fixed_polytope)

        return a list of subpolytopes B of container, such that V(B,B,fixed_polytope)<=m
        and V(B)<=m
    """    
    print("subpolytopes_using_sandwiches(container,fixed_polytope,m)") 
    print("list of subpolytopes B of container, such that V(B,B,fixed_polytope)<=m and V(B)<=m")
    
    if are_homothetic(fixed_polytope,container):
        # We could do only the else-case, but we do the case distinction to optimize the performance.
        # Think about the case fixed_polytope=Delta_3 and container is a multiple of Delta_3.
        # There can be many empty simplices inside the container. 
        # Since in the homothetic case, fixed_polytope and the container have the same GL_3(Z) symmetries up to translations.
        empty_simplices = empty_subsimplices_up_to_symmetry(container)
    else:
        empty_simplices = empty_subsimplices(container) 

    # one polytope is fixed (fixed_polytope), and 
    # the next one is to be fixed: it is sandwiched between the empty simplex and the container

    starting_simplices = []
    
    # test if the empty simplices themselves actually satisfy the prescribed bounds 
    # on volume and mixed volume. 
    for E in empty_simplices:
        if Volume(E)<=m and mixed_volume([E,E,fixed_polytope])<=m:
            starting_simplices.append(E)

    sf = create_translative_sandwich_factory(starting_simplices,container,m)    

    maxVolGap=max(sf.keys())
    
    while maxVolGap>0:
        print("Considering volume gap", maxVolGap)
        for sandwich in sf[maxVolGap]:
            A,B=sandwich
            
            for v in B.vertices(): # pick a vertex of B which is not in A
                if v not in A:
                    break
            
            # case distinction according to whether our second polytope to be fixed contains v
            # if it does we can enlarge A...
            blow_up_of_A=Polyhedron(list(A.vertices())+[v])
            # if it does not, we can reduce B
            reduction_of_B=Polyhedron([z for z in B.integral_points() if vector(z)!=vector(v)])
            
            if mixed_volume([blow_up_of_A,blow_up_of_A,fixed_polytope]) <= m: 
                append_sandwich_up_to_translation(sf,blow_up_of_A,iterative_pealing_off(blow_up_of_A,B,m))

            append_sandwich_up_to_translation(sf,A,reduction_of_B)
            
        del sf[maxVolGap]
        maxVolGap=max(sf.keys())
        
    result=[]
    
    for A,B in sf[0]:
        result.append(A)

    return result

def classify_equality_AB_pairs(mv):
    """
       mv - value of V(A,A,B)

       return all full-dimensional pairs A,B such that 
       V(A,A,B)=V(B,B,A)=mv holds. Note that this condition implies V(A)<=mv, V(B)<=mv
    """
    # NOTE: this could have been unified by introducing a class for storing stuff up to a normal form
    pairs_and_cayleys = {'tuples':[], 'cayleys':set([])}
    
    # we fix a possible volume for A (the first polytope)
    for m in range(1,mv+1):
        # we iterate over polytopes A of dimension 3 with the volume m
        for A in lattice_polytopes_with_given_dimension_and_volume(3,m):
            print("New A:")
            print(A.vertices())
            A = translative_normal_form(A)
            # around each B with V(A,A,B)=mv and dim(B)=3, we can circumscribe
            # a so-called bounding polytope bp with the property that 
            # each facet normal of bp is a facet normal of A.
            # we enumerate all such polytopes...
            bounding_polytopes = maximal_third_polytopes(A,A,mv,3)
            for bp in bounding_polytopes:
                print("New Bounding_Polytope:")
                print(bp.vertices())
                for B in subpolytopes_using_sandwiches(bp,A,mv):
                    if mixed_volume([A,A,B])==mv and mixed_volume([B,B,A])==mv:
                        insert_up_to_equivalence(pairs_and_cayleys,[A,B])
                
    return pairs_and_cayleys['tuples']

def classify_mv_inductive_approach(mv):
    """
        assuming that the files data/dim_3_mv_*.txt with *=1,...,mv-1
        contains all full-dimensional maximal triples of mixed volumes from 1 to mv-1, respectively,
        the function returns the list of all full-dimensional maximal triples A,B,C with V(A,B,C)=mv.
        
        By Proposition 5.2 there are two cases to be enumerated: 
        a) for some (A,B) in the triple t one has V(A,A,B)<mv
        b) for all (A,B) in the triple t one has V(A,A,B)=mv (and, as a consequence also V(A)<=m for each A in the triple t)
        The function is built on distinction of these two cases. 
        Case a) is the inductive part.
    """
    # Case a)
    
    fname='data/dim_3_mv_%d.txt'

    # we load the maximal full-dimensional triples of lower mixed volume from respective files
    maximal_triples = []
    for m in range(1,mv):
        maximal_triples = maximal_triples + polytope_tuples_from_file(fname % m)
    
    # we generate all full-dimensional triples out of the maximal ones that we've loaded
    lower_mv_triples = subtriples_up_to_equivalence(maximal_triples)
    
    
    # Pick the triples of the form (A,A,B) up to ordering and produce the list of pairs (A,B) out of them
    inductive_pairs = []
    
    for t in lower_mv_triples:
        p = get_pair(t)
        if p != None:
            inductive_pairs.append(p)
    
    # those are the triples that arise from case a)
    inductive_maximal_triples = classify_maximal_triples_from_pairs(inductive_pairs,mv,3)
    
    # non-inductive part ----------------------------------------------
    
    # first we classifiy (A,B) with V(A,A,B)=V(A,B,B)=mv
    new_pairs = classify_equality_AB_pairs(mv)

    # A and B are given, C is to be found so that (A,B,C) is maximal
    maximal_triples_from_new_pairs = classify_maximal_triples_from_pairs(new_pairs,mv,3)
    
    # Note that while V(A,A,B)=V(A,B,B)=mv is guaranteed we do not know if something like V(C,B,B)=mv holds. Such kind
    # of conditions are imposed by case b). So, we need to filter.
    new_maximal_triples = []

    for t in maximal_triples_from_new_pairs:
        # (A,B,C) is our triple, we already know that V(A,A,B)==V(B,A,A)==mv. 
        # But we don't know about V(C,C,A)==V(C,C,B)==V(C,A,A)==V(C,B,B)==mv. So we filter.
        if mixed_volume([t[2],t[2],t[0]])==mv and mixed_volume([t[2],t[2],t[1]])==mv and mixed_volume([t[2],t[0],t[0]])==mv and mixed_volume([t[2],t[1],t[1]])==mv:
            new_maximal_triples.append(t)

    # TODO: maybe, change the output format to have a list and not a nested list of two lists.
    return [inductive_maximal_triples,new_maximal_triples] 
    
