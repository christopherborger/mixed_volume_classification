load("polytope_tuples.sage")
load("volume_classification.sage")
from itertools import *


def irreducible_subtriples_up_to_equivalence(container_triples,max_nr_2d):
    """
        container_triples - list of irreducible triples
        max_nr_2d - maximal numbers of 2-dimensional polytopes allowed in the irreducible triple

        return a list of (up to equivalence) all irreducible
        subtriples of a given list of container_triples that contain at most max_nr_2d-many
        2-dimensional polytopes
    """
    
    print "irreducible_subtriples_up_to_equivalence"
    
    # we add this so that we do not have to bother about degenerate input
    if len(container_triples)==0:
        return [] 
    
    # Preprocessing: container_triples may contain triples with more than max_nr_2d two-dimensional polytopes. We remove those.
    # Since the elements of container triples are irreducible each polytope P in the triple is of dimension 2 or 3.
    # So, 3-c[0].dim() is 0 or 1. Which means that 9 - c[0].dim() - c[1].dim() - c[2].dim() is the number number of two-dim polytopes in the triple c.
    container_triples = [c for c in container_triples if 9 - (c[0].dim()+c[1].dim()+c[2].dim()) <= max_nr_2d]

    # We use total number of lattice points in a triple as the parameter defining complexity of the triple 
    # and generate subtriples by going through this parameter in the decreasing order.
    b = max([total_number_of_points_in_triple(t) for t in container_triples])
    
    # The key of the subtriples is the total number of lattice points
    subtriples={}
    # We store Cayley polytopes of triples to verify if the new generated triple is really new (up to equivalence).
    cayleys={}

    # initialization of subtriples and caleys (9 is the smallest possible total number of lattice points)
    for i in range(9,b+1):
        subtriples[i]=[]
        cayleys[i] = set([])

    # we add our container triples after initialization
    for t in container_triples:
        subtriples[total_number_of_points_in_triple(t)].append(t)
        cayleys[total_number_of_points_in_triple(t)].add(cayley_normal_form(t))
    
    # we go through possible total number of points downwards 
    for p in range(b,9,-1):
        print "total number of points in a triple equal to ",p
        for t in subtriples[p]:
            # choose a polytope in a triple
            for i in range(0,3):
                P = t[i]
                for v in P.vertices():
                    L= list(P.integral_points())
                    w = vector(v)
                    # The above conversion is based on the following observations: integral points are vectors, but vertices are an extra type
                    L.remove(w)
                    Q = Polyhedron([list(l) for l in L])
                    if Q.dim()>=2:
                        # make a copy of a triple
                        Q_triple=t[:] 
                        Q_triple[i]=Q
			if is_irreducible(Q_triple) and 9 - (Q_triple[0].dim()+Q_triple[1].dim()+Q_triple[2].dim()) <= max_nr_2d:
		                cayley = cayley_normal_form(Q_triple)
		                # the total number of lattice points for Q_triple is p-1 (went down by exactly one)
		                if not cayley in cayleys[p-1]:
		                    subtriples[p-1].append(Q_triple)
		                    cayleys[p-1].add(cayley)
    
    # putting all the computed subtriples into one big list and returning result
    
    subtriples_list = []
    for k in subtriples.keys():
        for t in subtriples[k]:
            subtriples_list.append(t)
    
    return subtriples_list


def create_boxes(A,m1,m2):
    """
        Input: A is two-dimensional in the two-dimensional space R^2 and is treated as A \times \{0 \} \subseteq R^3
        

        The function generates a list of bounding boxes for a lattice polytope B, derived from conditions V(A,A,B)=m1 and V(B,B,A)<=m2. 
        
        Without loss of generality it is assumed that B contains (0,0,0) and lies above the height 0 plane. Each computed bounding box 
        is supplied with a heighest point which is assumed to be in the polytope B. See Lemma 6.10 and Remark 6.11.
        
        The returned object is a list of pairs (Container for B, heighest point)
        
        Note that V(A,A,B) is fixing the height of B equal to be m1/ ( relative volume of A)
    """
    # initialize the list of boxes
    boxes = []
    # That's the prescribed height of B
    w = m1/(2*A.volume())
    # We make assumptions on B, first we assume that B contains (0,0,0)
    # Then, up to shearing, we can assume that at height w B has a point whose first and second components are between 0 and w-1
    possible_ps = [vector([x,y]) for x in range(w) for y in range(w)]
    # We iterate over all such points (p1,p2,w) of B at height w
    for p in possible_ps:
        # here, we'll collect points whose convex hull will give the bounding box for B
        points=[]
        # See section 6.5 for details
        bounding_slice = m2 * (A+(-A)).polar()
        trafo = matrix([[0,1/w],[-1/w,0]])
        #apply trafo to bounding_slice
        bounding_slice = Polyhedron([trafo*vector(v) for v in bounding_slice.vertices()])
        #go through the different possible heights
        for h in range(0,w+1):
            shift = 1/w*vector([p[0]*h,p[1]*h])
            # append more points (translate by shift vector, then embed into height h and append to the list of points)
            points=points+[list(vector(v)+shift)+[h] for v in bounding_slice.vertices_list()]
        container = Polyhedron(Polyhedron(points).integral_points())
        # add pair (container, heighest point of B)
        boxes.append([container,list(p)+[h]])
    return boxes

def pealed_off_container(A,B,volumeBound,mixedVolumeBound,fixed_polytope):
    """
        We've fixed one polytopes (fixed_polytope) and we want to fix another one, let us call it Q,
        so that V(fixed_polytopes,Q,Q) <= mixedVolumeBound is fulfilled. We work under the assumption that
        Q occurs in the sandwich (A,B) and the the relative volume of Q is at most volumeBound.
        
        This function modifies B to a smaller polytope by removing all vertices v of B such that 
        conv(A \cup \{v\}) is too large in the sense that Q = conv(A \cup \{v\}) does not satisfy 
        one of the conditions above. The returned object is such a smaller polytope. 
    """
    to_peal_off=[]
    for v in B.vertices():
        blow_up_of_A=Polyhedron(list(A.vertices())+[v])
        if relVolume(blow_up_of_A)>volumeBound or mixed_volume([blow_up_of_A,blow_up_of_A,fixed_polytope])>mixedVolumeBound:
            to_peal_off.append(vector(v))
    Z=[vector(z) for z in B.integral_points()]
    return Polyhedron([z for z in Z if z not in to_peal_off])

def iterative_pealing_off(A,B,volumeBound,mixedVolumeBound,fixed_polytope):
    pealed_off_B=pealed_off_container(A,B,volumeBound,mixedVolumeBound,fixed_polytope)
    while pealed_off_B != B:
        B=pealed_off_B
        pealed_off_B=pealed_off_container(A,B,volumeBound,mixedVolumeBound,fixed_polytope)
    return B


# ==================================
# Sandwich-based enumeration
# ==================================

# General remarks:
#        In comparison to the enumeration of single full-dimensional polytopes by volume,
#        we have the following issue with the sandwich-based approach. Relative volume 
#        cannot directly be used as the gap function, as it is not strictly
#        decreasing. One can use number of lattice points, instead

#        In contrast to the sandwiches in the volume-classification, here we consider
#        the tuples A,B up to common translations but not up to general common
#        affine transformations. The reason is that our sandwich factory is used 
#        for enumeration of the second polytopes in the triple when the first polytope is already fixed.

def append_sandwich_up_to_translation(sf,A,B):
    """
        sf - sandwich factory to append the sandwich A,B to
        A,B - sandwich to append to sf
        
    """
    gap=B.integral_points_count()-A.integral_points_count()
    
    if gap not in sf.keys():
        sf[gap] = set([])
    
    sf[gap].add(common_translative_normal_form([A,B]))

def create_sandwich_factory_for_subpolytopes(initial_polytopes,container,m_v,m,fixed_polytope):
    """
        initial_polytopes - list of initial values for the inner part of the sandwich
        container - common outer part of the sandwich

        return a sandwich factory containing tuples A,B where A is in 
        initial_polytopes and B is the respective reduced version of the container
    """
    sf = {}
    for S in initial_polytopes:
        B=iterative_pealing_off(S,container,m_v,m,fixed_polytope)
        append_sandwich_up_to_translation(sf,S,B)
    return sf                


def subpolytopes_of_boxes_using_sandwiches(container,point,fixed_polytope,m,m_v,min_d):
    """
        We've fixed one polytope and we provide choices for the second one in a given container. 
        
        Input: 
        
            fixed_polytope - polytope that has been fixed
            container - bounding box in which to search for the second polytopes Q
            point - we search for Q with (0,0,0) and point in Q
            m - upper bound for V(Q,Q,fixed_polytope)
            m_v - upper bound for the relative volume of Q
            min_d - lower bound on dim(Q)

        return a list of subpolytopes Q of container satisfying
            V(Q,Q,fixed_polytope)<=m
            V(Q)<=m_v 
            dim(Q)>=min_d
       
        Note that the (container,point) part of the input is assumed to be taken from the list returned by create_boxes 
    """
    
    print "subpolytopes_of_boxes_using_sandwiches"
    
    # Q is required to contain this segment
    starting_segment = Polyhedron([[0,0,0],point])

    sf = create_sandwich_factory_for_subpolytopes([starting_segment],container,m_v,m,fixed_polytope)    
    
    maxGap=max(sf.keys())
    
    # iterative through all occurring dimension gaps downwards (that is, from maxDimGap to 0)
    while maxGap > 0:
        
        print "Number of lattice points gap", maxGap
        
        for sandwich in sf[maxGap]:
            A,B=sandwich
                
            for v in B.vertices(): # pick a vertex of B which is not in A
                if v not in A:
                    break
                
            blow_up_of_A=Polyhedron(list(A.vertices())+[v])
            reduction_of_B=Polyhedron([z for z in B.integral_points() if vector(z)!=vector(v)])
                
            pealed_B = iterative_pealing_off(blow_up_of_A,B,m_v,m,fixed_polytope)
                
            if pealed_B.dim() >= min_d:
                append_sandwich_up_to_translation(sf,blow_up_of_A,pealed_B)
                
            if reduction_of_B.dim() >= min_d:
                append_sandwich_up_to_translation(sf,A,reduction_of_B)
               
        del sf[maxGap]
            
        # because of the lower bound on the dimension of Q, it may happen that 
        # the returned list of this search is going to be empty. This means
        # that the gap 0 entry of our sandwich factory will not exist.
        # If this is the case, at some point, sf will be an empty dictionary here. 
        # In this case, we terminate
        if len(sf.keys())==0:
            return []
            
        # we can restart the iteration by recomputing the maximum gap.
        maxGap=max(sf.keys())
    
    # if we reach this point, we have sandwiches with zero gap. We return A=B of these sandwiches.
    return [s[0] for s in sf[0]]
    

def subpolytopes_of_boxes_using_sandwiches_2d(container,point,fixed_polytope,m,m_v):
    """
        We've fixed one polytope and we provide choices for a two-dimensional second one in a given container.
        
        Input: 
        
            fixed_polytope - polytope that has been fixed
            container - bounding box in which to search for the second polytopes Q with dim(Q)=2
            point - we search for Q with (0,0,0) and point in Q
            m - upper bound for V(Q,Q,fixed_polytope)
            m_v - upper bound for the relative volume of Q

        return a list of subpolytopes Q of container satisfying
            V(Q,Q,fixed_polytope)<=m
            V(Q)<=m_v 
            dim(Q)==2
            
       
        Note that the (container,point) part of the input is assumed to be taken from the list returned by create_boxes 
    """
    print "subpolytopes_of_boxes_using_sandwiches_2d"
    
    starting_segment = Polyhedron([[0,0,0],point])

    sf = create_sandwich_factory_for_subpolytopes([starting_segment],container,m_v,m,fixed_polytope)    
    
    maxGap=max(sf.keys())
    
        
    while maxGap>0:
        print "Number of lattice points gap",maxGap
        for sandwich in sf[maxGap]:
            A,B=sandwich
                
            for v in B.vertices(): # pick a vertex of B which is not in A
                if v not in A:
                    break
                
            blow_up_of_A=Polyhedron(list(A.vertices())+[v])
            reduction_of_B=Polyhedron([z for z in B.integral_points() if vector(z)!=vector(v)])
                
            # note that by construction blow_up_of_A will always have dimension at least 2
            if B.dim()>2:
                B = Polyhedron(B.intersection(affine_hull(blow_up_of_A)).integral_points())
            
            # form the reduction of B with respect to blow_up_of_A
            pealed_B = iterative_pealing_off(blow_up_of_A,B,m_v,m,fixed_polytope)
            # and add the corresponding sandwich
            append_sandwich_up_to_translation(sf,blow_up_of_A,pealed_B)
            
            # we are only interested in polytopes of dim=2 and therefore don´t consider a sandwich (A,B)
            # if dim(B)<2
            if reduction_of_B.dim()>1:
                append_sandwich_up_to_translation(sf,A,reduction_of_B)
            
        del sf[maxGap]
               
        # NOTE: it might be that we don't need this conditional termination step, but it does not harm.
        if len(sf.keys())==0:
            return []
            
        maxGap=max(sf.keys())

    for s in sf[0]:
        if s[0].dim()==1:
            print s[0].vertices()
            print s[1].vertices()
            
    # if we reach this point, we have sandwiches with zero gap. We return A=B of these sandwiches.
    return [s[0] for s in sf[0]]

def classify_two_2d(mv):
    """
        Returns the list of maximal irreducible triples with at least two two-dimensional polytopes and mixed volume mv.

    The search is based on Proposition 5.2, case (2)
    The conditions V(P_1,P_1,P_2)<=mv and dim(P_1)=2 imply that the relative volume of P_1 is at most mv
    NOTE: the notation in the paper and the choice of names in the code do not quite match.
    A in the code is P_1 in the paper
    B in the code is P_2 in the paper
    mv in the code is m in the paper

    """
    
    print "classify_two_2d based on Proposition 5.2, case (2)"
    
    pairs_and_cayleys = {'tuples':[], 'cayleys':set([])}
    
    
    # we iterate over possible relative volumes of A
    for m in range(1,mv+1):
        print "Choice of first polytope, which is two-dimensional and has relative volume %d" % m
        possible_A = lattice_polytopes_with_given_dimension_and_volume(2,m)
        
        for A in possible_A:
            Ae = embed_polygon(A)
            
            # We are choosing the second polytope B.
            # if V_2(A)=m and V(A,A,B)<=mv, then the formula V(A,A,B)= V_2(A) * (height of B in the normal direction of A)
            # implies that V(A,A,B) is a multiple of m, which explains the following iteration
            # (m1 is a possible choice of V(A,A,B))
            for m1 in range(m,mv+1,m):
                # recall that the boxes list contains not only boxes. Each box is supplied
                # with a point such that (0,0,0),point is assumed to be in the second polytope B
                # This ensures that V(A,A,B)=m1 (equality, and not just upper bound!)
                # mv**2 is the upper bound on V(B,B,A); see Proposition 5.2, case (2)
                boxes = create_boxes(A,m1,mv**2)
                for box in boxes:
                    print "First polytope (two-dimensional) has vertices:", Ae.vertices()
                    print "Searching the second polytope inside the container with vertices:", box[0].vertices()
                    # box[0] is the container
                    # box[1] is a point at the top of the container to be contained in B.
                    # We have assumed V(B,B,A)<=mv**2 (so, the last but one parameter of subpolytopes_of_boxes_using_sandwiches_2d is chosen to be mv**2)
                    # We have V_2(B)<=V(B,B,A)<=mv**2 (so, the last parameter of subpolytopes_of_boxes_using_sandwiches_2d is chosen to be mv**2, too).
                    # The condition V_2(B)<=mv**2 is redundant for these concrete bounds.
                    for B in subpolytopes_of_boxes_using_sandwiches_2d(box[0],box[1],Ae,mv**2,mv**2):
                        insert_up_to_equivalence(pairs_and_cayleys,[Ae,B])
        
    triples = classify_maximal_triples_from_pairs(pairs_and_cayleys['tuples'],mv,2)
    return triples


def classify_one_2d(mv,only_b=False):
    """
        The search is based on Proposition 5.2, case (1)
        
        If only_b is set to true, only subcase 1b (the non-inductive subcase) is carried out
    """
    # NOTE: this function is contained of two disjoint codes for subcases (a) and (b)

    print "classify_one_2d, based on Proposition 5.2, case (1)"
    
    # ------------------ non-inductive part --------------------------------------------------
        
    print "subcase (b)"
    pairs_and_cayleys = {'tuples':[], 'cayleys':set([])}
    
    # we first choose the relative volume for the two-dimensional one (called A in the code and P_1 in the paper)
    for m in range(1,mv+1):
        # the parameter 2 is the dimension and the parameter m is the normalized volume
        possible_A = lattice_polytopes_with_given_dimension_and_volume(2,m)
        
        for A in possible_A:
            Ae = embed_polygon(A)
            # To enumerate B, we first choose a possible values m1 of V(A,A,B)
            # here m1 is a MULTIPLE of m satisfying m<=m1<=mv
            for m1 in range(m,mv+1,m):
                boxes = create_boxes(A,m1,mv)
                for box in boxes:
                    print "First polytope (two-dimensional) has vertices:", Ae.vertices()
                    print "Searching the second polytope inside the container with vertices:", box[0].vertices()
                
                    for B in subpolytopes_of_boxes_using_sandwiches(box[0],box[1],Ae,mv,mv**2,3):
                        if mixed_volume([B,B,Ae])==mv:
                            insert_up_to_equivalence(pairs_and_cayleys,[Ae,B])
    
    # if only_b is True, already form maximal triples here and return
    if only_b:
        triples = classify_maximal_triples_from_pairs(pairs_and_cayleys['tuples'],mv,2)
        
        return [t for t in triples if set([2,3])<= set([t[0].dim(), t[1].dim(), t[2].dim()])]
    
    # ------------------ inductive part ------------------------------------------------------
    # case 1: there is an irreducible triple (Pi,Pi,Pj) such that V(Pi,Pi,Pj)<m
    
    print "subcase (a)"

    pairs_and_cayleys_ind1 = {'tuples':[], 'cayleys':set([])}

    fname_full_d = 'data/dim_3_mv_%d.txt'
    fname_lower_d = 'data/dim_3_mv_%d_lower_d.txt'

    saturated_triples = []
    for m in range(1,mv):
        saturated_triples = saturated_triples + polytope_tuples_from_file(fname_full_d % m) + polytope_tuples_from_file(fname_lower_d % m)

    # the 1 is the maximal number of two-dimensional ones occurring in the triple, as we are interested in the triples where a three-dimensional
    # is repeated at least twice.
    triples = irreducible_subtriples_up_to_equivalence(saturated_triples,1)

    for t in triples:
        # if a three-dimensional is repeated twice, so we have the case (P,P,Q) with dim(P)=3 we return (P,Q)
        p = get_special_pair(t)
        if p != None:
            insert_up_to_equivalence(pairs_and_cayleys_ind1,p)

    # the results from both cases are united
    triples = classify_maximal_triples_from_pairs(pairs_and_cayleys['tuples']+pairs_and_cayleys_ind1['tuples'],mv,2)

    return [t for t in triples if set([2,3])<= set([t[0].dim(), t[1].dim(), t[2].dim()])]
 

 
