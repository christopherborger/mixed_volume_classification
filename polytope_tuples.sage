# ========================================================================
# Functions on tuples of polytopes 
# ========================================================================

load("polytopes.sage")

def get_pair(t):
    """
        t is a triple (A,B,C).
        If some two polytopes in the triple coincnide up to translation, we return the pair discarding the one which is repeated.
        Otherwise, we return None (nothing).
    """
    A = translative_normal_form(t[0])
    B = translative_normal_form(t[1])
    C = translative_normal_form(t[2])

    if A==B:
        return [A,C]
    if A==C:
        return [A,B]
    if B==C:
        return [A,B]

def get_special_pair(t):
    """
        If we have a triple with a repeated three-dimensional polytope, like (A,A,B), 
        we return the pair (A,B)
    """
    A = translative_normal_form(t[0])
    B = translative_normal_form(t[1])
    C = translative_normal_form(t[2])
    if A.dim()==3:
        if A==B:
            return [A,C]
        if A==C:
            return [A,B]
    if B.dim()==3:
        if B==C:
            return [B,A]

def common_translative_normal_form(p):
    """
        p is a pair of polytopes.
        This is used to check if two pairs coincide up to a common translation.
    """
    A_new = p[0] - vector(min(p[0].vertices()))
    B_new = p[1] - vector(min(p[0].vertices()))
    return tuple([A_new,B_new])

def mixed_volume(L):
    """
        mixed volume of a tuple
    """
    assert len(L) in [2,3]
    if len(L)==2:
        return mixed_volume_2(L[0],L[1])
    return mixed_volume_3(L[0],L[1],L[2])

def mixed_volume_3(P1,P2,P3):
    return (P1+P2+P3).volume()-(P1+P2).volume()-(P2+P3).volume()-(P1+P3).volume()+P1.volume()+P2.volume()+P3.volume()


def mixed_volume_2(P1,P2):
    return (P1+P2).volume()-P1.volume()-P2.volume()

def Cayley(P):
    """
        Cayley polytope of a list of polytopes.
        This function should replace cayley_polytope
    """
    k=len(P)
    A=map(tuple,[zero_vector(QQ,k-1)]+identity_matrix(k-1).rows())
    V=[]
    for i in xrange(k):
        for v in P[i].vertices():
            V.append(A[i]+tuple(v))
    return Polyhedron(V)

def cayley_normal_form(Tuple):
    """
        This is used to test if two tuples are equivalent. 
    """
    return affine_normal_form(Cayley([2 * P for P in Tuple]))

def is_irreducible(t):
    return t[0].dim()>=2 and t[1].dim()>=2 and t[2].dim()>=2 and (t[0]+t[1]).dim()==3 and (t[0]+t[2]).dim()==3 and (t[0]+t[1]).dim()==3

def is_maximal_in(t,i):
    # TODO: that's an old implementation, and we probably can remove it. Check this!
    j,k=list(set([0,1,2])-set([i]))
    mv = mixed_volume_3(t[0],t[1],t[2])
    Pi = translative_normal_form(t[i])
    Pj = translative_normal_form(t[j])
    Pk = translative_normal_form(t[k])
    return (Pi in maximal_third_polytopes(Pj,Pk,mv,2))

def total_number_of_points_in_triple(t):
    return len(t[0].integral_points())+len(t[1].integral_points())+len(t[2].integral_points())

def maximal_third_polytopes(A,B,mv,d_min):
    """
        Input: (A,B) irreducible pair (meaning dim(A)>=2, dim(B)>=2 and dim(A+B)>=3). A and B are first and second polytope, respectively.
        Output: all maximal C with V(A,B,C)=mv and dim(C)>=d_min. The list of C is the list of third polytopes.
    """
    measure=mixed_area_measure(A,B)
    bounding_polytopes=set([])
    # we assume zero vector is in C, wlog
    # we need to prescribe h(C,u) for u in measure.keys() (measure.keys() are primitive vectors in the support of the measure)
    # such that V(A,B,C)=mv. measure.values() are weights of the primitive vectors.
    for supFuncValuesVector in myWeightedIntegerVectors(mv, measure.values()):
        # inequalities are given by a coeff vector with coeff[0]  + coeff[1]* x[0] + ... + coeff[n] * x[n-1] >= 0 
        # so we convert from lhs*x <= rhs to rhs - lhs * x >=0.
        Cinequalities= [tuple([rhs]) + tuple(-vector(lhs))  for lhs, rhs in zip(measure.keys(),supFuncValuesVector)]
        Ccontainer=Polyhedron(ieqs=Cinequalities)
        # NOTE: since (A,B) is irreducible, Container is always compact (see Theorem~?? TODO: add theorem number)
        # so, we could actually omit the test Container.is_compact()
        if Ccontainer.is_compact():
            C = Polyhedron(Ccontainer.integral_points())
            if C.dim()>=d_min and mixed_volume([A,B,C])==mv:
                bounding_polytopes.add(translative_normal_form(C))
    return bounding_polytopes

def classify_maximal_triples_from_pairs(ab_pairs,mv,d_min):
    """
        returns a list of all maximal triples A,B,C with
        V(A,B,C)=mv and such that dim(C)>=d_min from a list of pairs A,B
    """
    print "classify_maximal_triples_from_pairs(ab_pairs,mv,d_min)"

    triples_and_cayleys = {'tuples':[], 'cayleys':set([])}

    cntr = 1

    for p in ab_pairs:
        print "Pair number", cntr
        cntr = cntr+1
        A = p[0]
        # NOTE: one would not really have to translate p[1], but we do this to make the output of the enumeration a bit nicer
        B = translative_normal_form(p[1])

        for C in maximal_third_polytopes(A,B,mv,d_min):
            # making translative_normal_form(C) is again for the puprose of making the output nicer.
            # we could have used just C.
            insert_up_to_equivalence(triples_and_cayleys,[A,B,translative_normal_form(C)])

    maximal_triples = []
    
    cntr = 1
    for t in triples_and_cayleys['tuples']:
        print "Checking triple ",cntr
        cntr += 1
        if is_maximal_triple(t):
            print "Found a new one"
            maximal_triples.append(t)

    return maximal_triples

def insert_up_to_equivalence(tuples_and_cayleys,new_tuple):
    """
        tuples_and_cayleys - dictionary containing a list
                             containing tuples of lattice polytope and a set 
                             containing affine normal forms of the cayley polytopes of the 
                             respective tuple
        new_tuple - the tuple to be inserted if there is not already an equivalent tuple
    """
    
    # construct normal of the tuple using Cayley polytopes
    anf_cayley = cayley_normal_form(new_tuple)

    # if tuple is not contained, yet..
    if not anf_cayley in tuples_and_cayleys['cayleys']:
        # ..add it
        tuples_and_cayleys['tuples'].append(new_tuple)
        tuples_and_cayleys['cayleys'].add(anf_cayley)
        
        
# ------------------------------------------


def maximal_second_polygons(A,mv,d_min):
    """
        compute all maximal B with mixed volume of (A,B) equal to mv
	of given minimal dimension d_min
    """
    measure=area_measure(A)
    bounding_polygons=set([])
    # we assume that the zero vector is in B, wlog
    # we need to prescribe h(B,u) for u in measure.keys() such that the mixed volume of (A,B) is equal to mv
    for supFuncValuesVector in myWeightedIntegerVectors(mv, measure.values()):
        Binequalities= [tuple([rhs]) + tuple(-vector(lhs))  for lhs, rhs in zip(measure.keys(),supFuncValuesVector)]
        Bcontainer=Polyhedron(ieqs=Binequalities)
        B = Polyhedron(Bcontainer.integral_points())
        if B.dim()>=d_min and mixed_volume_2(A,B)==mv:
            bounding_polygons.add(translative_normal_form(B))
    return bounding_polygons



def is_maximal_pair(p):
    """
        p - a pair of lattice polygons
        
        determine wether a given pair of lattice polygons is saturated,
        that is adding any point to any of the lattice polytopes 
        increases the mixed volume of the pair
    """
    A=p[0]
    B=p[1]
    mA=area_measure(A)
    mB=area_measure(B)
    Bplus=Polyhedron(ieqs=[ (support_function(B,u),)+ tuple(-vector(u)) for u in supp(mA)])
    Aplus=Polyhedron(ieqs=[ (support_function(A,u),)+ tuple(-vector(u)) for u in supp(mB)])
    for z in Bplus.integral_points():
        if z not in B:
            return False
    for z in Aplus.integral_points():
        if z not in A:
            return False
    return True

def is_maximal_triple(t):
    A=t[0]
    B=t[1]
    C=t[2]
    mAB=mixed_area_measure(A,B)
    mAC=mixed_area_measure(A,C)
    mBC=mixed_area_measure(B,C)
    Aplus=Polyhedron(ieqs=[ (support_function(A,u),)+ tuple(-vector(u)) for u in supp(mBC)])
    Bplus=Polyhedron(ieqs=[ (support_function(B,u),)+ tuple(-vector(u)) for u in supp(mAC)])
    Cplus=Polyhedron(ieqs=[ (support_function(C,u),)+ tuple(-vector(u)) for u in supp(mAB)])
    if not Aplus.is_compact() or not Bplus.is_compact() or not Cplus.is_compact():
        return False
    for z in Aplus.integral_points():
        if z not in A:
            return False
    for z in Bplus.integral_points():
        if z not in B:
            return False
    for z in Cplus.integral_points():
        if z not in C:
            return False
    return True

# -------------------------------------------

def mixed_area_measure(A,B):
    mAB=area_measure(A+B)
    mA=area_measure(A)
    mB=area_measure(B)
    measure=defaultdict(int)
    supp_container=set(mAB.keys()) | set(mA.keys()) | set(mB.keys())
    for u in supp_container:
        m=(mAB[u]-mA[u]-mB[u])/2
        if m != 0:
            measure[u]=m
    return measure

def R_saturated(A,B,C):
    s_A=supp(area_measure(A))
    s_B=supp(area_measure(B))
    s_C=supp(area_measure(C))
    s_AB=supp(mixed_area_measure(A,B))
    s_AC=supp(mixed_area_measure(A,C))
    s_BC=supp(mixed_area_measure(B,C))
    return (s_A<=s_BC) and (s_B<=s_AC) and (s_C<=s_AB)

# -------------------------------------------

# ==========================================
# Input/Output for polytope tuples
# ==========================================

def save_polytope_tuples(tuples,fname):
    f=open(fname,'w')
    print >>f, [[map(tuple,P.vertices()) for P in t] for t in tuples]
    f.close()

def polytope_tuples_from_file(fname):
    f=open(fname,'r')
    L=eval(f.read().replace('\n',' '))
    L=[map(Polyhedron,t) for t in L]
    return L

def recursive_map(f, it):
    return [recursive_map(f, x) if isinstance(x[0], list) else f(x) for x in it]

def minkowski_decompositions_from_file(fname):
    f=open(fname,'r')
    string=f.read()
    L=eval(string.replace('\n','').replace('\\',''))
    L=recursive_map(Polyhedron,L)
    return L
# ==========================================
# basic functions not directly related to polytopes or tuples of polytopes
# ==========================================

def myWeightedIntegerVectors(w,a):
    """
        NOTE: previously, we expirienced buggs in WeightedIntegerVectors. For example, the result for
         (6,[1,1,4,4]) contained vectors with non-integer components. Now, it seems to be okay, but we kept this version,
         which reimplements the function WeightedIntegerVectors.
        
        Interface: w - total weight of all items, a[i] is the weight of i-th item. 
        Output: yields all non-negative integer solutions x to the equation a[0]*x[0] + ... + a[n-1]*x[n-1] = w
    """
    n=len(a)
    # the set of all x with x>=0 and a[0]*x[0] + ... + a[n-1]*x[n-1] is an (n-1)-dimensional simplex P with vertices e[0]*w/a[0],...,e[n-1]*w/a[n-1] 
    # where e[0],...e[n-1] is a standard basis
    P=Polyhedron(vertices=[ w*e/Factor for e,Factor in zip(matrix.identity(n).columns(),a)])
    # we return the list of integral points of this simplex
    return P.integral_points()
