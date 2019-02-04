from collections import defaultdict

# ========================================================================
# Functions on polytopes which are not directly
# available as methods in the Polyhedron class
# ========================================================================

def embed_polygon(P):
    """
        Embed a polygon into height 0 in dimension three. 
        NOTE: this can be used for polytopes, too, but we use it only for polygons.
    """
    return Polyhedron([list(v) + [0] for v in P.vertices()])

def affine_hull(A):
    """
        affine hull of a polyhedron, represented as a Polyhedron object. 
    """
    return Polyhedron(ieqs=[vector(e) for e in A.equations()]+[-vector(e) for e in A.equations()])


def ineq_face(P,ineq):
    V=[v for v in P.vertices() if ineq.eval(v)==0]
    return Polyhedron(vertices=V)

def support_function(P,u):
    u=vector(u)
    return max([u.dot_product(v) for v in map(vector,P.vertices())])

def facets_and_outer_normals(P):
    """
        generator of pairs facet,outer_normal_of_facet
    """ 
    for ineq in P.Hrepresentation():
        F=ineq_face(P,ineq)
        normalF=-vector(ineq[1:])
        yield F,normalF
        
def interior_integral_points(P):
    return [z for z in P.integral_points() if P.interior_contains(z)]

def as_full_dim_polyhedron(P):
    P = P.lattice_polytope()
    
    L = P.lattice()
    V = P.vertices()
    
    if L.zero() not in V:
        V = [v-V[0] for v in V]
    
    S = L.span(V).saturation()
    
    return Polyhedron([S.coordinates(v) for v in V])


def Volume(P):
    """
        normalized volume (in the sense of the lattice polytope theory)
    """
    return factorial(P.dimension())*P.volume()

def relVolume(P):
    """
        relative normalized volume
    """
    return Volume(as_full_dim_polyhedron(P))

def affine_normal_form(P):
    V=[vector(v) for v in P.vertices()]
    s=sum(V)
    N=len(V)
    Q=N*P-s
    try:
        # in case pulp (hidden behind normal_form()) has some trouble to handle Q, we'll catch the exception and see what P was
        Q=Polyhedron(Q.lattice_polytope().normal_form())
    except Exception as e: 
        print "ERROR MESSAGE:"
        print "affine_normal_form(P) failed on the polytope P with the vertices"
        print map(tuple,P.vertices())
        print "Here is the exception:"
        print e
        print "END OF THE ERROR MESSAGE"
    Q=Q-vector(min(Q.vertices()))
    return Q/N

def translative_normal_form(P):
    """
        Return P minus the lexicographically minimal vertex of P.
        Properties of the returned polytope: 
            if P and Q coincide up to translations, then their translative normal form is the same.
            the zero is always a vertex of the returned polytope.
        
    """
    return P-vector(min(P.vertices()))

def are_homothetic(A,B):
    """
        Testing if A and B are homothetic by bringing the vertex sets of A and B to a kind of homothetic normal form.
    """
    A = translative_normal_form(A)
    B = translative_normal_form(B)
    
    A_factor = gcd([x for v in A.vertices() for x in list(v)])
    B_factor = gcd([x for v in B.vertices() for x in list(v)])

    return A/A_factor == B/B_factor


# =============================================
# Generating special polytopes
# =============================================

def std_simplex(d):
    return Polyhedron([d*(0,)]+identity_matrix(d).rows())

def empty_simplices(d,m):
    """
        returns a list of all d-dimensional empty simplices of volume at most m
        
        d must be 1,2 or 3 and m at least 1
    """
    # TODO: this function could have been implemented as a generator.
    assert d<=3
    assert m>=1
    emptySimplices=[std_simplex(d)]
    if d==3:
        for y in range(2,m+1):
            for x in range(y):
                if GCD(x,y)==1:
                    emptySimplices.append(Polyhedron([(0,0,0),(1,0,0),(0,0,1),(x,y,1)]))
    return emptySimplices

# ====================
# area measures
# ====================

def area_measure(P):
    assert P.ambient_dim()==P.dim() or (P.ambient_dim()==3 and P.dim()==2)
    if P.ambient_dim()==P.dim():
        return area_measure_full_dim(P)
    measure=defaultdict(int)
    for eq in P.equations():
        direction1=tuple(-vector(eq[1:]))
        direction2=tuple(vector(eq[1:]))
        area_of_P=relVolume(P)
        measure[direction1]=area_of_P
        measure[direction2]=area_of_P
    return measure

def area_measure_full_dim(P):
    assert P.ambient_dim()==P.dim()
    measure=defaultdict(int)
    for ineq in P.Hrepresentation():
        F=ineq_face(P,ineq)
        direction=tuple(-vector(ineq[1:]))
        measure[direction]=relVolume(F)
    return measure

def supp(measure):
    return set(measure.keys())


# ==========================================
# Input/Output for polytopes
# ==========================================

def save_polytopes(polytopes,fname):
    f=open(fname,'w')
    print >>f, [map(tuple,P.vertices()) for P in polytopes]
    f.close()

def polytopes_from_file(fname):
    f=open(fname,'r')
    L=eval(f.read().replace('\n',' '))
    L=map(Polyhedron,L)
    f.close()
    return L









