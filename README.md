# mixed_volume_classification
Sage code implementing the algorithms described in the paper [Classification of triples of lattice polytopes with a given mixed volume](https://arxiv.org/abs/1902.00891) as well as corresponding enumeration data. The folder python3 contains modified versions of the scripts, which are compatible with versions of SageMath running on top of Python3 (starting from version 9.0). This is joint work with Gennadiy Averkov and Ivan Soprunov.

## Content:

### Code:

* **polytopes.sage, polytope_tuples.sage**:   basic operations on single lattice polytopes and tuples of lattice polytopes respectively
* **full_dimensional_mvol_inductive.sage**:   classification of maximal full-dimensional triples of lattice polytopes of given 
                                        mixed volume in dimension 3 (following Algorithm 5.4 of the paper)
* **lower_dimensional_mvol_inductive.sage**:  classification of maximal irreducible triples of lattice polytopes of given mixed volume
                                        in dimension 3 (with at least one polytope not being full-dimensional, 
                                        following Algorithm 6.3 of the paper)
* **full_dimensional_mvol_polygons.sage**:    classification of maximal pairs of full-dimensional polygons of given mixed volume
* **volume_classification.sage**:             classification of lattice polytopes of given volume in dimensions 2 and 3

### Data:

This folder includes lists of maximal irreducible triples of mixed volume 1-4 as well as lists of pairs of full-dimensional
polygons of mixed volume 1-10. Where it makes sense, the files containing the triples are split on the one hand by whether 
all of the polytopes in the triple are full-dimensional or not and on the other hand by whether the triple has the additional
property of being R-maximal over the reals or not (see Definition 3.1) of the paper.
Both triples and tuples are decoded as lists of the vertices of the respective polytopes and can e.g. be read into Sagemath
using the function *polytope_tuples_from_file* from polytope_tuples.sage.

* **dim_3_mv_[1-2].txt**: maximal triples of full-dimensional lattice polytopes of mixed volume [1-2] in dimension 3
* **dim_3_mv_[3-4]_rmax.txt**: maximal triples of full-dimensional lattice polytopes of mixed volume [3-4] in dimension 3 that are R-maximal
* **dim_3_mv_[3-4]_zmax.txt**: maximal triples of full-dimensional lattice polytopes of mixed volume [3-4] in dimension 3 that are
                               not R-maximal
* **dim_3_mv_[2-4]_lower_d_rmax.txt**: maximal irreducible triples of lattice polytopes of mixed volume [2-4] in dimension 3 that
                                     contain at least one 2-dimensional polytope and are R-maximal
* **dim_3_mv_[2-4]_lower_d_zmax.txt**: maximal irreducible triples of lattice polytopes of mixed volume [2-4] in dimension 3 that
                                     contain at least one 2-dimensional polytope and are not R-maximal
* **dim_2_mv_[1-10].txt**: maximal pairs of full-dimensional polygons with mixed volume [1-10]
