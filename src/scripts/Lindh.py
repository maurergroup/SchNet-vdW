#!/usr/bin/python

"""
This script was provided by Dimitrii Maksimov.
Generate model Hessian matrix for a given geometry.  See
R. Lindh et al., Chem. Phys. Lett. 241, 423 (1995).

Usage examples:
  # Simply add Hessian to existing geometry.in:
  $ Lindh.py --append geometry.in
  # Convert an existing .xyz file:
  $ Lindh.py file.xyz geometry.in
  # Have a look at all the internal coordinates under the hood:
  $ Lindh.py -v file.xyz | less
  # Use --cutoff 5. if in hurry."""

import collections
import operator
import subprocess
import sys

VERSION_ERROR = """Incompatible Python version.

This script is supposed to be run with Python 2.x.
You might want to try "python2 Lindh.py" instead of "./Lindh.py"
to run it.
"""
if sys.version_info[0] != 2:
    sys.stderr.write(VERSION_ERROR)
    sys.exit(2)

import numpy as np
from numpy.linalg import norm

USAGE = "%prog [options] [infile [outfile]]\n" + __doc__

HUGE = 1e10

ABOHR = 0.52917721 # in AA
HARTREE = 27.211383 # in eV

K_BOND    = 0.450 * HARTREE / ABOHR**2
K_BENDING = 0.150 * HARTREE
K_TORSION = 0.005 * HARTREE

ALPHAS = np.array([[1.0000, 0.3949, 0.3949],
                   [0.3949, 0.2800, 0.2800],
                   [0.3949, 0.2800, 0.2800]]) * ABOHR**(-2)
REF_DS = np.array([[1.35, 2.10, 2.53],
                   [2.10, 2.87, 3.40],
                   [2.53, 3.40, 3.40]]) * ABOHR

COVRADS = dict(
    H=0.320, He=0.310,
    Li=1.630, Be=0.900, B=0.820, C=0.770,
    N=0.750, O=0.730, F=0.720, Ne=0.710,
    Na=1.540, Mg=1.360, Al=1.180, Si=1.110,
    P=1.060, S=1.020, Cl=0.990, Ar=0.980,
    K=2.030, Ca=1.740,
    Sc=1.440, Ti=1.320, V=1.220, Cr=1.180, Mn=1.170,
    Fe=1.170, Co=1.160, Ni=1.150, Cu=1.170, Zn=1.250,
    Ga=1.260, Ge=1.220, As=1.200, Se=1.160, Br=1.140, Kr=1.120,
    Rb=2.160, Sr=1.910,
    Y=1.620, Zr=1.450, Nb=1.340, Mo=1.300, Tc=1.270,
    Ru=1.250, Rh=1.250, Pd=1.280, Ag=1.340, Cd=1.480,
    In=1.440, Sn=1.410, Sb=1.400, Te=1.360, I=1.330, Xe=1.310,
    Cs=2.350, Ba=1.980,
    La=1.690, Ce=1.650, Pr=1.650, Nd=1.840, Pm=1.630, Sm=1.620,
    Eu=1.850, Gd=1.610, Tb=1.590, Dy=1.590, Ho=1.580, Er=1.570,
    Tm=1.560, Yb=2.000, Lu=1.560, Hf=1.440, Ta=1.340, W=1.300, Re=1.280,
    Os=1.260, Ir=1.270, Pt=1.300, Au=1.340, Hg=1.490, Tl=1.480, Pb=1.470,
    Bi=1.460, Po=1.460, At=2.000, Rn=2.000, Fr=2.000, Ra=2.000, Ac=2.000,
    Th=1.650, Pa=2.000, U=1.420)

ID = np.identity(3)
ZERO = np.zeros((3, 3))


def _acc_dict(key, d, val):
    """If key in dict, accumulate, otherwise set."""
    if key not in d:
        d[key] = np.zeros_like(val)
    d[key] += val


def isposvec(vec, eps=1e-10, noself=True):
    """Return if vector is in some special sense 'positive'.

    Positiveness is defined by the rightmost non-zero dimension:
    >>> isposvec(np.array([1., 0., -1.]))
    False
    >>> isposvec(np.array([1., 0., 0.]))
    True
    >>> isposvec(np.array([-1., 0., 1.]))
    True
    """
    for x in reversed(vec):
        if x > eps:
            return True
        elif x < -eps:
            return False
    if noself:
        raise ValueError("Self-referencer")
    else:
        return None
    


def canonize(atom_name):
    """Return canonical name of atom_name.

    The canonical name is the first capital of the string with an optional
    minuskel.

    Example:
    >>> print canonize("Ru"), canonize("-H3"), canonize("CT")
    Ru H C
    """
    name = atom_name
    while name and not name[0].isupper():
        name = name[1:]
    if len(name) > 1 and name[1].islower():
        name = name[0:2]
    else:
        name = name[0:1]
    return name

def name2row(atom_name):
    """Return row number of atom type (starting with 0, max 2)."""
    name = canonize(atom_name)
    if name in ('H', 'He'):
        return 0
    elif name in ('Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'):
        return 1
    else:
        return 2

class Damper(object):
    """Damper interface documentation

    A Damper is an object that can judge the importance of an atom pair from
    the interatomic distance and the atom numbers.
    """
    def exponent(self, AB, i, j):
        """Return exponent for distance AB and atom types i, j.

        What is actually done with AB, and in particular i, j
        is an implementation detail.  Most probably, the initializer
        should be passed some system related information, i.e. at
        least the atom types within the system.
        """
        raise NotImplementedError()
    def Rcut(self, max_exponent):
        """Return the maximum distance leading to max_exponent."""

class SimpleDamper(Damper):
    """Damper for maximum chain lenght (exponent==bondlength)."""
    def __init__(self, atom2any=None):
        self.atom2any = atom2any
    def exponent(self, AB, i, j):
        return norm(AB)
    def Rcut(self, max_exponent):
        """Return the maximum distance leading to max_exponent."""
        return max_exponent

class CovalenceDamper(Damper):
    """Damper class for covalent bonds (exponent in set([0, HUGE]))."""
    def __init__(self, atom2name, covrads=COVRADS):
        self.covrads = covrads
        self.atom2name = atom2name
    def exponent(self, AB, i, j):
        cr_i = self.covrads[canonize(self.atom2name[i])]
        cr_j = self.covrads[canonize(self.atom2name[j])]
        if norm(AB) < 1.3 * (cr_i + cr_j):
            return 0.
        else:
            return HUGE
    def Rcut(self, max_exponent):
        """Return the maximum distance leading to max_exponent."""
        return max(self.covrads[name] for name in self.atom2name.values())


class LindhExponent(Damper):
    """Class of the LINDH object which provides the exponent factors."""
    def __init__(self, atom2row, alphas=ALPHAS, ref_ds=REF_DS):
        self.alphas = alphas
        self.ref_ds = ref_ds
        self.atom2row = atom2row

    def exponent(self, AB, i_atom, j_atom):
        """Return the exponent for distance AB of given types."""
        i_row, j_row = self.atom2row[i_atom], self.atom2row[j_atom]
        alpha = self.alphas[i_row, j_row]
        ref_d = self.ref_ds[i_row, j_row]
        return alpha * (AB**2 - ref_d**2)

    def Rcut(self, max_exponent):
        """Return the maximum distance for given exponent."""
        lr_alpha = np.min(self.alphas)
        lr_ref_d = np.max(self.ref_ds)
        # max_exponent == lr_alpha * (Rcut**2 - lr_ref_d**2)
        # max_exponent / lr_alpha + ref_d**2 == Rcut**2
        Rcut = np.sqrt(max_exponent / lr_alpha + lr_ref_d**2)
        return Rcut


class Bravais(object):
    """Provide tools related to some given Bravais lattice.

    May be initialized by a list of one, two, or three Bravais vectors.
    Provides tools to fold a vector or several vectors into the central
    parallel epipede (into_pe) and to retrieve a list of lattice vectors
    within a given radius (all_within).
    """
    def __init__(self, lattice_vectors):
        """Initializes Bravais object."""
        if lattice_vectors is None:
            lattice_vectors = []
        else:
            lattice_vectors = list(lattice_vectors)
        self.n = len(lattice_vectors)
        if self.n > 0:
            self.bra = np.array(lattice_vectors)
            self.ibra = np.linalg.pinv(self.bra)
            self.rec = 2.*np.pi*np.transpose(self.ibra)
        else:
            self.bra = np.empty((0, 3))
            self.rec = np.empty((0, 3))

    def latvec(self, abc):
        """Return a lattice vector from integer Bravais indices."""
        vec = np.zeros(3)
        a, b, c = abc
        if a != 0: vec += a*self.bra[0,:]
        if b != 0: vec += b*self.bra[1,:]
        if c != 0: vec += c*self.bra[2,:]
        return vec

    def into_pe(self, vecs):
        """Fold vectors (last dimension 3) into parallel epipede.

        Examples:
        >>> lat = Bravais([[1,0,0], [0,2,0]])
        >>> np.allclose(lat.into_pe([3.2, 1.3, 0.5]), [0.2, -0.7, 0.5])
        True
        """
        vecs = np.asarray(vecs, dtype=float)
        shape = vecs.shape
        if shape[-1] != 3: raise ValueError("Last dim should be 3.")
        n_vec = np.product(shape[:-1])
        reslist = []
        for vec in vecs.reshape((n_vec, 3)):
            rcoeff = np.dot(vec, self.ibra)
            icoeff = np.around(rcoeff)
            res = vec - np.dot(icoeff, self.bra)
            reslist.append(res)
        return np.array(reslist).reshape(shape)

    def all_within(self, Rcut, add_base_PE=False):
        r"""Return a list of all lattice vector indices shorter than Rcut.

        If base_PE is True, add one parallel epipede (PE) to the region.

        Examples:
        >>> cos60, sin60 = np.cos(np.pi/3), np.sin(np.pi/3)
        >>> lat = Bravais([[1,0,0], [cos60, sin60, 0], [0,0,1.5]])
        >>> lat.all_within(0.5) == [(0, 0, 0)]
        True
        >>> (1, 0, 0) in lat.all_within(1.1)
        True
        >>> set1 = lat.all_within(1.2)
        >>> len(set1)
        7
        >>> len(lat.all_within(0.2, add_base_PE=True))
        27
        
        The resulting vectors are sorted:
        >>> lat.all_within(2.2)[0] == (0, 0, 0)
        True
        """
        len_rec = np.array([norm(v) for v in self.rec])
        ns = np.floor(len_rec * Rcut / (2*np.pi)) # cross check with diss
        n_cells = np.zeros(3, int)
        n_cells[:self.n] = ns
        abcs = set()
        for a in xrange(-n_cells[0], n_cells[0]+1):
            for b in xrange(-n_cells[1], n_cells[1]+1):
                for c in xrange(-n_cells[2], n_cells[2]+1):
                    vec = self.latvec((a, b, c))
                    if norm(vec) <= Rcut:
                        abcs.add((a, b, c))
        if add_base_PE:
            # add one in each direction
            def _around(d, n): return [-1, 0, 1] if d < n else [0]
            old_abcs = set(abcs)
            for i, j, k in old_abcs:
                for ii in _around(0, self.n):
                    for jj in _around(1, self.n):
                        for kk in _around(2, self.n):
                            abcs.add((i+ii, j+jj, k+kk))

        def _norm_of_abc(abc): return norm(self.latvec(abc))
        return sorted(abcs, key=_norm_of_abc)

def get_pairs(atoms1, atoms2, Rcut, use_scipy=True):
    if use_scipy:
        try:
            import scipy.spatial  # KDTree
            have_scipy = True
        except ImportError:
            have_scipy = False
    else:
        have_scipy = False

    if have_scipy:
        tree1 = scipy.spatial.KDTree(atoms1)
        tree2 = scipy.spatial.KDTree(atoms2)
        pairs = tree1.query_ball_tree(tree2, Rcut)
    else:
        sys.stderr.write("No scipy found; using fallback.\n")
        pairs = []
        for i, atom1 in enumerate(atoms1):
            this_pairs = []
            for j, atom2 in enumerate(atoms2):
                dist = norm(atom2 - atom1)
                if dist < Rcut:
                    this_pairs.append(j)
            pairs.append(this_pairs)
        sys.stderr.write("Searching done.\n")
    return pairs



class Pairs(object):
    """Find chains (pairs, triples, ...) of atoms.

    Example:
    # Chain of pairs, one at zero, one slightly distorted.
    >>> bra = [[1, 0, 0]]
    >>> atom = [[0, 0, 0], [0.01, 0.25, 0]]
    >>> pairs = Pairs(bra, atom, SimpleDamper(), 1.5)
    >>> bonds = list(pairs.chains(2, 1.01))   # all bonds up to length 1.
    >>> for damp, atlist in bonds:  # only intracell and its own images.
    ...     assert len(atlist) == 2
    ...     print "%6.4f %s" % (damp, atlist)
    0.2502 [(0, (0, 0, 0)), (1, (0, 0, 0))]
    1.0000 [(0, (0, 0, 0)), (0, (1, 0, 0))]
    1.0000 [(1, (0, 0, 0)), (1, (1, 0, 0))]
    >>> bendings = list(pairs.chains(3, 1.251))   # 1.251: one short 1->2 bond.
    >>> for damp, atlist in bendings:
    ...     assert len(atlist) == 3
    ...     print "%6.4f %s" % (damp, atlist)
    1.2502 [(0, (0, 0, 0)), (0, (-1, 0, 0)), (1, (-1, 0, 0))]
    1.2502 [(0, (0, 0, 0)), (0, (1, 0, 0)), (1, (1, 0, 0))]
    1.2502 [(0, (0, 0, 0)), (1, (0, 0, 0)), (1, (-1, 0, 0))]
    1.2502 [(0, (0, 0, 0)), (1, (0, 0, 0)), (1, (1, 0, 0))]
    """
    def __init__(self, bra, atom, damper, max_sing_thres):
        """Initialize Pairs object.

        Returns a Pairs object containing all pairs which the damper gives
        a value smaller than max_sing_thres
        """
        # save parameters
        self.lat = Bravais(bra)
        self.atom = np.array(atom, dtype=float)
        self.n_atom = self.atom.shape[0]
        if self.atom.shape != (self.n_atom, 3):
            raise ValueError("Invalid atom shape")
        self.damper = damper
        self.max_sing_thres = max_sing_thres
        Rcut = self.damper.Rcut(max_sing_thres)

        # get pairs
        self.abcs = self.lat.all_within(Rcut, add_base_PE=True)
        per_atom = []
        for abc in self.abcs:
            latvec = self.lat.latvec(abc)
            for atomvec in self.atom:
                per_atom.append(latvec + atomvec)
        self.per_atom = np.array(per_atom)
        pairs = get_pairs(self.atom, self.per_atom, Rcut)
        assert len(pairs) == self.n_atom

        # sort pairs
        self.pairs = []
        for i, partners in enumerate(pairs):
            proc_partners = []
            for jj in partners:
                if jj == i: continue
                Avec = self.atom[i]
                Bvec = self.per_atom[jj]
                vec = Bvec - Avec
                j = jj % self.n_atom    # original vector
                a = jj / self.n_atom
                abc = self.abcs[a]
                Rvec = self.lat.latvec(abc)
                assert np.allclose(vec, Rvec + self.atom[j] - self.atom[i])
                damp = damper.exponent(norm(vec), i, j)
                proc_partners.append((j, abc, damp))
            proc_partners.sort(key=operator.itemgetter(2)) # sort by damp
            self.pairs.append(proc_partners)

    def getvec(self, iabc):
        i, abc = iabc
        return self.lat.latvec(abc) + self.atom[i]

    def chains(self, n, thres, directed=True):
        """Return a list of (damp, [(i, abc), ...]) tuples of n-chains.

        This is the main workhorse and returns a weight-sorted list of
        bonds (n=2), bendings (n=3), or torsions (n=3).
        """
        res = []
        for i in xrange(self.n_atom):
            res.extend(self._chains_i(i, n, thres))
        if directed:
            final = []
            for damp, atlist in res:
                Avec = self.getvec(atlist[0])
                Dvec = self.getvec(atlist[-1])
                if isposvec(Dvec - Avec):
                    final.append((damp, atlist))
        final.sort()
        return final

    def _chains_i(self, i, n, thres):
        """Get all chains of length n from atom i."""
        if n < 2: raise ValueError("n should be at least 2.")
        res = []
        Avec = self.atom[i]
        for j, abc, damp in self.pairs[i]:
            if damp > thres: break    # they should be sorted
            if n == 2:
                # just pairs
                Bvec = self.lat.latvec(abc) + self.atom[j]
                tot_chain_damp = damp
                res.append((tot_chain_damp, [(i, (0,0,0)), (j, abc)]))
            else:
                # recursion
                rest_thres = thres - damp
                for chain_damp, atlist in self._chains_i(j, n-1, rest_thres):
                    shifted_atlist = [(i, (0,0,0))]
                    for (k, kabc) in atlist:
                        kabc = tuple([ai+ak for (ai, ak) in zip(abc, kabc)])
                        if i == k and kabc == (0, 0, 0):
                            break   # self reference
                        shifted_atlist.append((k, kabc))
                    else:
                        tot_chain_damp = damp + chain_damp
                        res.append((tot_chain_damp, shifted_atlist))
        return sorted(res)

class HessianBuilder(object):
    """Builder object for Hessians by rank-one additions.

    For a Hessian which is built successively as a sum of rank-one updates:
       H_{nu i, mu j} = \sum_k fac v_{nu i} v_{mu j}
    where each update vector is assumed to be sparse wrt the number of
    associated atoms.  The rank-one update is done by:
      HB.add_rank1(fac, vec)
    where vec == {nu: [x,y,z], ...}

    >>> HB = HessianBuilder(2, 0)
    >>> HB.add_rank1(1.0, {0: [1., 0., 0.]})
    >>> HB.add_rank1(0.5, {1: [0., 0., 1.]})
    >>> HD = np.zeros((2, 3, 2, 3))
    >>> HD[0,0,0,0] = 1.0
    >>> HD[1,2,1,2] = 0.5
    >>> np.allclose(HB.to_array(), HD)
    True
    """
    def __init__(self, n_atom, n_dyn_periodic):
        self.n_atom = n_atom
        self.n_dyn_periodic = n_dyn_periodic
        self.n_vec = n_atom + n_dyn_periodic
        self.Hdict = dict()

    def add_rank1(self, fac, vec):
        """Add rank-one term vec * vec * vec^T.

        Here, vec = {i_atom: np.array([xi, yi, zi]), ...}.
        """
        # Make sure that we have np.ndarrays
        for i_atom in vec:
            vec[i_atom] = np.asarray(vec[i_atom])
        # Perform dyadic product
        for i_atom, ivec in vec.iteritems():
            for j_atom, jvec in vec.iteritems():
                blk = fac * ivec[:,np.newaxis] * jvec[np.newaxis,:]
                _acc_dict((i_atom, j_atom), self.Hdict, blk)

    def add_rank1_from_atlist(self, fac, dq_datom, atlist):
        """Add rank-one term vec * vec * vec^T.

        Here, dq_atom = [np.array([xi, yi, zi]), ...], and
        atlist = [(i_tau, (a, b, c)), ...].
        """
        vecdict = dict()
        for (dqi, (i_atom, abc)) in zip(dq_datom, atlist):
            # Force term:
            _acc_dict(i_atom, vecdict, dqi)
            # Stress term:
            for i_bra, a_bra in enumerate(abc):
                i_vec = self.n_atom + i_bra
                if abs(a_bra) > 0 and i_bra < self.n_dyn_periodic:
                    _acc_dict(i_vec, vecdict, a_bra * dqi)
        self.add_rank1(fac, vecdict)

    def add_unity(self, fac):
        """Add multiple of unity."""
        blk = fac * np.identity(3)
        for i_vec in xrange(self.n_vec):
            _acc_dict((i_vec, i_vec), self.Hdict, blk)

    def to_array(self):
        """Construct full np.ndarray (only atomic coordinates, no stress)."""
        H = np.zeros((self.n_vec, 3, self.n_vec, 3))
        for i_atom, j_atom in self.Hdict:
            H[i_atom, :, j_atom, :] = self.Hdict[(i_atom, j_atom)]
        return H

def format_atlist(atlist, is_periodic):
    """Nicely format an atom list (atlist).

    Additionally adds 1 to atom numbers (-> start with 1):
    >>> print format_atlist([(0, (0, 0, 0)), (1, (0, -1, 0))], True)
      1( 0, 0, 0) --  2( 0,-1, 0)
    """
    if is_periodic:
        return " --".join(["%3i(%2i,%2i,%2i)" % ((i_atom+1,) + abc)
                           for (i_atom, abc) in atlist])
    else:
        return " --".join(["%3i" % (i_atom+1) for (i_atom, abc) in atlist])


def makeorthvec(orth):
    """Construct a (3 component) vector orthogonal to orth.

    >>> import numpy.random
    >>> vec = numpy.random.random(3)
    >>> assert np.dot(vec, makeorthvec(vec)) < 1e-12
    """
    orth /= norm(orth)
    vec = np.cross(orth, np.array([0., 0., 1.]))
    if (norm(vec) < 0.33):
        vec = np.cross(orth, np.array([1., 0., 0.]))
    return vec / norm(vec)


def model_matrix(bra, atom, builder, damper, thres, logfile=None):
    """Construct model Hessian.  Returns the HessianBuilder object.

    >>> builder = HessianBuilder(2, 0)
    >>> damper = LindhExponent(atom2row=[0, 0])    # Say, two hydrogens.
    >>> atom = np.array([[1., 0., 0.], [0., 0., 0.]])
    >>> HB = model_matrix(None, atom, builder, damper, 10.)
    >>> assert HB is builder
    >>> H = builder.to_array().reshape((3*2, 3*2))
    >>> assert np.allclose(H, H.T)
    >>> assert np.allclose(H[:,1:3], 0.)
    >>> assert np.allclose(H[:,4:6], 0.)
    >>> assert not np.allclose(H[0,0], 0.)
    >>> assert np.allclose(H[0,0], H[3,3])
    >>> assert np.allclose(H[0,3], -H[0,0])
    """
    thres_fac = np.exp(- thres) * K_BOND
    if logfile is not None:
        logfile.write("# Neglecting anything with a prefac < %8.3g "
                      "eV[/A^2]\n\n" % thres_fac)

    is_per = bra is not None and len(bra) > 0
    pairs = Pairs(bra, atom, damper, thres)

    # bonds:
    for damp, atlist in pairs.chains(2, thres):
        fac = np.exp(- damp) * K_BOND
        if fac < thres_fac:
            continue
        vecs = [pairs.getvec(at) for at in atlist]
        q, dq = q_bond(*vecs)
        if logfile is not None:
            logfile.write("# bond: %4.2f A %s  "
                          "[damp: %8.3g; prefac: %8.3g eV/A^2]\n" %
                          (q, format_atlist(atlist, is_per), damp, fac))
        builder.add_rank1_from_atlist(fac, dq, atlist)
    if logfile is not None:
        logfile.write("\n")

    # bendings:
    for damp, atlist in pairs.chains(3, thres):
        fac = np.exp(- damp) * K_BENDING
        if fac < thres_fac:
            continue
        vecs = [pairs.getvec(at) for at in atlist]
        q, dq = q_bending(*vecs)
        if logfile is not None:
            logfile.write("# angle: %4.0f deg %s  "
                          "[damp: %8.3g; prefac: %8.3g eV]\n" %
                          (np.rad2deg(q), format_atlist(atlist, is_per),
                           damp, fac))
        if 0.05*np.pi < q < 0.95*np.pi:
            builder.add_rank1_from_atlist(fac, dq, atlist)
        else:
            Avec, Bvec, Cvec = vecs
            wvec = makeorthvec(Cvec - Avec)
            q, dq = q_bending(Avec, Bvec, Cvec, direction=wvec)
            builder.add_rank1_from_atlist(fac, dq, atlist)
            wvec = np.cross(Bvec - Avec, wvec)
            wvec /= norm(wvec)
            q, dq = q_bending(Avec, Bvec, Cvec, direction=wvec)
            builder.add_rank1_from_atlist(fac, dq, atlist)
    if logfile is not None:
        logfile.write("\n")

    # torsions
    for damp, atlist in pairs.chains(4, thres):
        fac = np.exp(- damp) * K_TORSION
        if fac < thres_fac:
            continue
        vecs = [pairs.getvec(at) for at in atlist]
        try:
            q, dq = q_torsion(*vecs)
            if logfile is not None:
                logfile.write("# torsion: %4.0f deg %s  "
                              "[damp: %8.3g; prefac: %8.3g eV]\n" %
                              (np.rad2deg(q), format_atlist(atlist, is_per),
                               damp, fac))
            builder.add_rank1_from_atlist(fac, dq, atlist)
        except ValueError:
            if logfile is not None:
                logfile.write("# torsion: ---- deg %s "
                              "[damp: %8.3g; prefac: %8.3g eV]\n" %
                              (format_atlist(atlist, is_per), damp, fac))
    if logfile is not None:
        logfile.write("\n")

    return builder


##############################################################################
##################################### dd #####################################
##############################################################################

# The "dd" functions generally share the same interface.  The arguments are
# 2-tuples where the first object contains the actual parameter and the second
# object contains the derivatives of this parameter with respect to the
# original input parameters.  From this, the output value is calculated.
# Additionally, the derivatives of this value with respect to the original
# input parameters are evaluated by the chain rule.  The return value is a
# tuple of the result value and its derivative.

# Obviously, this is not the most efficient way to do it.  But at least,
# it works...

def _dd_matmat(val_shape, dval_di, di_dvar):
    val_rank = len(val_shape)
    assert val_shape == np.shape(dval_di)[:val_rank]
    di_shape = np.shape(dval_di)[val_rank:]
    di_rank = len(di_shape)
    assert di_shape == np.shape(di_dvar)[:di_rank]
    axes1 = range(val_rank, val_rank+di_rank)
    return np.tensordot(dval_di, di_dvar, (axes1, range(di_rank)))

def _dd_broadcast(val, dval):
    val_rank = len(np.shape(val))
    assert np.shape(val) == np.shape(dval)[:val_rank]
    dval_rank = len(np.shape(dval))
    newshape = np.shape(val) + (dval_rank-val_rank)*(1,)
    return np.reshape(val, newshape)

def dd_sum(*arg_ds):
    shape = np.shape(arg_ds[0][1])
    res = np.float(0.)
    dres = np.zeros(shape)
    for arg, darg in arg_ds:
        res += np.asarray(arg)
        dres += np.asarray(darg)
    return res, dres

def dd_mult(vec_d, fac):
    vec, dvec = vec_d
    return fac*vec, fac*dvec

def dd_prod(*arg_ds):
    shape = np.shape(arg_ds[0][1])
    res = np.float(1.)
    dres = np.zeros(shape)
    for arg, darg in arg_ds:
        dres *= arg                     # update previous derivs
        dres += np.asarray(darg) * res  # update with previous factors
        res *= arg                      # update value
    return res, dres

def dd_power(var_d, n):
    var, dvar = var_d
    val = var**n
    dval = n*(var**(n-1)) * dvar
    return val, dval

def dd_dot(vec1_d, vec2_d):
    vec1, dvec1 = vec1_d
    vec2, dvec2 = vec2_d
    res = np.dot(vec1, vec2)
    dres = (np.tensordot(vec1, dvec2, (-1, 0)) +
            np.tensordot(vec2, dvec1, (-1, 0)))
    return res, dres

def dd_cross(vec1_d, vec2_d):
    vec1, dvec1 = vec1_d
    vec2, dvec2 = vec2_d
    assert np.shape(vec1) == np.shape(vec2) == (3,)   # otherwise...
    res = np.cross(vec1, vec2)
    dres = - np.cross(vec2, dvec1, axisb=0).T + np.cross(vec1, dvec2, axisb=0).T
    return res, dres

def dd_norm(vec_d):
    return dd_power(dd_dot(vec_d, vec_d), 0.5)

def dd_normalized(vec_d):
    vec, dvec = vec_d
    fac, dfac = dd_power(dd_norm(vec_d), -1.)
    res = fac*vec
    dres = fac*dvec + vec[:,np.newaxis]*dfac[np.newaxis,:]
    return res, dres

def dd_cosv1v2(vec1_d, vec2_d):
    return dd_prod(dd_dot(vec1_d, vec2_d),
                   dd_power(dd_norm(vec1_d), -1.),
                   dd_power(dd_norm(vec2_d), -1.))

def dd_arccos(val_d):
    val, dval = val_d
    if 1. < abs(val) < 1. + 1e-10:
        val = np.sign(val)
    res = np.arccos(val)
    vval = _dd_broadcast(val, dval)
    dres = - 1. / np.sqrt(1. - vval**2) * dval
    return res, dres

def dd_arcsin(val_d):
    val, dval = val_d
    if 1. < abs(val) < 1. + 1e-10:
        val = np.sign(val)
    res = np.arcsin(val)
    vval = _dd_broadcast(val, dval)
    dres = 1. / np.sqrt(1. - vval**2) * dval
    return res, dres

def dd_angle(vec1_d, vec2_d):
    return dd_arccos(dd_cosv1v2(vec1_d, vec2_d))

def dd_bondlength(pos1_d, pos2_d):
    AB_d = dd_sum(pos2_d, dd_mult(pos1_d, -1.))
    return dd_norm(AB_d)

def dd_bondangle(pos1_d, pos2_d, pos3_d):
    BA_d = dd_sum(pos2_d, dd_mult(pos1_d, -1.))
    BC_d = dd_sum(pos2_d, dd_mult(pos3_d, -1.))
    return dd_angle(BA_d, BC_d)

def dd_bondangle_directed(pos1_d, pos2_d, pos3_d, dir_d):
    BA_d = dd_sum(pos2_d, dd_mult(pos1_d, -1.))
    BC_d = dd_sum(pos2_d, dd_mult(pos3_d, -1.))
    return dd_directed_angle(BA_d, BC_d, dir_d)

def dd_arctan2(y_d, x_d):
    y, dy = y_d
    x, dx = x_d
    phi = np.arctan2(y, x)
    tan, dtan = dd_prod(x_d, dd_power(y_d, -1.))
    tan = _dd_broadcast(tan, dtan)
    dphi = (1. + tan**2) * dtan
    return phi, dphi

def dd_directed_angle(vec1_d, vec2_d, dir_d):
    ndir_d = dd_normalized(dir_d)
    vv1_d = dd_cross(vec1_d, ndir_d)
    vv2_d = dd_cross(vec2_d, ndir_d)
    if (norm(vv1_d[0]) < 1e-7 or
        norm(vv2_d[0]) < 1e-7):
        return 0., np.zeros(np.shape(vec1_d[1])[1:])
    vv1_d = dd_normalized(vv1_d)
    vv2_d = dd_normalized(vv2_d)
    cosphi_d = dd_dot(vv1_d, vv2_d)
    vvv_d = dd_cross(vv1_d, vv2_d)
    sinphi_d = dd_dot(vvv_d, ndir_d)
    # phi_d = dd_arctan2(sinphi_d, cosphi_d)
    if (abs(cosphi_d[0]) < np.sqrt(0.5)):
        phi, dphi = dd_arccos(cosphi_d)
        if sinphi_d[0] < 0.:
            phi *= -1.
            dphi *= -1.
    else:
        phi, dphi = dd_arcsin(sinphi_d)
        if cosphi_d[0] < 0.:
            phi = - np.pi - phi
            if phi < np.pi: phi += 2*np.pi
            dphi *= -1.
    return phi, dphi

def dd_bondtorsion(pos1_d, pos2_d, pos3_d, pos4_d):
    BA_d = dd_sum(pos2_d, dd_mult(pos1_d, -1.))
    BC_d = dd_sum(pos2_d, dd_mult(pos3_d, -1.))
    CD_d = dd_sum(pos3_d, dd_mult(pos4_d, -1.))
    return dd_directed_angle(BA_d, CD_d, BC_d)

##############################################################################
###################################### q #####################################
##############################################################################

def q_bond(Avec, Bvec):
    """Bond length and derivative wrt vector AB.

    Test:
    >>> np.allclose(q_bond([0., 0., 0.], [1., 1., 1.])[0], np.sqrt(3.))
    True
    >>> assert _test_qgrad(q_bond, 2) < 1e-5
    """
    Avec_d = (np.asarray(Avec), np.c_[ID, ZERO])
    Bvec_d = (np.asarray(Bvec), np.c_[ZERO, ID])
    q, dq = dd_bondlength(Avec_d, Bvec_d)
    return q, dq.reshape((2, 3))

def q_bending(Avec, Bvec, Cvec, direction=None):
    """Bond angle and derivative wrt vectors AB and BC.

    Test:
    >>> A = np.array([ 1, 1, 1])
    >>> B = np.zeros(3)
    >>> C = np.array([-1,-1, 1])
    >>> print round(np.rad2deg(q_bending(A, B, C)[0]), 1)
    109.5
    >>> assert _test_qgrad(q_bending, 3) < 1e-5
    """
    Avec_d = (np.asarray(Avec), np.c_[ID, ZERO, ZERO])
    Bvec_d = (np.asarray(Bvec), np.c_[ZERO, ID, ZERO])
    Cvec_d = (np.asarray(Cvec), np.c_[ZERO, ZERO, ID])
    if direction is None:
        q, dq = dd_bondangle(Avec_d, Bvec_d, Cvec_d)
    else:
        dir_d = (direction, np.c_[ZERO, ZERO, ZERO])
        q, dq = dd_bondangle_directed(Avec_d, Bvec_d, Cvec_d, dir_d)
    return q, dq.reshape((3, 3))

def q_torsion(Avec, Bvec, Cvec, Dvec):
    """Bond torsion and derivative wrt vectors AB, BC, and CD.

    Test:
    >>> A = np.array([0., 0., 1.])
    >>> B = np.array([0., 0., 0.])
    >>> C = np.array([1., 0., 0.])
    >>> D = np.array([1., 1., 0.])
    >>> print round(np.rad2deg(q_torsion(A, B, C, D)[0]), 5)
    90.0
    >>> try:
    ...    assert _test_qgrad(q_torsion, 4) < 1e-5
    ... except ValueError:   # May happen with bad luck.
    ...    pass
    """
    ABvec = Bvec - Avec
    BCvec = Cvec - Bvec
    CDvec = Dvec - Cvec
    cosABC = np.dot(ABvec, BCvec) / (norm(ABvec) * norm(BCvec))
    cosBCD = np.dot(BCvec, CDvec) / (norm(BCvec) * norm(CDvec))
    if max(abs(cosABC), abs(cosBCD)) > 0.99:   # nearly linear angle
        raise ValueError("Nearly linear angle")
    else:
        Avec_d = (np.asarray(Avec), np.c_[ID, ZERO, ZERO, ZERO])
        Bvec_d = (np.asarray(Bvec), np.c_[ZERO, ID, ZERO, ZERO])
        Cvec_d = (np.asarray(Cvec), np.c_[ZERO, ZERO, ID, ZERO])
        Dvec_d = (np.asarray(Dvec), np.c_[ZERO, ZERO, ZERO, ID])
        q, dq = dd_bondtorsion(Avec_d, Bvec_d, Cvec_d, Dvec_d)
        return q, dq.reshape(4, 3)

##############################################################################
############################# unit test utilities ############################
##############################################################################

def _test_qgrad(q_func, n):
    import scipy.optimize
    import numpy.random
    x0 = np.random.standard_normal(3*n)
    def func(x):
        vecs = np.asarray(x).reshape((n,3))
        q, dq = q_func(*vecs)
        return q
    def grad(x):
        vecs = np.asarray(x).reshape((n,3))
        q, dq = q_func(*vecs)
        return np.reshape(dq, (3*n))
    return scipy.optimize.check_grad(func, grad, x0)

def testmod():
    import doctest
    doctest.testmod(raise_on_error=False)

##############################################################################
################################# read_.... ##################################
##############################################################################

def read_aims(in_):
    """Read FHI-aims geometry.in file.

    Only read 'atom' and 'lattice_vector' lines.
    >>> from StringIO import StringIO
    >>> f = StringIO("atom 0. 1. 2. He")
    >>> atom, atom_name, bra = read_aims(f)
    >>> np.allclose(atom, [[0., 1., 2.]])
    True
    >>> atom_name == ['He']
    True
    >>> bra.shape == (0, 3)
    True
    """
    atom = []
    bra = []
    atom_name = []
    for i, line in enumerate(in_):
        line, _, _ = line.partition("#")
        fields = line.split()
        if not fields: continue
        if fields[0] == "atom":
            if len(fields) != 5:
                raise ValueError("Invalid atom in line %i" % i)
            atom.append(map(float, fields[1:4]))
            atom_name.append(fields[4])
        elif fields[0] == "lattice_vector":
            if len(fields) != 4:
                raise ValueError("Invalid Bravais in line %i" % i)
            bra.append(map(float, fields[1:]))
    atom = np.array(atom)
    n_atom, three = atom.shape
    assert three == 3
    if bra:
        bra = np.array(bra)
    else:
        bra = np.empty((0, 3))
    n_periodic, three = np.shape(bra)
    assert three == 3
    if n_periodic > 3:
        raise ValueError("Too many lattice_vectors")
    return atom, atom_name, bra

def read_xyz(in_):
    r"""Read xmakemol xyz file.

    Example:
    >>> from StringIO import StringIO
    >>> f = StringIO("2\nComment\nH 0 0 0.3\nH 0 0 -0.3")
    >>> atom, atom_name, bra = read_xyz(f)
    >>> np.allclose(atom, [[0., 0., 0.3], [0., 0., -0.3]])
    True
    >>> atom_name == ['H', 'H']
    True
    >>> bra.shape == (0, 3)
    True
    """
    atom = []
    bra = []
    atom_name = []
    n_atom = int(in_.next().split()[0])
    comment = in_.next()

    for i, line in enumerate(in_):
        fields = line.split()
        atom_name.append(fields[0])
        atom.append([float(x) for x in fields[1:4]])
        rest = fields[4:]
        while "crystal_vector" in fields:
            pos = fields.index("crystal_vector")
            bra.append([float(x) for x in fields[pos+1:pos+4]])
            del fields[pos:pos+4]
    atom = np.array(atom)
    n_atom, three = atom.shape
    assert three == 3
    if bra:
        bra = np.array(bra)
    else:
        bra = np.empty((0, 3))
    n_periodic, three = np.shape(bra)
    assert three == 3
    if n_periodic > 3:
        raise ValueError("Too many lattice_vectors")
    return atom, atom_name, bra

def read_babel(in_, informat):
    r"""Read general file format using babel."""
    try:
        p = subprocess.Popen(["babel", "-i", informat, "-o", "xyz"],
                             stdin=in_, stdout=subprocess.PIPE)
    except OSError:
        sys.stderr.write(
            "Unknown input format and/or babel is missing.\n"
            "If the input file is in FHI-aims geometry.in format,\n"
            "try using '-i aims'.\n")
        sys.exit(2)
    atom, atom_name, bra = read_xyz(p.stdout)
    if p.wait() != 0:
        sys.exit(p.wait())
    return atom, atom_name, bra

def write_aims(out, atom, atom_name, bra):
    """Write out FHI-aims geometry.in file.

    Only the 'lattice_vector' and 'atom' lines are written, though.
    >>> write_aims(sys.stdout, np.array([[1., 1., 1.]]), ['H'], [])
    atom    1.0000000000000000    1.0000000000000000    1.0000000000000000  H
    <BLANKLINE>
    """
    for vec in bra:
        out.write("lattice_vector %f  %f  %f\n" % tuple(vec))
    if bra is not None and len(bra) > 0:
        out.write("\n")
    for vec, name in zip(atom, atom_name):
        out.write("atom  %20.16f  %20.16f  %20.16f %2s\n" %
                  (tuple(vec) + (name,)))
    out.write("\n")


##############################################################################
##################################### main ###################################
##############################################################################

def main():
    import optparse
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-i", "--informat", metavar="FORMAT",
                      action="store", default=None,
                      help="Input format extension "
                      "(default: taken from input suffix, 'in'=='aims')")
    parser.add_option("-t", "--test", action="store_true",
                      help="Perform doctests")
    parser.add_option("-v", "--verbose", action="store_true",
                      help="Add information about internal coordinats")
    parser.add_option("-c", "--cutoff", type=float, metavar="CUT",
                      default=15., help="Cut-off value in exponent [15.]")
    parser.add_option("-a", "--add-unity", metavar="MU",
                      type=float, default = 0.005,
                      help="Add MU (in eV/A^2) to "
                      "final Hessian [default: 0.005]")
    parser.add_option("-p", "--append", action="store_true",
                      help="Append instead of overwrite.")
    parser.add_option("--stress", action="store_true",
                      help="Also output (experimental) stress terms.")
    parser.add_option("--full", type="string",
                      help="Store full Hessian to file.")

    options, args = parser.parse_args()

    if options.test:
        testmod()
        sys.exit()

    inname = "-"
    outname = "-"

    if len(args) > 2:
        parser.error("Too many arguments")
    elif len(args) == 2:
        inname, outname = args
    elif len(args) == 1:
        inname, = args

    if options.append:
        if len(args) != 1:
            parser.error("Will only --append to the input file.")
        outname = inname

    if options.informat is not None:
        informat = options.informat
    elif "." in inname:
        _name, _dot, informat = inname.rpartition(".")
    else:
        informat = "in"

    if inname == "-":
        in_ = sys.stdin
    else:
        in_ = open(inname)

    # read geometry.in
    if informat in ("in", "aims"):
        atom, atom_name, bra = read_aims(in_)
        print atom_name
    elif informat == "xyz":
        atom, atom_name, bra = read_xyz(in_)
    else:
        atom, atom_name, bra = read_babel(in_, informat)
    if in_ != sys.stdin:
        in_.close()
    n_atom = atom.shape[0]
    assert len(atom_name) == n_atom

    if options.append or inname == outname != "-":
        out = open(outname, "a")
        out.write("\n")
    else:        
        if outname == "-":
            out = sys.stdout
        else:
            out = open(outname, "w")
            out.write("\n")
        write_aims(out, atom, atom_name, bra)

    if options.stress and bra is not None:
        n_dyn_periodic = len(bra)
    else:
        n_dyn_periodic = 0
    n_vec = n_atom + n_dyn_periodic
    builder = HessianBuilder(n_atom, n_dyn_periodic)
    damper = LindhExponent([name2row(name) for name in atom_name])
    logfile = out if options.verbose else None
    model_matrix(bra, atom, builder, damper, options.cutoff, logfile=logfile)
    if options.add_unity > 0.:
        builder.add_unity(options.add_unity)
    for i_atom, j_atom in sorted(builder.Hdict):
        if i_atom > j_atom: continue   # lower triangle
        blk = builder.Hdict[(i_atom, j_atom)]
        if np.all(np.abs(blk) < 5e-7): continue   # would be 0.000000s anyway
        blk_str = " ".join(["%11.6f" % num for num in blk.T.flatten()])
        if j_atom < n_atom:
            # Pure force term
            out.write("hessian_block %3i %3i   %s\n" % (
                    i_atom+1, j_atom+1, blk_str))
        elif options.stress:
            if i_atom < n_atom:
                # Mix term, transposed
                blk_str = " ".join(["%11.6f" % num for num in blk.flatten()])
                out.write("hessian_block_lv_atom %3i %3i   %s\n" % (
                        j_atom-n_atom+1, i_atom+1, blk_str))
            else:
                # Pure stress term
                out.write("hessian_block_lv %3i %3i   %s\n" % (
                        i_atom-n_atom+1, j_atom-n_atom+1, blk_str))


    if options.full is not None:
        H = builder.to_array().reshape((3*n_vec, 3*n_vec))
        np.savetxt(options.full, H)

    if options.verbose:
        H = builder.to_array().reshape((3*n_vec, 3*n_vec))
        # sys.stderr.write(str(H.round(2)) + "\n\n")
        w, v = np.linalg.eig(H)
        assert np.all(np.imag(w) < 1e-10)
        w = np.real(w)
        v = np.real(v)
        out.write("\n")
        out.write("# ev: " + " ".join(["%.3f" % num for num in sorted(w)]) +
                  "\n")
    
if __name__ == "__main__":
    main()
