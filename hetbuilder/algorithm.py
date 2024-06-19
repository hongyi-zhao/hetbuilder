import ase.io
from ase.atoms import Atoms
from ase.spacegroup import Spacegroup
from ase.geometry import permute_axes
from ase.geometry.analysis import Analysis

from itertools import combinations_with_replacement
from scipy.spatial import KDTree

from dataclasses import dataclass

import numpy as np

from scipy.linalg import polar

from hetbuilder.log import *
from hetbuilder.atom_checks import check_atoms

from ase.neighborlist import (
    NeighborList,
    natural_cutoffs,
    NewPrimitiveNeighborList,
    find_mic,
)

import sys

from hetbuilder.hetbuilder_backend import (
    double2dVector,
    double1dVector,
    int1dVector,
    int2dVector,
    CppAtomsClass,
    CppCoincidenceAlgorithmClass,
    CppInterfaceClass,
    get_number_of_omp_threads,
)


def ase_atoms_to_cpp_atoms(atoms: "ase.atoms.Atoms") -> "CppAtomsClass":
    """Converts :class:`~ase.atoms.Atoms` to the C++ CppAtomsClass."""
    lattice = atoms.cell.copy()
    positions = atoms.positions.copy()
    atomic_numbers = int1dVector([k for k in atoms.numbers])
    lattice = double2dVector([double1dVector(k) for k in lattice])
    positions = double2dVector([double1dVector(k) for k in positions])
    indices = int1dVector([k.index for k in atoms])
    magmoms = atoms.get_initial_magnetic_moments()
    magmoms = double1dVector([k for k in magmoms])
    return CppAtomsClass(lattice, positions, atomic_numbers, indices, magmoms)


def cpp_atoms_to_ase_atoms(cppatoms: "CppAtomsClass") -> "ase.atoms.Atoms":
    """Converts the C++ CppAtomsClass to :class:`~ase.atoms.Atoms`"""
    lattice = [[j for j in k] for k in cppatoms.lattice]
    positions = [[j for j in k] for k in cppatoms.positions]
    numbers = [i for i in cppatoms.atomic_numbers]
    magmoms = [i for i in cppatoms.magmoms]
    atoms = Atoms(
        numbers=numbers,
        positions=positions,
        cell=lattice,
        pbc=[True, True, True],
        magmoms=magmoms,
    )
    return atoms


def check_angles(
    angle_stepsize: float = 1, angle_limits: tuple = (0, 180), angles: list = []
) -> list:
    """ Helper function to assert correct input of angles."""
    if len(angles) == 0:
        a1 = angle_limits[0]
        a2 = angle_limits[1]
        assert a2 > a1, "Second angle must be larger than first one."
        assert angle_stepsize > 0, "Angle stepsize must be larger than zero."
        assert angle_stepsize < abs(
            a2 - a1
        ), "Angle stepsize must be larger then difference between angles."
        angles = list(np.arange(a1, a2, step=angle_stepsize)) + [a2]
        logger.info(
            "Searching {:d} angles between {:.1f} and {:.1f} degree with a stepsize of {:.1f} degree.".format(
                len(angles), a1, a2, angle_stepsize
            )
        )
        return angles
    elif angles != None:
        msg = ", ".join([str(k) for k in angles])
        logger.info("Calculating the following angles: {} in degree.".format(msg))
        return list(angles)
    else:
        logger.error("Angle specifications not recognized.")


def get_bond_data(atoms: "ase.atoms.Atoms", return_bonds: bool = True) -> tuple:
    """ Returns a tuple holding bond indices and average bond values.
    
    If the input structure is larger than 1000 atoms, the neighborlist is not computed and not all bonds are determined.
    Only a subsample is queried to determine the strain.
    """
    atoms = atoms.copy()
    symbs = set(atoms.get_chemical_symbols())
    bonds = None

    if len(atoms) > 1000 or return_bonds == False:
        # not computing neighborlists for all bonds here, too slow
        return_bonds = False
        s1 = set(atoms.get_chemical_symbols())
        s2 = set()
        i = 1
        p = atoms.positions
        tree = KDTree(atoms.positions)
        middle = np.mean(p, axis=0)
        middle[2] = np.min(p[:, 2])  # so we are not searching in vacuum
        while s1 != s2 and i < 5:  # 增加最大迭代次数
            sample = tree.query_ball_point(middle, r=i * 15)  # 增加采样半径
            subatoms = atoms[sample]
            s2 = set(subatoms.get_chemical_symbols())
#            print(f"Iteration {i}: Found symbols {s2}")
            i += 1
            if i == 5:
                raise Exception(
                    "Could not find all species in the subquery. This should not happen."
                )

        atoms = subatoms

    pairs = list(combinations_with_replacement(symbs, 2))
    cutoffs = np.array(natural_cutoffs(atoms)) * 1.25  # 增加截断半径
    nl = NeighborList(cutoffs, skin=0.0, primitive=NewPrimitiveNeighborList)
    nl.update(atoms)

    if return_bonds:
        nbonds = nl.nneighbors + nl.npbcneighbors
        bonds = []
        for a in range(len(atoms)):
            indices, offsets = nl.get_neighbors(a)
            for i, offset in zip(indices, offsets):
                startvector = atoms.positions[a]
                endvector = atoms.positions[i] + offset @ atoms.get_cell()
                if np.sum((endvector - startvector) ** 2) > 0:
                    bonds.append([a, i])

    unique_bonds = {}
    for p in pairs:
        bond_lengths = []
        for a in range(len(atoms)):
            indices, offsets = nl.get_neighbors(a)
            for i, offset in zip(indices, offsets):
                if {atoms[a].symbol, atoms[i].symbol} == set(p):
                    startvector = atoms.positions[a]
                    endvector = atoms.positions[i] + offset @ atoms.get_cell()
                    bond_length = np.linalg.norm(endvector - startvector)
                    if bond_length > 0:  # 过滤掉键长为 0 的键
                        bond_lengths.append(bond_length)
#                        print(f"Bond found: {atoms[a].symbol}-{atoms[i].symbol} with length {bond_length:.2f} Å")
#                    else:
#                        print(f"Bond length 0 Å found: {atoms[a].symbol}-{atoms[i].symbol}")
        if bond_lengths:
            avg = np.average(bond_lengths)
            unique_bonds[p] = avg
#        else:
#            print(f"No bonds found for pair: {p}")

    return bonds, unique_bonds


@dataclass
class Interface:
    """Exposes the C++ implementation of the CppInterfaceClass.
    
    Attributes:
        bottom (ase.atoms.Atoms): Lower layer as supercell.
        top (ase.atoms.Atoms): Upper layer as supercell.
        stack (ase.atoms.Atoms): Combined lower and upper layer as supercell.
        M (numpy.ndarray): Supercell matrix M.
        N (numpy.ndarray): Supercell matrix N.
        angle (float): Twist angle in degree.
        stress (float): Stress measure of the unit cell.
        strain (float): Strain measure of the bond lengths.

    """

    def __init__(
        self, interface: "CppInterfaceClass" = None, weight=0.5, **kwargs
    ) -> None:
        bottom = cpp_atoms_to_ase_atoms(interface.bottom)
        top = cpp_atoms_to_ase_atoms(interface.top)
        stack = cpp_atoms_to_ase_atoms(interface.stack)
        # First, ensure all atoms are within the unit cell boundaries.
        stack.wrap(pretty_translation=True)

        # Here, `recenter` was originally used, but the structure here is based on the results from atom_functions.cpp for 
        # subsequent operations and ultimately determines the output structure. 
        # The processing logic based on the `scale_cell_xy` method in atom_functions.cpp is already sufficient, 
        # and recenter cannot be used here again, otherwise, it will render the --vacuum parameter ineffective.
        self.bottom = bottom
        self.top = top
        self.stack = stack
        self.M = [[j for j in k] for k in interface.M]
        self.N = [[j for j in k] for k in interface.N]
        self.angle = interface.angle
        self._weight = weight
        self._stress = None
        self.bbl = kwargs.get("bottom_bond_lengths", None)
        self.tbl = kwargs.get("top_bond_lengths", None)
        self._bonds, self._bond_lengths = get_bond_data(self.stack)

    def __repr__(self):
        return "{}(M={}, N={}, angle={:.1f}, stress={:.1f})".format(
            self.__class__.__name__, self.M, self.N, self.angle, self.stress,
        )

    @property
    def stress(self) -> float:
        """Returns the stress measure."""
        return self.measure_stress()

    @property
    def strain(self) -> float:
        """Returns the strain measure."""
        return self.measure_strain()

    def measure_stress(self) -> float:
        """Measures the stress on both unit cells."""
        A = self.bottom.cell.copy()[:2, :2]
        B = self.top.cell.copy()[:2, :2]
        C = A + self._weight * (B - A)
        T1 = C @ np.linalg.inv(A)
        T2 = C @ np.linalg.inv(B)

        def measure(P):
            eps = P - np.identity(2)
            meps = np.sqrt(
                (
                    eps[0, 0] ** 2
                    + eps[1, 1] ** 2
                    + eps[0, 0] * eps[1, 1]
                    + eps[1, 0] ** 2
                )
                / 4
            )
            return meps

        U1, P1 = polar(T1)  # this one goes counterclockwise
        U2, P2 = polar(T2)  # this one goes clockwise
        # u is rotation, p is strain
        meps1 = measure(P1)
        meps2 = measure(P2)
        stress = meps1 + meps2
        # return (stress, P1 - np.identity(2), P2 - np.identity(2))
        return stress

    # def measure_strain(self) -> float:
        # """Measures the average strain on bond lengths on both substructures."""
        # bond_lengths = self.bond_lengths
        # bottom_strain = []
        # top_strain = []
        # for (k1, b1) in bond_lengths.items():
            # for k3, b3 in self.bbl.items():
                # if (k3 == k1) or (k3[::-1] == k1):
                    # d = np.abs((b3 - b1)) / b1 * 100
                    # bottom_strain.append(d)
            # for k3, b3 in self.tbl.items():
                # if (k3 == k1) or (k3[::-1] == k1):
                    # d = np.abs((b3 - b1)) / b1 * 100
                    # top_strain.append(d)
        # strain = np.average(bottom_strain) + np.average(top_strain)
        # return strain
    
    
    #To modify the measure_strain method according to the provided definition and
    #correctly handle potential issues (like divisions by zero or cases where bottom_strain
    #or top_strain lists are empty, which would lead to nan values), you can update the method to include checks for these cases:
    
    #In this updated version, the method now includes checks to avoid division by zero when calculating the strain differences.
    #Additionally, it handles cases where either bottom_strain or top_strain lists might be empty, which would previously result
    #in np.average returning nan. By providing default values (in this case, 0), you ensure that the method returns a valid numerical
    #value even when no strain data is available, thus avoiding the generation of nan values that could lead to errors in subsequent calculations.

    def measure_strain(self) -> float:
        """Measures the average strain on bond lengths on both substructures."""
        bond_lengths = self.bond_lengths
        bottom_strain = []
        top_strain = []
        for (k1, b1) in bond_lengths.items():
            for k3, b3 in self.bbl.items():
                if (k3 == k1) or (k3[::-1] == k1):
                    if b1 != 0:  # Check to avoid division by zero
                        d = np.abs((b3 - b1)) / b1 * 100
                        bottom_strain.append(d)
            for k3, b3 in self.tbl.items():
                if (k3 == k1) or (k3[::-1] == k1):
                    if b1 != 0:  # Check to avoid division by zero
                        d = np.abs((b3 - b1)) / b1 * 100
                        top_strain.append(d)
    
        # Handle cases where strain lists are empty to avoid 'nan' values
        if bottom_strain:
            avg_bottom_strain = np.average(bottom_strain)
        else:
            avg_bottom_strain = 0  # Default value or another suitable value
    
        if top_strain:
            avg_top_strain = np.average(top_strain)
        else:
            avg_top_strain = 0  # Default value or another suitable value
    
        # Calculate total average strain
        strain = avg_bottom_strain + avg_top_strain
        return strain
    

    @property
    def bonds(self):
        return self._bonds

    @property
    def bond_lengths(self):
        return self._bond_lengths


#相比于backend/coincidence_algorithm.cpp 中的原型函数：
#std::vector<Interface> CoincidenceAlgorithm::run(int Nmax,
#                                                 int Nmin,
#                                                 double1dvec_t angles,
#                                                 double tolerance,
#                                                 double weight,
#                                                 double distance,
#                                                 double vacuum,
#                                                 bool standardize,
#                                                 int no_idealize,
#                                                 double symprec,
#                                                 double angle_tolerance,
#                                                 int verbose)
#
# 此处进行了进一步的封装， 为了方便调用和控制，添加了下面的两个参数： 
#        angle_limits: tuple = (0, 90),
#        angle_stepsize: float = 1.0,
class CoincidenceAlgorithm:
    """Exposes the C++ implementation of the CppCoincidenceAlgorithmClass.
    
    Args:
        bottom (ase.atoms.Atoms): Lower layer, needs to be two-dimensional.
        top (ase.atoms.Atoms): Upper layer, needs to be two-dimensional.   
    """

    def __init__(self, bottom: "ase.atoms.Atoms", top: "ase.atoms.Atoms") -> None:
        # `check_atoms` here is used for necessary checks, confirmation, and preprocessing of the input structure, 
        # and this function calls `recenter`:
        self.bottom = check_atoms(bottom)
        self.top = check_atoms(top)
        _, self.bdl = get_bond_data(bottom)
        _, self.tbl = get_bond_data(top)

    def __repr__(self):
        return "{}(bottom={}, top={})".format(
            self.__class__.__name__, self.bottom, self.top
        )

    def run(
        self,
        Nmax: int = 10,
        Nmin: int = 0,
        angles: list = [],
        angle_limits: tuple = (0, 90),
        angle_stepsize: float = 1.0,
        tolerance: float = 0.1,
        weight: float = 0.5,
        distance: float = 3.5,
        vacuum: float = 15,
        standardize: bool = False,
        no_idealize: bool = False,
        symprec: float = 1e-5,
        angle_tolerance: float = 5,
        verbosity: int = 0,
    ) -> list:
        """Executes the coincidence lattice algorithm.

        Args:
            Nmax (int): Maximum number of translations. Defaults to 10.
            Nmin (int): Minimum number of translations. Defaults to 0.
            angles (list): List of angles in degree to search. Takes precedence over angle_limits and angle_stepsize.
            angle_limits (tuple): Lower and upper bound of angles to look through with given step size by angle_stepsize. Defaults to (0, 90) degree.
            angle_stepsize (float): Increment of angles to look through. Defaults to 1.0 degree.
            tolerance (float): Tolerance criterion to accept lattice match. Corresponds to a distance in Angström. Defaults to 0.1.
            weight (float): The coincidence unit cell is C = A + weight * (B-A). Defaults to 0.5.
            distance (float): Interlayer distance of the stacks. Defaults to 3.5 Angström.
            vacuum (float): Thickness of the vacuum layer of the stacks. Defaults to 15.0 Angström.
            standardize (bool): Perform spglib standardization. Defaults to False.
            no_idealize (bool): Does not idealize unit cell parameters in the spglib standardization routine. Defaults to False.
            symprec (float): Symmetry precision for spglib. Defaults to 1e-5 Angström.
            angle_tolerance (float): Angle tolerance for the spglib `spgat` routines. Defaults to 5.
            verbosity (int): Debug level for printout of Coincidence Algorithm. Defaults to 0.

        Returns:
            list : A list of :class:`~hetbuilder.algorithm.Interface`.

        """
        bottom = ase_atoms_to_cpp_atoms(self.bottom)
        top = ase_atoms_to_cpp_atoms(self.top)
        angles = check_angles(
            angle_limits=angle_limits, angle_stepsize=angle_stepsize, angles=angles
        )
        if (self.bottom == self.top) and (0 in angles):
            logger.warning("The bottom and top structure seem to be identical.")
            logger.warning(
                "Removing the angle 0° from the search because all lattice points would match."
            )
            angles = [k for k in angles if abs(k) > 1e-4]
            assert len(angles) > 0, "List of angles contains no values."

        assert Nmin < Nmax, "Nmin must be smaller than Nmax."
        assert Nmin >= 0, "Nmin must be larger than or equal 0."
        assert Nmax > 0, "Nmax must be larger than 0."
        assert distance > 0, "Interlayer distance must be larger than zero."
#        assert vacuum >= 10, "Thickness of the vacuum layer must be larger than or equal to 10 Anström."
        assert tolerance > 0, "Tolerance must be larger than zero."
        assert (
            angle_tolerance >= 0
        ), "Angle tolerance must be larger than or equal zero."
        assert (symprec) > 0, "Symmetry precision must be larger than zero."
        assert (weight >= 0) and (weight <= 1), "Weight factor must be between 0 and 1."
        assert verbosity in [0, 1, 2], "Verbose must be 0, 1, or 2."

        angles = double1dVector(angles)
        no_idealize = int(no_idealize)

        ncombinations = ((2 * (Nmax - Nmin)) ** 4) * len(angles)

        nthreads = get_number_of_omp_threads()
        logger.info("Using {:d} OpenMP threads.".format(nthreads))
        logger.info("Running through {:d} grid points...".format(ncombinations))
        alg = CppCoincidenceAlgorithmClass(bottom, top)
        results = alg.run(
            Nmax,
            Nmin,
            angles,
            tolerance,
            weight,
            distance,
            vacuum,
            standardize,
            no_idealize,
            symprec,
            angle_tolerance,
            verbosity,
        )
        if len(results) == 0:
            logger.error("Could not find any coincidence pairs for these parameters.")
            return None
        elif len(results) > 0:
            if len(results) == 1:
                logger.info("Found 1 result.")
            else:
                logger.info("Found {:d} results.".format(len(results)))

            interfaces = [
                Interface(
                    k,
                    weight=weight,
                    bottom_bond_lengths=self.bdl,
                    top_bond_lengths=self.tbl,
                )
                for k in results
            ]
            return interfaces
