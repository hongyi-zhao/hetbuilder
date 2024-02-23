#include <cmath>
#include "math_functions.h"
#include "logging_functions.h"
#include "atom_class.h"

#include "atom_functions.h"
using std::sin, std::cos, std::sqrt, std::pow, std::abs;

/**
 * Find all lattice points contained in a supercell.

    Adapted from pymatgen, which is available under MIT license:
    The MIT License (MIT) Copyright (c) 2011-2012 MIT & The Regents of the
    University of California, through Lawrence Berkeley National Laboratory
 */
double2dvec_t lattice_points_in_supercell(int2dvec_t &superCellMatrix)
{
    int2dvec_t diagonals = {{0, 0, 0},
                            {0, 0, 1},
                            {0, 1, 0},
                            {0, 1, 1},
                            {1, 0, 0},
                            {1, 0, 1},
                            {1, 1, 0},
                            {1, 1, 1}};
    int2dvec_t dpoints;
    int1dvec_t dotproduct;
    for (int row = 0; row < diagonals.size(); row++)
    {
        dotproduct = vec1x3_dot_3x3_matrix<int, int>(diagonals[row], superCellMatrix);
        dpoints.push_back(dotproduct);
    }

    int1dvec_t mins = {0, 0, 0};
    int1dvec_t maxes = {0, 0, 0};

    for (int row = 0; row < dpoints.size(); row++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (dpoints[row][j] < mins[j])
            {
                mins[j] = dpoints[row][j];
            }
            if (dpoints[row][j] > maxes[j])
            {
                maxes[j] = dpoints[row][j];
            }
        }
    }
    maxes = {maxes[0] + 1, maxes[1] + 1, maxes[2] + 1};

    int2dvec_t ar, br, cr;
    int1dvec_t subvec(3, 0);

    for (int a = mins[0]; a < maxes[0]; a++)
    {
        subvec = {a, 0, 0};
        ar.push_back(subvec);
    }
    for (int b = mins[1]; b < maxes[1]; b++)
    {
        subvec = {0, b, 0};
        br.push_back(subvec);
    }
    for (int c = mins[2]; c < maxes[2]; c++)
    {
        subvec = {0, 0, c};
        cr.push_back(subvec);
    }

    int2dvec_t allpoints;
    for (int i = 0; i < ar.size(); i++)
    {
        for (int j = 0; j < br.size(); j++)
        {
            for (int k = 0; k < cr.size(); k++)
            {
                subvec[0] = ar[i][0] + br[j][0] + cr[k][0];
                subvec[1] = ar[i][1] + br[j][1] + cr[k][1];
                subvec[2] = ar[i][2] + br[j][2] + cr[k][2];
                allpoints.push_back(subvec);
            }
        }
    }

    // convert integer matrix to doubles
    double2dvec_t allpoints_double;
    for (int row = 0; row < allpoints.size(); row++)
    {
        double1dvec_t doubleVec(allpoints[row].begin(), allpoints[row].end());
        allpoints_double.push_back(doubleVec);
    };

    double determinant = get_3x3_matrix_determinant<int>(superCellMatrix);
    double2dvec_t invSuperCellMatrix = invert_3x3_matrix<int>(superCellMatrix);
    double2dvec_t fracpoints;
    std::vector<double> dp;
    for (int row = 0; row < allpoints.size(); row++)
    {
        dp = vec1x3_dot_3x3_matrix<double, double>(allpoints_double[row], invSuperCellMatrix);
        fracpoints.push_back(dp);
    }

    double2dvec_t tvects;
    double fa, fb, fc;
    double1dvec_t fvec;
    for (int row = 0; row < fracpoints.size(); row++)
    {
        fa = fracpoints[row][0];
        fb = fracpoints[row][1];
        fc = fracpoints[row][2];
        if ((fa <= (1 - 1e-10) && (fa >= (-1e-10))) && (fb <= (1 - 1e-10) && (fb >= (-1e-10))) && (fc <= (1 - 1e-10) && (fc >= (-1e-10))))
        {
            fvec = {fa, fb, fc};
            tvects.push_back(fvec);
        }
    }
    try
    {
        int detsize = (int)determinant;
        if (detsize != tvects.size())
        {
            throw "Determinant of supercell does not match number of lattice points.";
        }
    }
    catch (const char *msg)
    {
        std::cout << msg << std::endl;
        tvects = {};
    }

    return tvects;
};

/**
 * Generate a supercell by applying a SuperCellMatrix to
    the input atomic configuration prim.

    Indices of the supercell atom map to the indices of the primitive cell for later use.
*/
Atoms make_supercell(Atoms &prim, int2dvec_t &superCellMatrix)
{
    double2dvec_t fracpoints = lattice_points_in_supercell(superCellMatrix);
    double2dvec_t cell = prim.lattice;
    double2dvec_t supercell = matrix3x3_dot_matrix3x3<int, double>(superCellMatrix, cell);

    double2dvec_t lattice_points;
    double1dvec_t dotproduct;
    for (int row = 0; row < fracpoints.size(); row++)
    {
        dotproduct = vec1x3_dot_3x3_matrix<double, double>(fracpoints[row], supercell);
        lattice_points.push_back(dotproduct);
    }

    double2dvec_t new_positions = prim.positions;

    int1dvec_t new_numbers = prim.atomic_numbers;
    int1dvec_t index_mapping = prim.indices;
    double1dvec_t new_magmoms = prim.magmoms;
    for (int i = 0; i < prim.numAtom; i++)
    {
        double1dvec_t atom_pos = prim.positions[i];
        int number = prim.atomic_numbers[i];
        int index = prim.indices[i];
        double magmom = prim.magmoms[i];
        for (int row = 0; row < lattice_points.size(); row++)
        {
            double1dvec_t lp = lattice_points[row];
            if ((std::abs(lp[0]) + std::abs(lp[1]) + std::abs(lp[2])) > 1e-6)
            {
                for (int k = 0; k < 3; k++)
                {
                    lp[k] += atom_pos[k];
                }
                new_positions.push_back(lp);
                new_numbers.push_back(number);
                index_mapping.push_back(index);
                new_magmoms.push_back(magmom);
            }
        }
    }
    Atoms superatoms = {supercell, new_positions, new_numbers, index_mapping, new_magmoms};
    return superatoms;
};

// Rotates Atoms object around z-axis for given angle theta in degrees.
Atoms rotate_atoms_around_z(Atoms &atoms, double &theta)
{
    double t = M_PI * theta / 180.0;
    double c = std::cos(t);
    double s = std::sin(t);
    double2dvec_t R = {{c, -s, 0}, {s, c, 0}, {0, 0, 1}};

    double2dvec_t positions = atoms.positions;
    double2dvec_t rotPositions;
    for (int row = 0; row < positions.size(); row++)
    {
        double1dvec_t vec = positions[row];
        double1dvec_t rotvec = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                rotvec[i] += (R[i][j] * vec[j]);
            }
        }
        rotPositions.push_back(rotvec);
    }

    double2dvec_t tCell = transpose_matrix3x3<double>(atoms.lattice);
    double2dvec_t rotCellTranspose = matrix3x3_dot_matrix3x3(R, tCell);
    double2dvec_t rotCell = transpose_matrix3x3<double>(rotCellTranspose);

    Atoms rotatoms(rotCell, rotPositions, atoms.atomic_numbers, atoms.indices, atoms.magmoms);
    return rotatoms;
};

// Translates Atoms object along z-direction for given shift.
void translate_atoms_z(Atoms &atoms, double shift)
{
    double2dvec_t pos1 = atoms.positions;
    double2dvec_t new_pos;
    for (int row = 0; row < atoms.numAtom; row++)
    {
        double1dvec_t subvec = {pos1[row][0],
                                pos1[row][1],
                                pos1[row][2] + shift};
        new_pos.push_back(subvec);
    };
    atoms.positions = new_pos;
};

// Helper function to get minimum and maximum z values of positions.
std::tuple<double, double> get_min_max_z(Atoms &atoms)
{
    double min_z = atoms.positions[0][2];
    double max_z = atoms.positions[atoms.numAtom - 1][2];
    for (int row = 0; row < atoms.numAtom; row++)
    {
        double z = atoms.positions[row][2];
        if (z < min_z)
        {
            min_z = z;
        }
        if (z > max_z)
        {
            max_z = z;
        }
    }
    return std::make_tuple(min_z, max_z);
};



// 进一步简化的考虑：  
// 当只考虑xy平面进行操作时，且希望采用类似 latticeA + weight * (latticeB - latticeA) 的简化表达式进行矩阵运算，我们可以通过将z轴方向的组件保持不变来实现这一点。这意味着你只对矩阵的前两行（代表xy平面的向量）进行操作，而第三行（代表z轴的向量）保持原样不变。

//在NumPy中的实现：

//如果你使用Python和NumPy，可以通过选择性地修改数组的一部分来实现这一点：

//import numpy as np

//# 假设 latticeA 和 latticeB 是两个3x3的NumPy数组
//latticeA = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
//latticeB = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]])

//weight = 0.5  # 插值权重

//# 只对前两行进行操作
//newcell = latticeA.copy()  # 创建一个副本以保留原始数据
//newcell[:2] = latticeA[:2] + weight * (latticeB[:2] - latticeA[:2])

//print(newcell)

// 
//在Eigen(C++)中的实现：

//在C++中，如果你使用Eigen库，可以通过对矩阵的特定行进行操作来实现相似的逻辑：

//#include <Eigen/Dense>
//#include <iostream>

//int main() {
//    Eigen::Matrix3d latticeA;
//    Eigen::Matrix3d latticeB;

//    // 填充latticeA和latticeB为示例数据
//    latticeA << 1, 2, 3, 4, 5, 6, 7, 8, 9;
//    latticeB << 9, 8, 7, 6, 5, 4, 3, 2, 1;

//    double weight = 0.5; // 插值权重

//    // 创建一个新的矩阵以存储结果，初始值为latticeA
//    Eigen::Matrix3d newcell = latticeA;

//    // 只修改前两行，对应于xy平面
//    newcell.block<2,3>(0,0) = latticeA.block<2,3>(0,0) + weight * (latticeB.block<2,3>(0,0) - latticeA.block<2,3>(0,0));

//    std::cout << "Newcell matrix:\n" << newcell << std::endl;

//    return 0;
//}

// 

//在这两种情况下，我们都通过仅对代表xy平面的矩阵部分进行操作，同时保留z轴方向（即第三行/列）的数据不变，来实现了只考虑xy平面的矩阵运算。这种方法允许你在保持z轴信息不变的同时，灵活地调整xy平面的参数。
// 

/**
 * Stacks two Atoms objects on top of each other with an interlayer distance given by distance.
 *
 * Returns a new Atoms object.
 *
 * The new unit cell is given by C = A + weight * (B - A).
 */
Atoms stack_atoms(Atoms bottom, Atoms top, double &weight, double &distance, double &vacuum)
{
 
// In 2D heterostructure interfaces, these scalings occur on the xy plane, thus necessitating the use of the logic implemented in the scale_cell_xy function to scale while keeping the z-direction Cartesian coordinates unchanged.
    
    // Adjusting z positions based on the thickness of the bottom and top layers
    auto [min_z1, max_z1] = get_min_max_z(bottom);
    auto [min_z2, max_z2] = get_min_max_z(top);
    translate_atoms_z(bottom, -min_z1);
    double bottom_thickness = max_z1 - min_z1;
    double top_thickness = max_z2 - min_z2;
    double shift = -min_z2 + bottom_thickness + distance;
    translate_atoms_z(top, shift);

    // Interpolating lattice vectors for the combined structure
    double2dvec_t latticeA = bottom.lattice;
    double2dvec_t latticeB = top.lattice;
    double2dvec_t newcell(3, std::vector<double>(3, 0));
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            newcell[i][j] = latticeA[i][j] + weight * (latticeB[i][j] - latticeA[i][j]);
        }
    }
    bottom.scale_cell_xy(newcell);
    top.scale_cell_xy(newcell);
    
    // Stacking the bottom and top layers
    Atoms stack = bottom + top;
    stack.lattice[2][2] = bottom_thickness + top_thickness + distance + vacuum;
    return stack;
};

