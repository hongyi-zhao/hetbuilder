#pragma once
#include "atom_class.h"

/**
 * Container class for an interface after the coincidence lattice search.
 * 
 * Contains the bottom layer and top layer in supercell form, as well as the
 * stack thereof.
 * 
 * Comparison operators are defined to allow for ordering and comparisons in a set.
 */
class Interface
{
public:
    Atoms bottomLayer;
    Atoms topLayer;
    Atoms stack;
    double angle;
    int2dvec_t M;
    int2dvec_t N;
    int spaceGroup;

    friend bool operator==(const Interface &c1, const Interface &c2);
    friend bool operator>(const Interface &c1, const Interface &c2);
    friend bool operator>=(const Interface &c1, const Interface &c2);
    friend bool operator<(const Interface &c1, const Interface &c2);
    friend bool operator<=(const Interface &c1, const Interface &c2);

    Interface(Atoms cBottomLayer,
              Atoms cTopLayer,
              Atoms cStack,
              double cAngle,
              int2dvec_t cMatrixM,
              int2dvec_t cMatrixN,
              int cSpaceGroup)
    {
        bottomLayer = cBottomLayer;
        topLayer = cTopLayer;
        stack = cStack;
        angle = cAngle;
        M = cMatrixM;
        N = cMatrixN;
        spaceGroup = cSpaceGroup;
    };
};

inline bool operator==(const Interface &c1, const Interface &c2)
{
    bool spgmatch = (c1.spaceGroup == c2.spaceGroup);
    bool nummatch = (c1.stack.numAtom == c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(std::abs(area1) - std::abs(area2)) < 1e-6;
    bool equals = (spgmatch && nummatch && areamatch);
    return equals;
}

inline bool operator>(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom > c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) > std::abs(area2);
    bool greater_than = (nummatch && areamatch);
    return greater_than;
}

inline bool operator>=(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom >= c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) >= std::abs(area2);
    bool greater_than = (nummatch && areamatch);
    return greater_than;
}

inline bool operator<(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom < c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) < std::abs(area2);
    bool less_than = (nummatch && areamatch);
    return less_than;
}

inline bool operator<=(const Interface &c1, const Interface &c2)
{
    bool nummatch = (c1.stack.numAtom <= c2.stack.numAtom);
    double area1 = c1.stack.lattice[0][0] * c1.stack.lattice[1][1] - c1.stack.lattice[0][1] * c1.stack.lattice[1][0];
    double area2 = c2.stack.lattice[0][0] * c2.stack.lattice[1][1] - c2.stack.lattice[0][1] * c2.stack.lattice[1][0];
    bool areamatch = std::abs(area1) <= std::abs(area2);
    bool less_than = (nummatch && areamatch);
    return less_than;
}