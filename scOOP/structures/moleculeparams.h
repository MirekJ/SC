#ifndef MOLECULEPARAMS_H
#define MOLECULEPARAMS_H

#include "macros.h"
#include <cassert>

/**
 * @brief Parameters for interactions of molecules(mainly inner interaction in chains), each molType has a MoleculeParams instance
 */
class MoleculeParams{
    class Bond{
        double bondEqDist;
        double bondEqConst;
    public:
        Bond() :
            bondEqDist( 0.0 ),
            bondEqConst( 0.0 ) // for bondEqConst == 0.0 basicly mean there is no bond at all
        {}

        Bond( double bondEqDist, double bondEqConst /*type enabling use of bondd bond1 etc.*/ ) :
            bondEqDist( bondEqDist ),
            bondEqConst( bondEqConst )
        {}

    };

public:
    char * name;        ///< \brief name of the molecule
    double bondsEq[ MAXCHL-1 ];
    double bondsEqConst[ MAXCHL-1 ];
    double bond1eq;     ///< \brief Equilibrium distance of harmonic bond between nearest neighbours
    double bond1c;      ///< \brief Spring constant for harmonic bond between nearest neighbours
    double bond2eq;     ///< \brief Equilibrium distance of harmonic bond between second nearest neighbours
    double bond2c;      ///< \brief Spring constant for harmonic bond between second nearest neighbours
    double bonddeq;     ///< \brief Equilibrium distance of directional harmonic bond between the nearest neighbours
    double bonddc;      ///< \brief Spring constant for directional harmonic bond between the nearest neighbours
    double angle1eq;    ///< \brief Equilibrium angle between two spherocylinders -neerest neighbours
    double angle1c;     ///< \brief Spring constant angle between two spherocylinders -nearest neighbours
    double angle2eq;    ///< \brief Equilibrium angle between two spherocylinder patches -nearest neighbours
    double angle2c;     ///< \brief Spring constant for angle between two spherocylinder patches -nearest neighbours

    double activity;                ///< \brief activity, specific to each molecule type
    double chemPot;                 ///< \brief mu + 3*ln(A), specific to each type
    std::vector<int> particleTypes;    ///< \brief 0..40 particle type
    std::vector<int> switchTypes;
    std::vector<double> deltaMu;
    int switchCount;        ///< \brief count of particles with defined switchtype

public:
    MoleculeParams() : name(NULL), bond1eq(-1.0), bond1c(-1.0), bond2eq(-1.0), bond2c(-1.0), /*bonddeq(-1),*/ bonddc(-1.0),
                       angle1eq(-1.0), angle1c(-1.0), angle2eq(-1.0), angle2c(-1.0),
                       activity(-1.0), chemPot(-1.0) {
        particleTypes.reserve(MAXCHL);
    }

    /**
     * @return True if Molecule of only one particle
     */
    bool isAtomic() {return (particleTypes.size() == 1);}

    bool isGrandCanonical() {return chemPot!=-1.0;}

    /**
     * @brief molSize
     * @return Number of particles in molecule
     */
    inline int molSize() {
        return particleTypes.size();
    }

    inline int switchTypeSeq(int seq) {
        for(int i=0; i<(int)switchTypes.size(); i++) {
            if(switchTypes[i] != -1)
                seq--;
            if(seq == -1)
                return i;
        }
        assert(false && "This should never happen!!!");
        return -1;
    }

    void info() {
        cout << name << ", activity: " << activity << ", atoms:\n";
        for(unsigned int i=0; i<particleTypes.size(); i++)
            cout << particleTypes[i] << ", ";
        cout << endl;
    }
};

#endif // MOLECULEPARAMS_H
