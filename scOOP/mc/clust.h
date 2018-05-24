#ifndef CLUST_H
#define CLUST_H

#include <algorithm>

#include "../structures/Conf.h"
#include "../structures/sim.h"
#include "totalenergycalculator.h"
#include "../structures/structures.h"

typedef struct{
    std::vector<int> particles;
    double
        energy = 0.0,
        nematic = 0.0,
        radiusOfGyration = 0.0;
    Vector nematicVector = {0,0,0};
} Cluster;


class ClusterSampler {
public:
    Conf* conf;
    Sim* sim;
    // actualy we just need options object

    TotalEnergyCalculator* calcEnergy;
    FileNames* files;

    // IDEA is to have vector of clusters where base on options in option file different
    // cluster definitions are used to construct clusters ... then diferent analysis
    // being done on cluster and then print out in different outup formats


    long * clusterstat;         ///< \brief Statistics about the size of cluster

    std::vector<Cluster> clusters; ///< \brief informations about the single clusters
    // what types of analysis are suppose to be calculated and included in output
    bool
        nematic = false,
        sizeDistribution = false,
        radiusOfGyration = false;
    // how to define cluster
    int clusterDefinition = 0; // if 0 energy based if 1 distance
    double clusterCutoff = -0.1; // cutoff used in cluster definition ither distance or energy


    long * clusterlist;         ///< \brief clusterlist[i] = cluster index of particle i
    long num_cluster;           ///< \brief number of single clusters
    unsigned long max_clust;             ///< \brief maximal clustersize

    ClusterSampler(Conf* conf, Sim* sim, TotalEnergyCalculator* calcEnergy, FileNames* files) : conf(conf),
    sim(sim), calcEnergy(calcEnergy), files(files) {

        // initialize statics base on options in options
        clusters.reserve(MAXN);
        if (clusterDefinition == 1) { // internaly we use only square of distance due to avoiding sqrt()
            clusterCutoff = clusterCutoff * clusterCutoff;
        }

        clusterstat = (long int*) malloc(sizeof(long) * max_clust);
        clusterlist = (long int*) malloc(sizeof(long) * MAXN);
        if(clusterlist == NULL){
            fprintf(stderr, "\nTOPOLOGY ERROR: Could not allocate memory for clusterlist!");
            exit(1);
        }

    }

    ~ClusterSampler() {
        if (clusterstat != NULL)
            free(clusterstat);

        if (clusterlist != NULL)
            free(clusterlist);
    }

    /**
     * @brief write_cluster write out all the cluster stat in files, if file name is given
     * @param cl_stat
     * @param cl
     * @param cl_list
     * @param decor
     * @param sweep
     * @return
     */
    int writeCluster(bool decor, long sweep);

private:

    /**
     * @brief same_cluster determines, wheter two particles are in the same cluster
     * @param fst
     * @param snd
     * @return
     */
    int sameCluster(long fst, long snd);

    /**
     * @brief gen_clusterlist generate the clusterlist
     * @return
     */
    int genClusterList();

    /**
     * @brief sort_clusterlist sort the clusterlist
     * @return
     */
    int sortClusterList();

    void clusterMovie();
    void computeEnergy();
    void computeClEnergy(Cluster *cl);

    void computeNematic();
    void computeClNematic(Cluster *cl);

    void computeRadiusOfGyration();
    void computeClRadiusOfGyration(Cluster *cl);

    int printClusterList(FILE *stream, bool decor);


    int printClusters(FILE *stream, bool decor);


    int printClusterStat(FILE *stream, bool decor);
};


#endif // CLUST_H
