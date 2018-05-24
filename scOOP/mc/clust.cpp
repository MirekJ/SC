#include "clust.h"


int ClusterSampler::writeCluster(bool decor, long sweep) {

    FILE *cl_stat = NULL;
    FILE *cl = NULL;

    cl_stat = fopen(files->clusterstatfile, "a");
    cl = fopen(files->clusterfile, "a");

    if (conf->pvec.empty()) { // well if we do not have particles in system (GC removed all) no need for sampling
        fprintf(cl_stat, "Sweep: %ld | Maximal size: %i\n", sweep, 0);
    } else {
        genClusterList();
        sortClusterList();
        computeEnergy();
        computeNematic();
        computeRadiusOfGyration();
        clusterMovie();
        if(cl_stat){
            // if no decor, this means usually into a file. Hence print info
            // about number of line per frame
            fprintf(cl_stat, "Sweep: %ld | Maximal size: %ld\n", sweep, max_clust);
            printClusterStat(cl_stat, decor);
        }
        if(cl){
            fprintf(cl, "Sweep: %ld | Number of clusters: %ld\n",
                sweep, num_cluster);
            printClusters(cl, decor);
        }
    }
    fclose(cl_stat);
    fclose(cl);

    return 0;
}


int ClusterSampler::sameCluster(long fst, long snd) {
    /*Here we add particle into cluster if particles fst and snd are from same molecule*/
    Molecule fstMol = conf->pvec.getMolOfPart(fst);
    if( std::find(fstMol.begin(), fstMol.end(), snd) != fstMol.end() ){
        return true;
    }

    switch (clusterDefinition) {
    case 0: /*cluster is made of attractively interacting particles*/
        // p2p ... include also conlist contribution ... which might not be the best for
        // describing some clusters ... should be without it???
        if(calcEnergy->p2p(fst, snd) < clusterCutoff ) {
            return true;
        } else {
            return false;
        }
        break;
    case 1: /*cluster is made of particles closer than some distance*/
        if ((conf->geo.image(&conf->pvec[fst].pos, &conf->pvec[snd].pos)).sizeSq() < clusterCutoff) {
            return true;
        } else {
            return false;
        }
        break;
    default:
            cerr << "Invalid cluster selection method! (" << clusterDefinition << ")" << endl;
            exit(1);
        break;
    }
}

int ClusterSampler::printClusterList(FILE *stream, bool decor) {
    if(decor){
        fprintf(stream, "\n"
                "-----------------------------------------------------\n"
                "  The Cluster List\n"
                "  (Index starts with 1)\n"
                "-----------------------------------------------------\n");
    }

    for(int i=0; i < (long)conf->pvec.size(); i++){
        fprintf(stream,"%3d %3ld %8.4lf %8.4f %8.4f", i + 1,
                clusterlist[i] + 1,
                conf->pvec[i].pos.x,
                conf->pvec[i].pos.y,
                conf->pvec[i].pos.z);
        fprintf(stream,"\n");
    }
    if(decor){
        fprintf(stream,"-----------------------------------------------------\n");
    }
    fflush(stream);
    return 0;
}

int ClusterSampler::printClusters(FILE *stream, bool decor) {
    if(decor){
        fprintf(stream, "\n"
                "-----------------------------------------------------\n"
                "  The Clusters\n"
                "  (Index starts with 1)\n"
                "-----------------------------------------------------\n");
    }
    for(int i = 0; i < num_cluster; i++){
        fprintf(stream, "%3d ( %+8.6f )[ %lu ]:", i + 1,clusters[i].energy,clusters[i].particles.size());
        for(unsigned int j = 0; j < clusters[i].particles.size(); j++){
            fprintf(stream, "%5i[%i]", clusters[i].particles[j] + 1, conf->pvec[clusters[i].particles[j]].molType);
        }
        fprintf(stream, "\t%lf\t%lf\n", clusters[i].nematic, clusters[i].radiusOfGyration);
    }
    if(decor){
        fprintf(stream,"---------------------------------------------------\n");
    }
    fflush(stream);
    return 0;
}

int ClusterSampler::printClusterStat(FILE *stream, bool decor) {
    if(decor){
        fprintf(stream, "\n"
                "-----------------------------------------------------\n"
                "   Cluster Distribution\n"
                "-----------------------------------------------------\n");
    }
    for(unsigned int i=0; i < max_clust; i++){
        fprintf(stream, "%5d\t%5ld\n", i + 1, clusterstat[i]);
    }
    if(decor){
        fprintf(stream, "--------------------------------------------------\n");
    }
    fflush(stream);
    return 0;
}

int ClusterSampler::genClusterList() {
    bool change = true; /* does it still change? */
    //long neighbour;
    long i, j, fst, snd, tmp, minnumber, maxnumber;

    // Set clusterindex to the corresponding index
    for( i = 0; i < (long)conf->pvec.size(); i++){
        clusterlist[i] = i;
    }

    // Start determining the cluster
    while(change){
        change = false;
        for(i = 0; i < (long)conf->pvec.size(); i++){
            /*If nore pairlist go over all pairs*/
            maxnumber = (long)conf->pvec.size();
            minnumber = i ;
            if (sim->pairlist_update) {
                maxnumber = conf->neighborList[i].neighborCount;
                //maxnumber = sim->pairlist[i].num_pairs; // del after
                minnumber=0;
            }
            /* Go over pairs to see if they are in the cluster */
            for(j = minnumber; j < maxnumber; j++){
                fst = i;
                snd = j;
                if (sim->pairlist_update) {
                    snd = conf->neighborList[i].neighborID[j];
                    //snd = sim->pairlist[i].pairs[j];
                }
                /*do cluster analysis only for spherocylinders*/
                /*
                if ( (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[0] < SP) && \
                   (topo.ia_params[conf->pvec[fst].type][conf->pvec[snd].type].geotype[1] < SP) ) {
                */
                /* if they are close to each other */
                if(sameCluster(fst, snd)){
                    if(fst > snd){
                        tmp = snd;
                        snd = fst;
                        fst = tmp;
                    }

                    if(clusterlist[fst] < clusterlist[snd]){
                        clusterlist[snd] = clusterlist[fst];
                        change = true;
                        break;
                        /* => will eventually start the i loop from new */
                    }
                    if(clusterlist[snd] < clusterlist[fst]){
                        clusterlist[fst] = clusterlist[snd];
                        change = true;
                        break;
                        /* => will eventually start the i loop from new */
                    }
                }
                //}
            }
            if(change){
                break;
            }
        }
    }

    return 0;
}

int ClusterSampler::sortClusterList() {
    long cluster_indices[(long)conf->pvec.size()];   /* holds the different cluster indices.
                        (currently too much memory) */
    long num_cluster = 0;                /* number of clusters, temporary needed */

    /* how many clusters are there? */
    long max_index = -1;
    for(int i = 0; i < (long)conf->pvec.size(); i++){
        if(max_index < clusterlist[i]){
            max_index = clusterlist[i];
            cluster_indices[num_cluster++] = max_index;
        }
    }

    clusters.clear();
    clusters.resize(num_cluster);

    /* fill in the particles belonging to one cluster */
    for(int i = 0; i < num_cluster; i++){
        for(int j = 0; j < (long)conf->pvec.size(); j++){
            if(clusterlist[j] == cluster_indices[i]){
                clusters[i].particles.push_back(j);
            }
        }
    }
    this->num_cluster = num_cluster;

    /* Find the biggest size */
    max_clust = 0;
    for(int i = 0; i < num_cluster; i++){
        if(clusters[i].particles.size() > max_clust){
            max_clust = clusters[i].particles.size();
        }
    }
    /* Set the statistics to zero */
    clusterstat = (long int*) realloc( clusterstat, sizeof(long) * max_clust); // OLD, MISTAKE? memmory dont have to be 0, no free
    memset(clusterstat, 0, sizeof(long) * max_clust);

    if (!clusterstat){
        fprintf(stderr, "Couldn't allocate any memory!\n");
        exit(1);
    }
    for(unsigned int i = 0; i < max_clust; i++) {
        clusterstat[i] = 0;
    }
    /* Do the statistics */
    for(int i = 0; i < num_cluster; i++){
        clusterstat[clusters[i].particles.size() - 1]++;
    }
    return 0;
}

void ClusterSampler::clusterMovie() {
    // we create array of files coresponding to clusters of different sizes
    std::vector<FILE *> clFiles = {NULL};
    char tmp[50];
    for (unsigned int i=1; i<=max_clust; i++) {
        sprintf(tmp, "%s_%i", files->clMovieFile, i);
        clFiles.push_back(fopen(tmp, "a"));
    }
    // now we draw clusters in each movie of cluster
    for (std::vector<Cluster>::iterator it = clusters.begin(); it != clusters.end(); ++it) {
        fprintf (clFiles[(*it).particles.size()],
                "%lu\nsweep %i; box %.10f %.10f %.10f\n",
                (*it).particles.size(),
                1,
                conf->geo.box.x,
                conf->geo.box.y,
                conf->geo.box.z
                );
        conf->draw(clFiles[(*it).particles.size()], (*it).particles);
    }
    for (unsigned int i=1; i<=max_clust; i++) {
        fclose(clFiles.back());
        clFiles.pop_back();
    }
}

void ClusterSampler::computeEnergy() {
    for (std::vector<Cluster>::iterator it = clusters.begin(); it != clusters.end(); ++it) {
        computeClEnergy(&(*it));
    }
}

void ClusterSampler::computeClEnergy(Cluster *cl) {
    cl->energy = 0.0;
    if (cl->particles.size() < 2)
        return;

    for (unsigned int i = 0; i < cl->particles.size();i++) {
        for (unsigned int j = i+1; j < cl->particles.size(); j++) {
            cl->energy += calcEnergy->p2p(cl->particles[i],cl->particles[j]);
        }
    }
}


void ClusterSampler::computeNematic() {
    for (std::vector<Cluster>::iterator it = clusters.begin(); it != clusters.end(); ++it) {
        computeClNematic(&(*it));
    }
}

void ClusterSampler::computeClNematic(Cluster *cl) {
    cl->nematic = 0.0;
    cl->nematicVector = {0,0,0};
    for (unsigned int i = 0; i < cl->particles.size(); i++){
        cl->nematicVector += conf->pvec[cl->particles[i]].dir;
    }
    cl->nematic = cl->nematicVector.size();
}

void ClusterSampler::computeRadiusOfGyration() {
    for (std::vector<Cluster>::iterator it = clusters.begin(); it != clusters.end(); ++it) {
        computeClRadiusOfGyration(&(*it));
    }
}

void ClusterSampler::computeClRadiusOfGyration(Cluster *cl) {
    cl->radiusOfGyration = 0;

    Vector
            clusterCM = conf->clusterCM(cl->particles),
            tmp(0.0, 0.0, 0.0);

    for (unsigned int i = 0; i < cl->particles.size(); i++){
        tmp = conf->pvec[cl->particles[i]].pos - clusterCM;
        cl->radiusOfGyration += tmp.size();
    }
}

