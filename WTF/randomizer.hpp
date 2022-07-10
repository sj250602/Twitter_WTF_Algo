#ifndef __RANDOMIZER
#define __RANDOMIZER

#include <random>
#include <iostream>
#include <bits/stdc++.h>
#include <thread>
#include <climits>
#include <ctime>

// It takes number_of_vertices and restart probability
class Randomizer{
    int restart_prob;
    int num_vertices;
    unsigned int *random_array, *random_offset, *reset_array;

    public:
        //Takes the seed, number of nodes, and restart probability
        Randomizer(int seed, int num_nodes, float prob);

        //returns -1 in case of restart
        //otherwise returns a random number in the range (0, UINT_MAX)
        //You can it with modulo to get a meaningful number.
        //user_id_in_process should be the node_id/vertex_id/user_id 
        //for whom you are running the RWR algorithm
        //i.e. if you're finding the recommandation for node 0, you should always pass 0
        //in this function. Even if you're at a random node(say 30) in the walk.
        //At each step of random walk this function should be called only once.
        //And there it will decide whether to restart or take next_step based on the value returned.
        int get_random_value(int user_id_in_process);
};

#endif