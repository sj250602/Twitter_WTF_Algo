#include <string>
#include <mpi.h>
#include <assert.h>
#include <bits/stdc++.h>
#include "randomizer.hpp"
#include<chrono>

using namespace std::chrono;
using namespace std;

bool compr(const pair<int,int> &a,const pair<int,int> &b){
    if(a.second>b.second){
        return true;
    }else if(a.second==b.second){
        if(a.first<b.first){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}

int main(int argc, char* argv[]){

    auto start = high_resolution_clock::now();
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL); std::cout.tie(NULL);
    assert(argc > 8);
    std::string graph_file = argv[1];
    int num_nodes = std::stoi(argv[2]);
    int num_edges = std::stoi(argv[3]);
    float restart_prob = std::stof(argv[4]);
    int num_steps = std::stoi(argv[5]);
    int num_walks = std::stoi(argv[6]);
    int num_rec = std::stoi(argv[7]);
    int seed = std::stoi(argv[8]);

    vector<int> adj[num_nodes];

    FILE *fp;
    unsigned char a[4];
    unsigned char b[4];
    int a1,b1;
    const char* str = graph_file.c_str();
    fp = fopen(str,"rb");
    if(fp==NULL){
        cout<<"error\n";
    }else{
        while(!feof(fp)){
            fread(&a,4,1,fp);
            fread(&b,4,1,fp);
            a1 = (int)a[3] | ((int)a[2]<<8) | ((int)a[1]<<16) | ((int)a[0]<<24);
            b1 = (int)b[3] | ((int)b[2]<<8) | ((int)b[1]<<16) | ((int)b[0]<<24);
            adj[a1].push_back(b1);
            // cout<<a1<<" "<<b1<<"\n";
        }
    }
    fclose(fp);
    adj[a1].pop_back();

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);

    vector<map<int,int>> circle(num_nodes);
    
    vector<int> rec_count(num_nodes);

    //Only one randomizer object should be used per MPI rank, and all should have same seed
    Randomizer random_generator(seed, num_nodes, restart_prob);
    int rank, size;
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<pair<int,int>> thread_time(size);
    int ed = num_nodes/size;
    for(int i=0;i<size;i++){
        if(i==size-1){
            thread_time[i] = {i*ed,num_nodes-1};
        }else{
            thread_time[i] = {i*ed,((i+1)*ed)-1};
        }
    }

    int strt = thread_time[rank].first,end = thread_time[rank].second;
    // for(int i=strt;i<=end;i++){
    //     for(int j=0;j<num_nodes;j++){
    //         circle[i][j] = {j,0};
    //     }
    // }

   // print_random(rank, num_nodes, random_generator,adj,num_steps,circle,num_walks,rec_count);
    for(int i=strt;i<=end;i++){
        vector<int> L = adj[i];
        for(int j=0;j<L.size();j++){ 
            for(int k=0;k<num_walks;k++){
                int node= L[j];
                for(int l=0;l<num_steps;l++){    
                    if(adj[node].size()==0){
                        node = L[j];
                    }else{
                        int w = random_generator.get_random_value(i);
                        if(w<0){
                            node = L[j];
                        }else{
                            w = w % adj[node].size();
                            node  = adj[node][w];
                            if(node!=i && find(L.begin(),L.end(),node)==L.end()){
                                circle[i][node]++;
                                if(circle[i][node]==1){
                                    rec_count[i]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // for(int i=strt;i<=end;i++){
    //     sort(circle[i].begin(),circle[i].end(),compr);
    // }

    if(rank!=0){
        for(int i=strt;i<=end;i++){

            vector<pair<int,int>> vect(circle[i].size());
            auto it = circle[i].begin();
            for(int k=0;k<vect.size();k++){
                vect[k] = *it;
                circle[i].erase(circle[i].begin());
                it = circle[i].begin();
            }

            sort(vect.begin(),vect.end(),compr);

           // cout<<vect.size()<<"\n";

            int send_buff[num_rec*2+1];
            send_buff[0] = rec_count[i];
            int l=1;
            if(vect.size()>=num_rec){
                for(int k=0;k<num_rec;k++){
                    send_buff[l++] = vect[k].first;
                    send_buff[l++] = vect[k].second;
                    //cout<<vect[k].first<<" "<<vect[k].second<<" ";
                }
                MPI_Send(&send_buff,2*num_rec+1,MPI_INT,0,rank,MPI_COMM_WORLD);
            }else{
                int k =0;
                for(;k<vect.size();k++){
                    send_buff[l++] = vect[k].first;
                    send_buff[l++] = vect[k].second;
                }
                for(;k<num_rec;k++){
                    send_buff[l++] = 0;
                    send_buff[l++] = 0;
                }
                MPI_Send(&send_buff,2*num_rec+1,MPI_INT,0,rank,MPI_COMM_WORLD);
            }
        }
    }else{

       vector<vector<pair<int,int>>> vect_adj(num_nodes,vector<pair<int,int>>(num_rec));
        for(int i=1;i<size;i++){
            int s = thread_time[i].first,e = thread_time[i].second;
            for(int j=s;j<=e;j++){
                int recv_buff[2*num_rec+1];
                MPI_Recv(&recv_buff,2*num_rec+1,MPI_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                rec_count[j] = recv_buff[0];
                int l = 1;
                for(int k=0;k<num_rec;k++){
                    vect_adj[j][k] = {recv_buff[l],recv_buff[l+1]};
                    l+=2;
                }
            }
        }

        for(int i=strt;i<=end;i++){
            vector<pair<int,int>> vect(circle[i].size());
            auto it = circle[i].begin();
            for(int k=0;k<vect.size();k++){
                vect[k] = *it;
                circle[i].erase(circle[i].begin());
                it = circle[i].begin();
            }

            sort(vect.begin(),vect.end(),compr);

            if(vect.size()>=num_rec){
                for(int k=0;k<num_rec;k++){
                    vect_adj[i][k] = {vect[k].first,vect[k].second};
                }
            }else{
                int k =0;
                for(;k<vect.size();k++){
                    vect_adj[i][k] = {vect[k].first,vect[k].second};
                }
                for(;k<num_rec;k++){
                    vect_adj[i][k] = {0,0};
                }
            }

        }

        FILE *fp1;
        unsigned char aw[4];
        unsigned char bw[4];
        string str_ = "output.dat";
        const char* str1= str_.c_str();
        fp1 = fopen(str1,"wb");
        if(fp1==NULL){
            cout<<"error\n";
        }else{
            for(int i=0;i<num_nodes;i++){
                int x = adj[i].size(),y= rec_count[i];
                aw[0] = ((x>>24) & 0xFF);aw[1] = ((x>>16) & 0xFF);aw[2] = ((x>>8) & 0xFF);aw[3] = ((x>>0) & 0xFF);
                fwrite(&aw,1,4,fp1);
                if(y<num_rec){
                    int j;
                    for(j=0;j<y;j++){
                        int z = vect_adj[i][j].first,v=vect_adj[i][j].second;
                        aw[0] = ((z>>24) & 0xFF);aw[1] = ((z>>16) & 0xFF);aw[2] = ((z>>8) & 0xFF);aw[3] = ((z>>0) & 0xFF);
                        fwrite(&aw,1,4,fp1);
                        aw[0] = ((v>>24) & 0xFF);aw[1] = ((v>>16) & 0xFF);aw[2] = ((v>>8) & 0xFF);aw[3] = ((v>>0) & 0xFF);
                        fwrite(&aw,1,4,fp1);
                    }
                    for(;j<num_rec;j++){
                        aw[0] = 'N';aw[1] = 'U';aw[2] = 'L';aw[3] = 'L';
                        fwrite(&aw,1,4,fp1);
                        fwrite(&aw,1,4,fp1);
                    }
                }else{
                    for(int j=0;j<num_rec;j++){
                        int z = vect_adj[i][j].first,v=vect_adj[i][j].second;
                        aw[0] = ((z>>24) & 0xFF);aw[1] = ((z>>16) & 0xFF);aw[2] = ((z>>8) & 0xFF);aw[3] = ((z>>0) & 0xFF);
                        fwrite(&aw,1,4,fp1);
                        aw[0] = ((v>>24) & 0xFF);aw[1] = ((v>>16) & 0xFF);aw[2] = ((v>>8) & 0xFF);aw[3] = ((v>>0) & 0xFF);
                        fwrite(&aw,1,4,fp1);
                    }
                }
            }
        }
        fclose(fp1);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop-start);
       // cout<<duration.count()<<" s\n";
    }

    MPI_Finalize();

}