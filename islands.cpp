#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <vector>
using namespace std;

void get_islands(double bmp[64][25], int kkk[6][64][25], int pos[6][64][25], int islands[64][25]){
    int island = 0;
    double visited[64][25];
    vector <int> vclust;
    vector <int> vpix;

    for (int count = 0; count < 25; count++)
        {
            for (int coun = 0; coun < 64; coun++)
            {
                visited[coun][count] = 0;
                islands[coun][count] = 0;
            }
        }

    for (int count = 0; count < 25; count++)
        {
            for (int coun = 0; coun < 64; coun+=2)
            {
                if(bmp[coun][count] == 0 || visited[coun][count] == 1){
                    continue;
                }
                visited[coun][count] = 1;
                island++;
                islands[coun][count] = island;
                for(int i = 0; i < 6; i++){
                    if(bmp[pos[i][coun][count]][kkk[i][coun][count]] == 0 || 
                    visited[pos[i][coun][count]][kkk[i][coun][count]] == 1 ||
                    (pos[i][coun][count] == 0 && kkk[i][coun][count] == 0)){
                        continue;
                    }
                    //cout << count << "\t" << coun << "\t" << kkk[i][coun][count] << "\t" << pos[i][coun][count] << endl;
                    visited[pos[i][coun][count]][kkk[i][coun][count]] = 1;
                    islands[pos[i][coun][count]][kkk[i][coun][count]] = island;
                    vclust.push_back(kkk[i][coun][count]);
                    vpix.push_back(pos[i][coun][count]);
                }
                
                while(vclust.size() > 0){
                    //cout << vclust.size() << "\t" << vclust[0] << "\t" << vpix[0] << endl;
                    for(int i = 0; i < 6; i++){
                        //cout << kkk[i][vpix[0]][vclust[0]] << "\t" << pos[i][vpix[0]][vclust[0]] << "\t" << bmp[pos[i][vpix[0]][vclust[0]]][kkk[i][vpix[0]][vclust[0]]] << endl;
                        if(bmp[pos[i][vpix[0]][vclust[0]]][kkk[i][vpix[0]][vclust[0]]] == 0 || 
                        visited[pos[i][vpix[0]][vclust[0]]][kkk[i][vpix[0]][vclust[0]]] == 1 ||
                        (pos[i][vpix[0]][vclust[0]] == 0 && kkk[i][vpix[0]][vclust[0]] == 0)){
                            continue;
                        }

                        vclust.push_back(kkk[i][vpix[0]][vclust[0]]);
                        vpix.push_back(pos[i][vpix[0]][vclust[0]]);
                    }
                    visited[vpix[0]][vclust[0]] = 1;
                    islands[vpix[0]][vclust[0]] = island;
                    vclust.erase(vclust.begin());
                    vpix.erase(vpix.begin());
                }
        }
    }
    islands[0][0] = island;
}

double get_bightest_island(double bmp[64][25], int islands[64][25]){
    int islands_numb = islands[0][0];
    double size_islands[islands_numb];
    for(int i = 0; i < islands_numb; i++){
        size_islands[i] = 0;
    }
    for(int coun = 0; coun < 25; coun++) {
		for(int count = 0; count < 64; count+=2) {
            if(bmp[count][coun] > 0){
                size_islands[islands[count][coun]-1]+=bmp[count][coun];
            }
        }
    }
    int brightest_island = 0;
    double brightest_island_size = 0;
	double sum_size = 0;
    for(int i = 0; i < islands_numb; i++){
        sum_size+=size_islands[i];
        if(size_islands[i] > brightest_island_size){
			brightest_island_size = size_islands[i];
            brightest_island = i+1;

        }
    }
    //cout << endl;
    for(int coun = 0; coun < 25; coun++) {
		for(int count = 0; count < 64; count+=2) {
            if(islands[count][coun] != brightest_island){
                bmp[count][coun] = 0;
            }
        }
    }
	//cout << brightest_island_size << "\t" << sum_size << endl;
	return brightest_island_size/sum_size;
}