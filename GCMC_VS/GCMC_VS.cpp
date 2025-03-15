// gcmc_vs.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <fstream>
#include <list>
#include <bitset>
#include <unordered_map>
#include <omp.h>
using namespace std;


//----------------------------------
//global params
//----------------------------------
//lattice length 
const int M = 64;
static const int particle_width = 1;
static const int particle_length = 8;


static void init_OCC(bitset<M> (& OCC)[M], int size) {
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++)
        {
            OCC[i][j] = 0;
        }
    }
   
}

//print the OCC lattice to the console
static void print_OCC(bitset<M>(&OCC)[M]) {
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {

            cout << OCC[j][i];
        }

        cout << "\n";
    }
}

static bool checkOccupancyHorizontal(bitset<M> (& OCC)[M], int particle_length ,vector<int> pos) {
    int x = pos[0];
    int y = pos[1];
   
    //option 1 
    
    for (int i = 0; i < particle_length; i++)
    {
        int current_x = x + i;
        //if overexceeding boundary jump back  
        if (current_x > M - 1)  current_x = x + i - M;
        if (OCC[current_x][y] == 1) {
            return true;
        }
    }
    
    return false;
   
    //option 2 is faster
    /*
    for (int i = 0; i < particle_length; i++)
    {
        int current_x = x + i;
        if (current_x > M - 1)  current_x = x + i - M;
        if ((OCC[current_x] >> y & bitset<M>((1 << particle_length) - 1)).count() > 0) {
            return true;
        }
    }
    return false;
    */
}
static bool checkOccupancyVertical(bitset<M>(&OCC)[M], int particle_length, vector<int> pos) {
    int x = pos[0];
    int y = pos[1];

    for (int i = 0; i < particle_length; i++)
    {
        int current_y = y + i;
        if (current_y > M - 1) current_y = y + i - M;
        if (OCC[x][current_y] == 1) {
            return true;
        }
    }
    return false;
}


// pos = {x,y}, adds a particle to the OCC lattice with pos being its defining point
static bool add_particle(bitset<M> (&OCC)[M], vector<int> pos, bool horizontal) {
    int x = pos[0];
    int y = pos[1];
    //if horizontal and can be placed
    if (horizontal && !checkOccupancyHorizontal(OCC, particle_length, pos))
    {
        //->place particle
        for (int i = 0; i < particle_length; i++)
        {
            int current_x = x + i;
            if (current_x > M - 1) current_x = x + i - M;
            //8-bit mask at y
            OCC[current_x][y] = 1;
        }
        return true;
    }
    //else if 
    else if(!horizontal && !checkOccupancyVertical(OCC, particle_length, pos)){
        for (int i = 0; i < particle_length; i++)
        {
            int current_y = y + i;
            if (current_y > M - 1) current_y = y + i - M;
            OCC[x][current_y] = 1; 
        }
        return true;
    }
    else {
        return false;
    }
    
}

//removes a particle from the OCC lattice, prop are the properties of the selected particle {x_pos, y_pos, orientation}
static void remove_particle(bitset<M> (&OCC)[M], vector<int> prop) {
    int x = prop[0];
    int y = prop[1];
  
    //if horizontal
    if (prop[2]) //if horizontal 
    {
        for (int i = 0; i < particle_length; i++)
        {
            int current_x = x + i;
            if (current_x > M - 1) current_x = x + i - M;
            OCC[current_x][y] = 0;

        }
    } else { //if vertical
        for (int i = 0; i < particle_length; i++)
        {
            int current_y = y + i;
            if (current_y > M - 1) current_y = y + i - M;
            OCC[x][current_y] = 0;

        }
    }
}


//z is activity, M is lattice length
static void GCMC(
    bitset<M> (& OCC)[M],
    vector<vector<int>> &rods,
    uniform_int_distribution<int> &unif_distr_pos,
    uniform_real_distribution<double> &unif_distr,
    mt19937 &gen, double z, int& N_v, int& N_h) {
    
    uniform_int_distribution<int> dis1(0,1);
    int N = N_v + N_h;
    //true -> insertion, false -> deletion
    bool insertion = dis1(gen);
    
    //true -> horizontal, false -> vertical
    bool horizontal = dis1(gen);
    
    // 
    //roll for next prob
    double r_alpha = unif_distr(gen);
    //prng for deletion or insertion probability
        //insertion
        if (insertion) {
            //inserting 
            //roll positions selected for insertion
            int r_x = unif_distr_pos(gen);
            int r_y = unif_distr_pos(gen);

            /*
            while (OCC[r_x][r_y] == 1) {
                 r_x = unif_distr_pos(gen);
                 r_y = unif_distr_pos(gen);
            }
           */

            //acceptance probability for insertion
            double alpha = min(2 * z * (pow(M, 2) / ((double)(N + 1))), 1.0);
           
           if (r_alpha <= alpha) { //try acceptance
               if (add_particle(OCC, { r_x, r_y }, horizontal)) {  //try adding particle, if successful add it to rods
                   rods.push_back({ r_x, r_y, horizontal });
                   horizontal ? N_h++ : N_v++; //increase vertical or horizontal counter
               }
           }else {
               return;
           }
            
        }
        //deletion
        else {
            
            //acceptance probability for  deletion
            double alpha = min(N/z * 1 / (double)(2 * (pow(M, 2))), 1.0);
            
            if (r_alpha <= alpha) { //try acceptance
                if (N == 0) return;
                //roll random index
                uniform_int_distribution<int> dis(0, N - 1);
                int rand = dis(gen);
                vector<int> props = rods[rand];
                remove_particle(OCC, props); //remove particle from OCC
                rods.erase(rods.begin() + rand); //remove particle froom list
                props[2] ? N_h-- : N_v--; //decrease vertical or horizontal counter
                
            }

        }




}


static void simulate(string filename, double z, int thermal_time, int steps_between_measurements, int number_of_measurements) {
    cout << "Thread: " << omp_get_thread_num() << "                  " << endl;
    //set seed for the RNG
    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_real_distribution<double> unif_distr(0.0, 1.0);
    uniform_int_distribution<int> unif_distr_pos(0, M - 1);

    const int dim_length = M;

    //data structure of OCC-lattice  
    bitset<M> OCC[M];
    init_OCC(OCC, M);

   
    int N_h = 0;
    int N_v = 0;
    //bookkeeping of rods {x_pos, y_pos, orientation}
    vector<vector<int>> rods = {};

    for (int i = 0; i < thermal_time; i++)
    {
        GCMC(OCC,rods ,unif_distr_pos, unif_distr, gen, z, N_v, N_h);
    }
   

    //write data to file
    ofstream data_file(filename);
    data_file << "it." << "\t" << "N" << "\t" << "N_h" << "\t" << "N_v" << "\n";
    for (unsigned int i = 0; i < number_of_measurements; i++) {
        random_device rand_dev;
        mt19937 gen(rand_dev());
        for (unsigned int j = 0; j < steps_between_measurements; j++)
        {
            GCMC(OCC, rods, unif_distr_pos, unif_distr, gen, z, N_v, N_h);
        }
        int N = N_h + N_v;
        data_file << i << "\t" << N << "\t" << N_h << "\t" << N_v << "\n";


    }

    //print_OCC(OCC);
   

    data_file.close();




}



int main()
{
    clock_t start = clock();
    string filename_1 = "particle_number.txt";
    //omp_set_num_threads(3);
    cout << "max av. threads" << omp_get_max_threads() << endl;
    #pragma omp parallel 
    {
        #pragma omp sections
        {
        #pragma omp section
        {
            simulate("particle_number_z=0.56.txt", 0.56, 1000, 1000, 10000);

        }
        
        #pragma omp section
        {
            simulate("particle_number_z=1.2.txt", 1.2, 1000, 1000, 10000);
        }
        
        #pragma omp section
        {
            simulate("particle_number_z=0.05.txt", 0.05, 1000, 1000, 10000);

        }
        }

    }
    
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    printf("execution time: %.3f sec", elapsed);


}







// Programm ausführen: STRG+F5 oder Menüeintrag "Debuggen" > "Starten ohne Debuggen starten"
// Programm debuggen: F5 oder "Debuggen" > Menü "Debuggen starten"

// Tipps für den Einstieg: 
//   1. Verwenden Sie das Projektmappen-Explorer-Fenster zum Hinzufügen/Verwalten von Dateien.
//   2. Verwenden Sie das Team Explorer-Fenster zum Herstellen einer Verbindung mit der Quellcodeverwaltung.
//   3. Verwenden Sie das Ausgabefenster, um die Buildausgabe und andere Nachrichten anzuzeigen.
//   4. Verwenden Sie das Fenster "Fehlerliste", um Fehler anzuzeigen.
//   5. Wechseln Sie zu "Projekt" > "Neues Element hinzufügen", um neue Codedateien zu erstellen, bzw. zu "Projekt" > "Vorhandenes Element hinzufügen", um dem Projekt vorhandene Codedateien hinzuzufügen.
//   6. Um dieses Projekt später erneut zu öffnen, wechseln Sie zu "Datei" > "Öffnen" > "Projekt", und wählen Sie die SLN-Datei aus.




