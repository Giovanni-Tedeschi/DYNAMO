#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <sstream>
#include <cfloat>
#include <unordered_map>

using namespace std;

class Atom{
    public:
        double Pos[3];
        double Vel[3];
        double Acc[3];
        double Mass;
        int ID;
        int MolID;

        double v(){
            return sqrt(pow(Vel[0],2) + pow(Vel[1],2) + pow(Vel[2],2));
        }

        Atom(){
            Pos[0] = 0; Pos[1] = 0; Pos[2] = 0;
            Vel[0] = 0; Vel[1] = 0; Vel[2] = 0;
            Acc[0] = 0; Acc[1] = 0; Acc[2] = 0;
            Mass = 1.;
        }
};

class Simulation{
    public:
        int snap_i = 0;
        int N_atoms;
        Atom *A; 
        double dt;
        double CFL;
        double R_coul;
        double f_coul;
        double R_harm;
        double f_harm;
        double R0_harm;
        double t_max;
        double dt_snap;
        double t = 0;
        std::string input_file;
        std::string output_path;

        void find_timestep(){
            /*dt = FLT_MAX;
            for(int i=0; i<N_atoms; i++){
                double dt_courant =  A[i].v() * CFL;
                dt = std::min(dt, dt_courant);
            } */
            dt = 0.0001;
        }

        void do_euler_step(){
            for(int i=0; i<N_atoms; i++){    
                for(int j=0; j<3; j++){
                    A[i].Pos[j] += A[i].Vel[j] * dt;
                }
            }
        }

        void first_kick(){
            for(int i=0; i<N_atoms; i++){    
                for(int j=0; j<3; j++){
                    A[i].Vel[j] += A[i].Acc[j] * 0.5 * dt;
                }
            }
        }

        void drift(){
            for(int i=0; i<N_atoms; i++){    
                for(int j=0; j<3; j++){
                    A[i].Pos[j] += A[i].Vel[j] * dt;
                }
            }
        }

        void second_kick(){
            for(int i=0; i<N_atoms; i++){    
                for(int j=0; j<3; j++){
                    A[i].Vel[j] += A[i].Acc[j] * 0.5 * dt;
                }
            }
        }

        double distance(Atom A, Atom B){
            return sqrt(pow(A.Pos[0]-B.Pos[0],2) + pow(A.Pos[1]-B.Pos[1],2) + pow(A.Pos[2]-B.Pos[2],2));
        }

        void compute_accelerations(){
            double r_ij = 0.;
            double a_ij = 0.;
            int i; int j; int k;
            for(i=0; i<N_atoms; i++){    
                for(k=0; k<3; k++) A[i].Acc[k] = 0.;
                for(j=0; j<N_atoms; j++){
                    r_ij = distance(A[i],A[j]);
                    a_ij = 0.;
                    if((r_ij > 0)){
                        // Coulomb interaction
                        if(r_ij < R_coul) a_ij += f_coul / pow(r_ij,2) / A[i].Mass;

                        // Harmonic interaction
                        if(r_ij < R_harm) a_ij += - f_harm * (r_ij - R0_harm) / A[i].Mass;

                        // Compute acceleration acting on particle i on all dimensions
                        for(k=0; k<3; k++){
                            A[i].Acc[k] += a_ij * (A[i].Pos[k] - A[j].Pos[k]) / r_ij;
                        }
                    }
                }
            }
        }

        void dump(){
            std::string output_file = output_path + std::to_string(snap_i) + ".txt";
            ofstream fp(output_file, std::ios::out);
            for(int i=0; i<N_atoms; i++){
                fp << A[i].Pos[0] << " " << A[i].Pos[1] << " " << A[i].Pos[2] << " ";
                fp << A[i].Vel[0] << " " << A[i].Vel[1] << " " << A[i].Vel[2] << " ";
                fp << A[i].Acc[0] << " " << A[i].Acc[1] << " " << A[i].Acc[2] << " ";
                fp << A[i].ID << " " << A[i].MolID;
                fp << "\n";
            }
            fp.close();
            snap_i ++;
        }

        void read_init(){

            // Read the initial condition file
            std::ifstream file;
            file.open(input_file);
            N_atoms = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
            A = static_cast<Atom*>(malloc(sizeof(Atom)*N_atoms));
            std::cout << N_atoms << "\n";           

            FILE* finp;
            finp = fopen(input_file.data(), "r");
            
            for(int i=0; i<N_atoms; i++){
                fscanf(finp, "%lf %lf %lf %lf %lf %lf %lf %d %d", &A[i].Pos[0], &A[i].Pos[1], &A[i].Pos[2], &A[i].Vel[0], &A[i].Vel[1], &A[i].Vel[2], &A[i].Mass, &A[i].ID, &A[i].MolID);
            }

            fclose(finp);
        }

        // Constructor for the Simulation class, reads input and sets up the A atoms array
        Simulation(std::string fname){

            // Read the param file
            std::unordered_map<std::string, std::string> configData;
            std::ifstream configFileStream(fname);
            for (std::string line{}; std::getline(configFileStream, line); ){
                std::istringstream iss{ line };
                if (std::string id{}, value{}; std::getline(std::getline(iss, id, ':') >> std::ws, value)) {
                    configData[id] = value;
                }
            }

            R_coul = stof(configData["R_coul"]);
            f_coul = stof(configData["f_coul"]);
            R_harm = stof(configData["R_harm"]);
            f_harm = stof(configData["f_harm"]);
            R0_harm = stof(configData["R0_harm"]);
            CFL     = stof(configData["CFL"]);
            t_max   = stof(configData["t_max"]);
            dt_snap   = stof(configData["dt_snap"]);
            input_file = configData["input_file"];
            input_file.erase(std::remove_if(input_file.begin(), input_file.end(), [](char c) { return c == '\n' || c == '\r';}), input_file.end());
            output_path = configData["output_path"];
            output_path.erase(std::remove_if(output_path.begin(), output_path.end(), [](char c) { return c == '\n' || c == '\r';}), output_path.end());

        }
};


int main(int argc, char* argv[]){
    // Read the param file name as terminal input
    std::string param_file = argv[1];

    // Initialize the simulation by reading and parsing the parameter file
    Simulation sim(param_file);

    // Read the initial condition file and fill the array of Atom ojects 
    sim.read_init();
    
    // Find the initial timestep
    sim.find_timestep();

    for(sim.t = 0; sim.t < sim.t_max; sim.t += sim.dt){
        sim.first_kick();
        sim.compute_accelerations();
        
        sim.drift();
        sim.second_kick();
        sim.find_timestep();

        if(sim.t > sim.snap_i*sim.dt_snap) sim.dump();
        std::cout << sim.t << " " << sim.dt << endl;
    }
    
}