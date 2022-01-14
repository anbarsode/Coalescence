//global timestep N-body simulation with contact potentials in 3D

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <iomanip>
#include <vector>
#include <sstream>
#include <chrono>

#include "Vector3D.cpp"
#include "Global_Constants.cpp"
#include "Particle.cpp"
#include "Utility_Functions.cpp"
#include "Input_Output.cpp"
#include "Simulator.cpp"

int main(int argc, char** argv)
{
    char* ic_name = argv[1];
    std::string out_base(""), out_dir("");
    double t0, dt;
    long log_freq = 1, out_freq = 1;
    long Nsteps;
    std::vector<particle> icstate;
    
    //icstate =  read_ic_file(ic_name, out_dir, out_base, t0, dt, Nsteps, log_freq, out_freq);
    icstate =  IC_generator(ic_name, out_dir, out_base, t0, dt, Nsteps, log_freq, out_freq);
    
    std::cout << "IC name: " << ic_name << std::endl;
    std::cout << "Output directory: " << out_dir << std::endl;
    std::cout << "Output base name: " << out_base << std::endl;
    std::cout << "Output frequency: " << out_freq << std::endl;
    std::cout << "Logging frequency: " << log_freq << std::endl;
    std::cout << "Total number of particles: " << icstate.size() << std::endl;
    std::cout << "Number of timesteps: " << Nsteps << std::endl;
    std::cout << "Starting time: " << t0 << std::endl;
    std::cout << "Timestep size: " << dt << std::endl;
    std::cout << "Ending time: " << (Nsteps*dt) << std::endl;
    
    GTCP3Dsimulation(icstate, t0, dt, Nsteps, out_dir + "/" + out_base, log_freq, out_freq);
}

