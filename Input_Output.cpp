void split(const std::string &s, char delim, std::vector<std::string> &elems) //splitting a line into substrings separated by a delimiter
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while(std::getline(ss, item, delim))
        elems.push_back(item);
}

// read ics and convert to a vector of particles (state)
std::vector<particle> read_ic_file(char* ic_name, std::string &out_dir, std::string &out_base, \
                                double &t0, double &dt, long &Nsteps, long &log_freq, long &out_freq)
{
    long Npart;
    std::string outlist_name("");
    std::ifstream f(ic_name, std::ios::in);
    assert(f.is_open());
    std::vector<particle> icstate;
    std::string line;
    bool reading_icstate = false;
    while(std::getline(f, line))
    {
        if((line[0] == '#') || (line[0] == 0)) continue;
        else
        {
            std::vector<std::string> data;
            split(line, '\t', data);
            if(!reading_icstate)
            {
                char keyword[100];
                strcpy(keyword, (data[0]).c_str());
                if(!strcmp(keyword, "output_directory")) out_dir = data[1];
                if(!strcmp(keyword, "output_file_base")) out_base = data[1];
                if(!strcmp(keyword, "output_frequency")) out_freq = stoi(data[1]);
                if(!strcmp(keyword, "logging_frequency")) log_freq = stoi(data[1]);
                if(!strcmp(keyword, "current_time")) t0 = stod(data[1]);
                if(!strcmp(keyword, "timestep_size")) dt = stod(data[1]);
                if(!strcmp(keyword, "number_of_timesteps")) Nsteps = stol(data[1]);
                if(!strcmp(keyword, "number_of_particles")) Npart = stol(data[1]);
                if(!strcmp(keyword, "STATE")) reading_icstate = true;
            }
            else
            {
                particle p;
                p.id = stol(data[0]);
                p.m = stod(data[1]);
                p.R = stod(data[2]);
                p.pos.x = stod(data[3]);
                p.pos.y = stod(data[4]);
                p.pos.z = stod(data[5]);
                p.vel.x = stod(data[6]);
                p.vel.y = stod(data[7]);
                p.vel.z = stod(data[8]);
                icstate.push_back(p);
            }
        }
    }
    f.close();
    
    if(Npart != icstate.size())
        std::cout << "Error: Particle number doesnt match" << std::endl;
    
    return icstate;
}

//write current state to a text file readable by ovito/python

void write_snapshot(std::string out_name, int timestepID, double time, std::vector<particle> &state)
{
    write_correlation(out_name + "_corr_" + std::to_string(timestepID) + ".txt", state);
    std::ofstream f(out_name + "_" + std::to_string(timestepID) + ".txt", std::ios::out);
    assert(f.is_open());
    f << "UnitLengthm\t" << UL << std::endl;
    f << "UnitMasskg\t" << UM << std::endl;
    f << "UnitTimes\t" << UT << std::endl;
    f << "current_time\t" << time << std::endl;
    f << "Lbox\t" << Lbox << std::endl;
    f << "STATE\t" << std::endl;
    f << "#ID\tm\tR\tx\ty\tz\tvx\tvy\tvz" << std::endl;
    for(long i=0; i<state.size(); i++)
        f << state[i].id << "\t" << state[i].m << "\t" << state[i].R << "\t"
        << state[i].pos.x << "\t" << state[i].pos.y << "\t" << state[i].pos.z << "\t"
        << state[i].vel.x << "\t" << state[i].vel.y << "\t" << state[i].vel.z << std::endl;
    f.close();
}

std::vector<particle> IC_generator(char* ic_name, std::string &out_dir, std::string &out_base, \
                                double &t0, double &dt, long &Nsteps, long &log_freq, long &out_freq)
{
    out_dir = output_directory;
    out_base = output_file_base;
    t0 = current_time;
    dt = timestep_size;
    Nsteps = number_of_timesteps;
    log_freq = logging_frequency;
    out_freq = output_frequency;
    
    std::vector<particle> state;
    for(long i=0; i<Nparticles; i++)
        state.push_back(create_random_particle());
    
    write_snapshot(ic_name, 0, 0, state);
    return state;
}
