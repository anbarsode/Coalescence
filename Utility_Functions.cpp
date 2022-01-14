//build lookup tables and interpolation methods for uext and merge_or_collide from text file data
bool merge_or_bounce(double costheta, double u21) // Make a lookup table here
{
    if(u21 > bounce_threshold_u12 && costheta < bounce_threshold_costheta)
        return false; // bounce
    else
        return true; // merge
}

vec3 u_ext_field_random(double t, vec3 pos) // Feed in the external velocity field as a function of spacetime (lookup table)
{
    srand(24);
    vec3 uext, ycap;
    uext.update(0,0,0,0);
    ycap.update(0,1,0,1);
    for(int i=0; i<10; i++)
    {
        vec3 center;
        center.update(Lbox*rand()/float(RAND_MAX), Lbox*0.5, Lbox*rand()/float(RAND_MAX), 0);
        double strength = 1.0, extent = 0.18 * Lbox;// + Lbox * rand() / float(RAND_MAX) * 0.1;
        //std::cout << i << " " << strength << " " << extent << std::endl;
        //center.print();
        double weight = strength * exp(-(pos - center).mag2() / extent / extent);
        //std::cout << weight << std::endl;
        uext = uext + ((pos - center) ^ (ycap * weight));// * (Lbox * pow((pos - center).mag(), -1) * rand() / RAND_MAX);
    }
    //uext.print();
    return uext * 1e11;
}

vec3 u_ext_field(double t, vec3 pos) // Taylor-Green vortex
{
    vec3 uext;
    uext.update(0,0,0,0);
    
    uext.x = sin(TG_kx * pos.x) * cos(TG_kz * pos.z) * sin(TG_w * t);
    uext.z = -cos(TG_kx * pos.x) * sin(TG_kz * pos.z) * sin(TG_w * t);
    
    return uext * TG_Strength;
}

void write_correlation(std::string corr_name, std::vector<particle> &state)
{
    std::vector<int> hist;
    int nbins = 100;
    double binsize = Lbox / nbins * sqrt(2);
    for(int i=0; i<nbins; i++)
        hist.push_back(0);
        
    std::ofstream f(corr_name, std::ios::out);
    assert(f.is_open());
    for(long i=0; i<state.size()-1; i++)
        for(long j=i+1; j<state.size(); j++)
        {
            double r12 = (state[i].pos - state[j].pos).mag();
            hist[int(r12 / binsize)]++;
        }
        
    for(int i=0; i<nbins; i++)
        f << hist[i] << std::endl;
            
    f.close();
}

particle create_random_particle(void)
{
    particle pnew;
    vec3 pos, vel;
    pos.update(double(rand())/RAND_MAX * Lbox, 0, double(rand())/RAND_MAX * Lbox);
    vel.update(double(rand())/RAND_MAX - 0.5, 0, double(rand())/RAND_MAX - 0.5);
    vel.mag();
    vel = vel * (birth_speed / vel.m);
    pnew.create(0, birth_mass, birth_radius, pos, vel);
    return pnew;
}

void PBC_precipitation(std::vector<particle> &state, double t)
{
    std::ofstream f(precipitation_file, std::ios::app);
    int ins = (Nparticles > state.size()) ? (Nparticles - state.size()) : 0;
    double precip = 0;
    
    for(long i=0; i<state.size(); i++) // PBC and birth-death
    {
        state[i].pos.x = (state[i].pos.x >= Lbox) ? (state[i].pos.x - Lbox) : (state[i].pos.x);
        state[i].pos.y = (state[i].pos.y >= Lbox) ? (state[i].pos.y - Lbox) : (state[i].pos.y);
        state[i].pos.z = (state[i].pos.z >= Lbox) ? (state[i].pos.z - Lbox) : (state[i].pos.z);
        
        state[i].pos.x = (state[i].pos.x < 0) ? (state[i].pos.x + Lbox) : (state[i].pos.x);
        state[i].pos.y = (state[i].pos.y < 0) ? (state[i].pos.y + Lbox) : (state[i].pos.y);
        state[i].pos.z = (state[i].pos.z < 0) ? (state[i].pos.z + Lbox) : (state[i].pos.z);
        
        if(state[i].R >= death_radius)
        {
            ins += int(pow(state[i].R / birth_radius, 3.0)) - 1;
            precip += state[i].m;
            state.erase(state.begin()+i);
            state.insert(state.begin()+i, create_random_particle());
        }
    }
    
    for(int i=0; i<ins; i++)
        state.push_back(create_random_particle());
    
    if(precip)
        f << t << "\t" << precip << std::endl;
    f.close();
}