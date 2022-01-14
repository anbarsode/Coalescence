vec3 acceleration(double t, vec3 pos, vec3 vel, double mass, double radius) // find acceleration based on gravity, stokes drag
{
    vec3 a;
    a.update(0, 0, 0, 0); 
    
    a.z -= g_eff; // Gravity - buoyancy
    
    vec3 u_ext = u_ext_field(t, pos); // ul/ut
    a = a - (vel - u_ext) * (six_pi_eta * radius / mass); // Stokes drag
    
    // Coriolis?
    //a.mag();
    //a.print();
    
    return a;
}

void move_particle(particle& p, double t, double dt) // RK4 integration
{
    vec3 xi = p.pos, vi = p.vel, xf, vf;
    vec3 c0, c1, c2, c3, d0, d1, d2, d3;
    
    c0 = vi * dt;
    d0 = acceleration(t, xi, vi, p.m, p.R) * dt;
    
    c1 = (vi + d0 * 0.5) * dt;
    d1 = acceleration(t + dt*0.5, xi + c0 * 0.5, vi + d0 * 0.5, p.m, p.R) * dt;
    
    c2 = (vi + d1 * 0.5) * dt;
    d2 = acceleration(t + dt*0.5, xi + c1 * 0.5, vi + d1 * 0.5, p.m, p.R) * dt;
    
    c3 = (vi + d2) * dt;
    d3 = acceleration(t + dt, xi + c2, vi + d2, p.m, p.R) * dt;
    
    xf = xi + (c0 + (c1 + c2) * 2 + c3) * (1/6.0);
    vf = vi + (d0 + (d1 + d2) * 2 + d3) * (1/6.0);
    
    p.update(xf, vf);
}

//https://www.geeksforgeeks.org/binary-insertion-sort/
long BinarySearch_xmin(std::vector<particle> &state, double item, long low, long high)
{
    while(low <= high)
    {
        long mid = low + (high - low) / 2;
        if(item == std::min(state[mid].pos.x, state[mid].Opos.x))
            return mid + 1;
        else if(item > std::min(state[mid].pos.x, state[mid].Opos.x))
            low = mid + 1;
        else
            high = mid - 1;
    }
    return low;
}

void InsertionSort_xmin(std::vector<particle> &state)
{
    long loc, j, k;
    particle p;
    for(long i = 1; i < state.size(); i++)
    {
        j = i - 1;
        p = state[i];
        loc = BinarySearch_xmin(state, std::min(p.pos.x, p.Opos.x), 0, j);
        state.erase(state.begin()+i);
        state.insert(state.begin()+loc, p);
    }
}

void particle_collisions(long loc, std::vector<long> &potential_colliders, std::vector<particle> &state, double t, double dt)
{
    particle p = state[loc];
    int ncol = 0;
    std::vector<long> del;
    
    for(long i=0; i<potential_colliders.size(); i++)
    {
        double frac = 0;
        particle q = state[potential_colliders[i]];
        if(p.detect_collision(q, frac))
        {
            ncol++;
            move_particle(p, t, dt*frac);
            move_particle(q, t, dt*frac);
            
            if(perform_collision(p, q))
            {
                particle pnew = p + q;
                move_particle(p, t, dt * (1-frac));
                state.erase(state.begin()+loc);
                state.insert(state.begin()+loc, pnew);
                del.push_back(potential_colliders[i]);
            }
            else
            {
                move_particle(p, t, dt * (1-frac));
                move_particle(q, t, dt * (1-frac));
                state.erase(state.begin()+loc);
                state.insert(state.begin()+loc, p);
                state.erase(state.begin()+potential_colliders[i]);
                state.insert(state.begin()+potential_colliders[i], q);
            }
        }
    }
    
    for(long i=del.size()-1; i>=0; i--)
    {
        //std::cout << "deleting " << std::endl;
        //state[del[i]].print();
        state.erase(state.begin()+del[i]);
    }
    
    /*if(del.size() != 0)
    {
        std::cout << "State after deletion " << std::endl;
        for(long i=state.size()-1; i>=0; i--)
        {
            state[i].print();
        }
    }

    if(ncol>=1)
        std::cout << "Warning: More than 1 collisions of a particle in a timestep" << std::endl;*/
}

void sweep_and_prune(std::vector<particle> &state, double t, double dt)
{
    if(deb) std::cout << "sweeping and pruning" << std::endl;
    InsertionSort_xmin(state);
    if(deb) std::cout << "sorted before pruning" << std::endl;
    
    double maxdisp = 0, disp;
    for(long i=0; i<state.size(); i++)
    {
        disp = fabs(state[i].pos.x - state[i].Opos.x);
        if(disp > maxdisp)
            maxdisp = disp;
    }
    if(deb) std::cout << "Maxdisp " << maxdisp << std::endl;
    
    for(long i=0; i<state.size()-1; i++)
    {
        std::vector<long> potential_colliders;
        double rlim, llim;
        particle p = state[i];
        
        if(p.pos.x > p.Opos.x)
        {
            rlim = p.pos.x;
            llim = p.Opos.x;
        }
        else
        {
            rlim = p.Opos.x;
            llim = p.pos.x;
        }
        
        //Forward search for potential colliders
        if(deb) std::cout << "Forward search" << std::endl;
        long j = i+1;
        while((std::min(state[j].pos.x, state[j].Opos.x) <= rlim) || (std::min(state[j].pos.x, state[j].Opos.x) <= rlim-Lbox))
        {
            double ysize = 0.5 * (fabs(p.pos.y - p.Opos.y) + fabs(state[j].pos.y - state[j].Opos.y));
            double ydist = fabs(0.5*((p.pos.y + p.Opos.y) - (state[j].pos.y + state[j].Opos.y)));
            if((ydist <= ysize) || (ydist >= (Lbox-ysize)))
            {
                double zsize = 0.5 * (fabs(p.pos.z - p.Opos.z) + fabs(state[j].pos.z - state[j].Opos.z));
                double zdist = fabs(0.5*((p.pos.z + p.Opos.z) - (state[j].pos.z + state[j].Opos.z)));
                if((zdist <= zsize) || (zdist >= (Lbox-zsize)))
                    potential_colliders.push_back(j);
            }
            
            j++;
            if(j >= state.size())
                j -= state.size();
            if(j == i)
            {
                std::cout << "Warning: Full loop" << std::endl;
                break;
            }
        }
        
        //Backward search for potential colliders
        if(deb) std::cout << "Backward search" << std::endl;
        j = i-1;
        while((std::min(state[j].pos.x, state[j].Opos.x) >= llim-maxdisp) || (std::min(state[j].pos.x, state[j].Opos.x) >= llim-maxdisp+Lbox))
        {
            if((std::max(state[j].pos.x, state[j].Opos.x) >= llim) || (std::max(state[j].pos.x, state[j].Opos.x) >= llim+Lbox))
            {
                double ysize = 0.5 * (fabs(p.pos.y - p.Opos.y) + fabs(state[j].pos.y - state[j].Opos.y));
                double ydist = fabs(0.5*((p.pos.y + p.Opos.y) - (state[j].pos.y + state[j].Opos.y)));
                if((ydist <= ysize) || (ydist >= (Lbox-ysize)))
                {
                    double zsize = 0.5 * (fabs(p.pos.z - p.Opos.z) + fabs(state[j].pos.z - state[j].Opos.z));
                    double zdist = fabs(0.5*((p.pos.z + p.Opos.z) - (state[j].pos.z + state[j].Opos.z)));
                    if((zdist <= zsize) || (zdist >= (Lbox-zsize)))
                        potential_colliders.push_back(j);
                }
            }
                
            j--;
            if(j < 0) j += state.size();
            if(j == i)
            {
                std::cout << "Warning: Full loop" << std::endl;
                break;
            }
            if(deb) std::cout << "i " << i << " backward search j " << j << std::endl;
        }
        if(deb) std::cout << "Search completed" << std::endl;
        particle_collisions(i, potential_colliders, state, t, dt);
    }
    if(deb) std::cout << "swept and pruned" << std::endl;
}

void check_intersections(std::vector<particle> &state, double t)
{
    for(long i=0; i<state.size()-1; i++)
    {
        particle p = state[i];
        for(long j=i+1; j<state.size(); j++) // PBC
        {
            particle q = state[j];
            if((p.pos - q.pos).mag() <= (p.R + q.R))
                std::cout << "Particles " << p.id << " and " << q.id << " are intersecting at time " << t << std::endl;
        }
    }
}

void move_all_particles(std::vector<particle> &state, double t, double dt)
{
    if(deb) std::cout << "Moving all particles for a timestep" << std::endl;
    for(long i=0; i<state.size(); i++)
        move_particle(state[i], t, dt); // Timestep evolution
        
    sweep_and_prune(state, t, dt); // Checking and performing all the collisions/mergers
    
    for(long i=0; i<state.size(); i++) // PBC
    {
        state[i].pos.x = (state[i].pos.x >= Lbox) ? (state[i].pos.x - Lbox) : (state[i].pos.x);
        state[i].pos.y = (state[i].pos.y >= Lbox) ? (state[i].pos.y - Lbox) : (state[i].pos.y);
        state[i].pos.z = (state[i].pos.z >= Lbox) ? (state[i].pos.z - Lbox) : (state[i].pos.z);
        
        state[i].pos.x = (state[i].pos.x < 0) ? (state[i].pos.x + Lbox) : (state[i].pos.x);
        state[i].pos.y = (state[i].pos.y < 0) ? (state[i].pos.y + Lbox) : (state[i].pos.y);
        state[i].pos.z = (state[i].pos.z < 0) ? (state[i].pos.z + Lbox) : (state[i].pos.z);
    }
    if(deb) std::cout << "Moved all particles for a timestep" << std::endl;
    check_intersections(state, t);
}

void brute_force_timestep(std::vector<particle> &state, double t, double dt)
{
    if(deb) std::cout << "Moving all particles for a timestep" << std::endl;
    for(long i=0; i<state.size(); i++)
        move_particle(state[i], t, dt); // Timestep evolution
        
    for(long i=0; i<state.size()-1; i++)
    {
        particle p = state[i];
        std::vector<long> potential_colliders;
        for(long j=i+1; j<state.size(); j++)
        {
            potential_colliders.push_back(j);
        }
        particle_collisions(i, potential_colliders, state, t, dt);
    }
    
    PBC_precipitation(state, t);
    
    /*for(long i=0; i<state.size(); i++) // PBC and birth-death
    {
        state[i].pos.x = (state[i].pos.x >= Lbox) ? (state[i].pos.x - Lbox) : (state[i].pos.x);
        state[i].pos.y = (state[i].pos.y >= Lbox) ? (state[i].pos.y - Lbox) : (state[i].pos.y);
        state[i].pos.z = (state[i].pos.z >= Lbox) ? (state[i].pos.z - Lbox) : (state[i].pos.z);
        
        state[i].pos.x = (state[i].pos.x < 0) ? (state[i].pos.x + Lbox) : (state[i].pos.x);
        state[i].pos.y = (state[i].pos.y < 0) ? (state[i].pos.y + Lbox) : (state[i].pos.y);
        state[i].pos.z = (state[i].pos.z < 0) ? (state[i].pos.z + Lbox) : (state[i].pos.z);
    }*/
    if(deb) std::cout << "Moved all particles for a timestep" << std::endl;
    // check_intersections(state, t);
}

void GTCP3Dsimulation(std::vector<particle> &icstate, double t0, double dt, long Nsteps, std::string out_name_base, long log_freq, long out_freq)
{
    std::vector<particle> state = icstate;
    long Npart = state.size();
    
    auto start = std::chrono::high_resolution_clock::now();
    for(long i=0; i<Nsteps; i++)
    {
        if(i % log_freq == 0)
        {
            auto elap = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
            std::cout << "Step " << i << " out of " << Nsteps << ", Npart: " << Npart << ", elapsed time(s) " << elap.count() / 1e6  << ";" << std::endl;
        }
        if(i % out_freq == 0)
        {
            write_snapshot(out_name_base, i, (t0 + dt * i), state);
            std::cout << "\tWritten snapshot at time(s): " << (t0 + dt * i)*UT << std::endl;
        }
        brute_force_timestep(state, (t0 + dt * i), dt);
        //move_all_particles(state, (t0 + dt * i), dt);
        Npart = state.size();
    }
    write_snapshot(out_name_base, Nsteps, (t0 + dt * Nsteps), state);
    auto elap = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "\nWritten final snapshot at time(s): " << (t0 + dt * Nsteps)*UT << ", elapsed time(s) " << elap.count() / 1e6  << ";" << std::endl;
}
