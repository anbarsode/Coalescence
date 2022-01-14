// Dimensionless class

class particle // Spherical particle moving in 3D Cartesian space
{
    public:
        long id;
        double m, R;
        vec3 pos, vel; // current state
        vec3 Opos; // previous state
        
        void create(long id_new, double m_new, double R_new, vec3& pos_new, vec3& vel_new)
        {
            id = id_new;
            m = m_new;
            R = R_new;
            pos = pos_new;
            vel = vel_new;
            Opos = pos_new;
        }
        
        void update(vec3& pos_new, vec3& vel_new)
        {
            Opos = pos;
            pos = pos_new;
            vel = vel_new;
        }
        
        void print() // Print the particle's data
        {
            std::cout << "ID: " << id << "\nMass: " << m << ", Radius: " << R << std::endl;
            std::cout << "Position: ";
            pos.print();
            std::cout << "Velocity: ";
            vel.print();
        }
        
        // Detect in-between collision with another particle p using linear trajectory approximation
        // and quadratic solution method in 3D
        bool detect_collision_quad(particle& p, double& frac)
        {
            if(deb) std::cout << "checking a collision" << std::endl;
            bool collision = false;
            
            double a,b,c,d;
            
            vec3 r12 = Opos - p.Opos;
            vec3 v12 = (pos - Opos) - (p.pos - p.Opos);
            
            a = v12.mag2();
            b = 2 * (r12 & v12) / a;
            c = (r12.mag2() - (R + p.R)*(R + p.R)) / a;
            if(c <= 0)
            {
                collision = true;
                frac = 0;
                return collision;
            }
            
            d = b*b - 4*c;
            
            if(d >= 0)
            {
                double tlo, thi, dr = sqrt(d), tol=1e-3;
                tlo = 0.5 * (-b - dr);
                if((tlo > tol) && (tlo < 1.0-tol))
                {
                    frac = tlo;
                    collision = true;
                    std::cout << "Particle " << id << " is expected to collide with particle " << p.id << " frac=" << frac << std::endl;
                }
                else
                {
                    thi = 0.5 * (-b + dr);
                    if((thi > tol) && (thi < 1.0-tol))
                    {
                        frac = thi;
                        collision = true;
                        std::cout << "Particle " << id << " is expected to collide with particle " << p.id << " frac=" << frac << std::endl;
                    }
                }
            }
            std::cout << frac << std::endl;
            if(deb) std::cout << "checked a collision" << std::endl;
            return collision;
        }
        
        // Detect in-between collision with another particle p using linear trajectory approximation
        // and particle-sphere intersection method in 3D
        bool detect_collision(particle& p, double& frac)
        {
            if(deb) std::cout << "checking a collision" << std::endl;
            bool collision = false;
            
            double Reff = R + p.R;
            
            vec3 disp_vec = pos - Opos + p.Opos - p.pos;
            vec3 isep_vec = p.Opos - Opos;
            
            double disp = disp_vec.mag();
            double isep = isep_vec.mag();
            //std::cout << disp << " " << isep << " " << Reff << std::endl;
            
            if(fabs(isep-disp) <= Reff)
            {
                double isep_costheta = (disp_vec & isep_vec) / disp;
                if(isep_costheta >= sqrt(isep*isep - Reff*Reff))
                {
                    frac = (isep_costheta - sqrt(isep_costheta*isep_costheta - isep*isep + Reff*Reff)) / disp;
                    collision = true;
                    //std::cout << "Particle " << id << " is expected to collide with particle " << p.id << std::endl;
                }
            }
            
            if(deb) std::cout << "checked a collision" << std::endl;
            return collision;
        }
        
        particle operator + (particle& p) // Merging with another particle
        {
            if(deb) std::cout << "merging two particles" << std::endl;
            particle q;
            q.id = id;
            double mf = m / (m + p.m);
            // Center of mass vectors
            q.pos = pos * mf + p.pos * (1-mf);
            q.vel = vel * mf + p.vel * (1-mf);
            
            q.R = R * pow((1.0 + p.m/m), 1.0/3.0); // Assuming constant density
            q.m = m + p.m;
            //std::cout << "merged two particles" << q.R << q.m << std::endl;
            return q;
            // Delete particle p after calling this operation
        }
};

bool perform_collision(particle &p, particle &q)
{
    //std::cout << "colliding/merging two particles" << std::endl;
    bool merge = false;
    
    vec3 r21_cap = q.pos - p.pos;
    r21_cap = r21_cap * (1 / r21_cap.mag());
    
    vec3 u21_cap = q.vel - p.vel;
    double u21 = u21_cap.mag();
    u21_cap = u21_cap * (1/u21);
    
    double costheta = - (u21_cap & r21_cap);
    if(costheta == 1 || costheta == -1) costheta *= (1 - 1e-8);
    
    if(merge_or_bounce(costheta, u21)) merge = true; //use a lookup table to decide based on impact param and rel vel
    
    else //bounce
    {
        double sintheta = sqrt(1 - costheta * costheta);
        double tantheta = sintheta / costheta;
        double mratio = q.m / p.m;
        double v21, v11, cosalpha, sinalpha;
        
        if(e != mratio)
        {
            double tanalpha = (1+mratio) * sintheta / (e-mratio) / costheta;
            cosalpha = 1 / sqrt(1+tanalpha*tanalpha);
            sinalpha = tanalpha * cosalpha;
            
            v21 = u21 * sintheta / sinalpha;
            v11 = (e * costheta - sintheta / tanalpha) * u21;
        }
        else
        {
            cosalpha = 0;
            sinalpha = 1;
            
            v21 = u21 * sintheta;
            v11 = e * u21 * costheta;
        }
        
        vec3 v2, v1;
        v1 = p.vel - r21_cap * v11;
        v2 = p.vel + r21_cap * (cosalpha * v21 + sinalpha * v21 * costheta / sintheta) + u21_cap * (sinalpha * v21 / sintheta);
        
        p.update(p.pos, v1);
        q.update(q.pos, v2);
    }
    
    if(deb) std::cout << "collided/merged two particles" << std::endl;
    return merge;
}