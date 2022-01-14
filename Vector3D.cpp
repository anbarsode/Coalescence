// Dimensionless class

class vec3 //3D Cartesian vector
{
    public:
        double x, y, z, m=0, m2=0;
        
        void update(double x_new, double y_new, double z_new, double m_new=0)
        {
            x = x_new;
            y = y_new;
            z = z_new;
            m = m_new;
        }
        
        void print() // Print the vector's components and magnitude
        {
            std::cout << "(" << x << ", " << y << ", " << z
            << "); length = " << m << std::endl;
        }
        
        double mag() // Calculate the magnitude of the vector
        {
            if(m == 0) m = sqrt(x*x + y*y + z*z);
            return m;
        }
        
        double mag2() // Calculate the magnitude of the vector
        {
            m2 = x*x + y*y + z*z;
            return m2;
        }
        
        vec3 operator - () // Negation
        {
            vec3 w;
            w.x = -x;
            w.y = -y;
            w.z = -z;
            w.m = m;
            return w;
        }  
            
        
        vec3 operator * (double s) // Scale up the vector by a constant factor s
        {
            vec3 w;
            w.x = x * s;
            w.y = y * s;
            w.z = z * s;
            w.m = m * fabs(s);
            return w;
        }    
        
        vec3 operator + (vec3 v) // Addition of two vectors
        {
            vec3 w;
            w.x = x + v.x;
            w.y = y + v.y;
            w.z = z + v.z;
            return w;
        }    
        
        vec3 operator - (vec3 v) // Subtraction of two vectors
        {
            vec3 w;
            w.x = x - v.x;
            w.y = y - v.y;
            w.z = z - v.z;
            return w;
        }
        
        double operator & (vec3 v) // Dot product of two vectors
        {
            double d;
            d = x*v.x + y*v.y + z*v.z;
            return d;
        }
        
        vec3 operator ^ (vec3 v) // Cross product of two vectors
        {
            vec3 w;
            w.x = y*v.z - z*v.y;
            w.y = z*v.x - x*v.z;
            w.z = x*v.y - y*v.z;
            return w;
        }
        
};

