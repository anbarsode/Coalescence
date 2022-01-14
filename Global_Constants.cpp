const double UL = 1.0e-5; // m
const double UT = 1.0e-3; // s
const double UM = 1.0e-12; // kg

const double Lbox = 0.1 / UL; // Box size in code units = SI / UL
const double g_eff = 9.8 * UT * UT / UL; // Effective acceleration due to gravity
                                         // after taking buoyancy into account
                                         // in code units = SI * UT * UT / UL

const double e = 1.0; // Coefficient of restitution
const double six_pi_eta = 3.342e-4 * UT * UL / UM; // Air viscosity * 6 * pi in code units = SI * UT * UL / UM

const double TG_Strength = 5.0e0 * UT / UL; // uext Taylor Green strength, m/s * UT / UL
const double TG_kx = 8 * 2 * 3.14159265 / Lbox; // uext Taylor Green kx
const double TG_kz = 8 * 2 * 3.14159265 / Lbox; // uext Taylor Green kz
const double TG_w = 2.5e1 * 2 * 3.14159265 * UT; // uext Taylor Green w

const double bounce_threshold_u12 = 1e2 * UT / UL;
const double bounce_threshold_costheta = 0.2;

const double birth_radius = 5.0e-5 / UL; // m / UL
const double death_radius = 1e-4 / UL; // Radius of droplet that will precipitate (be deleted), SI / UL
const double birth_mass = 4.0/3.0*3.14159265 * pow(birth_radius, 3.0); // Kg / UM
const double birth_speed = TG_Strength * pow((1.5e-5 / UL / birth_radius), 2.0); // 1.0e1 * UT / UL; // m/s * UT / UL

// Naming convention: base_<birthRadius>_<TGstrength>_<TGkx/2pi>
const std::string output_directory = "./snapshots_5.0_5.0_8";
const std::string output_file_base = "snapshot";
const std::string precipitation_file = "./precipitation_5.0_5.0_8.txt";
const int logging_frequency = 100;
const int output_frequency = 100;
const double current_time = 0.0;
const double timestep_size = 2.0e-5 / UT;
const long number_of_timesteps = 5e4; // 1 s of evolution
const long Nparticles = 2000; // Maintain these many particles

bool deb = false; // Debugging mode

bool merge_or_bounce(double costheta, double u21); // Make a lookup table here

vec3 u_ext_field(double t, vec3 pos); // Feed in the external velocity field as a function of spacetime (lookup table)

// R ~ 10-100 mu m, Lbox ~ 10 cm, ~1000 particles, ~1 ms timestep
