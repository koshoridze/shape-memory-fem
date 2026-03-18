#pragma once

enum PhaseType {
    DirectPhase,  // forward phase transformation (A -> M)
    ReversePhase  // reverse phase transformation (M -> A)
};

/* Material parameters */
constexpr double Ea = 70*1e9; // Young's modulus of austenite phase
constexpr double Em = 40*1e9; // Young's modulus of martensite phase
constexpr double alphaA = 1.2*1e-5; // Thermal expansion coefficient of austenite
constexpr double alphaM = 0.6*1e-6; // Thermal expansion coefficient of martensite
constexpr double nu = 0.3; // Poisson's ratio
constexpr double dS = 0.644*1e6; // Difference in volumetric entropy density between austenite and martensite
constexpr double Ms = 40; // Start temperature of forward transformation (A -> M)
constexpr double Mf = 20; // Finish temperature of forward transformation (A -> M)
constexpr double As = 50; // Start temperature of reverse transformation (M -> A)
constexpr double Af = 70; // Finish temperature of reverse transformation (M -> A)
constexpr double sigma_0 = 1e8; // Threshold stress for microstress distribution density in a representative volume
constexpr double rhoD = 0.14; // Maximum transformation strain intensity (forward transformation)
constexpr double eps0 = 0.001; // Linear strain associated with volumetric effect of forward martensitic transformation

constexpr double Fa = Ea / (1-2*nu);
constexpr double Fm = Em / (1-2*nu);
constexpr double Ga = Ea / (2+2*nu);
constexpr double Gm = Em / (2+2*nu);

constexpr double q0 = 0; // Initial martensite volume fraction
constexpr PhaseType phase0 = DirectPhase; // Initial phase
constexpr double Fz = 1e8;
constexpr double bz0 = 0;

constexpr double Time_end = 0.5;
constexpr int N1 = 201;
constexpr int N2 = 101;
constexpr int N3 = 101;
constexpr int NTime = N1+N2+N3;
constexpr double T0 = 100;
constexpr double T1 = 10;
constexpr double T2 = 50;
constexpr double T3 = 100;
