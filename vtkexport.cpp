#include <iostream>
#include <fstream>

#include "vtkexport.h"
#include "consts.h"
#include <fstream>
#include <filesystem>


double calcVonMisesStress(const Matrix3d& sigma)
{
    Matrix3d hydrst = (1./3.)*sigma.trace()*Matrix3d::Identity();
    Matrix3d der = sigma - hydrst;

    double vonmises = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            vonmises += pow(der(i,j),2);
    return sqrt(3.*vonmises/2.);
}

void ExportTimeStep(std::string meshfile,  const string outfile, const VectorXd T,
                    const VectorXd& U, const VectorXd& sigma, const VectorXd& omega, const VectorXd& eps,
                     const VectorXd& q, vector<PhaseType> Phase, int P)
{
    if (std::filesystem::exists(outfile))
            std::filesystem::remove(outfile);
        std::filesystem::copy(meshfile, outfile, std::filesystem::copy_options::update_existing);

    std::ofstream resFile(outfile, std::ios::app);

    resFile << "\nPOINT_DATA " << P << "\n";
    resFile << "Scalars T double\n";
    resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
        resFile << T(i) << "\n";
    resFile << "\n";

    resFile << "Vectors U double\n";
    //resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
        resFile << U(3*i) << " " << U(3*i+1) << " " << U(3*i+2) << "\n";
    resFile << "\n";

    resFile << "Tensors sigma double\n";
    //resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
    {
        resFile << sigma(6*i+0) << " " << sigma(6*i+3) << " " << sigma(6*i+4) << "\n"
                << sigma(6*i+3) << " " << sigma(6*i+1) << " " << sigma(6*i+5) << "\n"
                << sigma(6*i+4) << " " << sigma(6*i+5) << " " << sigma(6*i+2) << "\n";
    }
    resFile << "\n";

    resFile << "Scalars VM double\n";
    resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
    {
        Matrix3d stress;
        stress(0,0) = sigma(6*i+0);
        stress(1,1) = sigma(6*i+1);
        stress(2,2) = sigma(6*i+2);
        stress(0,1) = sigma(6*i+3);
        stress(0,2) = sigma(6*i+4);
        stress(1,2) = sigma(6*i+5);
        stress(1,0) = stress(0,1);
        stress(2,0) = stress(0,2);
        stress(2,1) = stress(1,2);
        resFile << calcVonMisesStress(stress) << "\n";
    }
    resFile << "\n";

    resFile << "Tensors omega double\n";
    //resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
    {
        resFile << omega(6*i+0)<< " " <<  omega(6*i+3) << " " << omega(6*i+4) << "\n"
                << omega(6*i+3) << " " << omega(6*i+1) << " " << omega(6*i+5) << "\n"
                << omega(6*i+4) << " " << omega(6*i+5) << " " << omega(6*i+2) << "\n";
    }
    resFile << "\n";

    resFile << "Tensors eps double\n";
    //resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
    {
        resFile << eps(6*i+0) << " " << eps(6*i+3) << " " << eps(6*i+4) << "\n"
                << eps(6*i+3) << " " << eps(6*i+1) << " " << eps(6*i+5) << "\n"
                << eps(6*i+4) << " " << eps(6*i+5) << " " << eps(6*i+2) << "\n";
    }
    resFile << "\n";

    resFile << "Scalars q double\n";
    resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
        resFile << q(i) << "\n";
    resFile << "\n";

    resFile << "Scalars phase double\n";
    resFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < P; i++)
        resFile << Phase[i] << "\n";
    resFile << "\n";

    resFile.close();
}
