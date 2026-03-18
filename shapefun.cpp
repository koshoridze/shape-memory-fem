#include "shapefun.h"

double xi_i(int i) {
    if ((i == 1) || (i == 2) || (i == 5) || (i == 6))
        return 1;
    else
        return -1;
}

double eta_i(int i) {
    if ((i == 2) || (i == 3) || (i == 6) || (i == 7))
        return 1;
    else
        return -1;
}

double zeta_i(int i) {
    if ((i == 4) || (i == 5) || (i == 6) || (i == 7))
        return 1;
    else
        return -1;
}

double N(int i, double xi, double eta, double zeta)
{
    return (1+xi_i(i)*xi) * (1+eta_i(i)*eta) * (1+zeta_i(i)*zeta) / 8;
}

double dxiN(int i, double /*xi*/, double eta, double zeta)
{
    return xi_i(i) * (1+eta_i(i)*eta) * (1+zeta_i(i)*zeta) / 8;
}

double detaN(int i, double xi, double /*eta*/, double zeta)
{
    return eta_i(i) * (1+xi_i(i)*xi) * (1+zeta_i(i)*zeta) / 8;
}

double dzetaN(int i, double xi, double eta, double /*zeta*/)
{
    return zeta_i(i) * (1+eta_i(i)*eta) * (1+xi_i(i)*xi) / 8;
}
