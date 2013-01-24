#include "neighbor.h"

namespace LAMMPS_NS {

int Neighbor::multi_levels()
{ return 0; }

void Neighbor::multi_levels(double &maxrad, double &minrad, int &nlevels)
{}

void Neighbor::stencil_gran_multi_3d_no_newton(NeighList *list,
                                               int sx, int sy, int sz)
{}

void Neighbor::granular_multi_no_newton(NeighList *list)
{}

}
