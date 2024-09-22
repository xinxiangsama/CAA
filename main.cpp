#include <iostream>
#include "DRP.h"

int main()
{   
    DRP_solution simulation(1.0, 1.0, 0.01);
    simulation.Runsim();
    std::cout <<"simulation done!"<<std::endl;
    std::getchar();
    return 0;
}