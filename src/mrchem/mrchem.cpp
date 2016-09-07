/** \mainpage The MRChem program
 *
 * \author Stig Rune Jensen
 *
 * \version 1.0
 *
 * \par Copyright:
 * GPLv4
 *
 */

#include "mrchem.h"
#include "MREnv.h"
#include "SCFDriver.h"

Getkw Input;

int main(int argc, char **argv) {
    Timer rolex;
    rolex.restart();

    MREnv::initializeMRCPP(argc, argv);

    SCFDriver driver(Input);
    driver.setup();
    driver.run();
    driver.clear();

    MREnv::finalizeMRCPP(rolex);

    return 0;
}

