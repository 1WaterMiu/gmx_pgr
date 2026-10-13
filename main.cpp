// SPDX-License-Identifier: LGPL-2.1-or-later
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/exceptions.h"
gmx::TrajectoryAnalysisModulePointer createPgr();
int main(int argc, char** argv) {
    bool pgr=argc>1 && std::strcmp(argv[1],"pgr")==0;
    bool help=argc>2 && std::strcmp(argv[1],"help")==0 && std::strcmp(argv[2],"pgr")==0;
    if(!pgr && !help) {
        execv(GMX_ORIGINAL,argv);
        std::perror("Cannot run original GROMACS");
        return 127;
    }
    auto& context=gmx::initForCommandLine(&argc,&argv);
    try {
        gmx::CommandLineModuleManager manager("gmx",&context);
        gmx::TrajectoryAnalysisCommandLineRunner::registerModule(&manager,"pgr","Volume-corrected proximal radial distribution function",&createPgr);
        int rc=manager.run(argc,argv);
        gmx::finalizeForCommandLine();
        return rc;
    } catch(const std::exception& e) {
        gmx::printFatalErrorMessage(stderr,e);
        return gmx::processExceptionAtExitForCommandLine(e);
    }
}
