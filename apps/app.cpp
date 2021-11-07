/*
 * This file is part of JKQ QCEC library which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum/ for more information.
 */

#include "CompilationFlowEquivalenceChecker.hpp"
#include "EquivalenceChecker.hpp"
#include "ImprovedDDEquivalenceChecker.hpp"
#include "SimulationBasedEquivalenceChecker.hpp"

#include <algorithm>
#include <iostream>
#include <locale>
#include <string>
#include <fstream> // exp
#include <csignal> // exp
#include <sys/time.h> //estimate time

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

/* usage
size_t currentSize = getCurrentRSS( );
size_t peakSize    = getPeakRSS( );
*/

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif


/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}





/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}



std::ofstream outFile;
double fid = 0;
bool isFid;

void signalHandler(int signum) 
{
    std::cout << "Interrupt signal (" << signum << ") received.\n";

    outFile << "TO/MO" << std::endl;
    outFile.close();

    // terminate program  
    exit(signum);  
}

void show_usage(const std::string& name) {
    std::cerr << "Usage: " << name << " <PATH_TO_FILE_1> <PATH_TO_FILE_2> (--method <method>)    " << std::endl;
    std::cerr << "Supported file formats:                                                        " << std::endl;
    std::cerr << "  .real                                                                        " << std::endl;
    std::cerr << "  .qasm                                                                        " << std::endl;
    std::cerr << "  .tfc                                                                         " << std::endl;
    std::cerr << "  .qc                                                                          " << std::endl;
    std::cerr << "Available methods:                                                             " << std::endl;
    std::cerr << "  reference                                                                    " << std::endl;
    std::cerr << "  naive                                                                        " << std::endl;
    std::cerr << "  proportional (default)                                                       " << std::endl;
    std::cerr << "  lookahead                                                                    " << std::endl;
    std::cerr << "  simulation (using 'classical', 'localquantum', or 'globalquantum' stimuli)   " << std::endl;
    std::cerr << "  compilationflow                                                              " << std::endl;
    std::cerr << "Result Options:                                                                                               " << std::endl;
    std::cerr << "  --ps:                                   Print statistics                                                    " << std::endl;
    std::cerr << "  --csv:                                  Print results as csv string                                         " << std::endl;
    std::cerr << "  --storeCEXinput:                        Store counterexample input state vector (for simulation method)     " << std::endl;
    std::cerr << "  --storeCEXoutput:                       Store resulting counterexample state vectors (for simulation method)" << std::endl;
    std::cerr << "Verification Parameters:                                                                          " << std::endl;
    std::cerr << "  --tol e (default 1e-13):                Numerical tolerance used during computation             " << std::endl;
    std::cerr << "  --nsims r (default 16):                 Number of simulations to conduct (for simulation method)" << std::endl;
    std::cerr << "  --fid F (default 0.999):                Fidelity limit for comparison (for simulation method)   " << std::endl;
    std::cerr << "  --stimuliType s (default 'classical'):  Type of stimuli to use (for simulation method)          " << std::endl;
    std::cerr << "Optimization Options:                                                                             " << std::endl;
    std::cerr << "  --swapReconstruction:                   reconstruct SWAP operations                             " << std::endl;
    std::cerr << "  --singleQubitGateFusion:                fuse consecutive single qubit gates                     " << std::endl;
    std::cerr << "  --removeDiagonalGatesBeforeMeasure:     remove diagonal gates before measurements               " << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        if (argc == 4) {
            std::string cmd = argv[1];
            std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return ::tolower(c); });
            if (cmd == "--help" || cmd == "-h")
                show_usage(argv[0]);
        } else {
            show_usage(argv[0]);
        }
        return 1;
    }

    // get filenames
    std::string file1 = argv[1];
    std::string file2 = argv[2];

    outFile.open(argv[3], std::ios::app);
    signal(SIGTERM, signalHandler);

    std::string isFid_str = argv[4];
    if (isFid_str == "-f") isFid = 1;
    else isFid = 0;

    ec::Configuration config{};

    // parse configuration options
    if (argc >= 6) {
        for (int i = 5; i < argc; ++i) {
            std::string cmd = argv[i];
            std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return ::tolower(c); });

            if (cmd == "--tol") {
                ++i;
                if (i >= argc) {
                    show_usage(argv[0]);
                    return 1;
                }
                cmd = argv[i];
                std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return ::tolower(c); });
                try {
                    config.tolerance = std::stod(cmd);
                } catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    show_usage(argv[0]);
                    return 1;
                }
            } else if (cmd == "--nsims") {
                ++i;
                if (i >= argc) {
                    show_usage(argv[0]);
                    return 1;
                }
                cmd = argv[i];
                std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return ::tolower(c); });
                try {
                    config.max_sims = std::stoull(cmd);
                } catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    show_usage(argv[0]);
                    return 1;
                }
            } else if (cmd == "--fid") {
                ++i;
                if (i >= argc) {
                    show_usage(argv[0]);
                    return 1;
                }
                cmd = argv[i];
                std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return ::tolower(c); });
                try {
                    config.fidelity_limit = std::stod(cmd);
                    if (config.fidelity_limit < 0. || config.fidelity_limit > 1.) {
                        std::cerr << "Fidelity should be between 0 and 1" << std::endl;
                        show_usage(argv[0]);
                        return 1;
                    }
                } catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    show_usage(argv[0]);
                    return 1;
                }
            } else if (cmd == "--method") {
                ++i;
                if (i >= argc) {
                    show_usage(argv[0]);
                    return 1;
                }
                cmd = argv[i];
                std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return ::tolower(c); });

                // try to extract method
                if (cmd == "reference") {
                    config.method = ec::Method::Reference;
                } else if (cmd == "naive") {
                    config.method   = ec::Method::G_I_Gp;
                    config.strategy = ec::Strategy::Naive;
                } else if (cmd == "proportional") {
                    config.method   = ec::Method::G_I_Gp;
                    config.strategy = ec::Strategy::Proportional;
                } else if (cmd == "lookahead") {
                    config.method   = ec::Method::G_I_Gp;
                    config.strategy = ec::Strategy::Lookahead;
                } else if (cmd == "compilationflow") {
                    config.method   = ec::Method::G_I_Gp;
                    config.strategy = ec::Strategy::CompilationFlow;
                } else if (cmd == "simulation") {
                    config.method = ec::Method::Simulation;
                } else {
                    show_usage(argv[0]);
                    return 1;
                }
            } else if (cmd == "--stimuliType") {
                ++i;
                if (i >= argc) {
                    show_usage(argv[0]);
                    return 1;
                }
                cmd = argv[i];
                std::transform(cmd.begin(), cmd.end(), cmd.begin(), [](unsigned char c) { return ::tolower(c); });

                if (cmd == "classical") {
                    config.stimuliType = ec::StimuliType::Classical;
                } else if (cmd == "localquantum") {
                    config.stimuliType = ec::StimuliType::LocalQuantum;
                } else if (cmd == "globalquantum") {
                    config.stimuliType = ec::StimuliType::GlobalQuantum;
                } else {
                    show_usage(argv[0]);
                    return 1;
                }
            } else if (cmd == "--storeCEXinput") {
                config.storeCEXinput = true;
            } else if (cmd == "--storeCEXoutput") {
                config.storeCEXoutput = true;
            } else if (cmd == "--swapReconstruction") {
                config.reconstructSWAPs = true;
            } else if (cmd == "--singleQubitGateFusion") {
                config.fuseSingleQubitGates = true;
            } else if (cmd == "--removeDiagonalGatesBeforeMeasure") {
                config.removeDiagonalGatesBeforeMeasure = true;
            } else {
                show_usage(argv[0]);
                return 1;
            }
        }
    }

    // read circuits
    qc::QuantumComputation qc1(file1);
    qc::QuantumComputation qc2(file2);

    struct timeval t1, t2;
    double elapsedTime;
    double runtime;

    // start timer
    gettimeofday(&t1, NULL);

    // perform equivalence check
    ec::EquivalenceCheckingResults results{};
    if (config.strategy == ec::Strategy::CompilationFlow) {
        ec::CompilationFlowEquivalenceChecker ec(qc1, qc2);
        results = ec.check(config);
    } else if (config.method == ec::Method::Simulation) {
        ec::SimulationBasedEquivalenceChecker ec(qc1, qc2);
        results = ec.check(config);
    } else {
        ec::ImprovedDDEquivalenceChecker ec(qc1, qc2);
        results = ec.check(config);
    }
    results.printJSON();

    //end timer
    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
    runtime = elapsedTime / 1000;

    if (isFid)  outFile << runtime << "," << getPeakRSS() << "," << fid << std::endl;
    else outFile << runtime << "," << getPeakRSS() << std::endl;
    outFile.close();

    return 0;
}
