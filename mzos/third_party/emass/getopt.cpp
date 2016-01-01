#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;
#include "getopt.h"

int optind = 1;             // start looking at argv[1]
int opterr = 1;             // print error messages
int optopt;
char *optarg;

static char error1[] = ": option requires an argument -- ";
static char error2[] = ": illegal option -- ";

int getopt(int argc, char *argv[], const char *opts) {

    static int sp = 1;      // initialize and preserve from call to call

    if (sp == 1) {
        if ( 
               optind >= argc              // end of argv array
               || argv[optind][0] != '-'   // not an option
               || argv[optind][1] == '\0'  // - followed by whitespace
            )
            return -1;                     // argument processing done
        else if (strcmp(argv[optind], "--") == 0) {
            optind++;       // -- is end of options flag, go to next argv
            return -1;
        }
    }

    char c = (unsigned char) argv[optind][sp];   // char after -
    optopt = c;
    const char *cp;
    if (c == ':' || (cp = strchr(opts, c)) == NULL) {
        if (opterr) {
            cerr << argv[0] << error2 << c << endl;
        }
        if (argv[optind][++sp] == '\0') {
            optind++;
            sp = 1;
        }
        return '?';         // conventional return value for unknown option
    }

    if (*++cp == ':') {     // option which requires a value to follow
        if (argv[optind][sp+1] != '\0')    // value follows without whitespace
            optarg = &argv[optind++][sp+1];
        else if (++optind >= argc) {
            if (opterr) {
                cerr << argv[0] << error1 << c << endl;
            }
            sp = 1;
            return '?';
        } else
            optarg = argv[optind++];    // set next argv as option value
        sp = 1;
    } else {                // option which requires no value
        if (argv[optind][++sp] == '\0') {
            sp = 1;
            optind++;
        }
        optarg = NULL;      // set null pointer as option value
    }
    return c;
}

