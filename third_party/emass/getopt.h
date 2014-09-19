#ifndef GETOPT_H
#define GETOPT_H

extern char *optarg;
extern int optind;
extern int opterr;
extern int optopt;

extern int getopt(int, char **, const char *);


#endif /* GETOPT_H */

