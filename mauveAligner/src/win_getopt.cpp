#ifdef _WIN32

#include "win_getopt.h"
#include <string.h>
#include <stdio.h>

char *optarg = NULL;
int optind = 1;
int opterr = 1;
int optopt = '?';

static char *nextchar = NULL;

int getopt(int argc, char *const argv[], const char *optstring) {
    char c;
    char *cp;

    if (optind >= argc || !argv[optind])
        return -1;

    if (argv[optind][0] != '-' || argv[optind][1] == '\0')
        return -1;

    if (argv[optind][0] == '-' && argv[optind][1] == '-' && argv[optind][2] == '\0') {
        optind++;
        return -1;
    }

    if (!nextchar || *nextchar == '\0') {
        nextchar = argv[optind] + 1;
    }

    c = *nextchar++;
    cp = strchr((char *)optstring, c);

    if (!cp || c == ':') {
        if (opterr)
            fprintf(stderr, "%s: illegal option -- %c\n", argv[0], c);
        if (*nextchar == '\0') {
            optind++;
            nextchar = NULL;
        }
        return '?';
    }

    if (cp[1] == ':') {
        if (*nextchar != '\0') {
            optarg = nextchar;
            optind++;
        } else if (optind + 1 < argc) {
            optarg = argv[++optind];
            optind++;
        } else {
            if (opterr)
                fprintf(stderr, "%s: option requires an argument -- %c\n", argv[0], c);
            if (optstring[0] == ':')
                c = ':';
            else
                c = '?';
        }
        nextchar = NULL;
    } else {
        if (*nextchar == '\0') {
            optind++;
            nextchar = NULL;
        }
        optarg = NULL;
    }

    return c;
}

int getopt_long(int argc, char *const argv[], const char *optstring,
                const struct option *longopts, int *longindex) {
    if (optind >= argc || !argv[optind]) return -1;
    
    if (argv[optind][0] == '-' && argv[optind][1] == '-') {
        char *cur_arg = argv[optind] + 2;
        char *arg_val = strchr(cur_arg, '=');
        int len = (arg_val) ? (int)(arg_val - cur_arg) : (int)strlen(cur_arg);

        for (int i = 0; longopts[i].name; i++) {
            if (strncmp(cur_arg, longopts[i].name, len) == 0 && len == strlen(longopts[i].name)) {
                if (longindex) *longindex = i;
                optind++;
                
                if (longopts[i].has_arg) {
                   if (arg_val) {
                       optarg = arg_val + 1;
                   } else if (optind < argc) {
                       optarg = argv[optind++];
                   } else {
                       if (opterr) fprintf(stderr, "%s: option requires an argument -- %s\n", argv[0], longopts[i].name);
                       return '?';
                   }
                }

                if (longopts[i].flag) {
                    *longopts[i].flag = longopts[i].val;
                    return 0;
                }
                return longopts[i].val;
            }
        }
        if (opterr) fprintf(stderr, "%s: illegal option -- %s\n", argv[0], argv[optind-1]);
        optind++; 
        return '?';
    }
    
    return getopt(argc, argv, optstring);
}

#endif // _WIN32
