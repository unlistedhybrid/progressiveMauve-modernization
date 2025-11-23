/*******************************************************************************
 * $Id: uniqueMerCount.cpp,v 1.1 2004/02/28 00:01:31 darling Exp $
 * Copyright 2002-2004 Aaron Darling. All rights reserved.
 * Please see the file COPYING for licensing and rights.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/DNAFileSML.h"
#include <iostream>
#include <memory>
#include <string>

namespace {
    void print_usage(const char* pname) {
        std::cerr << "Usage: " << pname << " <Sorted Mer List>\n";
    }
}

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        print_usage(argv[0]);
        return -1;
    }

    std::string sml_filename = argv[1];
    std::unique_ptr<genome::mems::DNAFileSML> file_sml = std::make_unique<genome::mems::DNAFileSML>();

    bool success = true;
    try {
        file_sml->LoadFile(sml_filename);
    } catch (const genome::gnException& gne) {
        success = false;
        std::cerr << gne << std::endl;
        return -1;
    }
    std::cout << std::endl << file_sml->UniqueMerCount() << std::endl;
    return 0;
}
