/*******************************************************************************
 * $Id: RepeatHashCat.cpp,v 1.2 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "RepeatHashCat.h"

namespace mems {

RepeatHashCat::RepeatHashCat()
	: RepeatHash()
{
}

RepeatHashCat::~RepeatHashCat()
{
}

RepeatHashCat::RepeatHashCat(const RepeatHashCat& rhc)
	: RepeatHash(rhc)
{
	concat_contig_start = rhc.concat_contig_start;
}

RepeatHashCat& RepeatHashCat::operator=(const RepeatHashCat& rhc)
{
	if(this != &rhc)
	{
		RepeatHash::operator=(rhc);
		concat_contig_start = rhc.concat_contig_start;
	}
	return *this;
}

RepeatHashCat* RepeatHashCat::Clone() const
{
	return new RepeatHashCat(*this);
}

} // namespace mems
