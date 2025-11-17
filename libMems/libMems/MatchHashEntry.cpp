/*******************************************************************************
 * $Id: Match.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MatchHashEntry.h"

using namespace std;
using namespace genome;
namespace mems {

MatchHashEntry::MatchHashEntry() : m_extended(false), m_mersize(0), m_offset(0)
{
}

MatchHashEntry::MatchHashEntry(const uint seq_count, const gnSeqI mersize, const MemType m_type)
	: m_mersize(mersize), m_offset(0)
{
	m_extended = (m_type == extended);
}

MatchHashEntry* MatchHashEntry::Clone() const
{
	return new MatchHashEntry(*this);
}

MatchHashEntry& MatchHashEntry::operator=(const MatchHashEntry& mhe)
{
	Match::operator=(mhe);
	m_extended = mhe.m_extended;
	m_mersize = mhe.m_mersize;
	m_offset = mhe.m_offset;
	return *this;
}

boolean MatchHashEntry::operator==(const MatchHashEntry& mhe) const
{
	if (!(Match::operator==(mhe)))
		return false;
	if (m_extended != mhe.m_extended)
		return false;
	if (m_mersize != mhe.m_mersize)
		return false;
	return true;
}

void MatchHashEntry::CalculateOffset()
{
	if (SeqCount() == 0)
	{
		m_offset = 0;
		return;
	}

	int64 first_start = LeftEnd(0);
	if (first_start == NO_MATCH)
		first_start = 0;
	
	m_offset = first_start;

	for (uint i = 1; i < SeqCount(); ++i)
	{
		int64 cur_start = LeftEnd(i);
		if (cur_start == NO_MATCH)
			cur_start = 0;
		
		int64 diff = cur_start - first_start;
		if (diff < 0)
			diff = -diff;
		
		m_offset += diff;
	}
}

boolean MatchHashEntry::offset_lessthan(const MatchHashEntry& a, const MatchHashEntry& b)
{
	return a.m_offset < b.m_offset;
}

boolean MatchHashEntry::start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b)
{
	if (a->FirstStart() < b->FirstStart())
		return true;
	if (a->FirstStart() > b->FirstStart())
		return false;
	return a->Offset() < b->Offset();
}

boolean MatchHashEntry::strict_start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b)
{
	for (uint i = 0; i < a->SeqCount() && i < b->SeqCount(); ++i)
	{
		int64 a_start = a->LeftEnd(i);
		int64 b_start = b->LeftEnd(i);
		
		if (a_start == NO_MATCH && b_start != NO_MATCH)
			return true;
		if (a_start != NO_MATCH && b_start == NO_MATCH)
			return false;
		if (a_start != NO_MATCH && b_start != NO_MATCH && a_start != b_start)
			return a_start < b_start;
	}
	return false;
}

int64 MatchHashEntry::end_to_start_compare(const MatchHashEntry& a, const MatchHashEntry& b)
{
	int64 a_end = a.RightEnd(0);
	int64 b_start = b.LeftEnd(0);
	
	if (a_end == NO_MATCH || b_start == NO_MATCH)
		return 0;
	
	return b_start - a_end;
}

int64 MatchHashEntry::start_compare(const MatchHashEntry& a, const MatchHashEntry& b)
{
	int64 a_start = a.LeftEnd(0);
	int64 b_start = b.LeftEnd(0);
	
	if (a_start == NO_MATCH || b_start == NO_MATCH)
		return 0;
	
	return a_start - b_start;
}

boolean MatchHashEntry::Contains(const MatchHashEntry& mhe) const
{
	if (Length() < mhe.Length())
		return false;

	for (uint i = 0; i < SeqCount(); ++i)
	{
		int64 this_start = LeftEnd(i);
		int64 mhe_start = mhe.LeftEnd(i);
		
		if (this_start == NO_MATCH && mhe_start != NO_MATCH)
			return false;
		if (this_start != NO_MATCH && mhe_start == NO_MATCH)
			continue;
		
		if (this_start != NO_MATCH && mhe_start != NO_MATCH)
		{
			int64 this_end = RightEnd(i);
			int64 mhe_end = mhe.RightEnd(i);
			
			if (mhe_start < this_start || mhe_end > this_end)
				return false;
		}
	}

	return true;
}

} // namespace mems
