/////////////////////////////////////////////////////////////////////////////
// File:            gnLocation.h
// Purpose:         Standard Location for Feature
// Description:     Feature location
// Author:          Aaron Darling 
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnLocation_h_
#define _gnLocation_h_

#include "libGenome/gnDefs.h"
#include <string>
#include <iostream>
#include "libGenome/gnClone.h"

namespace genome {

/**
 * This class is used to store sequence locations. gnBaseFeature uses it
 * to track feature locations. gnLocation is capable of representing any 
 * GenBank style location qualifier. gnLocation stores the start coordinate
 * and a startLength value which represents the length of an undetermined
 * region of sequence prior to the starting coordinate. A startLength of
 * GNSEQI_END represents an unknown start, whereas a startLength of 0
 * implies that the start coordinate is unambiguous. Any other value
 * indicates a range of coordinates for the starting position.
 * The range is bounded by the given start coordinate and extends
 * startLength characters upstream. Likewise, an endLength of GNSEQI_END
 * represents an unknown ambiguous end coordinate, and an endLength of 0
 * implies that the end coordinate is unambiguous. Any other value for
 * endLength denotes a range of ending coordinates starting at end and
 * continuing for endLength characters.
 */
class GNDLLEXPORT gnLocation : public gnClone
{
public:
    enum intersectRestriction {
        determinedRegions,
        undeterminedRegions,
        allRegions
    };

    enum gnLocationType {
        LT_Standard,
        LT_BetweenBases,
        LT_Complement,
        LT_Order,
        LT_Group,
        LT_OneOf,
        LT_Nothing
    };

    static constexpr gnSeqI Defined = 0;
    static constexpr gnSeqI Unknown = GNSEQI_END;

    // Constructors and assignment
    gnLocation();
    gnLocation(const gnLocation& s);
    gnLocation& operator=(const gnLocation& s) = default;
    gnLocation(gnSeqI start, gnSeqI end, gnLocationType type = LT_Standard, std::string contigName = "");
    gnLocation(gnSeqI start, gnSeqI startLength, gnSeqI end, gnSeqI endLength, gnLocationType type = LT_Standard, std::string contigName = "");

    gnLocation* Clone() const override;

    // Location operations
    void Clear();
    bool CropTo(const gnLocation& l);
    bool CropStart(gnSeqI start);
    bool CropEnd(gnSeqI end);
    bool Intersects(const gnLocation& l, intersectRestriction ir = determinedRegions) const;
    bool Contains(const gnLocation& l, intersectRestriction cr = determinedRegions) const;
    bool MovePositive(gnSeqI diff);
    bool MoveNegative(gnSeqI diff);
    bool MoveTo(int direction, gnSeqI diff);

    // Getters
    gnSeqI GetEnd() const;
    gnSeqI GetEndLength() const;
    gnSeqI GetLast() const;
    gnSeqI GetStart() const;
    gnSeqI GetStartLength() const;
    gnSeqI GetFirst() const;
    gnLocationType GetType() const;

    void GetBounds(gnSeqI& s, gnSeqI& sl, gnSeqI& e, gnSeqI& el) const;

    // Bound checks
    bool IsEndBoundLonger() const;
    bool IsStartBoundLonger() const;
    bool IsEndBoundShorter() const;
    bool IsStartBoundShorter() const;

    // Setters
    void SetEnd(gnSeqI end);
    void SetEnd(gnSeqI end, gnSeqI endLength);
    void SetEndLength(gnSeqI endLength);
    void SetStart(gnSeqI start);
    void SetStart(gnSeqI start, gnSeqI startLength);
    void SetStartLength(gnSeqI startLength);
    void SetType(gnLocationType lt);
    void SetBounds(gnSeqI start, gnSeqI startLength, gnSeqI end, gnSeqI endLength);
    void SetBounds(gnSeqI start, gnSeqI end);

    // Region algebra
    gnLocation GetUnion(const gnLocation& l) const;
    gnLocation GetIntersection(const gnLocation& l, intersectRestriction ir) const;

private:
    std::string m_name;
    gnSeqI m_start{0};
    gnSeqI m_startLength{0};
    gnSeqI m_end{0};
    gnSeqI m_endLength{0};
    gnLocationType m_type{LT_Standard};
};


// ==== Inline methods ====
inline gnSeqI gnLocation::GetEnd() const { return m_end; }
inline gnSeqI gnLocation::GetEndLength() const { return m_endLength; }
inline gnSeqI gnLocation::GetLast() const { return m_end + m_endLength; }
inline gnSeqI gnLocation::GetStart() const { return m_start; }
inline gnSeqI gnLocation::GetStartLength() const { return m_startLength; }
inline gnSeqI gnLocation::GetFirst() const { return m_start > m_startLength ? m_start - m_startLength : 0; }
inline gnLocation::gnLocationType gnLocation::GetType() const { return m_type; }
inline bool gnLocation::IsEndBoundLonger() const { return m_endLength > 0; }
inline bool gnLocation::IsStartBoundLonger() const { return m_startLength > 0; }
inline bool gnLocation::IsEndBoundShorter() const { return m_endLength == GNSEQI_END; }
inline bool gnLocation::IsStartBoundShorter() const { return m_startLength == GNSEQI_END; }
inline void gnLocation::SetEnd(gnSeqI end) { m_end = end; }
inline void gnLocation::SetEnd(gnSeqI end, gnSeqI endLength) { m_end = end; m_endLength = endLength; }
inline void gnLocation::SetEndLength(gnSeqI endLength) { m_endLength = endLength; }
inline void gnLocation::SetStart(gnSeqI start) { m_start = start; }
inline void gnLocation::SetStart(gnSeqI start, gnSeqI startLength) { m_start = start; m_startLength = startLength; }
inline void gnLocation::SetStartLength(gnSeqI startLength) { m_startLength = startLength; }
inline void gnLocation::SetType(gnLocationType lt) { m_type = lt; }

} // namespace genome

#endif // _gnLocation_h_
