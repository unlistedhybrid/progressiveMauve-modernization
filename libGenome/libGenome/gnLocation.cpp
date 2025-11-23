// gnLocation.cpp - Modernized for C++17
#include "libGenome/gnLocation.h"
#include "libGenome/gnDebug.h"
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace genome {

gnLocation::gnLocation() {
    Clear();
}

gnLocation::gnLocation(const gnLocation& s) {
    SetBounds(s.m_start, s.m_startLength, s.m_end, s.m_endLength);
    m_type = s.m_type;
    m_name = s.m_name;
}

gnLocation::gnLocation(gnSeqI start, gnSeqI startLength, gnSeqI end, gnSeqI endLength, gnLocationType type, std::string contigName)
    : m_start(start), m_startLength(startLength), m_end(end), m_endLength(endLength), m_type(type), m_name(std::move(contigName)) {}

gnLocation::gnLocation(gnSeqI start, gnSeqI end, gnLocationType type, std::string contigName)
    : m_start(start), m_startLength(0), m_end(end), m_endLength(0), m_type(type), m_name(std::move(contigName)) {}

gnLocation* gnLocation::Clone() const {
    return new gnLocation(*this);
}

void gnLocation::Clear() {
    m_start = 0;
    m_end = 0;
    m_startLength = 0;
    m_endLength = 0;
    m_type = LT_Nothing;
    m_name.clear();
}

void gnLocation::GetBounds(gnSeqI& s, gnSeqI& sl, gnSeqI& e, gnSeqI& el) const {
    s = m_start;
    sl = m_startLength;
    e = m_end;
    el = m_endLength;
}

void gnLocation::SetBounds(gnSeqI start, gnSeqI startLength, gnSeqI end, gnSeqI endLength) {
    SetStart(start, startLength);
    SetEnd(end, endLength);
}

void gnLocation::SetBounds(gnSeqI start, gnSeqI end) {
    m_start = start;
    m_end = end;
}

bool gnLocation::CropTo(const gnLocation& l) {
    gnSeqI start = l.GetStart();
    gnSeqI end = l.GetEnd();

    if (m_start < start) {
        gnSeqI tmp = std::min(start, m_end);
        m_startLength += tmp - m_start;
        m_start = tmp;
    }
    if (m_end < end) {
        gnSeqI tmp = std::max(end, m_start);
        m_endLength += m_end - tmp;
        m_end = tmp;
    }
    if (l.GetFirst() > GetFirst()) {
        if (l.GetFirst() <= m_end)
            m_startLength = m_start - l.GetFirst();
        else if (l.GetFirst() <= GetLast()) {
            m_end = l.GetFirst();
            m_start = l.GetFirst() + 1;
            m_startLength = 0;
        } else {
            Clear();
        }
    }
    if (l.GetLast() < GetLast()) {
        if (l.GetLast() >= m_start)
            m_endLength = l.GetLast() - m_end;
        else if (l.GetLast() >= GetFirst()) {
            m_start = l.GetLast();
            m_end = l.GetLast() - 1;
            m_endLength = 0;
        } else {
            Clear();
        }
    }
    return m_start != m_end;
}

bool gnLocation::CropStart(gnSeqI start) {
    if (m_start < start) {
        gnSeqI tmp = std::min(start, m_end);
        m_startLength += tmp - m_start;
        m_start = tmp;
    }
    return m_start != m_end;
}

bool gnLocation::CropEnd(gnSeqI end) {
    if (m_end > end) {
        gnSeqI tmp = std::max(end, m_start);
        m_endLength += m_end - tmp;
        m_end = tmp;
    }
    return m_start != m_end;
}

bool gnLocation::Intersects(const gnLocation& l, intersectRestriction ir) const {
    switch (ir) {
        case determinedRegions:
            return (l.GetStart() <= m_end) && (l.GetEnd() >= m_start);
        case undeterminedRegions:
            return ((l.GetFirst() <= m_start) && (l.GetLast() >= GetFirst())) ||
                   ((l.GetFirst() <= GetLast()) && (l.GetLast() >= m_end));
        case allRegions:
            return (l.GetFirst() <= GetLast()) && (l.GetLast() >= GetFirst());
        default:
            return false;
    }
}

bool gnLocation::Contains(const gnLocation& l, intersectRestriction ir) const {
    switch (ir) {
        case determinedRegions:
            return m_start <= l.GetStart() && l.GetEnd() <= m_end;
        case undeterminedRegions:
            return (GetFirst() <= l.GetFirst() && l.GetLast() < m_start) ||
                   (m_end < l.GetFirst() && l.GetLast() <= GetLast());
        case allRegions:
            return GetFirst() <= l.GetFirst() && l.GetLast() <= GetLast();
        default:
            return false;
    }
}

bool gnLocation::MovePositive(gnSeqI diff) {
    if (m_start > GNSEQI_END - diff || m_end > GNSEQI_END - diff)
        return false;
    m_start += diff;
    m_end += diff;
    return true;
}

bool gnLocation::MoveNegative(gnSeqI diff) {
    if (m_start < diff || m_end < diff)
        return false;
    m_start -= diff;
    m_end -= diff;
    return true;
}

bool gnLocation::MoveTo(int direction, gnSeqI diff) {
    return direction > 0 ? MovePositive(diff) : MoveNegative(diff);
}

// Not implemented: Unions of ambiguous regions require interval algebra.
gnLocation gnLocation::GetUnion(const gnLocation& l) const {
    throw std::logic_error("gnLocation::GetUnion -- not implemented");
}

gnLocation gnLocation::GetIntersection(const gnLocation& l, intersectRestriction ir) const {
    gnLocation inter_loc;
    if (ir == determinedRegions) {
        if ((l.GetStart() <= m_end) && (l.GetEnd() >= m_start)) {
            inter_loc.m_start = std::max(l.m_start, m_start);
            inter_loc.m_end = std::min(l.m_end, m_end);
        }
    } else {
        throw std::logic_error("gnLocation::GetIntersection -- only determinedRegions supported");
    }
    return inter_loc;
}

} // namespace genome
