// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements the DPCell for affine gap functions. It stores the score
// values for the three matrices: diagonal, vertical and horizontal.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_CELL_CONVEX_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_CELL_CONVEX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPCell                                                    [ConvexGaps]
// ----------------------------------------------------------------------------

// Specialization for the affine gap cost function.
// This class stores three values belonging to the matrix storing the current
// maximum, the matrix for horizontal gaps and the matrix for vertical gaps.
template <typename TScoreValue>
class DPCell_<TScoreValue, ConvexGaps>
{
public:
    TScoreValue _score               = DPCellDefaultInfinity<DPCell_>::VALUE;
    TScoreValue _horizontalScore     = DPCellDefaultInfinity<DPCell_>::VALUE;
    TScoreValue _verticalScore       = DPCellDefaultInfinity<DPCell_>::VALUE;
    TScoreValue _horizontalGapLength = 1;
    TScoreValue _verticalGapLength   = 1;

    DPCell_() = default;
    // Copy c'tor.
    DPCell_(DPCell_<TScoreValue, ConvexGaps> const & other) :
        _score(other._score),
        _horizontalScore(other._horizontalScore),
        _verticalScore(other._verticalScore),
        _horizontalGapLength(other._horizontalGapLength),
        _verticalGapLength(other._verticalGapLength)
    {}

    // Move c-tor
    DPCell_(DPCell_ && other) : DPCell_()
    {
        swap(*this, other);
    }

    // The assignment operator.
    DPCell_ &
    operator=(DPCell_ other)
    {
        swap(*this, other);
        return *this;
    }

    // Assign score to cell.
    DPCell_ &
    operator=(TScoreValue const & score)
    {
        _score = score;
        return *this;
    }

    ~DPCell_() = default;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

// Needed for banded chain alignment for the set.
template <typename TStream, typename TScore>
inline TStream& operator<<(TStream & stream,
                           DPCell_<TScore, ConvexGaps> const & dpCell)
{
    stream << "<S = " << dpCell._score << " H = " << dpCell._horizontalScore
           << " V = " << dpCell._verticalScore
           << " HL = " << dpCell._horizontalGapLength
           << " VL = " << dpCell._verticalGapLength  << ">";
    return stream;
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

// Needed for banded chain alignment for the set.
template <typename TScoreValueLeft, typename TScoreValueRight>
inline bool operator<(DPCell_<TScoreValueLeft, ConvexGaps> const & left,
                      DPCell_<TScoreValueRight, ConvexGaps> const & right)
{
    return left._score < right._score && left._horizontalScore < right._horizontalScore &&
           left._verticalScore < right._verticalScore;
}

// ----------------------------------------------------------------------------
// Function _verticalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> >::Type
_verticalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell)
{
    return dpCell._verticalScore;
}

template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> const>::Type
_verticalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> const & dpCell)
{
    return dpCell._verticalScore;
}

// ----------------------------------------------------------------------------
// Function _verticalGapLengthOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> >::Type
_verticalGapLengthOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell)
{
    return dpCell._verticalGapLength;
}

template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> const>::Type
_verticalGapLengthOfCell(DPCell_<TScoreValue, ConvexGaps> const & dpCell)
{
    return dpCell._verticalGapLength;
}

// ----------------------------------------------------------------------------
// Function _setVerticalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >,void)
_setVerticalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell, TScoreValue const & newVerticalScore)
{
    dpCell._verticalScore = newVerticalScore;
}

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >,void)
_setVerticalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell, TScoreValue const & newVerticalScore)
{
    dpCell._verticalScore = newVerticalScore;
}

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >,void)
_setVerticalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell, TScoreValue const & newVerticalScore, TScoreValue const & mask)
{
    dpCell._verticalScore = blend(dpCell._verticalScore, newVerticalScore, mask);
}

// ----------------------------------------------------------------------------
// Function _horizontalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for horizontal-gaps of the given cell.
template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> >::Type
_horizontalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell)
{
    return dpCell._horizontalScore;
}

template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> const>::Type
_horizontalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> const & dpCell)
{
    return dpCell._horizontalScore;
}

// ----------------------------------------------------------------------------
// Function _horizontalGapLengthOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for horizontal-gaps of the given cell.
template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> >::Type
_horizontalGapLengthOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell)
{
    return dpCell._horizontalGapLength;
}

template <typename TScoreValue>
inline typename  Reference<DPCell_<TScoreValue, ConvexGaps> const>::Type
_horizontalGapLengthOfCell(DPCell_<TScoreValue, ConvexGaps> const & dpCell)
{
    return dpCell._horizontalGapLength;
}

// ----------------------------------------------------------------------------
// Function _setHorizontalScoreOfCell()
// ----------------------------------------------------------------------------

// Returns the score of the matrix for vertical-gaps of the given cell.
template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Not<Is<SimdVectorConcept<TScoreValue> > >,void)
_setHorizontalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell, TScoreValue const & newHorizontalScore)
{
    dpCell._horizontalScore = newHorizontalScore;
}

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >,void)
_setHorizontalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell, TScoreValue const & newHorizontalScore)
{
    dpCell._horizontalScore = newHorizontalScore;
}

template <typename TScoreValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TScoreValue> >,void)
_setHorizontalScoreOfCell(DPCell_<TScoreValue, ConvexGaps> & dpCell, TScoreValue const & newHorizontalScore, TScoreValue const & mask)
{
    dpCell._horizontalScore = blend(dpCell._horizontalScore, newHorizontalScore, mask);
}

template <typename TScoreValue>
inline void
swap(DPCell_<TScoreValue, ConvexGaps> & lhs,
     DPCell_<TScoreValue, ConvexGaps> & rhs)
{
    std::swap(lhs._score, rhs._score);
    std::swap(lhs._horizontalScore, rhs._horizontalScore);
    std::swap(lhs._verticalScore, rhs._verticalScore);
    std::swap(lhs._horizontalGapLength, rhs._horizontalGapLength);
    std::swap(lhs._verticalGapLength, rhs._verticalGapLength);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_CELL_AFFINE_H_
