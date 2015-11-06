// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Defines the methods to compute the score when using linear gap costs.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _internalComputeScore()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL,
                      TTraceValueR,
                      TracebackOff const &)
{
//    TScoreValue _MASK = -static_cast<TScoreValue>(activeCell._score < rightCompare);
//    activeCell._score = (rightCompare & _MASK) | (activeCell._score & ~_MASK);
    if (activeCell._score < rightCompare)
        activeCell._score = rightCompare;
    return TraceBitMap_::NONE;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<SingleTrace, TGapsPlacement> > const &)
{
//    const TScoreValue MASK = -static_cast<TScoreValue>(activeCell._score < rightCompare);
//    activeCell._score = (rightCompare & MASK) | (activeCell._score & ~MASK) ;
//    return (rightTrace & MASK) | (leftTrace & ~MASK);
    if (activeCell._score < rightCompare)
    {
        activeCell._score = rightCompare;
        return rightTrace;
    }
    return leftTrace;
}

template <typename TScoreValue, typename TTraceValueL, typename TTraceValueR, typename TGapsPlacement>
inline typename TraceBitMap_::TTraceValue
_internalComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                      TScoreValue const & rightCompare,
                      TTraceValueL leftTrace,
                      TTraceValueR rightTrace,
                      TracebackOn<TracebackConfig_<CompleteTrace, TGapsPlacement> > const &)
{
//    TScoreValue MASK = -static_cast<TScoreValue>(activeCell._score < rightCompare);
//    activeCell._score = (rightCompare & MASK) | (activeCell._score & ~MASK);
//    static typename TraceBitMap_::TTraceValue traceTmp = (rightTrace & MASK) | (leftTrace & ~MASK);
//    MASK = -static_cast<TScoreValue>(activeCell._score == rightCompare);
//    return ((traceTmp | rightTrace) & MASK) | (traceTmp & ~MASK);
    if (activeCell._score <= rightCompare)
    {
        if (activeCell._score == rightCompare)
            return leftTrace | rightTrace;
        activeCell._score = rightCompare;
        return rightTrace;
    }
    return leftTrace;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                 [RecursionDirectionAll, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
                DPCell_<TScoreValue, LinearGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionAll const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_::TTraceValue TTraceValue;
    activeCell._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);

    TScoreValue tmpScore = _scoreOfCell(previousVertical) + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);

    TTraceValue tv = _internalComputeScore(activeCell, tmpScore, TraceBitMap_::DIAGONAL, TraceBitMap_::VERTICAL | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX, TTracebackConfig());
    tmpScore = _scoreOfCell(previousHorizontal) + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    return _internalComputeScore(activeCell, tmpScore, tv, TraceBitMap_::HORIZONTAL | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX, TTracebackConfig());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionDirectionUpperDiagonal, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
                DPCell_<TScoreValue, LinearGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionUpperDiagonal const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    activeCell._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    TScoreValue tmpScore = _scoreOfCell(previousHorizontal)
                          + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);

    return _internalComputeScore(activeCell, tmpScore, TraceBitMap_::DIAGONAL, TraceBitMap_::HORIZONTAL | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX, TTracebackConfig());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore       [RecursionDirectionLowerDiagonal, LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & previousDiagonal,
                DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, LinearGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionLowerDiagonal const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    activeCell._score = _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal);
    TScoreValue tmpScore = _scoreOfCell(previousVertical)
                           + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);

    return _internalComputeScore(activeCell, tmpScore, TraceBitMap_::DIAGONAL, TraceBitMap_::VERTICAL | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX, TTracebackConfig());
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                      [RecursionDirectionHorizontal]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, LinearGaps> const & previousHorizontal,
                DPCell_<TScoreValue, LinearGaps> const & /*previousVertical*/,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionHorizontal const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    activeCell._score = _scoreOfCell(previousHorizontal)
                        + scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);

    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
        return TraceBitMap_::NONE;

    return TraceBitMap_::HORIZONTAL | TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
}

// ----------------------------------------------------------------------------
// Function _doComputeScore                        [RecursionDirectionVertical]
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_::TTraceValue
_doComputeScore(DPCell_<TScoreValue, LinearGaps> & activeCell,
                DPCell_<TScoreValue, LinearGaps> const & /*previousDiagonal*/,
                DPCell_<TScoreValue, LinearGaps> const & /*previousHorizontal*/,
                DPCell_<TScoreValue, LinearGaps> const & previousVertical,
                TSequenceHValue const & seqHVal,
                TSequenceVValue const & seqVVal,
                TScoringScheme const & scoringScheme,
                RecursionDirectionVertical const &,
                DPProfile_<TAlgorithm, LinearGaps, TTracebackConfig, TExecPolicy> const &)
{
    activeCell._score = _scoreOfCell(previousVertical)
                        + scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal);

    if (!IsTracebackEnabled_<TTracebackConfig>::VALUE)
        return TraceBitMap_::NONE;

    return TraceBitMap_::VERTICAL | TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_LINEAR_H_
