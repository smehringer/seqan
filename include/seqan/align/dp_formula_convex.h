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
// Implements the affine gap cost functions.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_CONVEX_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_CONVEX_H_

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
// Function _computeScore                 [RecursionAllDirection, ConvexGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, ConvexGaps> & current,
              DPCell_<TScoreValue, ConvexGaps> & previousDiagonal,
              DPCell_<TScoreValue, ConvexGaps> const & previousHorizontal,
              DPCell_<TScoreValue, ConvexGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionAll const &,
              DPProfile_<TAlgorithm, ConvexGaps, TTracebackConfig, TExecPolicy> const &)
{
    SEQAN_ASSERT(_horizontalGapLengthOfCell(current) != 0);
    // Compute intermediate diagonal result.
    TScoreValue intermediate = static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                                        score(scoringScheme, seqHVal, seqVVal));
    // Cache previous Diagonal
    _scoreOfCell(previousDiagonal) = _scoreOfCell(previousHorizontal);
    // Compute horizontal direction.
    auto tmp  = _maxScore(_horizontalScoreOfCell(current),
                          _horizontalScoreOfCell(previousHorizontal) +
                              std::min(scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal),
                                       scoreGapOpenHorizontal(scoringScheme, seqHVal, seqHVal)/_horizontalGapLengthOfCell(previousHorizontal)),
                          _scoreOfCell(previousHorizontal) +
                              scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal),
                          TraceBitMap_<TScoreValue>::HORIZONTAL,
                          TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN,
                          TTracebackConfig{});

    if (tmp & TraceBitMap_<TScoreValue>::HORIZONTAL)
        _horizontalGapLengthOfCell(current) = _horizontalGapLengthOfCell(previousHorizontal) + 1;
    else
        _horizontalGapLengthOfCell(current) = 1;

    // Compute vertical direction.
    tmp |= _maxScore(_verticalScoreOfCell(previousVertical),
                     _verticalScoreOfCell(previousVertical) +
                              std::min(scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal),
                                       scoreGapOpenVertical(scoringScheme, seqHVal, seqHVal)/_verticalGapLengthOfCell(previousVertical)),
                     _scoreOfCell(previousVertical) +
                         scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal),
                     TraceBitMap_<TScoreValue>::VERTICAL,
                     TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                     TTracebackConfig{});


    if (tmp & TraceBitMap_<TScoreValue>::VERTICAL)
        _verticalGapLengthOfCell(current) = _verticalGapLengthOfCell(previousVertical) + 1;
    else
        _verticalGapLengthOfCell(current) = 1;


    auto tmp2 = _maxScore(_scoreOfCell(current),
                    _verticalScoreOfCell(previousVertical),
                    _horizontalScoreOfCell(current),
                    TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                    TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                    TTracebackConfig{});
    tmp = _maxScore(_scoreOfCell(current),
                    intermediate,
                    _scoreOfCell(current),
                    TraceBitMap_<TScoreValue>::DIAGONAL | tmp,
                    tmp2 | tmp,
                    TTracebackConfig{});
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tmp = _maxScore(_scoreOfCell(current),
                        TraceBitMap_<TScoreValue>::NONE,
                        _scoreOfCell(current),
                        TraceBitMap_<TScoreValue>::NONE,
                        tmp,
                        TTracebackConfig{});
    }
    // Cache score for previous vertical.
    _scoreOfCell(previousVertical) = _scoreOfCell(current);
    _verticalGapLengthOfCell(previousVertical) = _verticalGapLengthOfCell(current);
    return tmp;
}

// ----------------------------------------------------------------------------
// Function _computeScore       [RecursionUpperDiagonalDirection, ConvexGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, ConvexGaps> & current,
              DPCell_<TScoreValue, ConvexGaps> & previousDiagonal,
              DPCell_<TScoreValue, ConvexGaps> const & previousHorizontal,
              DPCell_<TScoreValue, ConvexGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionUpperDiagonal const &,
              DPProfile_<TAlgorithm, ConvexGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;
    SEQAN_ASSERT(_horizontalGapLengthOfCell(current) != 0);

    // Compute intermediate diagonal result.
    TScoreValue intermediate = static_cast<TScoreValue>(_scoreOfCell(previousDiagonal) +
                                                        score(scoringScheme, seqHVal, seqVVal));
    // Cache previous Diagonal
    _scoreOfCell(previousDiagonal) = _scoreOfCell(previousHorizontal);
    // Compute horiztonal direction.
    // _horizontalScoreOfCell(current) = _horizontalScoreOfCell(previousHorizontal) +
    //                                      scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal);
    TTraceValue tv = _maxScore(_horizontalScoreOfCell(current),
                               _horizontalScoreOfCell(previousHorizontal) +
                                    std::min(scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal),
                                             scoreGapOpenHorizontal(scoringScheme, seqHVal, seqHVal)/_horizontalGapLengthOfCell(previousHorizontal)),
                               _scoreOfCell(previousHorizontal) +
                                    scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal),
                               TraceBitMap_<TScoreValue>::HORIZONTAL,
                               TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN,
                               TTracebackConfig());

    if (tv & TraceBitMap_<TScoreValue>::HORIZONTAL)
        _horizontalGapLengthOfCell(current) = _horizontalGapLengthOfCell(previousHorizontal) + 1;
    else
        _horizontalGapLengthOfCell(current) = 1;

    // Ignore vertical direction in upper diagonal.
    _verticalScoreOfCell(previousVertical) = DPCellDefaultInfinity<DPCell_<TScoreValue, ConvexGaps> >::VALUE;
    // Compute diagonal direction and compare with horizontal.
    // _scoreOfCell(current) = _horizontalScoreOfCell(current);
    tv = _maxScore(_scoreOfCell(current),
                   intermediate,
                   _horizontalScoreOfCell(current),
                   tv | TraceBitMap_<TScoreValue>::DIAGONAL,
                   tv | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX,
                   TTracebackConfig());
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tv = _maxScore(_scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TTracebackConfig{});
    }
    _scoreOfCell(previousVertical) = _scoreOfCell(current);
    _verticalGapLengthOfCell(previousVertical) = _verticalGapLengthOfCell(current);
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore       [RecursionDirectionLowerDiagonal, ConvexGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, ConvexGaps> & current,
              DPCell_<TScoreValue, ConvexGaps> const & previousDiagonal,
              DPCell_<TScoreValue, ConvexGaps> const & /*previousHorizontal*/,
              DPCell_<TScoreValue, ConvexGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionLowerDiagonal const &,
              DPProfile_<TAlgorithm, ConvexGaps, TTracebackConfig, TExecPolicy> const &)
{
    typedef typename TraceBitMap_<TScoreValue>::Type TTraceValue;
    SEQAN_ASSERT(_horizontalGapLengthOfCell(current) != 0);

    // Compute vertical direction.
    TTraceValue tv = _maxScore(_verticalScoreOfCell(previousVertical),
                               _verticalScoreOfCell(previousVertical) +
                              std::min(scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal),
                                       scoreGapOpenVertical(scoringScheme, seqHVal, seqHVal)/_verticalGapLengthOfCell(previousVertical)),
                               _scoreOfCell(previousVertical) +
                                    scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal),
                              TraceBitMap_<TScoreValue>::VERTICAL,
                              TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                              TTracebackConfig());


    if (tv & TraceBitMap_<TScoreValue>::VERTICAL)
        _verticalGapLengthOfCell(current) = _verticalGapLengthOfCell(previousVertical) + 1;
    else
        _verticalGapLengthOfCell(current) = 1;


    // Ignore horizontal direction in lower diagonal.
    _horizontalScoreOfCell(current) = DPCellDefaultInfinity<DPCell_<TScoreValue, ConvexGaps> >::VALUE;

    // Compute diagonal direction and compare with vertical.
    tv = _maxScore(_scoreOfCell(current),
                   _scoreOfCell(previousDiagonal) + score(scoringScheme, seqHVal, seqVVal),
                   _verticalScoreOfCell(previousVertical),
                   tv | TraceBitMap_<TScoreValue>::DIAGONAL,
                   tv | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX,
                   TTracebackConfig());
    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        tv = _maxScore(_scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       _scoreOfCell(current),
                       TraceBitMap_<TScoreValue>::NONE,
                       tv,
                       TTracebackConfig{});
    }
    return tv;
}

// ----------------------------------------------------------------------------
// Function _computeScore                      [RecursionHorizontalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, ConvexGaps> & current,
              DPCell_<TScoreValue, ConvexGaps> & previousDiagonal,
              DPCell_<TScoreValue, ConvexGaps> const & previousHorizontal,
              DPCell_<TScoreValue, ConvexGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionHorizontal const &,
              DPProfile_<TAlgorithm, ConvexGaps, TTracebackConfig, TExecPolicy> const &)
{
    SEQAN_ASSERT(_horizontalGapLengthOfCell(current) != 0);
    // Cache previous diagonal value.
    _scoreOfCell(previousDiagonal) = _scoreOfCell(previousHorizontal);
    // Compute horizontal direction.
    auto traceDir = _maxScore(_horizontalScoreOfCell(current),
                              _horizontalScoreOfCell(previousHorizontal) +
                                std::min(scoreGapExtendHorizontal(scoringScheme, seqHVal, seqVVal),
                                         scoreGapOpenHorizontal(scoringScheme, seqHVal, seqHVal)/_horizontalGapLengthOfCell(previousHorizontal)),
                              _scoreOfCell(previousHorizontal) +
                                scoreGapOpenHorizontal(scoringScheme, seqHVal, seqVVal),
                              TraceBitMap_<TScoreValue>::HORIZONTAL,
                              TraceBitMap_<TScoreValue>::HORIZONTAL_OPEN,
                              TTracebackConfig()) | TraceBitMap_<TScoreValue>::MAX_FROM_HORIZONTAL_MATRIX;

    if (traceDir & TraceBitMap_<TScoreValue>::HORIZONTAL)
        _horizontalGapLengthOfCell(current) = _horizontalGapLengthOfCell(previousHorizontal) + 1;
    else
        _horizontalGapLengthOfCell(current) = 1;

    // Ignore vertical direction.
    _verticalScoreOfCell(previousVertical) = DPCellDefaultInfinity<DPCell_<TScoreValue, ConvexGaps> >::VALUE;
    _scoreOfCell(current) = _horizontalScoreOfCell(current);

    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        traceDir = _maxScore(_scoreOfCell(current),
                             TraceBitMap_<TScoreValue>::NONE,
                             _scoreOfCell(current),
                             TraceBitMap_<TScoreValue>::NONE,
                             traceDir,
                             TTracebackConfig{});
    }
    _scoreOfCell(previousVertical) = _scoreOfCell(current);
    _verticalGapLengthOfCell(previousVertical) = _verticalGapLengthOfCell(current);
    return traceDir;
}

// ----------------------------------------------------------------------------
// Function _computeScore                        [RecursionVerticalDirection]
// ----------------------------------------------------------------------------

template <typename TScoreValue,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme,
          typename TAlgorithm, typename TTracebackConfig, typename TExecPolicy>
inline typename TraceBitMap_<TScoreValue>::Type
_computeScore(DPCell_<TScoreValue, ConvexGaps> & current,
              DPCell_<TScoreValue, ConvexGaps> const & /*previousDiagonal*/,
              DPCell_<TScoreValue, ConvexGaps> const & /*previousHorizontal*/,
              DPCell_<TScoreValue, ConvexGaps> & previousVertical,
              TSequenceHValue const & seqHVal,
              TSequenceVValue const & seqVVal,
              TScoringScheme const & scoringScheme,
              RecursionDirectionVertical const &,
              DPProfile_<TAlgorithm, ConvexGaps, TTracebackConfig, TExecPolicy> const &)
{
    SEQAN_ASSERT(_horizontalGapLengthOfCell(current) != 0);
    // Compute vertical direction.
    auto traceDir = _maxScore(_verticalScoreOfCell(previousVertical),
                              _verticalScoreOfCell(previousVertical) +
                              std::min(scoreGapExtendVertical(scoringScheme, seqHVal, seqVVal),
                                       scoreGapOpenVertical(scoringScheme, seqHVal, seqHVal)/_verticalGapLengthOfCell(previousVertical)),
                              _scoreOfCell(previousVertical) +
                                scoreGapOpenVertical(scoringScheme, seqHVal, seqVVal),
                              TraceBitMap_<TScoreValue>::VERTICAL,
                              TraceBitMap_<TScoreValue>::VERTICAL_OPEN,
                              TTracebackConfig()) | TraceBitMap_<TScoreValue>::MAX_FROM_VERTICAL_MATRIX;


    if (traceDir & TraceBitMap_<TScoreValue>::VERTICAL)
        _verticalGapLengthOfCell(current) = _verticalGapLengthOfCell(previousVertical) + 1;
    else
        _verticalGapLengthOfCell(current) = 1;


    // Ignore horizontal direction.
    _horizontalScoreOfCell(current) = DPCellDefaultInfinity<DPCell_<TScoreValue, ConvexGaps> >::VALUE;
    _scoreOfCell(current) = _verticalScoreOfCell(previousVertical);

    if (IsLocalAlignment_<TAlgorithm>::VALUE)
    {
        traceDir = _maxScore(_scoreOfCell(current),
                             TraceBitMap_<TScoreValue>::NONE,
                             _scoreOfCell(current),
                             TraceBitMap_<TScoreValue>::NONE,
                             traceDir,
                             TTracebackConfig{});
    }
    _scoreOfCell(previousVertical) = _scoreOfCell(current);
    _verticalGapLengthOfCell(previousVertical) = _verticalGapLengthOfCell(current);
    return traceDir;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_FORMULA_AFFINE_H_
