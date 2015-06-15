// ==========================================================================
//                               laganAlignment
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Svenja Mehringer <svenja.mehringer@fu-berlin.de>
// ==========================================================================

#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/index.h>
#include "laganAlignment_impl.h"

using namespace seqan;

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
// Function laganAlignment()                                            [Align]
// ----------------------------------------------------------------------------
/*!
 * @fn seeding
 * @brief Computes a pairwise global Alignment given two sequences using the LAGAN algorithm.
 *
 * @signature int laganAlignment(alignment, lagan_parameter, scoreSchemeAnchor,
 *                 scoreSchemeGap, alignConfig, bandExtension);
 *
 * @tparam alignment            An <a href="seqan:Align">Align</a> object that stores the alignment. The number of rows must be 2 and the sequences must have already been set. row(align, 0) is the horizontal sequence in the alignment matrix, row(align, 1) is the vertical sequence.
 * @tparam lagan_parameter      A <a href="seqan:StringSet">StringSet</a> containg a variable number of String(unsigned) container. Each container keeps three lagan paramters: length of seeds, maximum distance between seeds, maximum bandwidth between seeds.
 * @tparam scoreSchemeAnchor    The scoring scheme used for the alignment. If scoringSchemeGap is specified, then scoringSchemeAnchor is used for the regions around the seeds and scoringSchemeGap for the gap regions between two consecutive seeds. Types: Score
 * @tparam scoreSchemeGap       The optional scoring scheme for the gap regions between two anchors. Types: Score
 * @tparam alignConfig          The <a href="seqan:AlignConfig">AlignConfig</a> to use for the alignment.
 * @tparam bandExtension        Optional extension of the band around the seeds. At the moment only band extensions greater or equal 1 are allowed. Type: nolink:int. Default: 15.
 * 
 * The <tt>laganAlignment()</tt> function computes a pairwise global alignment
 * given two sequences of the same type <tt>TSequence</tt>. First an <tt>QGramIndex</tt>
 * is build over one sequenced and then the second is queried for common seeds. The 
 * retrieved seeds are then chained globally and an alignment using <tt>bandedChainAlignment</tt>
 * is computed.
 * 
 * @see QGramIndex
 * @see chainSeedsGlobally
 * @see bandedChainAlignment
 */

// given only one scoring scheme
template<typename TSequence, typename TAlignSpec, typename TScoreValue, 
         typename TScoreSpecAnchor>
int laganAlignment(Align<TSequence, TAlignSpec> & alignment,
                   String<unsigned> & lagan_parameter,
                   Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor)
{
    int result = computelaganAlignment(alignment, lagan_parameter, scoreSchemeAnchor);
    return result;
}
