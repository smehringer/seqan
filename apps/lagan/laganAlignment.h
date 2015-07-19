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
#include "laganAlignment_impl3.h"
#include "processEvents2.h"

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
// ==================================TJournaledString journal;
//    setHost(journal, host);==========================================
template<typename TJournal, typename TDeltaEvents>
int transformBack(TJournal & journal, unsigned i, TDeltaEvents & records)
{
	int vp = 0; // stores virtual position when altering the journaled string

	for (unsigned r = 0; r < length(records); ++r)
	{
		DeltaEvent & e = records[r];
		if (e.seqs[i] == true)
		{
			erase(journal, vp+e.pos, vp+endPos(e));
			if (length(e.ins)!=0)
				insert(journal, vp+e.pos, e.ins);
			vp += length(e.ins);
			vp -= e.del;
		}
	}
	return 0;
}

template<typename TDeltaEvents>
unsigned countRecords(TDeltaEvents & records)
{
	unsigned c = 0;
	for (unsigned i = 0; i < length(records); ++i)
	{
		if ((records[i].type == 0) | (records[i].type == 1))
			c += size(records[i].seqs);
		else
			c += 2*size(records[i].seqs);

	}
	return c;
}

template<typename TDeltaEvents>
int eraseZeros(TDeltaEvents & records)
{
	TDeltaEvents tmp_recs;
	for (unsigned i = 0; i < length(records); ++i)
		if (size(records[i].seqs)!=0)
			appendValue(tmp_recs, records[i]);
	records = tmp_recs;
	return 0;
}

int printEvent(DeltaEvent & e)
{
	if (e.type == 0)
		std::cout << "SNP at " << e.pos  << " = " << e.ins << " seqs: "
		<< e.seqs[0]<< e.seqs[1]<< e.seqs[2] << e.seqs[3] << std::endl;
	else if (e.type == 1)
		std::cout << "DEL at " << e.pos  << " = " << e.del << " seqs: "
		<< e.seqs[0]<< e.seqs[1]<< e.seqs[2] << e.seqs[3] << std::endl;
	else if (e.type == 2)
		std::cout << "INS at " << e.pos  << " = " << e.ins << " seqs: "
		<< e.seqs[0]<< e.seqs[1]<< e.seqs[2] << e.seqs[3] << std::endl;
	else
		std::cout << "SV at " << e.pos  << " = " << e.del << " " << e.ins << " seqs: "
		<< e.seqs[0]<< e.seqs[1]<< e.seqs[2] << e.seqs[3] << std::endl;
	return 0;
}
// ----------------------------------------------------------------------------
// Function laganAlignment()                                            [Align]
// ----------------------------------------------------------------------------

// given only one scoring scheme
template<typename TSequence /*, typename TScoreValue, typename TScoreSpecAnchor*/>
int laganAlignment(TSequence & ref, String<TSequence> & seqs,
                   String<unsigned> & lagan_parameter
                   /*, Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor*/)
{
	typedef Index<TSequence, IndexQGram<SimpleShape, OpenAddressing> > TIndex;

	std::cout << "## Building up index over reference sequence...\n";
	TIndex index(ref);
	resize(indexShape(index), lagan_parameter[0]);
	indexRequire(index, QGramSADir());
	std::cout << "## Done Building up index.\n\n";

	String<DeltaEvent> records;

	for (unsigned i = 0 ; i < length(seqs); ++i)
	{
		std::cout << "## SEQUENCE " << i <<"\n";
        getDeltaEvents(records, ref, seqs[i], i, length(seqs), index, lagan_parameter/*, scoreSchemeAnchor*/);
        std::cout << "\n";
	}

	unsigned count = countRecords(records);
	std::cout << "--------------------------------------------------------------\n";
	std::cout << "Before compression:" << count <<  " Journal Entries\n";
	std::cout << "--------------------------------------------------------------\n\n";
	sort(records, CompareByPosAndTypeLessThan_()); // todo:: needed here?
//	printEvent(records[0]);
//	for (unsigned i = 1; i < length(records); ++i)
//	{
//		if (not (isEqual(records[i],records[i-1])))
//			printEvent(records[i]);
//		else
//			std::cout << "+ 1\n";
//	}

	//std::cout << ref << std::endl;
	processDeltaEventsOnReference(records, ref);
	eraseZeros(records);

	count = countRecords(records);
	std::cout << "\n--------------------------------------------------------------\n";
	std::cout << "After compression:" << count <<  " Journal Entries\n";
	std::cout << "--------------------------------------------------------------\n\n";
//	for (unsigned i = 0; i < length(records); ++i)
//		printEvent(records[i]);

	std::cout << ref << std::endl;

	sort(records, CompareByPosAndTypeLessThan_());

	typedef typename Value<TSequence>::Type TSeqValue;
	typedef String<TSeqValue, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournaledString;

	for (unsigned i = 0; i < length(seqs); ++i)
	{
		TJournaledString journal;
		setHost(journal, ref);
		transformBack(journal, i, records);
		std::cout << journal << std::endl;
		std::cout << "Successfully coded and decoded Seq" << i << " ?:" << (journal == seqs[i]) << std::endl;
		//SEQAN_ASSERT(journal == seqs[i]);
	}
    return 0;
}
