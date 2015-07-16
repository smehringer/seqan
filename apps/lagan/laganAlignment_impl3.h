// ==========================================================================
//                             laganAlignment_imp
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
// This File implements the LAGAN Algorithm.
// The laganAlignment() function computes a pairwise global alignment
// given two sequences of the same type. First an QGramIndex
// is build over one sequenced and then the second is queried for common seeds.
// The retrieved seeds are then chained globally and an alignment using
// bandedChainAlignment is computed.
// ==========================================================================

#ifndef SEQAN_DEMOS_MINILAGAN_LAGANALIGNMENT_IMPL2_H
#define SEQAN_DEMOS_MINILAGAN_LAGANALIGNMENT_IMPL2_H

#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/index.h>
#include <seqan/sequence_journaled.h>
#include "processEvents.h"

#include <time.h>
#include <omp.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


// ----------------------------------------------------------------------------
// Struct SearchField
// ----------------------------------------------------------------------------
/*
 * Stores the boundaries of a search field, starting with the end position of a
 * seed to the start position of the next in both sequences.
 */
struct SearchField
{
    // positions of search field
    unsigned beginH, endH, beginV, endV;

    SearchField(unsigned bH, unsigned bV, unsigned eH, unsigned eV) :
        beginH(bH), endH(eH), beginV(bV), endV(eV)
    {}
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function timestamp()
// ----------------------------------------------------------------------------
// Only to inform time measurements on std::cout
void timestamp()
{
    time_t ltime;
    ltime=time(NULL);
    std::cout << "\t\t\t\t" << asctime(localtime(&ltime));
}

// ----------------------------------------------------------------------------
// Function seeding()
// ----------------------------------------------------------------------------
/*
 * Given two sequence (fragments as infixes) the functions builds an index over
 * the reference and afterwards occurences of common subsequences of size q are added
 * to a seedSet.
 */
template<typename TSeed, typename TSeq, typename TInfix, typename TIndexTag>
int nativeSeeding(SeedSet<TSeed> & seedSet, TSeq & ref,
            TInfix & seq, unsigned q, TIndexTag /*tag*/)
{
    typedef Index<TSeq, IndexQGram<SimpleShape, TIndexTag> > TIndex;
    typedef typename Iterator<TInfix>::Type TIterator;
    typedef String<typename SAValue<Index<TSeq, TIndex> >::Type> TOccurrences;

    unsigned distance = 0;

    TIndex index(ref);
    resize(indexShape(index),q);

    hashInit(indexShape(index), begin(seq));
    for(TIterator it = begin(seq); it != (end(seq)-q+1); goNext(it))
    {
    	unsigned repeat_limit = 0;
        hashNext(indexShape(index), it);
        TOccurrences occs = getOccurrences(index, indexShape(index));
        for (unsigned i = 0; i < length(occs); ++i)
        {
        	if ((beginPosition(ref) + occs[i] > repeat_limit) | (beginPosition(ref)+occs[i] == 0))
        	{
				TSeed seed = TSeed(beginPosition(ref) + occs[i], beginPosition(seq) + position(it, seq), q);
				extendSeed(seed, ref, seq, EXTEND_RIGHT, MatchExtend());
				if (endPositionV(seed)>endPosition(seq))
					setEndPositionV(seed, endPosition(seq));
				if (endPositionH(seed)>endPosition(ref))
					setEndPositionH(seed, endPosition(ref));
				if (!addSeed(seedSet, seed, distance, Merge()))
					addSeed(seedSet, seed, Single());
        	}
        	repeat_limit = beginPosition(ref)+occs[i]+q; // q = length of initial seed
        }
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function parallelSeeding()
// ----------------------------------------------------------------------------
/*
 * Given two sequence (fragments as infixes) the functions builds an index over
 * the reference and afterwards occurrences of common subsequences of size q are added
 * to a seedSet.
 */
template<typename TSeed, typename TSeq, typename TInfix, typename TIndexTag>
int parallelSeeding(SeedSet<TSeed> & seedSet, TSeq & ref, TInfix & seq,
                    unsigned q, unsigned num_threads, TIndexTag /*tag*/)
{
    typedef Index<TSeq, IndexQGram<SimpleShape, TIndexTag> > TIndex;
    typedef typename Iterator<TInfix>::Type TIterator;
    typedef typename Iterator<SeedSet<TSeed> >::Type TIter;
    typedef String<typename SAValue<Index<TSeq, TIndex> >::Type> TOccurrences;

    unsigned distance = 0;

    if (int(length(seq)/num_threads) < q)
    	num_threads = 1;

    //set up for parallelism
    String<SeedSet<TSeed> > tmp_sets;
    resize(tmp_sets, num_threads);
    omp_set_num_threads(num_threads);
    TIndex index(ref);
    resize(indexShape(index), q);

    // each thread fills a temporary seedSet
    #pragma omp parallel for
    for (unsigned t = 0; t < num_threads; ++t)
    {
        typename Fibre<TIndex, QGramShape>::Type shape = indexShape(index);

        TIterator begin_it = begin(seq)+ t*(int)(length(seq)/num_threads);
        TIterator end_it;
        if (t != num_threads-1)
            end_it = begin_it+ (int)(length(seq)/num_threads);
        else
            end_it = end(seq) - q+1;

        //hashInit(shape, begin_it);
        for(TIterator it = begin_it; it != end_it; it += q)
        {
            hash(shape, it);
            TOccurrences occs = getOccurrences(index, shape);
            for (unsigned i = 0; i < length(occs); ++i)
            {
                TSeed seed = TSeed(beginPosition(ref) + occs[i], beginPosition(seq) + position(it, seq), q);
                if (!addSeed(tmp_sets[t], seed, distance, Merge()))
                    addSeed(tmp_sets[t], seed, Single());
            }
        }
    }

    // combine seedSets
    for (unsigned t = 0; t < num_threads; ++t)
    {
        //std::cout << length(tmp_sets[t]) << std::endl;
        for (TIter it = begin(tmp_sets[t], Standard()); it != end(tmp_sets[t], Standard()); ++it)
        {

            TSeed seed = *it;
            // if begin position of seed lies whithin an overlay area it must be merged
            if (beginPositionV(seed) < (t * (int)(length(seq)/num_threads) + q + distance))
            {
                if (!addSeed(seedSet, seed, distance, Merge()))
                     addSeed(seedSet, seed, Single());
            }
            else
            {
                addSeed(seedSet, seed, Single());
            }
        }

    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function scoreSeed()
// ----------------------------------------------------------------------------
template<typename TSeed>
long scoreSeed(TSeed & seed)
{
	long n = (long)seedSize(seed);
	SEQAN_ASSERT_EQ(upperDiagonal(seed), lowerDiagonal(seed));

	long score = n - (long)abs(lowerDiagonal(seed));
	return (score);
}

// ----------------------------------------------------------------------------
// Function fastFirstSeeding()
// ----------------------------------------------------------------------------
template<typename TSeed, typename TIndex, typename TSeq, typename TInfix>
int parallelFastFirstSeeding(SeedSet<TSeed> & seedSet, TIndex & index, TSeq & ref,
		             TInfix & seq, unsigned q, unsigned num_threads)
{
	typedef typename Iterator<TInfix>::Type TIterator;
    typedef typename Iterator<SeedSet<TSeed> >::Type TIter;
	typedef String<typename SAValue<Index<TSeq, TIndex> >::Type> TOccurrences;

	num_threads += num_threads;
	unsigned kband = (length(seq)/3 + length(ref)/3)/2;

	if (int(length(seq)/num_threads) < q)
		num_threads = 1;

	//set up for parallelism
	String<SeedSet<TSeed> > tmp_sets;
	resize(tmp_sets, num_threads);
	omp_set_num_threads(num_threads/2); // set to actual number of threads (without doubling)

	 // each thread fills a temporary seedSet
	#pragma omp parallel for
	for (unsigned t = 0; t < (num_threads); ++t)
	{
		typename Fibre<TIndex, QGramShape>::Type shape = indexShape(index);

		TIterator begin_it = begin(seq)+ t*(int)(length(seq)/num_threads);
		TIterator end_it;
		if (t != num_threads-1)
			end_it = begin_it + (int)(length(seq)/num_threads);
		else
			end_it = end(seq) - q+1;

		while (position(begin_it, seq) < position(end_it, seq))
		{
			unsigned offset = 1;
			//std::cout << "\r" << "\t" << t << "\t" << position(begin_it,seq) << std::endl;
			//if (not(infix(seq, position(it, seq), position(it, seq)+q) == "NNNNNNNNNNNNNNNNNNNN"))
			//{
				hash(shape, begin_it);
				TOccurrences occs = getOccurrences(index, shape);
				long max_score = -maxValue<long>();
				unsigned repeat_limit = 0;
				// if there are no hits the offset is 1 otherwise its the end position
				// of the best scoring seed.
				// This way long runs will speed up seeding and result in less memory
				// usage.
				//unsigned len = length(occs);
				for (unsigned i = 0; i < length(occs); ++i)
				{
					// avoid repeats:
					// a repeat can be identified if the subsequent seed is found within
					// the preceding seed (endPos(pre_seed) > beginPos(post_seed)
					// the next if clause therefore is avoiding repeats
					if (((beginPosition(ref)+occs[i] > repeat_limit) || (beginPosition(ref)+occs[i] == 0)) && (abs(position(begin_it,seq)-occs[i]) < kband))
					{
						TSeed seed = TSeed(beginPosition(ref) + occs[i], beginPosition(seq) + position(begin_it, seq), q);
						extendSeed(seed, ref, seq, EXTEND_BOTH, MatchExtend());
						
						if (seedSize(seed)>=100)
						{
						
							addSeed(tmp_sets[t], seed, Single());

							long score = scoreSeed(seed);

							if (score > max_score)
							{
								max_score = score;
								offset = endPositionV(seed) - position(begin_it, seq); // cannot take length of seed because it extendeds to BOTH ends
							}
						}
					}
					repeat_limit = beginPosition(ref)+occs[i]+q; // q = length of initial seed
				}
			//}
				begin_it += offset;
		}
	std::cout << "Done " << t << "\tFound " << length(tmp_sets[t]) << " seeds." << std::endl;
	}

	// combine seedSets
	for (unsigned t = 0; t < num_threads; ++t)
	{
		//std::cout << length(tmp_sets[t]) << std::endl;
		for (TIter it = begin(tmp_sets[t], Standard()); it != end(tmp_sets[t], Standard()); ++it)
			addSeed(seedSet, *it, Single());
	}
    return 0;
}

// ----------------------------------------------------------------------------
// Function fastFirstSeeding()
// ----------------------------------------------------------------------------
template<typename TSeed, typename TIndex, typename TSeq, typename TInfix>
int fastFirstSeeding(SeedSet<TSeed> & seedSet, TIndex & index, TSeq & ref,
		             TInfix & seq, unsigned q)
{
	//typedef Index<TSeq, IndexQGram<SimpleShape, OpenAddressing> > TIndex;

	typedef typename Iterator<TInfix>::Type TIterator;
	typedef String<typename SAValue<Index<TSeq, TIndex> >::Type> TOccurrences;

	unsigned distance = 0;

	//TIndex index(ref);
	//resize(indexShape(index),q);
	typename Fibre<TIndex, QGramShape>::Type shape = indexShape(index);

	TIterator it = begin(seq);
	while (position(it, seq) < length(seq)-q+1)
	{
		unsigned offset = 1;
		//std::cout << position(it,seq) << std::endl;
		//if (not(infix(seq, position(it, seq), position(it, seq)+q) == "NNNNNNNNNNNNNNNNNNNN"))
		//{
			hash(shape, it);
			TOccurrences occs = getOccurrences(index, shape);
			long max_score = -maxValue<long>();
			unsigned repeat_limit = 0;
			// if there are no hits the offset is 1 otherwise its the end position
			// of the best scoring seed.
			// This way long runs will speed up seeding and result in less memory
			// usage.
			//unsigned len = length(occs);
			for (unsigned i = 0; i < length(occs); ++i)
			{
				// avoid repeats:
				// a repeat can be identified if the subsequent seed is found within
				// the preceding seed (endPos(pre_seed) > beginPos(post_seed)
				// the next if clause therefore is avoiding repeats
				if ((beginPosition(ref)+occs[i] > repeat_limit) | (beginPosition(ref)+occs[i] == 0))
				{
					TSeed seed = TSeed(beginPosition(ref) + occs[i], beginPosition(seq) + position(it, seq), q);
					extendSeed(seed, ref, seq, EXTEND_BOTH, MatchExtend());

					if (!addSeed(seedSet, seed, distance, Merge()))
						addSeed(seedSet, seed, Single());

					long score = scoreSeed(seed);

					if (score > max_score)
					{
						max_score = score;
						offset = endPositionV(seed) - position(it, seq); // cannot take length of seed because it extendeds to BOTH ends
					}
				}
				repeat_limit = beginPosition(ref)+occs[i]+q; // q = length of initial seed
			}
		//}
		it += offset;
	}
    return 0;
}

// ----------------------------------------------------------------------------
// Function computeSearchFields()
// ----------------------------------------------------------------------------
template<typename TPair, typename TSeq>
int computeSearchFields(StringSet<SearchField> & fields, TPair & posV, TPair & posH,
						TSeq & seq, TSeq & ref, unsigned qV, unsigned qH)
{
	// posV and posH store the begin and end positions of seeds
	// when seeds are chained globally, sorting the positions will lead
	// to correct search field borders from left to right.
	sort(posV);
	sort(posH);
	unsigned bPH = 0; // initial begin Position H
	unsigned bPV = 0; // initial begin Position V
	// compute search fields
	for (unsigned j = 0; j < length(posV); ++j)
	{
		// if clause must avoid that a search field can have 0
		// characters if qV and qH are 0
	    if ((posH[j].i1 - bPH > qH and posV[j].i1 - bPV >= qV) or
	    	(posH[j].i1 - bPH >= qH and posV[j].i1 - bPV > qV))
	    {
	        SearchField field(bPH, bPV, posH[j].i1, posV[j].i1);
	        SEQAN_ASSERT_LEQ(bPH, posH[j].i1);
	        SEQAN_ASSERT_LEQ(bPV, posV[j].i1);
	        appendValue(fields, field);
	    }
	    bPH = posH[j].i2;
	    bPV = posV[j].i2;
	}
	if ((length(ref)-bPH > qH and length(seq)-bPV >= qV) or
		(length(ref)-bPH >= qH and length(seq)-bPV > qV))
	{
	    SearchField last(bPH, bPV, length(ref), length(seq));
	    appendValue(fields, last);
	}
	return 0;
}
// ----------------------------------------------------------------------------
// Function IterativeSeeding()
// ----------------------------------------------------------------------------
/*
 * This function enables an iterative loop when finding seeds. At first, search
 * fields are calculated between previously found seeds. The search field
 * boundaries are used to create sequence infixes that are each given to the
 * seeding() function to process.
 */
template<typename TSeed, typename TPairSet, typename TSeq>
int iterativeSeeding(SeedSet<TSeed> & seedSet, TPairSet & posV, TPairSet & posH,
		             TSeq & ref, TSeq & seq, String<unsigned> & lagan_parameter)
{
    typedef SeedSet<TSeed> TSeedSet;
	typedef Pair<unsigned, unsigned> TPair;
    typedef typename Infix<TSeq>::Type TInfix;
    typedef typename Value<TSeq>::Type TAlphabet;

    unsigned alphSize = ValueSize<TAlphabet>::VALUE;
    int closedAdressingLimit = (int)(log(4000000000)/log(alphSize)); // replace 4000000000 ?

    StringSet<SearchField> fields;

    for(unsigned i = 0; i < length(lagan_parameter); ++i)
    {
        int q = lagan_parameter[i];
        std::cout << "q = " << q << std::endl;

        clear(fields);
        computeSearchFields(fields, posV, posH, seq, ref, q, q);

        for (unsigned sf = 0; sf < length(fields); ++sf)
        {
//        	std::cout << fields[sf].beginV << " "<< fields[sf].endV << " "
//        			  << fields[sf].beginH << " " << fields[sf].endH << std::endl;
            TSeedSet tmp_seedSet;
            String<TSeed> tmp_seedChain;
            TInfix seq_fragment = infix(seq, fields[sf].beginV, fields[sf].endV);
            TInfix ref_fragment = infix(ref, fields[sf].beginH, fields[sf].endH);

            //generate seeds
            if(q < closedAdressingLimit)
            {
            	nativeSeeding(tmp_seedSet, ref_fragment, seq_fragment, lagan_parameter[i], Default());
            }
            else
            {
            	nativeSeeding(tmp_seedSet, ref_fragment, seq_fragment, lagan_parameter[i], OpenAddressing());
            }

            //chain temporary seeds globally
            if(length(tmp_seedSet) != 0)
                chainSeedsGlobally(tmp_seedChain, tmp_seedSet, SparseChaining());

//            std::cout << "after chaining \n ";
//            for (unsigned j = 0; j < length(tmp_seedChain); ++j)
//                std::cout << tmp_seedChain[j] << std::endl;

            //add remaining seeds to global seedSet
            for (unsigned j = 0; j < length(tmp_seedChain); ++j)
            {
                TSeed & seed = tmp_seedChain[j];
                TPair pairH(beginPositionH(seed), endPositionH(seed));
                TPair pairV(beginPositionV(seed), endPositionV(seed));
                appendValue(posH, pairH);
                appendValue(posV, pairV);
                addSeed(seedSet, seed, Single());
            }
        }
    //std::cout << "Found all in all " << length(seedSet) << " seeds" << std::endl;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function transformIntoJournal()
// ----------------------------------------------------------------------------
template<typename TSeed, typename TSeq>
int transformIntoJournal(String<DeltaEvent> & records, String<TSeed> & seedChain,
		                 TSeq & host, TSeq & sequence, unsigned s, unsigned n)
{
    typedef Pair<unsigned, unsigned > TPair;
    typedef StringSet<TPair> TPairSet;

    StringSet<SearchField> fields;
    TPairSet posV;
    TPairSet posH;

    for (unsigned i = 0; i < length(seedChain); ++i)
    {
    	TSeed & seed = seedChain[i];
		TPair pairH(beginPositionH(seed), endPositionH(seed));
		TPair pairV(beginPositionV(seed), endPositionV(seed));
		appendValue(posH, pairH);
		appendValue(posV, pairV);
    }

    computeSearchFields(fields, posV, posH, sequence, host, 0, 0);
    //computeDeletions(fields, posV, posH, sequence, host);

    // transform search fields into delta events in a journaled String
//    unsigned vp = 0; // stores virtual position when altering the journaled string
    for (unsigned sf = 0; sf < length(fields); ++sf)
	{
    	unsigned eH = fields[sf].endH;
    	unsigned bH = fields[sf].beginH;
    	unsigned eV = fields[sf].endV;
    	unsigned bV = fields[sf].beginV;
    	if(eV-bV == 1 && eH-bH == 1)
    	{
//    		std::cout << "SNP at " << bH << " " << host[bH] << "->" << sequence[bV] << std::endl;
		String<Dna5> ins;
		appendValue(ins, sequence[bV];
    		DeltaEvent rec = DeltaEvent(bH, s, n, 0, ins, 1);
    		appendValue(records, rec);
    	}
    	else
    	{
    		if ((eH-bH >= 1) & (eV-bV >= 1))
    		{
//    			std::cout << "Structural Variation: del at " << bH << "-" << eH << " " <<  " ins at " << bV << "-" << eV << std::endl;

    			Infix<Dna5String>::Type infix(sequence, bV, eV);
			String<Dna5> ins;
			getString(ins, infix);
    			DeltaEvent rec = DeltaEvent(bH, s, n, 3, ins, (eH-bH));
    			appendValue(records, rec);
    		}
    		else if (eH-bH >= 1)
			{
//				std::cout << "Deletion at " << bH << "-" << eH << " " << infix(host, bH, eH) << std::endl;

				String<Dna5> ins;
				DeltaEvent rec = DeltaEvent(bH, s, n, 1, ins, (eH-bH));
				appendValue(records, rec);
			}
    		else if (eV-bV >= 1)
    		{
//    			std::cout << "Insertion at " << bH << "-" << eH << " " << infix(sequence, bV, eV) << std::endl;

    			Infix<Dna5String>::Type infix(sequence, bV, eV);
    			String<Dna5> ins;
    			getString(ins, infix);
    			DeltaEvent rec = DeltaEvent(bH, s, n, 2, ins, 0);
    			appendValue(records, rec);
    		}
    	}
	}


    return 0;
}

// ----------------------------------------------------------------------------
// Function seeding()
// ----------------------------------------------------------------------------
template<typename TSeed, typename TIndex, typename TSeq>
int seeding(SeedSet<TSeed> & seedSet, TIndex & index, TSeq & ref, TSeq & seq,
					 String<unsigned> & lagan_parameter)
{
	typedef SeedSet<TSeed> TSeedSet;
	typedef String<TSeed> TSeedChain;
	typedef Pair<unsigned, unsigned> TPair;
	typedef StringSet<TPair> TPairSet;

	TSeedSet tmp_seedSet;
	TSeedChain tmp_seedChain;
	TPairSet posV;
	TPairSet posH;

	parallelFastFirstSeeding(tmp_seedSet, index, ref, seq, lagan_parameter[0], 14);

//	std::cout << "\n before chaining \n ";
//	typedef typename Iterator<TSeedSet>::Type TIter;
//	for (TIter it = begin(tmp_seedSet, Standard()); it != end(tmp_seedSet, Standard()); ++it)
//		std::cout << *it << std::endl;

	std::cout << "Done Seeding. Start Chaining...\n";
	if(length(tmp_seedSet) != 0)
	    chainSeedsGlobally(tmp_seedChain, tmp_seedSet, SparseChaining());

//	std::cout << "\n after chaining \n ";
	for (unsigned j = 0; j < length(tmp_seedChain); ++j)
	{
		TSeed & seed = tmp_seedChain[j];
//		std::cout << seed << std::endl;
		TPair pairH(beginPositionH(seed), endPositionH(seed));
		TPair pairV(beginPositionV(seed), endPositionV(seed));
		appendValue(posH, pairH);
		appendValue(posV, pairV);
		addSeed(seedSet, seed, Single());
	}

	iterativeSeeding(seedSet, posV, posH, ref, seq, lagan_parameter);
	return 0;
}
// ----------------------------------------------------------------------------
// Function computelaganAlignment()                                     [Align]
// ----------------------------------------------------------------------------
//main
template<typename TSequence, typename TIndex /*, typename TScoreValue, typename TScoreSpecAnchor*/>
int getDeltaEvents(String<DeltaEvent> & records, TSequence & ref, TSequence & seq,
		           unsigned s, unsigned n, TIndex & index, String<unsigned> & lagan_parameter
                          /*, Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor*/)
{
    typedef Seed<Simple> TSeed;
    typedef SeedSet<TSeed> TSeedSet;
    TSeedSet seedSet;
    String<TSeed> seedChain;

    // -----------------------------------------------------------------------
    // Scan query for seeds
    // -----------------------------------------------------------------------
    std::cout << "### Seeding...";
    timestamp();
    seeding(seedSet, index, ref, seq, lagan_parameter);
    SEQAN_ASSERT(length(seedSet)!=0);
    
    std::cout << "### Chaining...";
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
//    for (unsigned i = 0; i < length(seedChain); ++i)
//        std::cout << seedChain[i] << std::endl;
    std::cout << "### " << length(seedChain) << " seeds were found.";
    timestamp();
    // -----------------------------------------------------------------------
    // Variant Calling
    // -----------------------------------------------------------------------
    std::cout << "### Variant Calling...\n";
    transformIntoJournal(records, seedChain, ref, seq, s, n);

    //timestamp();
    return 0;
}

#endif //SEQAN_DEMOS_MINILAGAN_LAGANALIGNMENT_IMPL2_eeds.

