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
#include <seqan/sequence.h>

#include <time.h>
#include <omp.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
/*
struct CompareByPosAndTypeLessThan_
{
	template <typename TRecord>
	inline bool operator()(TRecord const &lhs, TRecord const &rhs)
	{
		if (lhs.pos == rhs.pos)
			return lhs.type < rhs.type;
		return lhs.pos < rhs.pos;
	}
};
*/
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
 * seqH and afterwards occurences of common subsequences of size q are added
 * to a seedSet.
 */
template<typename TSeed, typename TSeq, typename TInfix, typename TIndexTag>
int nativeSeeding(SeedSet<TSeed> & seedSet, TSeq & seqH,
            TInfix & seqV, unsigned q, TIndexTag /*tag*/)
{
    typedef Index<TSeq, IndexQGram<SimpleShape, TIndexTag> > TIndex;
    typedef typename Iterator<TInfix>::Type TIterator;
    typedef String<typename SAValue<Index<TSeq, TIndex> >::Type> TOccurrences;

    unsigned distance = 0;

    TIndex index(seqH);
    resize(indexShape(index),q);

    hashInit(indexShape(index), begin(seqV));
    for(TIterator it = begin(seqV); it != (end(seqV)-q+1); goNext(it))
    {
    	unsigned repeat_limit = 0;
        hashNext(indexShape(index), it);
        TOccurrences occs = getOccurrences(index, indexShape(index));
        for (unsigned i = 0; i < length(occs); ++i)
        {
        	if (beginPosition(seqH) + occs[i] > repeat_limit)
        	{
				TSeed seed = TSeed(beginPosition(seqH) + occs[i], beginPosition(seqV) + position(it, seqV), q);
				if (!addSeed(seedSet, seed, distance, Merge()))
					addSeed(seedSet, seed, Single());
        	}
        	repeat_limit = beginPosition(seqH)+occs[i]+q; // q = length of initial seed
        }
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function parallelSeeding()
// ----------------------------------------------------------------------------
/*
 * Given two sequence (fragments as infixes) the functions builds an index over
 * seqH and afterwards occurences of common subsequences of size q are added
 * to a seedSet.
 */
template<typename TSeed, typename TSeq, typename TInfix, typename TIndexTag>
int parallelSeeding(SeedSet<TSeed> & seedSet, TSeq & seqH, TInfix & seqV,
                    unsigned q, unsigned num_threads, TIndexTag /*tag*/)
{
    typedef Index<TSeq, IndexQGram<SimpleShape, TIndexTag> > TIndex;
    typedef typename Iterator<TInfix>::Type TIterator;
    typedef typename Iterator<SeedSet<TSeed> >::Type TIter;
    typedef String<typename SAValue<Index<TSeq, TIndex> >::Type> TOccurrences;

    unsigned distance = 0;

    if (int(length(seqV)/num_threads) < q)
    	num_threads = 1;

    //set up for parallelism
    String<SeedSet<TSeed> > tmp_sets;
    resize(tmp_sets, num_threads);
    omp_set_num_threads(num_threads);
    TIndex index(seqH);
    resize(indexShape(index), q);

    // each thread fills a temporary seedSet
    #pragma omp parallel for
    for (unsigned t = 0; t < num_threads; ++t)
    {
        typename Fibre<TIndex, QGramShape>::Type shape = indexShape(index);

        TIterator begin_it = begin(seqV)+ t*(int)(length(seqV)/num_threads);
        TIterator end_it;
        if (t != num_threads-1)
            end_it = begin_it+ (int)(length(seqV)/num_threads);
        else
            end_it = end(seqV) - q+1;

        //hashInit(shape, begin_it);
        for(TIterator it = begin_it; it != end_it; it += q)
        {
            hash(shape, it);
            TOccurrences occs = getOccurrences(index, shape);
            for (unsigned i = 0; i < length(occs); ++i)
            {
                TSeed seed = TSeed(beginPosition(seqH) + occs[i], beginPosition(seqV) + position(it, seqV), q);
                if (!addSeed(tmp_sets[t], seed, distance, Merge()))
                    addSeed(tmp_sets[t], seed, Single());
            }
        }
    }

    // combine seedSets
    for (unsigned t = 0; t < num_threads; ++t)
    {
        std::cout << length(tmp_sets[t]) << std::endl;
        for (TIter it = begin(tmp_sets[t], Standard()); it != end(tmp_sets[t], Standard()); ++it)
        {

            TSeed seed = *it;
            // if begin position of seed lies whithin an overlay area it must be merged
            if (beginPositionV(seed) < (t * (int)(length(seqV)/num_threads) + q + distance))
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
int scoreSeed(TSeed & seed)
{
	unsigned n = seedSize(seed);
	SEQAN_ASSERT_EQ(upperDiagonal(seed), lowerDiagonal(seed));

	int score = n - n * abs(lowerDiagonal(seed));
	return (score);
}

// ----------------------------------------------------------------------------
// Function fastFirstSeeding()
// ----------------------------------------------------------------------------
template<typename TSeed, typename TIndex, typename TSeq, typename TInfix>
int fastFirstSeeding(SeedSet<TSeed> & seedSet, TIndex & index, TSeq & seqH,
		             TInfix & seqV, unsigned q)
{
	//typedef Index<TSeq, IndexQGram<SimpleShape, OpenAddressing> > TIndex;

	typedef typename Iterator<TInfix>::Type TIterator;
	typedef String<typename SAValue<Index<TSeq, TIndex> >::Type> TOccurrences;

	unsigned distance = 0;

	//TIndex index(seqH);
	//resize(indexShape(index),q);
	typename Fibre<TIndex, QGramShape>::Type shape = indexShape(index);

	TIterator it = begin(seqV);
	while (position(it, seqV) < length(seqV)-q+1)
	{
		hash(shape, it);
		TOccurrences occs = getOccurrences(index, shape);
		int max_score = -maxValue<int>();
		unsigned offset = 1;
		unsigned repeat_limit = 0;

		// if there are no hits the offset is 1 otherwise its the end position
		// of the best scoring seed.
		// This way long runs will speed up seeding and result in less memory
		// usage.
		for (unsigned i = 0; i < length(occs); ++i)
		{
			// avoid repeats:
			// a repeat can be identified if the subsequent seed is found within
			// the preceding seed (endPos(pre_seed) > beginPos(post_seed)
			// the next if clause therefore is avoiding repeats
			if (beginPosition(seqH)+occs[i] > repeat_limit)
			{
				TSeed seed = TSeed(beginPosition(seqH) + occs[i], beginPosition(seqV) + position(it, seqV), q);
				extendSeed(seed, seqH, seqV, EXTEND_BOTH, MatchExtend());

				if (!addSeed(seedSet, seed, distance, Merge()))
					addSeed(seedSet, seed, Single());

				int score = scoreSeed(seed);

				if (score > max_score)
				{
					max_score = score;
					offset = endPositionV(seed) - position(it, seqV); // cannot take length of seed because it extendeds to BOTH ends
				}
			}
			repeat_limit = beginPosition(seqH)+occs[i]+q; // q = length of initial seed
		}
		it += offset;
	}
    return 0;
}

// ----------------------------------------------------------------------------
// Function computeSearchFields()
// ----------------------------------------------------------------------------
template<typename TPair, typename TSeq>
int computeSearchFields(StringSet<SearchField> & fields, TPair & posV, TPair & posH,
						TSeq & seqV, TSeq & seqH, unsigned qV, unsigned qH)
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
	if ((length(seqH)-bPH > qH and length(seqV)-bPV >= qV) or
		(length(seqH)-bPH >= qH and length(seqV)-bPV > qV))
	{
	    SearchField last(bPH, bPV, length(seqH), length(seqV));
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
		             TSeq & seqH, TSeq & seqV, String<unsigned> & lagan_parameter)
{
    typedef SeedSet<TSeed> TSeedSet;
	typedef Pair<unsigned, unsigned> TPair;
    typedef typename Infix<TSeq>::Type TInfix;
    typedef typename Value<TSeq>::Type TAlphabet;

    unsigned alphSize = ValueSize<TAlphabet>::VALUE;
    unsigned closedAdressingLimit = (int)(log(4000000000)/log(alphSize)); // replace 4000000000 ?

    StringSet<SearchField> fields;

    for(unsigned i = 0; i < length(lagan_parameter); ++i)
    {
        int q = lagan_parameter[i];
        std::cout << "q = " << q << std::endl;

        clear(fields);
        computeSearchFields(fields, posV, posH, seqV, seqH, q, q);

        for (unsigned sf = 0; sf < length(fields); ++sf)
		{
        	std::cout << fields[sf].beginV << " "<< fields[sf].endV << " "
        			<< fields[sf].beginH << " " << fields[sf].endH << std::endl;
		}
        // search seeds within search fields
        for (unsigned sf = 0; sf < length(fields); ++sf)
        {

            TSeedSet tmp_seedSet;
            String<TSeed> tmp_seedChain;
            TInfix seqV_fragment = infix(seqV, fields[sf].beginV, fields[sf].endV);
            TInfix seqH_fragment = infix(seqH, fields[sf].beginH, fields[sf].endH);

            //generate seeds
            if(q < closedAdressingLimit)
            {
            	nativeSeeding(tmp_seedSet, seqH_fragment, seqV_fragment, lagan_parameter[i], Default());
            }
            else
            {
            	nativeSeeding(tmp_seedSet, seqH_fragment, seqV_fragment, lagan_parameter[i], OpenAddressing());
            }

/*            //generate seeds
            if(q < closedAdressingLimit)
            	parallelSeeding(tmp_seedSet, seqH_fragment, seqV_fragment, lagan_parameter[i], 1, Default());
            else
                parallelSeeding(tmp_seedSet, seqH_fragment, seqV_fragment, lagan_parameter[i], 1, OpenAddressing());
*/
            //chain temporary seeds globally
            if(length(tmp_seedSet) != 0)
                chainSeedsGlobally(tmp_seedChain, tmp_seedSet, SparseChaining());

            for (unsigned j = 0; j < length(tmp_seedChain); ++j)
                std::cout << tmp_seedChain[j] << std::endl;

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
    std::cout << "Found all in all " << length(seedSet) << " seeds" << std::endl;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function transformIntoJournal()
// ----------------------------------------------------------------------------
template<typename TSeed, typename TSeq>
int transformIntoJournal(String<TSeed> & seedChain, TSeq & host, TSeq & sequence)
{
	typedef typename Value<TSeq>::Type TSeqValue;
    typedef String<TSeqValue, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournaledString;
    typedef Pair<unsigned, unsigned> TPair;
    typedef StringSet<TPair> TPairSet;

    TJournaledString journal;
    setHost(journal, host);

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
    unsigned vp = 0; // stores virtual position when altering the journaled string
    for (unsigned sf = 0; sf < length(fields); ++sf)
	{
    	unsigned eH = fields[sf].endH;
    	unsigned bH = fields[sf].beginH;
    	unsigned eV = fields[sf].endV;
    	unsigned bV = fields[sf].beginV;
    	if(eV-bV == 1 && eH-bH == 1)
    	{
    		std::cout << "SNP at " << bH << " " << host[bH] << "->" << sequence[bV] << std::endl;
    		//if (bH != bV)
    			//std::cout << "oh oh " << std::endl;
    		//SEQAN_ASSERT_EQ(bH, bV); does not work if not on main diagonal
    		//SEQAN_ASSERT_EQ(eH, eV);
    		assignValue(journal, vp+bH, sequence[bV]);
//    		erase(journal, vp+bH, vp+bH+1);
//    		insert(journal, vp+bH, sequence[bV]);
    		// no update of the virtual position needed
    	}
    	else
    	{
    		if (eH-bH > 1)
			{
				std::cout << "Deletion at " << bH << " " << infix(host, bH, eH) << std::endl;
				erase(journal, vp+bH, vp+eH);
				vp -= (eH-bH);
			}
    		if (eV-bV > 1)
    		{
    			std::cout << "Insertion at " << bV << " " << infix(sequence, bV, eV) << std::endl;
    			insert(journal, vp+bV, vp+eV);
    			vp += (eV-bV);
    		}
    	}
	}

    //String<int> records;
    //sort(records, CompareByPosAndTypeLessThan_());
    //std::cout << "Host: " << host << std::endl;
    //std::cout << "Journal: " << journal << std::endl;
    SEQAN_ASSERT(journal == sequence);
    //std::cout << "Nodes: " << journal._journalEntries << std::endl;
    //std::cout << std::endl;

    return 0;
}

// ----------------------------------------------------------------------------
// Function seeding()
// ----------------------------------------------------------------------------
template<typename TSeed, typename TIndex, typename TSeq>
int seeding(SeedSet<TSeed> & seedSet, TIndex & index, TSeq & seqH, TSeq & seqV,
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

	fastFirstSeeding(tmp_seedSet, index, seqH, seqV, lagan_parameter[0]);

	if(length(tmp_seedSet) != 0)
	    chainSeedsGlobally(tmp_seedChain, tmp_seedSet, SparseChaining());

	for (unsigned j = 0; j < length(tmp_seedChain); ++j)
	{
		TSeed & seed = tmp_seedChain[j];
		TPair pairH(beginPositionH(seed), endPositionH(seed));
		TPair pairV(beginPositionV(seed), endPositionV(seed));
		appendValue(posH, pairH);
		appendValue(posV, pairV);
		addSeed(seedSet, seed, Single());
	}

	iterativeSeeding(seedSet, posV, posH, seqH, seqV, lagan_parameter);
	return 0;
}
// ----------------------------------------------------------------------------
// Function computelaganAlignment()                                     [Align]
// ----------------------------------------------------------------------------
//main
template<typename TSequence, typename TIndex /*, typename TScoreValue, typename TScoreSpecAnchor*/>
int computelaganAlignment(TSequence & ref, TSequence & seq, TIndex index,
						  String<unsigned> & lagan_parameter
                          /*, Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor*/)
{
    typedef Seed<Simple> TSeed;
    typedef SeedSet<TSeed> TSeedSet;
    TSeedSet seedSet;

    // -----------------------------------------------------------------------
    // Step 1: Scan query for seeds
    // -----------------------------------------------------------------------
    std::cout << "# # Step 1: Iterative Seeding and Chaining...\n";
    timestamp();
    seeding(seedSet, index, ref, seq, lagan_parameter);

    // -----------------------------------------------------------------------
    // Step 2: Compute global chain
    // -----------------------------------------------------------------------
    std::cout << "# # Step 2: Chaining...\n";
    timestamp();
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    //for (unsigned i = 0; i < length(seedChain); ++i)
      //  std::cout << seedChain[i] << std::endl;
    std::cout << "# # # " << length(seedChain) << " seeds remained\n";

    // -----------------------------------------------------------------------
    // Step 3: Transform to journaled string
    // -----------------------------------------------------------------------
    std::cout << "# # Step 3: Transforming...\n";
    timestamp();
    transformIntoJournal(seedChain, ref, seq);

    //timestamp();
    return 0;
}

#endif //SEQAN_DEMOS_MINILAGAN_LAGANALIGNMENT_IMPL2_H
