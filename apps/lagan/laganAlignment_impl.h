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

#ifndef SEQAN_DEMOS_MINILAGAN_LAGANALIGNMENT_IMPL_H
#define SEQAN_DEMOS_MINILAGAN_LAGANALIGNMENT_IMPL_H

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
int seeding(SeedSet<TSeed> & seedSet, TSeq & seqH,
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
        hashNext(indexShape(index), it);
        TOccurrences occs = getOccurrences(index, indexShape(index));
        for (unsigned i = 0; i < length(occs); ++i)
        {
            TSeed seed = TSeed(beginPosition(seqH) + occs[i], beginPosition(seqV) + position(it, seqV), q);
            if (!addSeed(seedSet, seed, distance, Merge()))
                addSeed(seedSet, seed, Single());
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
	    if (posH[j].i1 - bPH >= qH and posV[j].i1 - bPV >= qV)
	    {
	        SearchField field(bPH, bPV, posH[j].i1, posV[j].i1);
	        SEQAN_ASSERT_LEQ(bPH, posH[j].i1);
	        SEQAN_ASSERT_LEQ(bPV, posV[j].i1);
	        appendValue(fields, field);
	    }
	    bPH = posH[j].i2;
	    bPV = posV[j].i2;
	}
	if (length(seqH)-1-bPH >= qH and length(seqV)-1-bPV >= qV)
	{
	    SearchField last(bPH, bPV, length(seqH)-1, length(seqV)-1);
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
template<typename TSeed, typename TSeq>
int iterativeSeeding(SeedSet<TSeed> & seedSet, TSeq & seqH, TSeq & seqV,
					 String<unsigned> & lagan_parameter)
{
    typedef SeedSet<TSeed> TSeedSet;
    typedef Pair<unsigned, unsigned> TPair;
    typedef StringSet<TPair> TPairSet;
    typedef typename Infix<TSeq>::Type TInfix;
    typedef typename Value<TSeq>::Type TAlphabet;

    unsigned alphSize = ValueSize<TAlphabet>::VALUE;
    unsigned closedAdressingLimit = (int)(log(4000000000)/log(alphSize)); // replace 4000000000 ?

    TPairSet posV;
    TPairSet posH;
    StringSet<SearchField> fields;

    for(unsigned i = 0; i < length(lagan_parameter); ++i)
    {
        int q = lagan_parameter[i];
        std::cout << "q = " << q << std::endl;

        clear(fields);
        sort(posV);
		sort(posH);
		unsigned bPH = 0; // initial begin Position H
		unsigned bPV = 0; // initial begin Position V
		// compute search fields
		for (unsigned j = 0; j < length(posV); ++j)
		{
			if (posH[j].i1 - bPH >= q and posV[j].i1 - bPV >= q)
			{
				SearchField field(bPH, bPV, posH[j].i1, posV[j].i1);
				appendValue(fields, field);
			}
			bPH = posH[j].i2;
			bPV = posV[j].i2;
		}
		if (length(seqH)-1-bPH >= q and length(seqV)-1-bPV >= q)
		{
			SearchField last(bPH, bPV, length(seqH)-1, length(seqV)-1);
			appendValue(fields, last);
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
                seeding(tmp_seedSet, seqH_fragment, seqV_fragment, lagan_parameter[i], Default());
            else
            	seeding(tmp_seedSet, seqH_fragment, seqV_fragment, lagan_parameter[i], OpenAddressing());

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
    typedef String<TSeq, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournaledString;
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

    computeSearchFields(fields, posV, posH, sequence, host, 1, 0);
    computeDeletions(fields, posV, posH, sequence, host);

    // transform search fields into delta events in a journaled String
    for (unsigned sf = 0; sf < length(fields); ++sf)
	{
    	unsigned eH = fields[sf].endH;
    	unsigned bH = fields[sf].beginH;
    	unsigned eV = fields[sf].endV;
    	unsigned bV = fields[sf].beginV;
    	if(eV-bV == 1 && eH-bH == 1)
    	{
    		std::cout << "SNP at " << bH << " " << host[bH] << "->" << sequence[bV] << std::endl;
    		SEQAN_ASSERT_EQ(bH, bV);
    		SEQAN_ASSERT_EQ(eH, eV);
    	}
    	else
    	{
    		if (eV-bV > 1)
    			std::cout << "Insertion at " << bV << " " << infix(sequence, bV, eV) << std::endl;
    		if (eH-bH > 1)
    		{
    			std::cout << "Deletion at " << bH << " " << infix(host, bH, eH) << std::endl;
    		}
    	}
	}

    return 0;
}

// ----------------------------------------------------------------------------
// Function computeSearchFields()
// ----------------------------------------------------------------------------
template<typename TPair, typename TSeq>
int computeDeletions(StringSet<SearchField> & fields, TPair & posV, TPair & posH,
						TSeq & seqV, TSeq & seqH)
{
	// posV and posH are already sorted from before!
	unsigned bPH = 0; // initial begin Position H
	unsigned bPV = 0; // initial begin Position V
	// compute search fields
	for (unsigned j = 0; j < length(posV); ++j)
	{
	    if (posH[j].i1 - bPH >= 1 and posV[j].i1 - bPV == 0)
	    {
	        SearchField field(bPH, bPV, posH[j].i1, posV[j].i1);
	        appendValue(fields, field);
	    }
	    bPH = posH[j].i2;
	    bPV = posV[j].i2;
	}
	if (length(seqH)-1-bPH >= 1 and length(seqV)-1-bPV == 0)
	{
	    SearchField last(bPH, bPV, length(seqH)-1, length(seqV)-1);
	    appendValue(fields, last);
	}
	return 0;
}

// ----------------------------------------------------------------------------
// Function computelaganAlignment()                                     [Align]
// ----------------------------------------------------------------------------
//main
template<typename TSequence, typename TAlignSpec, typename TScoreValue, 
         typename TScoreSpecAnchor, typename TScoreSpecGap>
int computelaganAlignment(Align<TSequence, TAlignSpec> & alignment,
						  String<unsigned> & lagan_parameter,
                          Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                          Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                          unsigned bandExtension)
{
    typedef Seed<Simple> TSeed;

    typedef SeedSet<TSeed> TSeedSet;

    StringSet<TSequence> seqs = stringSet(alignment);
    TSequence seqH = seqs[0];
    TSequence seqV = seqs[1];
    TSeedSet seedSet;

    // -----------------------------------------------------------------------
    // Step 1: Build index and scan query for seeds
    // -----------------------------------------------------------------------
    std::cout << "# # Step 1: Iterative Seeding and Chaining...\n";
    timestamp();
    iterativeSeeding(seedSet, seqH, seqV, lagan_parameter);

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
    transformIntoJournal(seedChain, seqH, seqV);

    // -----------------------------------------------------------------------
    // Step 3: Compute Alignment
    // -----------------------------------------------------------------------
    std::cout << "# # Step 3: Aligning...\n";
    timestamp();
    AlignConfig<false, false, false, false> alignConfig;
    int result = bandedChainAlignment(alignment, seedChain, scoreSchemeAnchor,
                                      scoreSchemeGap, alignConfig, bandExtension);

    //timestamp();
    return result;
}

#endif //SEQAN_DEMOS_MINILAGAN_LAGANALIGNMENT_IMPL_H
