// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Basic header file.
// ==========================================================================

#ifndef EXTRAS_APPS_GDF_TOOLS_GDF_TOOLS_H_
#define EXTRAS_APPS_GDF_TOOLS_GDF_TOOLS_H_

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>
#include <seqan/journaled_set.h>
#include <seqan/journaled_string_tree.h>
#include <seqan/parallel.h>
#include <seqan/seq_io.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct ReadVariantInformation_;
typedef Tag<ReadVariantInformation_> WithoutLoadingJournalData;

struct ReadVariantInformationAndJournalData_;
typedef Tag<ReadVariantInformationAndJournalData_> WithLoadingJournalData;

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// The general app options.
struct AppOptions
{
    unsigned verbosity;  // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    // TODO(rmaerker): referenceFile, jseqIndexFile, bjseq file, djseq file
    unsigned numThreads;

    CharString referenceFile;        // We need the reference file (should be an external file!)! It's mandatory!
    CharString jseqFile;
    CharString cglIndexFile;     // This is the input index file which we might want to generate.

    AppOptions() : verbosity(1), numThreads(1)
    {}
};

// --------------------------------------------------------------------------
// Class ConverterOptions
// --------------------------------------------------------------------------

struct ConverterOptions
{
    unsigned verbosity;  // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    unsigned numThreads;

    CharString inputFile;
    CharString outputFile;

    unsigned refAlphabet;
    unsigned varAlphabet;

    // VCF Options:
    CharString vcfReferenceFile;
    bool readGenotype;
    bool includeReference;
    String<unsigned> haplotypes;


    unsigned numIndividuals;

    bool suppressSVs;
    bool selftest;
    CharString compareFile;

    ConverterOptions() :
        verbosity(1),
        numThreads(1),
        readGenotype(false),
        includeReference(false),
        numIndividuals(1),
        suppressSVs(false),
        selftest(false)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _printJournalStat()
// ----------------------------------------------------------------------------

template <typename TJournalString>
inline void _printJournalStat(TJournalString const & journal)
{
    typedef typename JournalType<TJournalString const>::Type TJournalEntries;
    typedef typename Value<TJournalEntries>::Type TJournalNode;

    std::cerr << "Number of journal entries: " << length(journal._journalEntries._journalNodes)  << std::endl;
    std::cerr << "Lenght of Insertionbuffer: " << length(journal._insertionBuffer) << std::endl;
    std::cerr << "sizeof(TJournal()) " << sizeof(TJournalNode) << std::endl;
    unsigned memUsage = sizeof(TJournalNode) * length(journal._journalEntries._journalNodes) + length(journal._insertionBuffer);
    std::cerr << "Memory usage: " <<  memUsage + sizeof(TJournalString) << " Bytes." << std::endl;
    std::cerr << "Memory ref: " <<  length(host(journal)) << " Bytes." << std::endl;
}

// ----------------------------------------------------------------------------
// Function _loadSequenceFasta()
// ----------------------------------------------------------------------------

template <typename TId, typename TSequence>
inline int _loadSequenceFasta(TId & idString,
                              TSequence & sequence,
                              CharString file)
{
    std::ifstream fileStream(toCString(file), std::ios_base::in);
    if (!fileStream.good())
    {
        std::cerr << "Can't read file <" << file << "!"<< std::endl;
        return -1;
    }

    RecordReader<std::ifstream, SinglePass<> > recReader(fileStream);
    int res = readRecord(idString, sequence, recReader, Fasta());
    if (res != 0)
        return -1;
    return 0;
}


}

#endif  // EXTRAS_APPS_GDF_TOOLS_GDF_TOOLS_H_
