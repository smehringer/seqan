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
// Selftest functionality for the gdf converter.
// ==========================================================================

#ifndef EXTRAS_APPS_GDF_TOOLS_GDF_CONVERTER_SELFTEST_H_
#define EXTRAS_APPS_GDF_TOOLS_GDF_CONVERTER_SELFTEST_H_

#include "gdf_tools.h"
#include "gdf_converter_vcf.h"


namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

//template <typename TReference>
//struct VcfRecordTranslator<StringSet<TReference> >
//{
//    typedef StringSet<TReference> TSet;
//
//    TSet * setPtr;
//
//    VcfRecordTranslator(TSet & _set) : setPtr(&_set)
//    {}
//
//    tempalte
//    operator(
//};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TRefAlphabet, typename TVarAlphabet>
inline bool _runSelftest(ConverterOptions const & options, TRefAlphabet /*ref*/, TVarAlphabet /*var*/)
{
    typedef DeltaMap<unsigned, TVarAlphabet> TDeltaMap;
    typedef String<TRefAlphabet, Alloc<> > TReference;

    typedef JournaledStringTree<TDeltaMap> TJst;
//    typedef StringSet<TReference> TSet;

    TJst jst;
    GdfHeader header;

    SEQAN_TRY
    {
        open(jst, header, options.outputFile);
    }
    SEQAN_CATCH(GdfIOException e)
    {
        std::cerr << e.what() << std::endl;
        return false;
    }

    // Create full journaled string tree.
    create(jst, 1);

    // Read the reference file.
    std::ifstream refStream;
    refStream.open(toCString(options.vcfReferenceFile), std::ios_base::in);
    if (!refStream.good())
    {
        std::cerr << "[LOG] Error while opening " << options.vcfReferenceFile << "!\n[LOG] Abort ..." << std::endl;
        return false;
    }
    RecordReader<std::ifstream, SinglePass<> > refReader(refStream);

    CharString refId;
    TReference ref;
    while (!atEnd(refReader))
    {
        if (readRecord(refId, ref, refReader, Fasta()))
        {
            std::cerr << "[LOG] Error while reading  " << options.vcfReferenceFile << "!\n[LOG] Abort ..." << std::endl;
            return false;
        }
    }

    // B) Open vcf and generate fasta based on converter.
    std::ifstream compStream;
    compStream.open(toCString(options.compareFile), std::ios_base::in);
    if (!compStream.good())
    {
        std::cerr << "[LOG] Error while opening " << options.compareFile << "!\n[LOG] Abort ..." << std::endl;

        return false;
    }
    // Now we parse the vcf file and check if the delta is present at the given position.
    RecordReader<std::ifstream, SinglePass<> > reader(compStream);

    bool result = true;
    unsigned counter = 0;
    while (!atEnd(reader))
    {
        CharString seqId;
        TReference seq;
        if (readRecord(seqId, seq, reader, Fasta()))
        {
            std::cerr << "[LOG] Error while reading  " << options.compareFile << "!\n[LOG] Abort ..." << std::endl;
            return false;
        }

        if (seqId != header.nameStore[counter])
        {
            std::cerr << "[LOG] Sequence Ids don't match. Gdf-ID at [" << counter << "]: " << header.nameStore[counter] << " and Fasta-ID: " << seqId << "!" << std::endl;
            result = false;
        }
        if (seq != stringSet(jst)[counter])
        {
            std::cerr << "[LOG] Sequences don't match. Gdf-seq at [" << counter << "]: " << stringSet(jst)[counter] << " and Fasta-seq: " << seq << "!" << std::endl;
//            for (unsigned i = 0; i < _min(length(stringSet(jst)[counter]), length(seq)); ++i)
//            {
//                if (stringSet(jst)[counter][i] != seq[i])
//                {
//                    std::cerr << "[LOG] Seq ID = " << seqId << std::endl;
//                    std::cerr << "[LOG] Mismatch at " << i << ": " << stringSet(jst)[counter][i] << " != " << seq[i] << std::endl;
//                    std::cerr << "[LOG] hostPos = " << virtualToHostPosition((stringSet(jst)[counter])._journalEntries, i) << std::endl;
//                    std::cerr << "[LOG] entry = " << *findInJournalEntries((stringSet(jst)[counter])._journalEntries, i) << std::endl;
//                    std::cerr << "//______COMPARE_____________________________________________" << std::endl;
//                    std::cerr << "[LOG] hostPos = " << virtualToHostPosition((stringSet(jst)[counter])._journalEntries, i - 1) << std::endl;
//                    std::cerr << "[LOG] entry = " << *findInJournalEntries((stringSet(jst)[counter])._journalEntries, i - 1) << std::endl;
//                    std::cerr << "//______COMPARE_____________________________________________" << std::endl;
//                    std::cerr << "[LOG]                        Pos: " << "..|....|...*|....|...." << std::endl;
//                    std::cerr << "[LOG] infix(jrn, 185923, 185945): " << infix(stringSet(jst)[counter], 185923, 185945) << std::endl;
//                    std::cerr << "[LOG] infix(seq, 185923, 185945): " << infix(seq, 185923, 185945) << std::endl;
//                    std::cerr << "[LOG] infix(ref, 185734, 185756): " << infix(ref, 185755, 185766) << infix(ref, 186597, 186608) << std::endl;
//
//                    break;
//                }
//            }
            result = false;
        }
        ++counter;
    }
    return result;
}


template <typename TRefAlphabet>
inline bool _runSelftest(ConverterOptions const & options, TRefAlphabet /*alph*/)
{
    switch (options.varAlphabet)
    {
        case 0: return _runSelftest(options, TRefAlphabet(), Dna());
        case 1: return _runSelftest(options, TRefAlphabet(), Dna5());
        default: return false;
    }
}

inline bool _runSelftest(ConverterOptions const & options)
{
    switch (options.refAlphabet)
    {
        case 0: return _runSelftest(options, Dna());
        case 1: return _runSelftest(options, Dna5());
        default: return false;
    }
}

}
#endif  // EXTRAS_APPS_GDF_TOOLS_GDF_CONVERTER_SELFTEST_H_
