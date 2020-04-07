/*==========================================================================
               SeqAn - The Library for Sequence Analysis
                         http://www.seqan.de
============================================================================
Copyright (C) 2007

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
==========================================================================*/

//#define SEQAN_TCOFFEE_DEBUG

#include <cstdlib>

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

using namespace seqan;

struct minimizer
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor. Attention: all operations on a solely default constructed decorator,
    //!       except assigning a new range, are UB.
    constexpr minimizer() = default;
    //!\brief Copy constructor.
    constexpr minimizer(minimizer const &) = default;
    //!\brief Copy construction via assignment.
    constexpr minimizer & operator=(minimizer const &) = default;
    //!\brief Move constructor.
    constexpr minimizer(minimizer && rhs) = default;
    //!\brief Move assignment.
    constexpr minimizer & operator=(minimizer && rhs) = default;
    //!\brief Use default deconstructor.
    ~minimizer() = default;

    //!\brief Copy constructor from uint64_t.
    constexpr minimizer(uint64_t const v) : value{v} {};
    //!\brief Copy construction via assignment from uint64_t.
    constexpr minimizer & operator=(uint64_t const v) { value = v; return *this; };
    //!\brief Move constructor from uint64_t.
    constexpr minimizer(uint64_t && v) : value{v} {};
    //!\brief Move assignment from uint64_t.
    constexpr minimizer & operator=(uint64_t && v) { value = v; return *this; };

    constexpr minimizer(uint64_t const v, uint64_t const p) : value{v}, position{p} {};

    operator uint64_t() const
    {
        return value;
    }

    // operator unsigned() const
    // {
    //     return value;
    // }
    //!\}

    uint64_t value{};
    uint64_t position{};

    constexpr friend bool operator==(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value == rhs.value;
    }

    constexpr friend bool operator!=(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value != rhs.value;
    }

    constexpr friend bool operator<(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value < rhs.value;
    }
    constexpr friend bool operator<=(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value <= rhs.value;
    }
    constexpr friend bool operator>(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value > rhs.value;
    }
    constexpr friend bool operator>=(minimizer const & lhs, minimizer const & rhs) noexcept
    {
        return lhs.value >= rhs.value;
    }
};

struct Minimizer
{
public:

    // Random, but static value for xor for hashes. Counteracts consecutive minimizers.
    // E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
    uint64_t const seed{0x8F3F73B5CF1C9ADE};
    // Shape for forward hashes
    Shape<Dna, SimpleShape> kmerShape;
    // Shape for hashes on reverse complement
    Shape<Dna, SimpleShape> revCompShape;
    // k-mer size
    uint16_t k{20};
    // window size
    uint32_t w{20};
    // start positions of minimizers
    std::vector<uint64_t> minBegin;
    // end positions of minimizers
    std::vector<uint64_t> minEnd;

    template<typename TIt>
    inline void hashInit(TIt it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIt>
    inline auto hashNext(TIt it)
    {
        return seqan::hashNext(kmerShape, it);
    }

    template<typename TIt>
    inline void revHashInit(TIt it)
    {
        seqan::hashInit(revCompShape, it);
    }

    template<typename TIt>
    inline auto revHashNext(TIt it)
    {
        return seqan::hashNext(revCompShape, it);
    }

    inline auto length()
    {
        return seqan::length(kmerShape);
    }

    inline void resize(uint16_t newKmerSize, uint32_t neww)
    {
        k = newKmerSize;
        w = neww;
        seqan::resize(kmerShape, k);
        seqan::resize(revCompShape, k);
    }

    String<minimizer>  getMinimizer(String<Dna> const & text)
    {
        if (k > seqan::length(text))
            return String<minimizer>{};

        uint64_t possible = seqan::length(text) > w ? seqan::length(text) - w + 1 : 1;
        uint32_t windowKmers = w - k + 1;

        String<minimizer> kmerHashes{};
        reserve(kmerHashes, possible); // maybe rather reserve to expected?

        // Stores hash, begin and end for all k-mers in the window
        std::deque<uint64_t> windowValues;

        auto it = begin(text);
        hashInit(it);

        // Initialisation. We need to compute all hashes for the first window.
        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = hashNext(it) ^ seed;

            windowValues.push_back(kmerHash);

            ++it;
        }

        auto min = std::min_element(std::begin(windowValues), std::end(windowValues));
        appendValue(kmerHashes, minimizer{*min, static_cast<uint64_t>(std::distance(std::begin(windowValues), min))});
        // appendValue(kmerHashPoss, );

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        uint64_t current_pos{kmerHashes[0].position};
        for (uint64_t pos = 1; pos < possible; ++pos)
        {
            windowValues.pop_front();

            uint64_t kmerHash = hashNext(it) ^ seed;
            windowValues.push_back(kmerHash);
            ++it;

            min = std::min_element(std::begin(windowValues), std::end(windowValues));

            if (current_pos != pos + std::distance(std::begin(windowValues), min))
            {
                current_pos = pos + std::distance(std::begin(windowValues), min);
                appendValue(kmerHashes, minimizer{*min, current_pos});
            }
        }

        return kmerHashes;
    }
};


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

inline void
_addVersion(ArgumentParser & parser)
{
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TSeqSet, typename TNameSet>
bool
_loadSequences(TSeqSet& sequences, TNameSet& fastaIDs, const char *fileName)
{
    SeqFileIn inFile;
    if (!open(inFile, fileName))
    {
        std::cerr << "Could not open " << fileName << "for reading!" << std::endl;
        return false;
    }
    readRecords(fastaIDs, sequences, inFile);
    return (length(fastaIDs) > 0u);
}

template <typename TNameSet>
bool
_loadSequences(StringSet<String<minimizer>, Owner<>> & minimizer_sequences,
               TNameSet& fastaIDs,
               const char *fileName)
{
    SeqFileIn inFile;
    if (!open(inFile, fileName))
    {
        std::cerr << "Could not open " << fileName << "for reading!" << std::endl;
        return false;
    }

    StringSet<String<Dna>, Owner<>> sequences;
    {
        StringSet<String<Iupac>, Owner<>> iupac_sequences;
        readRecords(fastaIDs, iupac_sequences, inFile);
        resize(sequences, length(iupac_sequences));

        for (size_t idx = 0; idx < length(iupac_sequences); ++idx)
        {
            for (size_t jdx = 0; jdx < length(iupac_sequences[idx]); ++jdx)
            {
                appendValue(sequences[idx], iupac_sequences[idx][jdx]);
            }
        }
    }

    // compute minimizers per sequence and store the corresponding chain in minimizer_sequences
    resize(minimizer_sequences, length(sequences));
    Minimizer mini;
    mini.resize(12, 20);
    for (size_t idx = 0; idx < length(sequences); ++idx)
    {
        minimizer_sequences[idx] = mini.getMinimizer(sequences[idx]);
        std::cout << length(minimizer_sequences[idx]) << std::endl;
    }

    return (length(fastaIDs) > 0u);
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TScore>
inline void
customizedMsaAlignment(MsaOptions<TAlphabet, TScore> const& msaOpt)
{
    typedef String<TAlphabet> TSequence;
    StringSet<TSequence, Owner<> > sequenceSet;
    StringSet<String<char> > sequenceNames;
    typedef VirtualStream<char, Output> TOutStream;

    // Alignment of the sequences
    typedef Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > TGraph;
    TGraph gAlign;
    typedef typename Size<StringSet<TSequence, Dependent<> >>::Type TSize;

    segment_generation_config<TSize> seg_gen_config;
    seg_gen_config.seqfiles = msaOpt.seqfiles;

    if (length(msaOpt.seqfiles) == 1)
    {
        // compute MSA on all pairs of pairs in multi fasta file
        _loadSequences(sequenceSet, sequenceNames, msaOpt.seqfiles[0].c_str());
        selectPairs(sequenceSet, seg_gen_config.global_alignment_pairs); // all-to-all
    }
    else // multi multi-record fasta files
    {
        assert(length(msaOpt.seqfiles) > 1);
        // TODO
    }

    // MSA
    try
    {
        globalMsaAlignment(gAlign, sequenceSet, sequenceNames, seg_gen_config, msaOpt);
    }
    catch (const std::bad_alloc & exception)
    {
        std::cerr << "Allocation for globalAlignment failed. Use smaller data or try a seeded alignment. \n"
                  << exception.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Alignment output
    TOutStream outStream;
    if (!open(outStream, toCString(msaOpt.outfile)))
    {
        std::cerr << "Can't open " << msaOpt.outfile << " for writing!" << std::endl;
        return;
    }
    if (guessFormatFromFilename(msaOpt.outfile, Fasta()))
        write(outStream, gAlign, sequenceNames, FastaFormat());
    else
        write(outStream, gAlign, sequenceNames, MsfFormat());

    String<String<char> > nodeMap;
    String<String<char> > edgeMap;
    _createEdgeAttributes(gAlign, edgeMap);

    // create node attributes
    typedef typename Id<TGraph>::Type TIdType;
    resizeVertexMap(nodeMap, gAlign);

    typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
    TConstIter it(gAlign);
    for(;!atEnd(it);++it) {
        TIdType id = sequenceId(gAlign, *it);
        std::ostringstream outs;
        outs << "label = \"";
        outs << "[";
        auto regStart = fragmentBegin(gAlign, *it);
        if (regStart == 0)
            outs << "0"; // if it is the very first minimizer, include beginning of the sequence
        else
            outs << sequenceSet[id][regStart].position;
        outs << ",";
        auto regEnd = fragmentBegin(gAlign, *it) + fragmentLength(gAlign, *it);
        if (regEnd >= length(sequenceSet[id]))
            outs << "end";
        else
            outs << sequenceSet[id][regEnd].position;
        outs << ")";
        outs << "\", group = ";
        outs << id;
        append(property(nodeMap, *it), outs.str().c_str());
        //std::cout << property(nodeMap, *it) << std::endl;
    }

    std::ofstream dotFile("graph.dot");
    writeRecords(dotFile, gAlign, nodeMap, edgeMap,DotDrawing());
    dotFile.close();
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TScore, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, TScore>&, TSc)
{
    // No operation
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TScore, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, TScore>&, TSc)
{
    // No operation
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc msc)
{
    msaOpt.sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc mmsc)
{
    msaOpt.sc.data_mismatch = mmsc;
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TScore>
inline void
_initMsaParams(ArgumentParser& parser, TScore& scMat)
{

    // Msa configuration
    MsaOptions<TAlphabet, TScore> msaOpt;

    // Set main options
    getOptionValue(msaOpt.outfile, parser, "outfile");

    std::string tmpVal;
    for (unsigned int optNo = 0; optNo < getOptionValueCount(parser, "seq"); ++optNo)
    {
        getOptionValue(tmpVal, parser, "seq");
        appendValue(msaOpt.seqfiles, tmpVal);
    }

    String<char> optionVal;

 /*   getOptionValue(optionVal, parser, "format");
    if (optionVal == "fasta")
        msaOpt.outputFormat = 0;
    else if (optionVal == "msf")
        msaOpt.outputFormat = 1;
*/
    // *********************************************

    // Set segment match generation options
    for (unsigned int optNo = 0; optNo < getOptionValueCount(parser, "method"); ++optNo)
    {
        getOptionValue(tmpVal, parser, "method", optNo);
        if (tmpVal == "global")
        {
            appendValue(msaOpt.method, 0);
        }
        else if (tmpVal == "local")
        {
            appendValue(msaOpt.method, 1);
        }
        else if (tmpVal == "overlap")
        {
            appendValue(msaOpt.method, 2);
        }
        else if (tmpVal == "lcs")
        {
            appendValue(msaOpt.method, 3);
        }
        else if (tmpVal == "none")
        {
            clear(msaOpt.method);
        }
    }

    msaOpt.isDefaultPairwiseAlignment = !isSet(parser, "pairwise-alignment");
    if (!msaOpt.isDefaultPairwiseAlignment)
    {
        getOptionValue(tmpVal, parser, "pairwise-alignment");
        if (tmpVal == "unbanded")
        {
            msaOpt.pairwiseAlignmentMethod = 1;
        }
        else if (tmpVal == "banded")
        {
            msaOpt.pairwiseAlignmentMethod = 2;
        }
    }
    if (isSet(parser, "band-width")) {
        if (msaOpt.pairwiseAlignmentMethod == 1)
        {
            std::cerr << "Ambiguous pairwise alignment method. Band width cannot be specified for an unbanded method" << std::endl;
            std::exit(0);
        }
        msaOpt.isDefaultPairwiseAlignment = false;
        msaOpt.pairwiseAlignmentMethod = 2;
    }
    getOptionValue(msaOpt.bandWidth, parser, "band-width");

    for (unsigned int optNo = 0; optNo < getOptionValueCount(parser, "libraries"); ++optNo)
    {
        getOptionValue(tmpVal, parser, "libraries", optNo);
        if(endsWith(tmpVal,".blast"))
            appendValue(msaOpt.blastfiles, tmpVal);
        else if(endsWith(tmpVal,".mummer"))
                appendValue(msaOpt.mummerfiles, tmpVal);
            else if(endsWith(tmpVal,".lib"))
                    appendValue(msaOpt.libfiles, tmpVal);
                else if(endsWith(tmpVal,".aln"))
                        appendValue(msaOpt.alnfiles, tmpVal);
    }

/*
    for (unsigned int optNo = 0; optNo < getOptionValueCount(parser, "blast"); ++optNo)
    {
        getOptionValue(tmpVal, parser, "blast", optNo);
        appendValue(msaOpt.blastfiles, tmpVal);
    }

    for (unsigned int optNo = 0; optNo < getOptionValueCount(parser, "mummer"); ++optNo)
    {
        getOptionValue(tmpVal, parser, "mummer", optNo);
        appendValue(msaOpt.mummerfiles, tmpVal);
    }

    for (unsigned int optNo = 0; optNo < getOptionValueCount(parser, "lib"); ++optNo)
    {
        getOptionValue(tmpVal, parser, "lib", optNo);
        appendValue(msaOpt.libfiles, tmpVal);
    }

    for (unsigned int optNo = 0; optNo < getOptionValueCount(parser, "aln"); ++optNo)
    {
        getOptionValue(tmpVal, parser, "aln", optNo);
        appendValue(msaOpt.alnfiles, tmpVal);
    }
*/

// Set scoring options
    msaOpt.sc = scMat;
    getOptionValue(msaOpt.sc.data_gap_open, parser, "gop");
    getOptionValue(msaOpt.sc.data_gap_extend, parser, "gex");
    int msc = 0;
    getOptionValue(msc, parser, "msc");
    _setMatchScore(msaOpt, msc);
    int mmsc = 0;
    getOptionValue(mmsc, parser, "mmsc");
    _setMismatchScore(msaOpt, mmsc);

    // Set guide tree options
    getOptionValue(msaOpt.treefile, parser, "usetree");
    getOptionValue(optionVal, parser, "build");
    if (optionVal == "nj")
        msaOpt.build = 0;
    else if (optionVal == "min")
        msaOpt.build = 1;
    else if (optionVal == "max")
        msaOpt.build = 2;
    else if (optionVal == "avg")
        msaOpt.build = 3;
    else if (optionVal == "wavg")
        msaOpt.build = 4;

    // Set alignment evaluation	options
    getOptionValue(msaOpt.infile, parser, "infile");

    // Check if any segment-match generation procedure is selected, otherwise set the default
    if ((empty(msaOpt.blastfiles)) && (empty(msaOpt.mummerfiles)) && (empty(msaOpt.libfiles))
        && (empty(msaOpt.alnfiles)) && (empty(msaOpt.method)))
    {
        appendValue(msaOpt.method, 0);
        appendValue(msaOpt.method, 1);
    }

    // Evaluation mode?
    if (isSet(parser, "infile"))
    {
        throw "evaluation mode not available.";
    }
    else
    { // or alignment mode?
        if (!isSet(parser, "seq"))
        {
            printShortHelp(parser, std::cerr);	// print short help and exit
            exit(0);
        }
        customizedMsaAlignment(msaOpt);
    }
}

//////////////////////////////////////////////////////////////////////////////////
inline void
_initScoreMatrix(ArgumentParser& parser, minimizer const)
{
    // String<char> matrix;
    // getOptionValue(matrix, parser, "matrix");
    // if (isSet(parser, "matrix"))
    // {
    //     Score<int, ScoreMatrix<> > sc;
    //     loadScoreMatrix(sc, toCString(matrix));
    //     _initMsaParams<Dna5>(parser, sc);
    // }
    // else
    // {
    std::cout << "i am here" << std::endl;
        Score<int> sc;
        _initMsaParams<minimizer>(parser, sc);
    // }
}

void
_setUpArgumentParser(ArgumentParser & parser)
{
    _addVersion(parser);
    setDate(parser, SEQAN_DATE);
    setAppName(parser,"seqan_tcoffee");
    setCategory(parser, "Sequence Alignment");

    setShortDescription(parser, "Multiple sequence alignment");

    addUsageLine(parser, "-s <\\fIFASTA FILE\\fP> [\\fIOPTIONS\\fP]");
    addDescription(parser, "SeqAn::T-Coffee is a multiple sequence alignment tool.");
    addDescription(parser, "(c) Copyright 2009 by Tobias Rausch");

    addSection(parser, "Main Options:");
    addOption(parser, ArgParseOption("s", "seq", "Name of multi-fasta input file.", ArgParseArgument::INPUT_FILE, "foo", true));
    // setValidValues(parser, "seq", getFileExtensions(Fasta()));  // allow only fasta files as input

    addOption(parser, ArgParseOption("a", "alphabet", "The used sequence alphabet.", ArgParseArgument::STRING));
    setValidValues(parser, "alphabet", "protein dna rna iupac minimizer");
    setDefaultValue(parser, "alphabet", "protein");

    addOption(parser, ArgParseOption("o", "outfile", "Name of the output file.", ArgParseArgument::OUTPUT_FILE));
    setDefaultValue(parser, "outfile", "out.fasta");
    std::vector<std::string> outputFormats = getFileExtensions(Fasta());
    outputFormats.push_back(".msf");
    setValidValues(parser, "outfile", outputFormats);


  //  addOption(parser, ArgParseOption("f", "format", "Format of the output.", ArgParseArgument::STRING));
 //   setValidValues(parser, "format", "fasta msf");
 //   setDefaultValue(parser, "format", "fasta");

    addSection(parser, "Segment Match Generation Options:");
    addOption(parser,
              ArgParseOption("m", "method", "Defines the generation method for matches. "
                             "To select multiple generation methods recall this option with different arguments. ",
                             ArgParseArgument::STRING, "", true));
    setValidValues(parser, "method", "global local overlap lcs none");
    setDefaultValue(parser, "method", "global");
    addDefaultValue(parser, "method", "local");

    addOption(parser,
              ArgParseOption("l", "libraries", "Name of match file. "
                             "To select multiple files recall this option with different arguments.",
                             ArgParseArgument::INPUT_FILE, "", true));

    setValidValues(parser, "l", "blast mums aln lib");  // allow blast, mummer aln and tcoffee lib files

    addOption(parser,
              ArgParseOption("pa", "pairwise-alignment", "Pairwise alignment method. "
                             "Default: \\fIunbanded\\fP for usual alignments (< 50 sequences), "
                             "\\fIbanded\\fP for deep alignments (>= 50 sequences)",
                             ArgParseArgument::STRING));
    setValidValues(parser, "pa", "unbanded banded");

    addOption(parser,
              ArgParseOption("bw", "band-width", "Band width. "
                             "This option automatically select \\fIbanded\\fP pairwise alignment "
                             "(see \\fBpa\\fP for details)",
                             ArgParseArgument::INTEGER));
    setDefaultValue(parser, "bw", 60);
    setMinValue(parser, "bw", "2");

    // code before KNIME adaption
    /*   addOption(parser,
              ArgParseOption("bl", "blast", "Name of \\fIBLAST\\fP match file. "
                             "To select multiple \\fIBLAST\\fP files recall this option with different arguments.",
      addOption(parser,
              ArgParseOption("mu", "mummer", "Name of \\fIMUMmer\\fP match file. "
                             "To select multiple \\fIMUMmer\\fP files recall this option with different arguments.",
                             ArgParseArgument::INPUT_FILE, "", true));
    addOption(parser,
              ArgParseOption("al", "aln", "Name of \\fIFASTA\\fP align file."
                             "To select multiple \\fIFASTA\\fP files recall this option with different arguments.",
                             ArgParseArgument::INPUT_FILE, "", true));
    addOption(
            parser,
            ArgParseOption("li", "lib", "Name of \\fIT-Coffee\\fP library. "
                           "To select multiple \\fIT-Coffee\\fP libraries recall this option with different arguments.",
                           ArgParseArgument::INPUT_FILE, "", true));
     */

    addSection(parser, "Scoring Options:");
    addOption(parser, ArgParseOption("g", "gop", "gap open penalty", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "gop", -13);
    addOption(parser, ArgParseOption("e", "gex", "gap extension penalty", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "gex", -1);
    addOption(parser, ArgParseOption("ma", "matrix", "score matrix", ArgParseArgument::STRING));
    setDefaultValue(parser, "matrix", "Blosum62");
    addOption(parser, ArgParseOption("ms", "msc", "match score", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "msc", 5);
    addOption(parser, ArgParseOption("mm", "mmsc", "mismatch penalty", ArgParseArgument::INTEGER));
    setDefaultValue(parser, "mmsc", -4);

    addSection(parser, "Guide Tree Options:");
    addOption(
            parser,
            ArgParseOption("u", "usetree", "Name of the file containing the \\fINewick Guide Tree\\fP.",
                           ArgParseArgument::STRING));
    addOption(
            parser,
            ArgParseOption(
                    "b",
                    "build",
                    "Method to build the tree. "
                    "Following methods are provided: \\fINeighbor-Joining\\fP (\\fBnj\\fP), \\fIUPGMA single linkage\\fP "
                    "(\\fBmin\\fP), \\fIUPGMA complete linkage\\fP (\\fBmax\\fP), \\fIUPGMA average linkage\\fP "
                    "(\\fBavg\\fP), \\fIUPGMA weighted average linkage\\fP (\\fBwavg\\fP). "
                    "\\fINeighbor-Joining\\fP creates an unrooted tree, which we root at the last joined pair.",
                    ArgParseArgument::STRING));
    setDefaultValue(parser, "build", "nj");
    setValidValues(parser, "build", "nj min max avg wavg");

    addSection(parser, "Alignment Evaluation Options:");
    addOption(
            parser,
            ArgParseOption("i", "infile", "Name of the alignment file <\\fIFASTA FILE\\fP>",
                           ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "infile", getFileExtensions(Fasta()));  // allow only fasta files as input

}

//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv [])
{
#ifdef SEQAN_TCOFFEE_DEBUG
    double totalStartTime = sysTime();
    std::cout << std::fixed << std::setprecision(5);
#endif

    // Command line parsing
    ArgumentParser parser("seqan_tcoffee");
    _setUpArgumentParser(parser);

    if (argc == 1)
    {
        printShortHelp(parser, std::cerr);	// print short help and exit
        return 0;
    }

    ArgumentParser::ParseResult res = parse(parser, argc, argv, std::cout, std::cerr);

    if (res != ArgumentParser::PARSE_OK)
    {
        return res == ArgumentParser::PARSE_ERROR;
    }
    if (isSet(parser, "help") || isSet(parser, "version"))
        return 0;	// print help or version and exit

    // Basic command line options
    String<char> alphabet;
    getOptionValue(alphabet, parser, "alphabet");

    // Initialize scoring matrices
    if  (alphabet == "minimizer")
        _initScoreMatrix(parser, minimizer());
    else
        throw "you need minimizer!";

#ifdef SEQAN_TCOFFEE_DEBUG
    std::cout << std::setw(30) << std::left << "Total time:" << std::setw(10) << std::right << sysTime() - totalStartTime << "  s" << std::endl;
#endif

    return 0;
}
