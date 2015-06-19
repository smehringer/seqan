// ==========================================================================
//                                   lagan app
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

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include "laganAlignment.h"

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class LaganOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
struct LaganOptions
{
    typedef seqan::String<unsigned> Tuns;
    // arguments
    Tuns lagan_parameter;
    // allocate space for 3 parameter-sets of 3 unsigned

    unsigned q; // lagan parameter
    unsigned bandExtension;
    unsigned scoreM, scoreMM, scoreG; // scoringsheme

    seqan::CharString filename;
    seqan::CharString filenameOUT;

    LaganOptions() :
        q(4),
        bandExtension(2),
        scoreM(0), scoreMM(-1), scoreG(-1) // edit distance
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(LaganOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("my_app");
    // Set short description, version, and date.
    setShortDescription(parser, "This app can compute the LAGAN algorithm on two input sequences");
    setVersion(parser, "0.1");
    setDate(parser, "April 2015");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP]");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    //  options
    addOption(parser, seqan::ArgParseOption(
        "i", "input-file", "Path to the input file.",
        seqan::ArgParseArgument::INPUT_FILE, "FASTA"));
    setValidValues(parser, "input-file", "fa");
    setRequired(parser, "i");

    addOption(parser, seqan::ArgParseOption(
        "q", "q_gram",
        "neccessary lagan parameter: "
        "q - Size of the q-gram(-index) when finding seeds. Define multiple times for iterative search ",
        seqan::ArgParseArgument::INTEGER, "INT",true,1));
    setRequired(parser, "q");

    addOption(parser, seqan::ArgParseOption(
        "s", "scoringsheme", 
        "Scoring Sheme: (Match, Mismatch, Gap)",
        seqan::ArgParseArgument::INTEGER, "INT INT INT",false,3));

    addOption(parser, seqan::ArgParseOption(
        "be", "bandExtension", 
        "The allowed extension around seeds that will be used for bandedChainAlignment()",
        seqan::ArgParseArgument::INTEGER, "INT"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.filename, parser, "input-file");

    getOptionValue(options.q, parser, "q_gram",0);
    resize(options.lagan_parameter, getOptionValueCount(parser, "q_gram"));

    for (unsigned i = 0; i < getOptionValueCount(parser, "q_gram"); ++i)
    {
        getOptionValue(options.lagan_parameter[i], parser, "q_gram",i);
    }

    getOptionValue(options.scoreM, parser, "scoringsheme",0);
    getOptionValue(options.scoreMM, parser, "scoringsheme",1);
    getOptionValue(options.scoreG, parser, "scoringsheme",2);

    getOptionValue(options.bandExtension, parser, "bandExtension");

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
using namespace seqan;

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    LaganOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "LAGAN ALGORITHM" 
              << "\n=====================================================\n";

    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn;

    std::cout << "Reading in Sequences...\n";
    if (!open(seqFileIn, toCString(options.filename)))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    try{
        readRecords(ids, seqs, seqFileIn);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    DnaString seqH;
    DnaString seqV;
    seqH = seqs[0];
    seqV = seqs[1];
    std::cout << "# Length of Sequences. SeqH: " << length(seqH) 
              << "  SeqV: " << length(seqV) <<"\n";

    Align<DnaString, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seqH);
    assignSource(row(alignment, 1), seqV);

    Score<int, Simple> scoreScheme(options.scoreM, options.scoreMM, options.scoreG);

    // -----------------------------------------------------------------------
    // Compute Alignment
    // -----------------------------------------------------------------------
    std::cout << "# Computing Alignment...\n";
    int result = laganAlignment(seqH, seqV, options.lagan_parameter /*, scoreScheme*/ );
    std::cout << "# Alignment Score: " << result << std::endl;

    std::cout << "=====================================================\n" << "DONE\n";
    return 0;
}
