#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char* argv[])
{

    CharString idH;
    DnaString seqH;
    CharString idV;
    DnaString seqV;


    try {
        SeqFileIn seqFile1(argv[1]);
        readRecord(idH, seqH, seqFile1);

        SeqFileIn seqFile2(argv[2]);
        readRecord(idV, seqV, seqFile2);
    } catch (...) {
        std::cout << "Caught exception!" << std::endl;
        return -1;

    }

    int score = 0;
    for (unsigned i = 0; i < 100; ++i)
        score = globalAlignmentScore(seqH, seqV, Score<int, Simple>(4, -5, -2, -6));

    std::cout << "Score: " << score << std::endl;
}