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

    std::cout << "|s1| = " << length(seqH) << std::endl;
    std::cout << "|s2| = " << length(seqV) << std::endl;
    std::cout << "|M|  = " << (length(seqH) + 1) * (length(seqV) + 1) << std::endl;

    double time = 0.0;
    int score = 0;
    for (unsigned i = 0; i < 100; ++i)
    {
        double start = sysTime();
        score = globalAlignmentScore(seqH, seqV, Score<int, Simple>(4, -5, -2, -6));
        time += sysTime() - start;
    }

    std::cout << "Score: " << score << std::endl;
    std::cout << "Time: " << (time / 100.0) << std::endl;
    return 0;
}