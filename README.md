# SeqAn - The Library for Sequence Analysis

[![build status][1]][2]
[![license][3]][4]
[![latest release][5]][6]
[![platforms][7]][8]
[![start][9]][10]
[![twitter][11]][12]

[1]: https://img.shields.io/github/workflow/status/seqan/seqan/CI%20on%20Linux/develop?style=flat&logo=github&label=SeqAn%20CI "Open GitHub actions page"
[2]: https://github.com/seqan/seqan/actions?query=branch%3Adevelop
[3]: https://img.shields.io/badge/license-BSD-green.svg "Open license file"
[4]: https://github.com/seqan/seqan/blob/develop/LICENSE
[5]: https://img.shields.io/github/release/seqan/seqan.svg "Get the latest release"
[6]: https://github.com/seqan/seqan/releases/latest
[7]: https://img.shields.io/badge/platform-linux%20%7C%20bsd%20%7C%20osx%20%7C%20win-informational.svg "Open our API documentation"
[8]: https://docs.seqan.de/seqan/develop/
[9]: https://img.shields.io/github/stars/seqan/seqan.svg?style=social "See who starred us"
[10]: https://github.com/seqan/seqan/stargazers
[11]: https://img.shields.io/twitter/follow/SeqAnLib.svg?label=follow&style=social "Follow us on Twitter"
[12]: https://twitter.com/seqanlib

| **NOTE <br> [SeqAn3 is out and hosted in a different repository](https://github.com/seqan/seqan3)**                        |
|:--------------------------------------------------------------------------------------------------------------------------:|
| All new applications will be based on SeqAn3. <br> Existing applications that receive updates will successively be ported. |

## What Is SeqAn?

SeqAn is an open source C++ library of efficient algorithms and data structures for the analysis of sequences with the focus on biological data.
Our library applies a unique generic design that guarantees high performance, generality, extensibility, and integration with other libraries.
SeqAn is easy to use and simplifies the development of new software tools with a minimal loss of performance.

## License

The SeqAn library itself, the tests and demos are licensed under the very permissive 3-clause BSD License.
The licenses for the applications themselves can be found in the LICENSE files.

## Prerequisites

Older compiler versions might work but are not supported.

### Linux, macOS, FreeBSD
  * GCC ≥ 7
  * Clang/LLVM ≥ 9
  * Intel Compiler ≥ 2021.3.0 (Intel OneAPI)

### Windows
  * Visual C++ ≥ 16.0 / Visual Studio ≥ 2019

### Architecture support
  * Intel/AMD platforms, including optimisations for modern instruction sets (`POPCNT`, `SSE4`, `AVX2`, `AVX512`)
  * All Debian release architectures supported, including most ARM and all PowerPC platforms.

### Build system
  * To build tests, demos, and official SeqAn applications you also need CMake ≥ 3.12.

Some official applications might have additional requirements or only work on a subset of platforms.

## Documentation Resources

* [Getting Started](https://seqan.readthedocs.io/en/main/Tutorial/GettingStarted)
* [Manual](https://seqan.readthedocs.io/en/main)
* [Tutorial](https://seqan.readthedocs.io/en/main/index.html#tutorials)
* [How-Tos](https://seqan.readthedocs.io/en/main/Tutorial/HowTo)
* [API Documentation (stable)](https://docs.seqan.de/seqan/main/)

## Contact

* [Mailing List](https://lists.fu-berlin.de/listinfo/seqan-dev#subscribe)
* [GitHub Project (issues, source code)](https://github.com/seqan/seqan)
