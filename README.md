[![Build Status][travis-image]][travis-url]
[![Coverage Status][codecov-image]][codecov-url]
[travis-image]: https://img.shields.io/travis/chihming/proNet-core.svg
[travis-url]: https://travis-ci.org/chihming/proNet-core
[codecov-image]: https://img.shields.io/codecov/c/github/chihming/proNet-core.svg
[codecov-url]: https://codecov.io/gh/chihming/proNet-core

# proNet-core
This is a C++ framework for variant weighted network embeddings. We currently release the command line interface for following models:
- [DeepWalk](http://dl.acm.org/citation.cfm?id=2623732)
- [Walklets](https://arxiv.org/abs/1605.02115)
- [LINE](http://dl.acm.org/citation.cfm?id=2741093)
- [HPE](http://dl.acm.org/citation.cfm?id=2959169)

In the near future, we will redesign the framework making some solid APIs for fast development on different network embedding techniques.

# Developed Environment
- g++ > 4.9 (In macOS, it needs OpenMP-enabled compilers. Try ```brew install gcc --without-multilib```)

# Related Work
You can find related work from [awesome-network-embedding](https://github.com/chihming/awesome-network-embedding).

# Compilation
```
$ git clone https://github.com/chihming/proNet-core
$ cd proNet-core
$ make
```

# Task
Given a network input:
```txt
userA itemA 3
userA itemC 5
userB itemA 1
userB itemB 5
userC itemA 4
```
The model learns the representations of each vertex:
```
6 5
userA 0.0815412 0.0205459 0.288714 0.296497 0.394043
itemA -0.207083 -0.258583 0.233185 0.0959801 0.258183
itemC 0.0185886 0.138003 0.213609 0.276383 0.45732
userB -0.0137994 -0.227462 0.103224 -0.456051 0.389858
itemB -0.317921 -0.163652 0.103891 -0.449869 0.318225
userC -0.156576 -0.3505 0.213454 0.10476 0.259673
```

# Command Line Interface
Directly call cli to see the comment usage:
```
./cli
```
An example comment:
```
./cli -model DeepWalk -train net.txt -output rep.txt -dimensions 64 -window_size 5 -negative_samples 5 -alpha 0.025 -threads 4
```

# Example Script
This shell script will help obtain the representations of the Youtube links in [Youtube-links](http://socialnetworks.mpi-sws.mpg.de/data/youtube-links.txt.gz) dataset.
```sh
cd example
sh train_youtube.sh
```
