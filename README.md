[![Build Status][travis-image]][travis-url]
[![license][MIT-licence]](https://raw.githubusercontent.com/cnclabs/smore/master/LICENSE)
[![Gitter chat for developers at https://gitter.im/dmlc/xgboost][gitter-image]](https://gitter.im/proNet-core/Lobby)

[gitter-image]: https://badges.gitter.im/Join%20Chat.svg
[travis-image]: https://api.travis-ci.org/cnclabs/smore.svg?branch=master
[travis-url]: https://travis-ci.org/cnclabs/smore
[MIT-licence]: https://img.shields.io/badge/License-MIT-yellow.svg

## [RecSys'19 Slides Here](https://www.slideshare.net/changecandy/recsys19-smore)

# SMORe

This is a C++ framework for variant **weighted** network embedding techniques. We currently release the command line interface for following models:
- DeepWalk
  - [DeepWalk: online learning of social representations](http://dl.acm.org/citation.cfm?id=2623732)
- Walklets
  - [Don't Walk, Skip! Online Learning of Multi-scale Network Embeddings](https://arxiv.org/abs/1605.02115)
- LINE(**L**arge-scale **I**nformation **N**etwork **E**mbedding)
  - [LINE: Large-scale Information Network Embedding](http://dl.acm.org/citation.cfm?id=2741093) 
- HPE (**H**eterogeneous **P**reference **E**mbedding)
  - [Query-based Music Recommendations via Preference Embedding](http://dl.acm.org/citation.cfm?id=2959169)
- APP (**A**symmetric **P**roximity **P**reserving graph embedding)
  - [Scalable Graph Embedding for Asymmetric Proximity](https://aaai.org/ocs/index.php/AAAI/AAAI17/paper/view/14696)
- MF (**M**atrix **F**actorization)
- BPR (**B**ayesian **P**ersonalized **R**anking)
  - [BPR: Bayesian personalized ranking from implicit feedback](https://dl.acm.org/citation.cfm?id=1795167)
- WARP-like
  - [WSABIE: Scaling Up To Large Vocabulary Image Annotation](https://dl.acm.org/citation.cfm?id=2283856)
  - [Learning to Rank Recommendations with the k-Order Statistic Loss](https://dl.acm.org/citation.cfm?id=2507157.2507210)
- HOP-REC
  - [HOP-Rec: High-Order Proximity for Implicit Recommendation](https://dl.acm.org/citation.cfm?id=3240381)
- CSE (named nemf & nerank in cli)
  - [Collaborative Similarity Embedding for Recommender Systems](https://arxiv.org/abs/1902.06188)

In the near future, we will redesign the framework making some solid APIs for fast development on different network embedding techniques.

# Developed Environment
- g++ > 4.9 (In macOS, it needs OpenMP-enabled compilers. e.g. brew reinstall gcc6 --without-multilib)

# Compilation
```
$ git clone https://github.com/cnclabs/smore
$ cd smore
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
Directly call the execution file to see the usage like:
```
./cli/deepwalk
./cli/walklets
./cli/line
./cli/hpe
./cli/app
./cli/mf
./cli/bpr
./cli/warp
./cli/hoprec
```
then you will see the options description like:
```
Options Description:
        -train <string>
                Train the Network data
        -save <string>
                Save the representation data
        -dimensions <int>
                Dimension of vertex representation; default is 64
        -undirected <int>
                Whether the edge is undirected; default is 1
        -negative_samples <int>
                Number of negative examples; default is 5
        -window_size <int>
                Size of skip-gram window; default is 5
        -walk_times <int>
                Times of being staring vertex; default is 10
        -walk_steps <int>
                Step of random walk; default is 40
        -threads <int>
                Number of training threads; default is 1
        -alpha <float>
                Init learning rate; default is 0.025
Usage:
./deepwalk -train net.txt -save rep.txt -undirected 1 -dimensions 64 -walk_times 10 -walk_steps 40 -window_size 5 -negative_samples 5 -alpha 0.025 -threads 1
```

# Example Script
This shell script will help obtain the representations of the Youtube links in [Youtube-links](http://socialnetworks.mpi-sws.mpg.de/data/youtube-links.txt.gz) dataset.
```sh
cd example
sh train_youtube.sh
```
Changing the number of threads in *train_youtube.sh* could speedup the process.

# Related Work
You can find related work from [awesome-network-embedding](https://github.com/chihming/awesome-network-embedding).

# Citation
```
@inproceedings{smore,
author = {Chen, Chih-Ming and Wang, Ting-Hsiang and Wang, Chuan-Ju and Tsai, Ming-Feng},
title = {SMORe: Modularize Graph Embedding for Recommendation},
year = {2019},
booktitle = {Proceedings of the 13th ACM Conference on Recommender Systems},
series = {RecSys ’19}
}
```
```
@article{pronet2017,
  title={Vertex-Context Sampling for Weighted Network Embedding},
  author={Chih-Ming Chen and Yi-Hsuan Yang and Yian Chen and Ming-Feng Tsai},
  journal={arXiv preprint arXiv:{1711.00227}},
  year={2017}
}
```

# Note
for HOP-REC & CSE, it is required to assign the field of each vertex in "vertex field" form:
```
userA u
userB u
userC u
itemA i
itemB i
itemC i
itemD i
```
by ``-field`` argument.
