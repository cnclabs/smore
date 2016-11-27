# proNet-core
pair-wise representations optimization Network, a c++ framework for varient network embedding models.

Developed Environment:
- g++ 4.9

# Installation
```
$ git clone https://github.com/chihming/proNet-core
$ cd proNet-core
$ make
```

# Task Explaination
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
userA 0.0815412 0.0205459 0.288714 0.296497 0.394043
itemA -0.207083 -0.258583 0.233185 0.0959801 0.258183
itemC 0.0185886 0.138003 0.213609 0.276383 0.45732
userB -0.0137994 -0.227462 0.103224 -0.456051 0.389858
itemB -0.317921 -0.163652 0.103891 -0.449869 0.318225
userC -0.156576 -0.3505 0.213454 0.10476 0.259673
```

# Command Line Interface
Example Usage:
```
./cli -model DeepWalk -train net.txt -output rep.txt -window_size 5 -negative_samples 5 -alpha 0.025 -threads 4
```

