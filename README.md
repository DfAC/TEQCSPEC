# TEQCSPEC

This is an updated version of [TEQC multipath metric](http://www.mathworks.com/matlabcentral/fileexchange/12886-teqcspec/content/Teqcspec/teqcspec.m) originally written by Clement Ogaja.
Description of functions can be found in [his paper][paper].

Introduced changes allows to run script in new version of Matlab (past 2012) and to use RINEX 2.11. There have been small changes in layout allowing easier use for the students of H24VLP MSc Module at the University of Nottingham.

This code is used during [24VSP Satellite Based Positioning Practical Project 1 (Multipath)](https://github.com/DfAC/TeachingSlides/tree/master/H24VLP_P1_MP).

## Input

Software require COMPACT2 qc files produced by [teqc](http://bit.ly/1KfxvZM) (*teqc +qc +plot data.16o*). You need to use old binaries as the latest ones (2015+) use COMPACT3 format. For more information check [this paper][paper].


Note that sky plot will be based on selected file (*.sn* or *.mp*) yet they are only meaningful for *.mp* ones.

## Contributions

* Clement Ogaja
* Lei Yang
* Lukasz K Bonenberg



[paper]: http://link.springer.com/article/10.1007%2Fs10291-006-0052-6 "Clement Ogaja 2007 paper"
