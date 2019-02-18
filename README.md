# TEQCSPEC

This is an updated version of [TEQC multipath metric](http://www.mathworks.com/matlabcentral/fileexchange/12886-teqcspec/content/Teqcspec/teqcspec.m) originally written by Clement Ogaja.
Description of functions can be found in [his paper][paper]. This code was used during [24VSP Satellite Based Positioning Practical Project 1 (Multipath)](https://github.com/DfAC/TeachingSlides/tree/master/H24VLP_P1_MP).

Introduced changes allows to run script in new version of Matlab (past 2012) and to use RINEX 2.11. There also have been small changes in layout allowing easier use for the students of H24VLP MSc Module at the University of Nottingham.


## Input

Software **require [COMPACT2 format][COMPACT2_format] produced by [teqc](http://bit.ly/1KfxvZM) (*teqc +qc +plot data.16o*). You need to use old binaries as the latest ones (2015+) use COMPACT3 format**. For more information check [this paper][paper] and [the issue description](https://github.com/DfAC/TEQCSPEC/issues/2), which also have link to old Windows binaries.

The sky plot will be based on selected file (*.sn* or *.mp*) yet they are only meaningful for *.mp* files.


Note that command for COMPACT3 qc have changed  `teqc +qc -nav *.nav,*.gnav,*.lnav +plot -max_rx_SVs 40 your_filename`. It will still produce COMPACT3 files only!



## Contributions

* Clement Ogaja
* Lei Yang
* Lukasz K Bonenberg



[paper]: http://link.springer.com/article/10.1007%2Fs10291-006-0052-6 "Clement Ogaja 2007 paper"
[COMPACT2_format]: https://postal.unavco.org//pipermail/teqc/2009/000827.html
