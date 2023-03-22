<p align="center">

# `WISP`

</p>

**CURRENTLY WIP**

<img align="right" src="devwks/logo.png" alt="wisp logo" width="300"/>

This Python program is meant to determine to which taxa a bacteria is belonginig to, from long reads (>10.000bp), solely based upon alignment-free methods. As of now, the main focus is upon kmers proportions. It aims to do binning over a collection of samples, giving a probable class to each read.

The core functionnalities relies on a class probabiliy attribution to discriminate reads that might not be good indicators for our specie to be determined. As many other options, you can choose the ratio and the selection function to suit best your biological context.

**WISP is research software**. If you want to use it, please source the code. 