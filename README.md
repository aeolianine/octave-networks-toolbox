# octave-networks-toolbox: A set of graph/networks analysis functions in Octave, 2012-2023

## Quick description

This is a repository of functions relevant to network/graph analysis, organized by functionality. These routines are useful for someone who wants to start hands-on work with networks fairly quickly, explore simple graph statistics, distributions, simple visualization and compute common network theory metrics.

## History

The original (2006-2011) version of these routines was written in Matlab, and is still hosted by strategic.mit.edu (http://strategic.mit.edu/downloads.php?page=matlab_networks). The octave-networks-toolbox inherits the original BSD open source license and copyright, provided at the end of this file. Many of the routines might still be compatible with Matlab. For Octave/Matlab differences, see http://en.wikibooks.org/wiki/MATLAB_Programming/Differences_between_Octave_and_MATLAB.

## Installation

The code currently runs on GNU Octave version 6.3.0. No specific library installation necessary. Dependencies between functions are documented in the function headers. The routines can be called directly from the Octave prompt, either in the same directory or from anywhere if the toolbox folder is added to the path. For example:

```matlab
octave:1> % running numNodes.m
octave:1> numNodes([0 1 1; 1 0 1; 1 1 0])
ans =  3
```

To run tests:
```matlab
octave:10> % test a single function
octave:10> test numEdges  
PASSES 1 out of 1 test
octave:11>
octave:11> % test all functions
octave:11> oruntests
```

## Matlab compatibility

With newer versions of Matlab, the Octave branch may not always be Matlab-compatible, for example due to syntax changes. Consider exploring forks that focus on Matlab compatibility. There is currently no plan for the Octave original version to be synchronized with Matlab.

## Authorship

This code was originally written and posted by Gergana Bounova. It is undergoing continuous expansion and development. Collaborators are very welcome. Thank you for the many comments and bug reports so far! Contributions via email are usually mentioned in the function header. Please use github for comments, questions, suggestions, issues, bugs or simply fork.

## Organization

The functions are organized in 11 categories: basic routines, diagnostic routines, representation routines, centralities, distances, simple motif routines, linear algebra functions, modularity routines, graph construction models, visualization and auxiliary. These categories reflect roles/functionality and topics in the literature, but they are arbitrary, and mostly used for documentation purposes.

## Documentation

Documentation is available in the [Functions Manual](https://github.com/aeolianine/octave-networks-toolbox/blob/master/functions_manual.pdf). The manual contains general background information, function headers, code examples, and references. For some functions, additional background, definitions or derivations are included. 

## Citation

If you want to cite this code, you can use DOI: [10.5281/zenodo.22398](http://dx.doi.org/10.5281/zenodo.22398). This citation refers to the second release, from August 2 2015.

## License/Copyright

Copyright (c) 2013-2023, Massachusetts Institute of Technology.
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the Massachusetts Institute of Technology nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.