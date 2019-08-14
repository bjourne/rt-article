# Performance Survey of Modern Ray/Triangle Intersection Algorithms

This repository contains raw data and information on how to reproduce
benchmark figures from my unpublished article "Performance Survey of
Modern Ray/Triangle Intersection Algorithms."

## Reproduction instructions

The benchmark figures comes from measuring the performance of the
implementations of the algorithms in Embree. Consequently my fork of
Embree needs to be checked out:

    $ git clone https://github.com/bjourne/embree

All the work is done on the branch "my-hacksK" (fix this). Build my
version of Embree on that branch using the following:

    $ mkdir build
    $ cd build
    $ ccmake -DCMAKE_INSTALL_PREFIX=/tmp/root/ \
        -DEMBREE_ISPC_EXECUTABLE=/path/here \
        -DEMBREE_GEOMETRY_QUAD=OFF -DEMBREE_GEOMETRY_POINT=OFF \
        -DEMBREE_GEOMETRY_CURVE=OFF -DEMBREE_GEOMETRY_GRID=OFF \
        -DEMBREE_GEOMETRY_INSTANCE=OFF -DEMBREE_GEOMETRY_SUBDIVISION=OFF \
        -DEMBREE_GEOMETRY_USER=OFF ..
    $ make

Refer to Embree's official repository for general compilation
instructions.

When running the programs viewer, viewer_stream, viewer_ispc,
pathtracer and pathtracer_ispc, Embree will print out the message

    Running with triangle intersector <something>

indicating which triangle intersection algorithm is active. To change
intersection algorithm, find the line

    #define ISECT_METHOD ISECT_<something>

in `kernels/geometry/triangle.h`, modify it accordingly and recompile
Embree.


## Algorithm implementations

For ease of access, implementations of the benchmarked algorithms can
also be found in the `src` directory in this repository. However,
there may be slight differences in them if I forget to update them in
tandem with the Embree fork.
