#!/bin/bash
set -eux

########################################################################
# Please configure these paths!
########################################################################
model_dir=/tmp/models
crown=$model_dir/crown/crown.ecs
sanm=$model_dir/sanm/sanm.obj
lucy=$model_dir/Alucy.obj
asian_dragon=$model_dir/asian_dragon/asian_dragon.ecs
embree=/tmp/sources/embree
########################################################################

embree_build=${embree}/build
triangle_h=${embree}/kernels/geometry/triangle.h

# Viewpoints
sanm_args1="--vp 8.854205132 1.856649399 7.360687256 \
         --vi 13.88171768 2.228513956 -2.522052765 \
         --vu 0 1 0 --fov 90 --righthanded"
sanm_args2="--vp 22.80886841 1.800413132 2.112316132 \
        --vi 19.92428589 1.44120419 2.852566957 \
        --vu 0 1 0 --fov 90 --righthanded \
        --size 1024 1024"
lucy_args1="--vp 0.1194648743 16.96469498 3.015714645 \
        --vi 0.1176664829 19.96464157 3.024817467 \
        --vu 0 1 0 --fov 90 --righthanded"
lucy_args2="--vp -283.4347229 645.4871216 184.485611 \
        --vi -280.9531555 644.796936 182.9482117 \
        --vu 0 1 0 --fov 90 --righthanded"

bench_args="--benchmark 400 10000 --threads 1"

single=${embree_build}/viewer
ispc=${embree_build}/viewer_ispc
stream=${embree_build}/viewer_stream

if [ $# = 0 ]; then
    echo "Please specify model: lucy1, lucy2, sanm1, sanm2, crown or dragon."
    exit 1
fi
case $1 in
    lucy1)
    bench_args="-i $lucy $lucy_args1 $bench_args"
    ;;
    lucy2)
    bench_args="-i $lucy $lucy_args2 $bench_args"
    ;;
    sanm1)
    bench_args="-i $sanm $sanm_args1 $bench_args"
    ;;
    sanm2)
    bench_args="-i $sanm $sanm_args2 $bench_args"
    ;;
    crown)
    bench_args="-c $crown $bench_args"
    ;;
    dragon)
    bench-args="-c $asian_dragon $bench_args"
    ;;
    *)
    "Echo wrong model choice."
    exit 1
esac
isects=(EMBREE HH HH2 SF01 MT BW12 BW9 SHEV DS EMBREE2)
for isect in ${isects[@]}; do
    sed_repl="s@#define ISECT_METHOD .*@#define ISECT_METHOD ISECT_${isect}@g"
    sed -i "${sed_repl}" $triangle_h
    cd $embree_build
    make > /dev/null
    time $single $bench_args
    time $ispc $bench_args
done
