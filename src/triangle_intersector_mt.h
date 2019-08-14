// Copyright (C) 2019 Björn Lindqvist <bjourne@gmail.com>
#pragma once

template<int K>
static __forceinline bool
isectAlgo(const Vec3vf<K>& o, const Vec3vf<K>& d,
          const vfloat<K>& tn, const vfloat<K>& tf,
          const Vec3vf<K>& v0,
          const Vec3vf<K>& e1,
          const Vec3vf<K>& e2,
          MTHit<K>& hit, const vbool<K>& valid0) {
  Vec3vf<K> pvec = cross(d, e2);
  vfloat<K> det = dot(e1, pvec);

  vfloat<K> rdet = rcp(det);
  Vec3vf<K> tvec = o - v0;

  vfloat<K> u = dot(tvec, pvec) * rdet;
  Vec3vf<K> qvec = cross(tvec, e1);
  vfloat<K> v = dot(d, qvec) * rdet;

  vfloat<K> t = dot(e2, qvec) * rdet;
  const vbool<K> valid = (u >= 0.0f) & (v >= 0.0f) & (u + v <= 1.0f)
    & (tn < t) & (t <= tf);
  if (likely(none(valid)))
    return false;

  new (&hit) MTHit<K>(valid, u, v, t, cross(e2, e1));
  return true;
}

template<int M, int K>
static __forceinline bool
intersectKRaysMTris(const Vec3vf<K>& o,
                    const Vec3vf<K>& d,
                    const vfloat<K>& tn, const vfloat<K>& tf,
                    size_t i,
                    const vbool<K>& valid0,
                    MTHit<K>& hit,
                    const TriangleM<M>& tri)
{
  Vec3vf<K> v0 = broadcast<vfloat<K>>(tri.v0, i);
  Vec3vf<K> e1 = broadcast<vfloat<K>>(tri.e1, i);
  Vec3vf<K> e2 = broadcast<vfloat<K>>(tri.e2, i);
  return isectAlgo(o, d, tn, tf, v0, e1, e2, hit, valid0);
}

template<int M>
static __forceinline bool
intersect1RayMTris(const Vec3vf<M>& o, const Vec3vf<M>& d,
                   const vfloat<M>& tn, const vfloat<M>& tf,
                   const TriangleM<M>& tri, MTHit<M>& hit)
{
  return isectAlgo<M>(o, d,
                      tn, tf,
                      tri.v0, tri.e1, tri.e2,
                      hit, true);
}
