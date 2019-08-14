// Copyright (C) 2019 Björn Lindqvist <bjourne@gmail.com>
#pragma once

template<int K>
static __forceinline bool
isect(const Vec3vf<K>& o, const Vec3vf<K>& d,
      const vfloat<K>& tn, const vfloat<K>& tf,
      const Vec3vf<K>& n0, const Vec3vf<K>& n1, const Vec3vf<K>& n2,
      const vfloat<K>& d0, const vfloat<K>& d1, const vfloat<K>& d2,
      MTHit<K>& hit, const vbool<K>& valid0)
{
  vfloat<K> t_o = dot(o, n2) + d2;
  vfloat<K> t_d = dot(d, n2);
  vfloat<K> t = -t_o * rcp(t_d);
  Vec3vf<K> wr = o + d * t;
  vfloat<K> u = dot(wr, n0) + d0;
  vfloat<K> v = dot(wr, n1) + d1;
  const vbool<K> valid = valid0
    & (u >= 0.0f) & (v >= 0.0f) & (u + v <= 1.0f)
    & (tn < t) & (t <= tf);
  if (likely(none(valid)))
    return false;
  new (&hit) MTHit<K>(valid, u, v, t, n2);
  return true;
}

template<int M, int K>
static __forceinline bool
intersectKRaysMTris(const Vec3vf<K>& o, const Vec3vf<K>& d,
                    const vfloat<K>& tn, const vfloat<K>& tf,
                    size_t i,
                    const vbool<K>& valid0,
                    MTHit<K>& hit,
                    const TriangleM<M>& tri)
{
  Vec3vf<K> n0 = broadcast<vfloat<K>>(tri.n0, i);
  Vec3vf<K> n1 = broadcast<vfloat<K>>(tri.n1, i);
  Vec3vf<K> n2 = broadcast<vfloat<K>>(tri.n2, i);
  vfloat<K> d0 = vfloat<K>(tri.d0[i]);
  vfloat<K> d1 = vfloat<K>(tri.d1[i]);
  vfloat<K> d2 = vfloat<K>(tri.d2[i]);
  return isect<K>(o, d, tn, tf,
                  n0, n1, n2, d0, d1, d2,
                  hit, valid0);
}

template<int M>
static __forceinline bool
intersect1RayMTris(const Vec3vf<M>& o, const Vec3vf<M>& d,
                   const vfloat<M>& tn, const vfloat<M>& tf,
                   const TriangleM<M>& tri, MTHit<M>& hit)
{
  return isect<M>(o, d, tn, tf,
                  tri.n0, tri.n1, tri.n2,
                  tri.d0, tri.d1, tri.d2,
                  hit, true);
}
