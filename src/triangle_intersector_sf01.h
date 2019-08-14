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
  Vec3vf<K> v1 = v0 - e1;
  Vec3vf<K> v2 = e2 + v0;

  Vec3vf<K> v0o = v0 - o;
  Vec3vf<K> v1o = v1 - o;
  Vec3vf<K> v2o = v2 - o;

  vfloat<K> w2 = dot(d, cross(v1o, v0o));
  vfloat<K> w0 = dot(d, cross(v2o, v1o));
  vbool<K> s2 = w2 >= vfloat<K>(zero);
  vbool<K> s0 = w0 >= vfloat<K>(zero);

  vfloat<K> w1 = dot(d, cross(v0o, v2o));
  vbool<K> s1 = w1 >= vfloat<K>(zero);
  vbool<K> valid = valid0 & (s2 == s1) & (s0 == s2);

  Vec3vf<K> n = cross(e2, e1);
  vfloat<K> den = dot(n, d);
  vfloat<K> t = dot(n, v0o) * rcp(den);
  valid &= (tn < t) & (t <= tf);
  if (likely(none(valid)))
    return false;

  vfloat<K> rcpSum = rcp(w0 + w1 + w2);
  vfloat<K> u = w1 * rcpSum;
  vfloat<K> v = w2 * rcpSum;

  new (&hit) MTHit<K>(valid, u, v, t, n);
  return true;
}

template<int M, int K>
static __forceinline bool
intersectKRaysMTris(Vec3vf<K> o, Vec3vf<K> d,
                    vfloat<K> tn, vfloat<K> tf,
                    size_t i,
                    const vbool<K>& valid0,
                    MTHit<K>& hit,
                    const TriangleM<M>& tri)
{
  Vec3vf<K> v0 = broadcast<vfloat<K>>(tri.v0, i);
  Vec3vf<K> e1 = broadcast<vfloat<K>>(tri.e1, i);
  Vec3vf<K> e2 = broadcast<vfloat<K>>(tri.e2, i);
  return isectAlgo(o, d, tn, tf,
                   v0, e1, e2,
                   hit, valid0);
}

template<int M>
static __forceinline bool
intersect1RayMTris(const Vec3vf<M>& o, Vec3vf<M> d,
                   vfloat<M> tn, vfloat<M> tf,
                   const TriangleM<M>& tri, MTHit<M>& hit)
{
  return isectAlgo<M>(o, d, tn, tf,
                      tri.v0, tri.e1, tri.e2,
                      hit, true);
}
