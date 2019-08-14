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

  Vec3vf<K> n = cross(e2, e1);
  vfloat<K> a = dot(n, v0 - o);
  vfloat<K> b = dot(n, d);
  vfloat<K> t = a * rcp(b);
  vbool<K> valid = valid0 & (tn < t) & (t <= tf);
  if (likely(none(valid))) {
    return false;
  }

  vfloat<K> u_u = dot(e1, e1);
  vfloat<K> u_v = dot(e1, e2);
  vfloat<K> v_v = dot(e2, e2);
  Vec3vf<K> w = o + t * d - v0;
  vfloat<K> w_u = dot(w, e1);
  vfloat<K> w_v = dot(w, e2);
  vfloat<K> det = u_v * u_v - u_u * v_v;
  vfloat<K> rdet = rcp(det);

  vfloat<K> u = (u_v * w_v - v_v * w_u) * rdet;
  vfloat<K> v = (u_v * w_u - u_u * w_v) * rdet;
  valid &= (u >= 0.0f) & (v >= 0.0f) & (u + v <= 1.0f);
  if (likely(none(valid)))
    return false;

  new (&hit) MTHit<K>(valid, u, v, t, n);
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
