// Copyright (C) 2019 Björn Lindqvist <bjourne@gmail.com>
#pragma once

template<int K>
static __forceinline bool
isect(const Vec3vf<K>& o,
      const Vec3vf<K>& d,
      const vfloat<K>& tn, const vfloat<K>& tf,
      const Vec3vf<K>& n0, const Vec3vf<K>& n1, const Vec3vf<K>& n2,
      const vfloat<K>& d0, const vfloat<K>& d1, const vfloat<K>& d2,
      MTHit<K>& hit, const vbool<K>& valid0)
{
  vfloat<K> det = dot(n0, d);
  vbool<K> valid = det != vfloat<K>(zero);
  if (unlikely(none(valid)))
    return false;
  vfloat<K> t = d0 - dot(o, n0);
  Vec3vf<K> wr = o * det + d * t;
  vfloat<K> u = dot(wr, n1) + det * d1;
  vfloat<K> v = dot(wr, n2) + det * d2;
  vfloat<K> p0 = ((det - u - v) ^ u) | (u ^ v)
    | ((tf * det - t) ^ t);
  valid &= asInt(signmsk(p0)) == vint<K>(zero);
  if (likely(none(valid)))
    return false;
  vfloat<K> rdet = rcp(det);
  new (&hit) MTHit<K>(valid, u * rdet, v * rdet, t * rdet, n0);
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
  return isect<K>(o, d,
                  tn, tf,
                  n0, n1, n2,
                  d0, d1, d2,
                  hit, valid0);
}

/*! Intersect 1 ray with one of M triangles. */
template<int M>
static __forceinline bool
intersect1RayMTris(Vec3vf<M> o, Vec3vf<M> d,
                   vfloat<M> tn, vfloat<M> tf,
                   const TriangleM<M>& tri, MTHit<M>& hit)
{
  return isect<M>(o, d, tn, tf,
                  tri.n0, tri.n1, tri.n2,
                  tri.d0, tri.d1, tri.d2,
                  hit, true);
}
