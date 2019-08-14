// ======================================================================== //
// Copyright 2009-2018 Intel Corporation                                    //
// Copyright 2019 Björn Lindqvist <bjourne@gmail.com>                       //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //
//
// This file contains Embree's default medium build quality triangle
// intersection algorithm.
#pragma once

template<int M, int K>
static __forceinline bool
intersectKRaysMTris(Vec3vf<K> o, Vec3vf<K> d,
                    vfloat<K> tn, vfloat<K> tf,
                    size_t i,
                    const vbool<K>& valid0,
                    MTHit<K>& hit,
                    const TriangleM<M>& tri)
{
  // Do broadcasting and cross prod.
  Vec3<vfloat<K>> v0 = broadcast<vfloat<K>>(tri.v0, i);
  Vec3<vfloat<K>> e1 = broadcast<vfloat<K>>(tri.e1, i);
  Vec3<vfloat<K>> e2 = broadcast<vfloat<K>>(tri.e2, i);
  Vec3<vfloat<K>> ng = cross(e2, e1);

  /* calculate denominator */
  const Vec3vf<K> c = v0 - o;
  const Vec3vf<K> r = cross(c, d);
  const vfloat<K> den = dot(ng, d);
  const vfloat<K> absDen = abs(den);
  const vfloat<K> sgnDen = signmsk(den);

  /* test against edge e2 v0 */
  const vfloat<K> u = dot(e2, r) ^ sgnDen;
  vbool<K> valid = valid0 & (u >= 0.0f);
  if (likely(none(valid))) {
    return false;
  }

  /* test against edge v0 e1 */
  const vfloat<K> v = dot(e1, r) ^ sgnDen;
  valid &= v >= 0.0f;
  if (likely(none(valid))) {
    return false;
  }

  /* test against edge e1 e2 */
  const vfloat<K> w = absDen-u-v;
  valid &= w >= 0.0f;
  if (likely(none(valid))) {
    return false;
  }

  /* perform depth test */
  const vfloat<K> t = dot(ng, c) ^ sgnDen;
  valid &= (absDen * tn < t) & (t <= absDen * tf);
  if (unlikely(none(valid))) {
    return false;
  }

  valid &= den != vfloat<K>(zero);
  if (unlikely(none(valid))) {
    return false;
  }

  /* calculate hit information */
  const vfloat<K> rcpAbsDen = rcp(absDen);
  new (&hit) MTHit<K>(valid,
                      u * rcpAbsDen,
                      v * rcpAbsDen,
                      t * rcpAbsDen,
                      ng);
  return true;
}

/*! Intersect 1 ray with one of M triangles. */
template<int M>
static __forceinline bool
intersect1RayMTris(const Vec3vf<M>& o, const Vec3vf<M>& d,
                   const vfloat<M>& tn, const vfloat<M>& tf,
                   const TriangleM<M>& tri, MTHit<M>& hit)
{
  const Vec3vf<M> c = tri.v0 - o;
  const Vec3vf<M> r = cross(c,d);
  const Vec3vf<M> ng = cross(tri.e2, tri.e1);
  const vfloat<M> den = dot(ng, d);
  const vfloat<M> absDen = abs(den);
  const vfloat<M> sgnDen = signmsk(den);

  /* perform edge tests */
  const vfloat<M> u = dot(r, tri.e2) ^ sgnDen;
  const vfloat<M> v = dot(r, tri.e1) ^ sgnDen;

  vbool<M> valid = (den != vfloat<M>(zero)) &
    (u >= 0.0f) & (v >= 0.0f) & (u + v <= absDen);
  if (likely(none(valid)))
    return false;

  /* perform depth test */
  const vfloat<M> t = dot(ng, c) ^ sgnDen;
  valid &= ((absDen * tn) < t) & (t <= absDen * tf);
  if (likely(none(valid)))
    return false;

  /* update hit information */
  const vfloat<M> rcpAbsDen = rcp(absDen);
  new (&hit) MTHit<M>(valid,
                      u * rcpAbsDen,
                      v * rcpAbsDen,
                      t * rcpAbsDen,
                      ng);
  return true;
}
