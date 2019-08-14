// Copyright (C) 2019 Björn Lindqvist <bjourne@gmail.com>
#pragma once

typedef union {
    unsigned int u;
    float f;
} uint_or_float;

template<int M>
static __forceinline bool
isect1(float ox, float oy, float oz,
       float dx, float dy, float dz,
       float tn, float tf,
       size_t i,
       const TriangleM<M>& tri,
       float *u, float *v, float *t,
       float *nx, float *ny, float *nz) {
  float np = tri.np[i];
  float nu = tri.nu[i];
  float nv = tri.nv[i];
  float pu = tri.pu[i];
  float pv = tri.pv[i];
  float e1u = tri.e1u[i];
  float e1v = tri.e1v[i];
  float e2u = tri.e2u[i];
  float e2v = tri.e2v[i];
  int ci = tri.ci[i];
  float dett, det, du, dv;
  if (ci == 0) {
    dett = np - (oy * nu + oz * nv + ox);
    det = dy * nu + dz * nv + dx;
    du = dy * dett - (pu - oy) * det;
    dv = dz * dett - (pv - oz) * det;
  } else if (ci == 1) {
    dett = np - (ox * nu + oz * nv + oy);
    det = dx * nu + dz * nv + dy;
    du = dx * dett - (pu - ox) * det;
    dv = dz * dett - (pv - oz) * det;
  } else {
    dett = np - (ox * nu + oy * nv + oz);
    det = dx * nu + dy * nv + dz;
    du = dx * dett - (pu - ox) * det;
    dv = dy * dett - (pv - oy) * det;
  }
  *u = e2v * du - e2u * dv;
  *v = e1u * dv - e1v * du;
  uint_or_float pdetu = { .f = *u };
  uint_or_float pdetv = { .f = *v };
  uint_or_float pdet0 = { .f = det - *u - *v };
  if (((pdet0.u ^ pdetu.u) | (pdetu.u ^ pdetv.u)) & 0x80000000) {
    return false;
  }
  float rdet = 1 / det;
  *t = dett * rdet;
  if (*t <= tn || *t > tf)
    return false;
  *u = *u * rdet;
  *v = *v * rdet;
  if (ci == 0)
    *nx = 1, *ny = nu, *nz = nv;
  else if (ci == 1)
    *nx = nu, *ny = 1, *nz = nv;
  else
    *nx = nu, *ny = nv, *nz = 1;
  return true;
}

template<int K, int M>
static __forceinline bool
isectK(const Vec3vf<K>& o,
       const Vec3vf<K>& d,
       const vfloat<K>& tn, const vfloat<K>& tf,
       size_t i,
       const vbool<K>& valid0,
       MTHit<K>& hit,
       const TriangleM<M>& tri)
{
  vfloat<K> np = vfloat<K>(tri.np[i]);
  vfloat<K> nu = vfloat<K>(tri.nu[i]);
  vfloat<K> nv = vfloat<K>(tri.nv[i]);
  vfloat<K> pu = vfloat<K>(tri.pu[i]);
  vfloat<K> pv = vfloat<K>(tri.pv[i]);
  vfloat<K> e1u = vfloat<K>(tri.e1u[i]);
  vfloat<K> e1v = vfloat<K>(tri.e1v[i]);
  vfloat<K> e2u = vfloat<K>(tri.e2u[i]);
  vfloat<K> e2v = vfloat<K>(tri.e2v[i]);
  vfloat<K> dett, det, du, dv;
  int ci = tri.ci[i];
  if (ci == 0) {
    dett = np - (o.y * nu + o.z * nv + o.x);
    det = d.y * nu + d.z * nv + d.x;
    du = d.y * dett - (pu - o.y) * det;
    dv = d.z * dett - (pv - o.z) * det;
  } else if (ci == 1) {
    dett = np - (o.x * nu + o.z * nv + o.y);
    det = d.x * nu + d.z * nv + d.y;
    du = d.x * dett - (pu - o.x) * det;
    dv = d.z * dett - (pv - o.z) * det;
  } else {
    dett = np - (o.x * nu + o.y * nv + o.z);
    det = d.x * nu + d.y * nv + d.z;
    du = d.x * dett - (pu - o.x) * det;
    dv = d.y * dett - (pv - o.y) * det;
  }
  vfloat<K> u = e2v * du - e2u * dv;
  vfloat<K> v = e1u * dv - e1v * du;
  vfloat<K> p0 = ((det - u - v) ^ u) | (u ^ v);
  vbool<K> valid = valid0 & asInt(signmsk(p0)) == vint<K>(zero);
  if (likely(none(valid)))
    return false;
  vfloat<K> rdet = rcp(det);
  vfloat<K> t = dett * rdet;
  valid &= (tn < t) & (t <= tf);
  if (likely(none(valid)))
    return false;
  vfloat<K> x, y, z;
  if (ci == 0)
    x = 1, y = nu, z = nv;
  else if (ci == 1)
    x = nu, y = 1, z = nv;
  else
    x = nu, y = nv, z = 1;
  new (&hit) MTHit<K>(valid, u * rdet, v * rdet, t,
                      Vec3vf<K>(x, y, z));
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
  return isectK(o, d, tn, tf, i, valid0, hit, tri);
}

template<int M>
static __forceinline bool
intersect1RayMTris(const Vec3vf<M>& o, const Vec3vf<M>& d,
                   const vfloat<M>& tn, const vfloat<M>& tf,
                   const TriangleM<M>& tri, MTHit<M>& hit)
{
  float u, v, t, nx, ny, nz;
  float ox = o.x[0], oy = o.y[0], oz = o.z[0];
  float dx = d.x[0], dy = d.y[0], dz = d.z[0];
  float tn0 = tn[0], tf0 = tf[0];
  bool any_hit = false;
  for (size_t i = 0; i < M; i++) {
    if (tri.geomIDs[i] == -1) {
      break;
    }
    if (isect1<M>(ox, oy, oz,
                  dx, dy, dz,
                  tn0, tf0,
                  i,
                  tri,
                  &u, &v, &t, &nx, &ny, &nz)) {
      if (!any_hit) {
        new (&hit) MTHit<M>();
        any_hit = true;
      }
      set(hit.valid, i);
      hit.u[i] = u;
      hit.v[i] = v;
      hit.t[i] = t;
      hit.ng.x[i] = nx;
      hit.ng.y[i] = ny;
      hit.ng.z[i] = nz;
    }
  }
  return any_hit;
}
