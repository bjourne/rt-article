// Copyright (C) 2019 Björn Lindqvist <bjourne@gmail.com>
#pragma once

template<int K, int M>
static __forceinline bool
isectK(const Vec3vf<K>& o, const Vec3vf<K>& d,
      const vfloat<K>& tn, const vfloat<K>& tf,
      size_t i,
      const vbool<K>& valid0,
      MTHit<K>& hit,
      const TriangleM<M>& tri) {

  vfloat<K> t0 = vfloat<K>(tri.T[i][0]);
  vfloat<K> t1 = vfloat<K>(tri.T[i][1]);
  vfloat<K> t2 = vfloat<K>(tri.T[i][2]);
  vfloat<K> t3 = vfloat<K>(tri.T[i][3]);
  vfloat<K> t4 = vfloat<K>(tri.T[i][4]);
  vfloat<K> t5 = vfloat<K>(tri.T[i][5]);
  vfloat<K> t6 = vfloat<K>(tri.T[i][6]);
  vfloat<K> t7 = vfloat<K>(tri.T[i][7]);
  vfloat<K> t8 = vfloat<K>(tri.T[i][8]);
  int ci = tri.ci[i];
  vfloat<K> u, v, t;
  if (ci == 0) {
    vfloat<K> t_o = o.x + t6 * o.y + t7 * o.z + t8;
    vfloat<K> t_d = d.x + t6 * d.y + t7 * d.z;
    t = -t_o * rcp(t_d);
    Vec3vf<K> wr = o + d * t;
    u = t0 * wr.y + t1 * wr.z + t2;
    v = t3 * wr.y + t4 * wr.z + t5;
  } else if (ci == 1) {
    vfloat<K> t_o = t6 * o.x + o.y + t7 * o.z + t8;
    vfloat<K> t_d = t6 * d.x + d.y + t7 * d.z;
    t = -t_o * rcp(t_d);
    Vec3vf<K> wr = o + d * t;
    u = t0 * wr.x + t1 * wr.z + t2;
    v = t3 * wr.x + t4 * wr.z + t5;
  } else {
    vfloat<K> t_o = t6 * o.x + t7 * o.y + o.z + t8;
    vfloat<K> t_d = t6 * d.x + t7 * d.y + d.z;
    t = -t_o * rcp(t_d);
    Vec3vf<K> wr = o + d * t;
    u = t0 * wr.x + t1 * wr.y + t2;
    v = t3 * wr.x + t4 * wr.y + t5;
  }
  const vbool<K> valid = valid0
    & (u >= 0.0f) & (v >= 0.0f) & (u + v <= 1.0f)
    & (tn < t) & (t <= tf);
  if (likely(none(valid)))
    return false;
  vfloat<K> x, y, z;
  if (ci == 0)
    x = 1, y = t6, z = t7;
  else if (ci == 1)
    x = t6, y = 1, z = t7;
  else
    x = t6, y = t7, z = 1;
  new (&hit) MTHit<K>(valid, u, v, t, Vec3vf<K>(x, y, z));
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
intersect1RayMTris(Vec3vf<M> o, Vec3vf<M> d,
                   vfloat<M> tn, vfloat<M> tf,
                   const TriangleM<M>& tri, MTHit<M>& hit)
{
  float ox = o.x[0];
  float oy = o.y[0];
  float oz = o.z[0];
  float dx = d.x[0];
  float dy = d.y[0];
  float dz = d.z[0];
  float tn0 = tn[0];
  float tf0 = tf[0];
  bool any_hit = false;
  for (size_t i = 0; i < M; i++) {
    // Test if necessary
    if (tri.geomIDs[i] == -1)
      return any_hit;

    float t0 = tri.T[i][0];
    float t1 = tri.T[i][1];
    float t2 = tri.T[i][2];
    float t3 = tri.T[i][3];
    float t4 = tri.T[i][4];
    float t5 = tri.T[i][5];
    float t6 = tri.T[i][6];
    float t7 = tri.T[i][7];
    float t8 = tri.T[i][8];
    int ci = tri.ci[i];
    float t, to, td, wrx, wry, wrz, u, v;
    if (ci == 0) {
      to = ox + t6 * oy + t7 * oz + t8;
      td = dx + t6 * dy + t7 * dz;
      t = -to / td;
      wry = oy + dy * t;
      wrz = oz + dz * t;
      u = t0 * wry + t1 * wrz + t2;
      v = t3 * wry + t4 * wrz + t5;
    } else if (ci == 1) {
      to = t6 * ox + oy + t7 * oz + t8;
      td = t6 * dx + dy + t7 * dz;
      t = -to / td;
      wrx = ox + dx * t;
      wrz = oz + dz * t;
      u = t0 * wrx + t1 * wrz + t2;
      v = t3 * wrx + t4 * wrz + t5;
    } else {
      to = t6 * ox + t7 * oy + oz + t8;
      td = t6 * dx + t7 * dy + dz;
      t = -to / td;
      wrx = ox + dx * t;
      wry = oy + dy * t;
      u = t0 * wrx + t1 * wry + t2;
      v = t3 * wrx + t4 * wry + t5;
    }
    if (u >= 0 && v >= 0 && u + v <= 1 && tn0 < t && t <= tf0) {
      if (!any_hit) {
        new (&hit) MTHit<M>();
        any_hit = true;
      }
      float nx, ny, nz;
      if (ci == 0)
        nx = 1, ny = t6, nz = t7;
      else if (ci == 1)
        nx = t6, ny = 1, nz = t7;
      else
        nx = t6, ny = t7, nz = 1;
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
