#pragma once

#include "defs.h"
#include "floatn.h"

typedef struct {
  union {
    float3 p;
    struct {
      float x, y, z;
    };
  };
} obj_v_t;

typedef struct {
  union {
    uint32_t is[3];
    struct {
      uint32_t i0, i1, i2;
    };
  };
  obj_v_t vs[3];
  float area;
  float3 n;
  float3 c;
} obj_f_t;

typedef struct {
  char file[256];
  uint32_t num_vs;
  uint32_t num_fs;
  obj_v_t* vs;
  obj_f_t* fs;
} obj_t;

void obj_term(obj_t* obj);
int obj_load(obj_t* obj, const char objfile[], int term);
int obj_write(const obj_t* obj, const char objfile[]);
