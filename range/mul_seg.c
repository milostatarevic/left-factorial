/*
  This file is part of the left-factorial repository.
  Copyright (C) 2019 Milos Tatarevic

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, write to the Free Software Foundation,
  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "defs.h"

typedef struct {
  ulong h, i, j, l, phase;
  inverse_t *inverse;
} mul_thd_t;

ulong segment_size(ulong h, ulong segment)
{
  ulong segment_start;

  if (h >= state.h_switch)
    return 1;

  segment_start = segment * state.segment_len[h];
  if (segment_start + state.segment_len[h] <= state.level_len[h])
    return state.segment_len[h];

  return state.level_len[h] - segment_start;
}

void reduce_segments(ulong h, ulong segment)
{
  ulong i, pi, len;
  mpz_t *list0, *list1, *list_p;

  len = segment_size(h, segment);
  if (len < 2)
    return;

  list0 = (mpz_t *) malloc(len * sizeof(mpz_t));
  list1 = (mpz_t *) malloc(len * sizeof(mpz_t));
  list_p = (mpz_t *) malloc(len * sizeof(mpz_t));

  for(i = 0; i < len; i ++)
  {
    mpz_init(list0[i]);
    mpz_init(list1[i]);
    mpz_init(list_p[i]);
  }

  load_segment(list_p, PHASE_P_TREE, len, h, segment, 0);
  load_segment(list0, PHASE_N_TREE, len, h, segment, 0);
  load_segment(list1, PHASE_N_TREE, len, h, segment, 1);

  for(i = 0; i < len; i ++)
  {
    pi = i + 1;
    if (pi >= len) pi = i;

    if ((i & 1) == 0)
    {
      fast_mod(&list0[i], &list_p[pi]);
      fast_mod(&list1[i], &list_p[pi]);
    }
    else
    {
      mpz_set_ui(list0[i], 0);
      mpz_set_ui(list1[i], 0);
    }
  }

  save_segment(list0, PHASE_N_TREE, len, h, segment, 0);
  save_segment(list1, PHASE_N_TREE, len, h, segment, 1);

  for(i = 0; i < len; i ++)
  {
    mpz_clear(list0[i]);
    mpz_clear(list1[i]);
    mpz_clear(list_p[i]);
  }

  free(list0);
  free(list1);
  free(list_p);
}

void mul_values(mpz_t *dest, mpz_t *src_0, mpz_t *src_1,
                ulong src_end, ulong len, inverse_t *inverse)
{
  if (len <= src_end)
  {
    mpz_set(*dest, *src_0);
    return;
  }

  m_mul(dest, src_0, src_1);
  reduce_mod(dest, inverse);
}

void mul_values_double(mpz_t *dest0, mpz_t *dest1,
                       mpz_t *src0_0, mpz_t *src0_1,
                       mpz_t *src1_0, mpz_t *src1_1,
                       ulong src_end, ulong len,
                       inverse_t *inverse)
{
  if (len <= src_end)
  {
    mpz_set(*dest0, *src0_0);
    mpz_set(*dest1, *src0_1);
    return;
  }

  mpz_t t;

  mpz_init(t);
  m_mul(&t, src0_0, src1_1);
  mpz_add(*dest1, t, *src0_1);
  mpz_clear(t);

  m_mul(dest0, src0_0, src1_0);

  reduce_mod(dest0, inverse);
  reduce_mod(dest1, inverse);
}

ulong list_mul(mpz_t *r, mpz_t *list, ulong len, ulong phase, ulong segment)
{
  ulong h, i, j, d;
  mpz_t tmp;

  h = 0;
  if (len == 1)
  {
    mpz_set(*r, list[0]);
  }
  else
  {
    for (i = len; i > 1; i = (i + 1) / 2, h ++)
    {
      d = (i + 1) / 2;
      for(j = 0; j < d; j ++)
      {
        mpz_init_set_ui(tmp, 1);
        mul_values(&tmp, &list[2 * j], &list[2 * j + 1], 2 * j + 1, i, NULL);
        mpz_set(list[j], tmp);
        mpz_clear(tmp);
      }

      for(j = d; j < i; j ++)
        mpz_clear(list[j]);

      if (d >= 2)
        save_segment(list, phase, d, h + 1, segment, 0);
    }
    mpz_set(*r, list[0]);
  }

  mpz_clear(list[0]);
  free(list);

  return h;
}

ulong list_mul_double(mpz_t *r0, mpz_t *r1,
                      mpz_t *list0, mpz_t *list1, ulong len,
                      ulong phase, ulong segment)
{
  ulong h, i, j, d;
  mpz_t tmp0, tmp1;

  h = 0;
  if (len == 1)
  {
    mpz_set(*r0, list0[0]);
    mpz_set_ui(*r1, 1);
  }
  else
  {
    for (i = len; i > 1; i = (i + 1) / 2, h ++)
    {
      d = (i + 1) / 2;
      for(j = 0; j < d; j ++)
      {
        mpz_init_set_ui(tmp0, 1);
        mpz_init_set_ui(tmp1, 1);

        mul_values_double(&tmp0, &tmp1,
                          &list0[2 * j], &list1[2 * j],
                          &list0[2 * j + 1], &list1[2 * j + 1],
                          2 * j + 1, i, NULL);

        mpz_set(list0[j], tmp0);
        mpz_set(list1[j], tmp1);
        mpz_clear(tmp0);
        mpz_clear(tmp1);
      }

      for(j = d; j < i; j ++)
      {
        mpz_clear(list0[j]);
        mpz_clear(list1[j]);
      }

      if (phase == PHASE_N_TREE)
      {
        if (d >= 2 && h >= 1)
        {
          save_segment(list0, phase, d, h + 1, segment, 0);
          save_segment(list1, phase, d, h + 1, segment, 1);

          if (h >= 2)
            reduce_segments(h, segment);
        }
      }
    }
    mpz_set(*r0, list0[0]);
    mpz_set(*r1, list1[0]);

    if (phase == PHASE_N_TREE && h - 1 >= 2)
      reduce_segments(h - 1, segment);
  }

  mpz_clear(list0[0]);
  mpz_clear(list1[0]);

  free(list0);
  free(list1);

  return h;
}

void list_to_double(mpz_t **list0, mpz_t **list1, ulong n, ulong len)
{
  ulong i, v;

  *list0 = (mpz_t *) malloc(len * sizeof(mpz_t));
  *list1 = (mpz_t *) malloc(len * sizeof(mpz_t));

  for(i = 0, v = n; i < len; i ++, v ++)
  {
    mpz_init_set_ui((*list0)[i], v);
    mpz_init_set_ui((*list1)[i], 1);
  }
}

void *mul_thread(void *thread_data)
{
  ulong ch, h, i, j, l, phase, g;
  mpz_t *values0, *values1, r0, r1;
  inverse_t *inverse;
  mul_thd_t *data;

  data = (mul_thd_t *) thread_data;
  h = data->h;
  i = data->i;
  j = data->j;
  l = data->l;
  phase = data->phase;
  inverse = data->inverse;

  mpz_init(r0);
  mpz_init(r1);

  if (phase == PHASE_P_TREE)
  {
    values0 = (mpz_t *) malloc(l * sizeof(mpz_t));

    for(g = 0; g < l; g ++)
      mpz_init_set_ui(values0[g], 1);
    load_segment(values0, phase, l, 0, j, 0);

    ch = list_mul(&r0, values0, l, phase, j);
  }
  else
  {
    list_to_double(&values0, &values1, i + 1, l);
    ch = list_mul_double(&r0, &r1, values0, values1, l, phase, j);
  }

  if (phase == PHASE_P_TREE || phase == PHASE_N_TREE)
  {
    if (ch < 2 && phase == PHASE_N_TREE) ch = 2;
    if (ch == 0 && phase == PHASE_P_TREE) ch = 1;

    for(; ch < h; ch ++)
    {
      save_segment(&r0, phase, 1, ch, j, 0);
      if (phase == PHASE_N_TREE)
        save_segment(&r1, phase, 1, ch, j, 1);
    }
  }

  save_element(&r0, phase, h, j, 0, inverse);
  if (phase != PHASE_P_TREE)
    save_element(&r1, phase, h, j, 1, inverse);

  mpz_clear(r0);
  mpz_clear(r1);

  pthread_exit(NULL);
}

ulong mul_all_segments(ulong phase, ulong start, ulong end, ulong step, ulong h, inverse_t *inverse)
{
  int ts, t;
  ulong i, j, l;
  void *status;
  pthread_t *threads;
  pthread_attr_t attr;
  mul_thd_t *mul_thd;

  threads = (pthread_t *) malloc(state.num_threads * sizeof(pthread_t));
  mul_thd = (mul_thd_t *) malloc(state.num_threads * sizeof(mul_thd_t));

  for (j = 0, i = start; i < end;)
  {
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(t = 0; t < state.num_threads && i < end; t++)
    {
      printf("thread %lu - %d\n", j, t);

      l = (i + step >= end) ? (end - i) : step;

      mul_thd[t] = (mul_thd_t) {
        .h = h, .i = i, .j = j, .l = l, .phase = phase, .inverse = inverse
      };
      ts = pthread_create(&threads[t], &attr, mul_thread, (void *) &mul_thd[t]);
      if (ts)
        error("thread creation error: %d", ts);

      i += step;
      j ++;
    }

    pthread_attr_destroy(&attr);

    for(t--; t >= 0; t--)
    {
     ts = pthread_join(threads[t], &status);
     if (ts)
       error("error on joining: %d", ts);
    }
  }
  free(threads);
  free(mul_thd);

  return j;
}
