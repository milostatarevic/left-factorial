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
  ulong h, j, segment, r_phase, rc_phase;
} remainder_thd_t;

void m_simple(mpz_t *r, ulong h, ulong i, ulong k)
{
  ulong v;

  if (h >= 2)
    return;

  if (h == 0)
  {
    if (k == 0)
      mpz_set_ui(*r, state.p_start + i);
    else if (k == 1)
      mpz_set_ui(*r, 1);
  }
  else if (h == 1 && 2 * i >= state.p_size - 1)
  {
    if (k == 0)
      mpz_set_ui(*r, state.p_start + 2 * i);
    else if (k == 1)
      mpz_set_ui(*r, 1);
  }
  else if (h == 1)
  {
    v = state.p_start + 2 * i;

    if (k == 0)
    {
      // no need to cover 2i == state.p_size
      mpz_set_ui(*r, v);
      mpz_mul_ui(*r, *r, v + 1);
    }
    else if (k == 1)
      mpz_set_ui(*r, v + 1);
  }
}

void *remainder_thread(void *thread_data)
{
  ulong h, j, s, segment, r_phase, rc_phase, len, len_1;
  mpz_t *rc_list_0, *rc_list_1, *n_list_0, *n_list_1,
        *p_list, *save_list_0, *save_list_1;
  mpz_t *_p, *_rc0, *_rc1;
  mpz_t a0, a1;
  remainder_thd_t *data;

  data = (remainder_thd_t *) thread_data;
  h = data->h;
  j = data->j;
  segment = data->segment;
  r_phase = data->r_phase;
  rc_phase = data->rc_phase;

  len = segment_size(h, segment);
  len_1 = segment_size(h + 1, segment);

  mpz_init(a0);
  mpz_init(a1);

  rc_list_0 = mpz_list_init(len_1);
  rc_list_1 = mpz_list_init(len_1);

  n_list_0 = mpz_list_init(len);
  n_list_1 = mpz_list_init(len);
  p_list = mpz_list_init(len);

  save_list_0 = (mpz_t *) malloc(len * sizeof(mpz_t));
  save_list_1 = (mpz_t *) malloc(len * sizeof(mpz_t));


  if (h + 1 == state.h_switch)
  {
    load_element(&rc_list_0[0], r_phase, h + 1, j / 2, 0);
    load_element(&rc_list_1[0], r_phase, h + 1, j / 2, 1);
  }
  else
  {
    load_segment(rc_list_0, r_phase, len_1, h + 1, segment, 0);
    load_segment(rc_list_1, r_phase, len_1, h + 1, segment, 1);
  }

  if (h >= 2)
  {
    load_segment(n_list_0, PHASE_N_TREE, len, h, segment, 0);
    load_segment(n_list_1, PHASE_N_TREE, len, h, segment, 1);
  }

  load_segment(p_list, PHASE_P_TREE, len, h, segment, 0);

  //

  for(s = 0; s < len && j < state.level_len[h]; j += 2)
  {
    _rc0 = &rc_list_0[s / 2];
    _rc1 = &rc_list_1[s / 2];

    if (j + 1 < state.level_len[h])
    {
      if (h >= 2)
      {
        mpz_mul(a0, *_rc0, n_list_0[s]);
        mpz_mul(a1, *_rc0, n_list_1[s]);
      }
      else
      {
        m_simple(&a0, h, j, 0);
        m_simple(&a1, h, j, 1);

        mpz_mul(a0, *_rc0, a0);
        mpz_mul(a1, *_rc0, a1);
      }

      mpz_add(a1, a1, *_rc1);
    }

    _p = &p_list[s];

    mpz_mod(*_rc0, *_rc0, *_p);
    mpz_mod(*_rc1, *_rc1, *_p);

    mpz_init_set(save_list_0[s], *_rc0);
    mpz_init_set(save_list_1[s], *_rc1);

    s ++;

    if (j + 1 >= state.level_len[h])
      break;

    _p = &p_list[s];

    mpz_mod(a0, a0, *_p);
    mpz_mod(a1, a1, *_p);

    mpz_init_set(save_list_0[s], a0);
    mpz_init_set(save_list_1[s], a1);

    s ++;
  }

  save_segment(save_list_0, rc_phase, s, h, segment, 0);
  save_segment(save_list_1, rc_phase, s, h, segment, 1);

  mpz_list_free(rc_list_0, len_1);
  mpz_list_free(rc_list_1, len_1);

  mpz_list_free(n_list_0, len);
  mpz_list_free(n_list_1, len);
  mpz_list_free(p_list, len);

  mpz_list_free(save_list_0, s);
  mpz_list_free(save_list_1, s);

  mpz_clear(a0);
  mpz_clear(a1);

  pthread_exit(0);
}

ulong remainders(mpz_t *r0, mpz_t *r1)
{
  slong h;
  ulong j, tmp, r_phase, rc_phase, segment;
  mpz_t p, n0, n1, rc0, rc1, a0, a1;

  pthread_attr_t attr;
  pthread_t *threads;
  remainder_thd_t *remainder_thd;
  int ts, t;
  void *status;

  r_phase = MAX_PHASE;
  rc_phase = r_phase + 1;

  save_element(r0, r_phase, state.height - 1, 0, 0, NULL);
  save_element(r1, r_phase, state.height - 1, 0, 1, NULL);

  mpz_init(rc0);
  mpz_init(rc1);
  mpz_init(n0);
  mpz_init(n1);
  mpz_init(a0);
  mpz_init(a1);
  mpz_init(p);

  state.use_ntt = state.use_ntt_mod = 1;
  for (h = state.height - 2; h >= state.h_switch; h --)
  {
    printf("h: %lu, use ntt\n", h);

    for(j = 0; j < state.level_len[h]; j += 2)
    {
      load_element(&rc0, r_phase, h + 1, j / 2, 0);
      load_element(&rc1, r_phase, h + 1, j / 2, 1);

      if (j + 1 < state.level_len[h])
      {
        load_element(&n0, PHASE_N_TREE, h, j, 0);
        load_element(&n1, PHASE_N_TREE, h, j, 1);

        m_mul(&a0, &rc0, &n0);
        m_mul(&a1, &rc0, &n1);
        mpz_add(a1, a1, rc1);
      }

      load_element(&p, PHASE_P_TREE, h, j, 0);

      fast_mod(&rc0, &p);
      fast_mod(&rc1, &p);

      save_element(&rc0, rc_phase, h, j, 0, NULL);
      save_element(&rc1, rc_phase, h, j, 1, NULL);

      if (j + 1 >= state.level_len[h])
        continue;

      load_element(&p, PHASE_P_TREE, h, j + 1, 0);

      fast_mod(&a0, &p);
      fast_mod(&a1, &p);

      save_element(&a0, rc_phase, h, j + 1, 0, NULL);
      save_element(&a1, rc_phase, h, j + 1, 1, NULL);
    }

    // switch phases
    tmp = r_phase;
    r_phase = rc_phase;
    rc_phase = tmp;

    // delete unused files
    exec("rm %s/c%lu-l%lu-*", folders.tree, rc_phase, h + 1);
  }

  mpz_clear(rc0);
  mpz_clear(rc1);
  mpz_clear(n0);
  mpz_clear(n1);
  mpz_clear(a0);
  mpz_clear(a1);
  mpz_clear(p);

  ////
  threads = (pthread_t *) malloc(state.num_threads * sizeof(pthread_t));
  remainder_thd = (remainder_thd_t *) malloc(state.num_threads * sizeof(remainder_thd_t));

  for (h = state.h_switch - 1; h >= 0; h --)
  {
    printf("h: %lu, len: %lu, segment size: %lu\n",
           h, state.level_len[h], state.segment_len[h]);

    j = 0;
    for(segment = 0; j < state.level_len[h]; )
    {
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

      for(t = 0; t < state.num_threads && j < state.level_len[h]; t++)
      {
        remainder_thd[t] = (remainder_thd_t) {
          .h = h, .j = j, .segment = segment, .r_phase = r_phase, .rc_phase = rc_phase,
        };

        ts = pthread_create(&threads[t], &attr, remainder_thread, (void *) &remainder_thd[t]);
        if (ts)
          error("thread creation error: %d", ts);

        segment ++;
        j += state.segment_len[h];
      }

      pthread_attr_destroy(&attr);

      for(t--; t >= 0; t--)
      {
         ts = pthread_join(threads[t], &status);
         if (ts)
           error("error on joining: %d", ts);
      }
    }

    // switch phases
    tmp = r_phase;
    r_phase = rc_phase;
    rc_phase = tmp;

    // delete unused files
    exec("rm %s/c%lu-l%lu-*", folders.tree, rc_phase, h + 1);
  }
  free(threads);
  free(remainder_thd);

  return r_phase;
}

void collect_remainders(ulong r_phase)
{
  ulong j, s, len, p, r_p, segment;
  mpz_t *p_list, *r_list;
  FILE *f_out;

  f_out = open_file("w", "out.remainders.%lu-%lu", state.p_start, state.p_end);

  j = 0;
  for (segment = 0; j < state.level_len[0]; segment ++)
  {
    len = segment_size(0, segment);

    p_list = mpz_list_init(len);
    r_list = mpz_list_init(len);

    load_segment(p_list, PHASE_P_TREE, len, 0, segment, 0);
    load_segment(r_list, r_phase, len, 0, segment, 1);

    for (s = 0; s < len && j < state.level_len[0]; s ++, j ++)
    {
      if (mpz_cmp_ui(p_list[s], 1))
      {
        p = mpz_get_ui(p_list[s]);
        r_p = mpz_get_ui(r_list[s]);

        r_p -= 1;
        if (r_p > p)
          r_p += p;

        fprintf(f_out, "%lu\t%lu\n", p, r_p);
      }
    }

    mpz_list_free(p_list, len);
    mpz_list_free(r_list, len);
  }

  fclose(f_out);
}
