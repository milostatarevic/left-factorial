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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <time.h>

#include <gmp.h>
#include <flint.h>
#include <ulong_extras.h>

#include "ntt/ntt-internal.h"
#include "gmp-ntt/gmp-ntt.h"

#define MAX_HEIGHT       64
#define MUL_STEP_DEPTH   3
#define STEP_BITS        24
#define STEP_SIZE        (1 << STEP_BITS)

enum {
  PHASE_P_TREE,
  PHASE_MUL_DIFFERENCE,
  PHASE_MUL,
  PHASE_N_TREE,
  MAX_PHASE,
};

typedef struct {
  ulong use_ntt, use_ntt_mod, num_threads,
        p_start, p_end, p_saved_at, p_size, height, h_switch,
        level_len[MAX_HEIGHT], segment_len[MAX_HEIGHT], last_h[MAX_PHASE];
} state_t;

typedef struct {
  char *tree, *saved, *backup, *inverse;
} folders_t;

typedef struct {
  ulong bits, h_end, blank, value_sum, inverse_sum;
  mpz_t value, inverse;
} inverse_t;

// globals
extern state_t state;
extern folders_t folders;

// inverse.c
void init_inverse(inverse_t *, ulong);
void reduce_mod(mpz_t *, inverse_t *);

// mpz.c
mpz_t *mpz_list_init(ulong);
void  mpz_list_free(mpz_t *, ulong);
void  mpz_write(FILE *, mpz_t *);
void  mpz_read(FILE *, mpz_t *);
void  fast_mod(mpz_t *, mpz_t *);
void  m_mul(mpz_t *, mpz_t *, mpz_t *);

// mul_el.c
ulong mul_all_elements(ulong, ulong, ulong, inverse_t *);

// mul_seg.c
ulong segment_size(ulong, ulong);
ulong mul_all_segments(ulong, ulong, ulong, ulong, ulong, inverse_t *);

// primes.c
void generate_primes(void);

// remainders.c
ulong remainders(mpz_t *, mpz_t *);
void collect_remainders(ulong);

// store.c
void load_segment(mpz_t *, ulong, ulong, ulong, ulong, ulong);
void save_segment(mpz_t *, ulong, ulong, ulong, ulong, ulong);
void load_element(mpz_t *, ulong, ulong, ulong, ulong);
void save_element(mpz_t *, ulong, ulong, ulong, ulong, inverse_t *);
slong last_stored_index(slong, int);
void move_to_saved(ulong, int);

// util.c
void exec(const char*, ...);
void error(const char*, ...);
FILE *open_file(const char *, const char*, ...);
void start_time(struct timespec *);
void print_time(struct timespec *, struct timespec *);
