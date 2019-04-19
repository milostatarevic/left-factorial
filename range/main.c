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

state_t state;
folders_t folders;

void init_state(ulong p_start, ulong p_end, ulong p_saved_at, int num_threads)
{
  ulong i, h, p_size;

  state.use_ntt = state.use_ntt_mod = 0;

  state.p_start = p_start;
  state.p_end = p_end;
  state.p_saved_at = p_saved_at;
  state.num_threads = num_threads;

  p_size = p_end - p_start + 1;
  state.p_size = p_size;

  state.height = FLINT_CLOG2(p_size) + 1;
  for(i = 0, h = p_size; i < state.height; i ++, h = (h + 1) / 2)
    state.level_len[i] = h;

  state.h_switch = STEP_BITS;
  if (state.h_switch > state.height || state.h_switch == 0)
    state.h_switch = 1;

  printf("height: %lu, h-switch: %lu\n", state.height, state.h_switch);

  for(h = 0; h < state.h_switch; h ++)
  {
    state.segment_len[h] = 1 << (state.h_switch - h);
    printf("h: %lu, len: %lu, segment size: %lu\n",
           h, state.level_len[h], state.segment_len[h]);
  }

  for(i = 0; i < MAX_PHASE; i ++)
    state.last_h[i] = 0;
}

void init_folders(char *tree_folder, char *saved_folder,
                  char *backup_folder, char *inverse_folder)
{
  folders.tree = tree_folder;
  folders.saved = saved_folder;
  folders.backup = backup_folder;
  folders.inverse = inverse_folder;

  exec("mkdir -p %s %s %s %s",
        folders.tree, folders.saved, folders.backup, folders.inverse);
  exec("rm -f %s/* %s/*", folders.tree, folders.inverse);
}

void run()
{
  ulong h, h_end, step, phase, start, end, len, r_phase;
  mpz_t r0, r1;
  inverse_t inverse;
  struct timespec t_start, t_finish;

  printf("range: [%lu..%lu], saved_at: %lu\n",
        state.p_start, state.p_end, state.p_saved_at);

  // generate the list of primes

  printf("generating the list of primes...\n");
  start_time(&t_start);

  generate_primes();
  print_time(&t_start, &t_finish);

  // main loop

  mpz_init(r0);
  mpz_init(r1);

  inverse.blank = 1;
  for(phase = 0; phase < MAX_PHASE; phase ++)
  {
    start_time(&t_start);
    state.use_ntt = state.use_ntt_mod = 0;

    switch(phase)
    {
      case PHASE_P_TREE:
      case PHASE_N_TREE:
        h = state.h_switch;
        step = STEP_SIZE;

        start = state.p_start - 1;
        end = state.p_end;

        break;
      case PHASE_MUL_DIFFERENCE:
        h = 0;
        step = STEP_SIZE;

        start = state.p_saved_at;
        end = state.p_start - 1;

        inverse.blank = 1; // don't reduce

        break;
      case PHASE_MUL:
        end = start = step = h = 0;
        inverse.blank = 0;

        break;
    }

    switch(phase)
    {
      case PHASE_P_TREE:
        printf("building p-tree...\n");
        break;
      case PHASE_N_TREE:
        printf("building n-tree...\n");
        break;
      case PHASE_MUL_DIFFERENCE:
        printf("multiply all [%lu, %lu)\n", start, end);
        break;
    }

    printf("phase: %lu, interval size: %lu, step: %lu\n", phase, end - start, step);

    len = mul_all_segments(phase, start, end, step, h, &inverse);
    if (phase == PHASE_MUL)
      len = last_stored_index(-1, 0) + 1;

    printf("multiplying elements...\n");

    state.use_ntt = state.use_ntt_mod = 1;
    h_end = mul_all_elements(phase, len, h, &inverse);

    switch(phase)
    {
      case PHASE_P_TREE:
        init_inverse(&inverse, h_end);
        break;
      case PHASE_MUL_DIFFERENCE:
        move_to_saved(phase, 0);
        break;
      case PHASE_MUL:
        load_element(&r0, phase, h_end, 0, 0);
        load_element(&r1, phase, h_end, 0, 1);
        break;
    }

    print_time(&t_start, &t_finish);
  }

  // copy the remaining list

  printf("copying...\n");
  start_time(&t_start);

  move_to_saved(PHASE_N_TREE, 1);
  print_time(&t_start, &t_finish);

  // we don't need inverse any more
  mpz_clear(inverse.value);

  // collect remainders

  printf("collecting remainders...\n");
  start_time(&t_start);

  r_phase = remainders(&r0, &r1);
  collect_remainders(r_phase);

  print_time(&t_start, &t_finish);
}

int main(int argc, char ** argv)
{
  int num_threads;
  ulong ui_p_start, ui_p_end, ui_p_saved_at;
  mpz_t p_start, p_end, p_saved_at;
  struct timespec t_start, t_finish;

  if (argc != 9)
  {
    printf("usage: %s p_start p_end p_saved_at num_threads "
           "tree_folder saved_folder backup_folder inverse_folder\n", argv[0]);
    return 1;
  }

  start_time(&t_start);
  setbuf(stdout, NULL);

  ui_p_start = atoll(argv[1]);
  ui_p_end = atoll(argv[2]);
  ui_p_saved_at = atoll(argv[3]);
  num_threads = atoll(argv[4]);

  mpz_init_set_ui(p_start, ui_p_start);
  mpz_init_set_ui(p_end, ui_p_end);

  mpz_nextprime(p_start, p_start);
  mpz_nextprime(p_end, p_end);

  ui_p_start = mpz_get_ui(p_start);
  ui_p_end = mpz_get_ui(p_end);

  if (ui_p_saved_at > 0)
  {
    mpz_init_set_ui(p_saved_at, ui_p_saved_at - 1);
    mpz_nextprime(p_saved_at, p_saved_at);
    ui_p_saved_at = mpz_get_ui(p_saved_at);

    if (ui_p_saved_at == ui_p_start)
      ui_p_start ++;
  }

  ntt_init();
  set_gmp_ntt_threads(num_threads);

  init_state(ui_p_start, ui_p_end, ui_p_saved_at, num_threads);
  init_folders(argv[5], argv[6], argv[7], argv[8]);

  run();

  printf("done.\n");
  print_time(&t_start, &t_finish);

  return 0;
}
