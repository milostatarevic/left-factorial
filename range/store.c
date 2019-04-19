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

void load_segment(mpz_t *list, ulong phase, ulong len, ulong h, ulong segment, ulong t)
{
  ulong i, n, buffer_len, write_len;
  char *write_buffer, marker;
  FILE *f;

  f = open_file("rb", "%s/c%lu-l%lu-s%lu-t%lu",
                folders.tree, phase, h, segment, t);
  if (h <= 1)
  {
    fread(&write_len, sizeof(ulong), 1, f);

    write_buffer = (char *) malloc(write_len + 16);
    fread(write_buffer, write_len, 1, f);

    buffer_len = 0;
    for (i = 0; i < len; i ++)
    {
      marker = write_buffer[buffer_len];
      buffer_len ++;

      if (marker <= 1)
        mpz_set_ui(list[i], marker);
      else
      {
        memcpy(&n, write_buffer + buffer_len, sizeof(ulong));
        buffer_len += sizeof(ulong);

        mpz_set_ui(list[i], n);
      }

      if (buffer_len > write_len)
        error("overflow? phase: %lu, h: %lu, segment: %lu, t: %lu "
              "len: %lu, i: %lu, buffer_len: %lu, write len: %lu",
              phase, h, segment, t, len, i, buffer_len, write_len);
    }
    free(write_buffer);
  }
  else
  {
    for(i = 0; i < len; i ++)
      mpz_read(f, &list[i]);
  }

  fclose(f);
}

void save_segment(mpz_t *list, ulong phase, ulong len, ulong h, ulong segment, ulong t)
{

  ulong i, n, buffer_len;
  char *write_buffer, marker;
  FILE *f;

  f = open_file("wb", "%s/c%lu-l%lu-s%lu-t%lu",
                folders.tree, phase, h, segment, t);
  if (h <= 1)
  {
    write_buffer = (char *) malloc(len * (sizeof(ulong) + 1));

    buffer_len = 0;
    for (i = 0; i < len; i ++)
    {
      if (list[i][0]._mp_size > 1)
        error("unexpected mp_size > 1");

      n = mpz_get_ui(list[i]);
      marker = (n <= 1) ? n : sizeof(ulong);

      write_buffer[buffer_len] = marker;
      buffer_len ++;

      if (marker > 1)
      {
        memcpy(write_buffer + buffer_len, &n, sizeof(ulong));
        buffer_len += sizeof(ulong);
      }
    }

    fwrite(&buffer_len, sizeof(ulong), 1, f);
    fwrite(write_buffer, buffer_len, 1, f);

    free(write_buffer);
  }
  else
  {
    for(i = 0; i < len; i ++)
      mpz_write(f, &list[i]);
  }


  fclose(f);
}

void load_element(mpz_t *r, ulong phase, ulong h, ulong i, ulong t)
{
  FILE *f;
  char buf[128];

  if (phase == PHASE_MUL && h == 0)
    sprintf(buf, "%s/s-e%lu-t%lu", folders.saved, i, t);
  else
    sprintf(buf, "%s/c%lu-l%lu-e%lu-t%lu", folders.tree, phase, h, i, t);

  f = fopen(buf, "rb");
  if (f)
  {
    mpz_read(f, r);
    fclose(f);
  }
  else
    mpz_set_ui(*r, 1);
}

void save_element(mpz_t *r, ulong phase, ulong h, ulong i, ulong t, inverse_t *inverse)
{
  FILE *f;

  // backup before reducing
  if (phase == PHASE_N_TREE && h == state.last_h[phase] && inverse)
  {
    f = open_file("wb", "%s/n-e%lu-t%lu", folders.backup, i, t);
    mpz_write(f, r);
    fclose(f);
  }

  f = open_file("wb", "%s/c%lu-l%lu-e%lu-t%lu", folders.tree, phase, h, i, t);
  reduce_mod(r, inverse);
  mpz_write(f, r);
  fclose(f);
}

slong last_stored_index(slong phase, int mode)
{
  slong index, i;
  FILE *pf;
  char buf0[128], buf1[128];

  if (phase == PHASE_N_TREE && mode == 1)
    sprintf(buf0, "%s", folders.backup);
  else if (phase > 0)
    sprintf(buf0, "%s/c%lu-l%lu-*", folders.tree, phase, state.last_h[phase]);
  else
    sprintf(buf0, "%s/s-*", folders.saved);

  sprintf(buf1, "ls -1tr %s | awk -F\"-e\" '{print $2}' | awk -F\"-\" '{print $1}'", buf0);
  pf = popen(buf1, "r");

  index = -1;
  while(fgets(buf0, 128, pf))
    if ((i = atoll(buf0)) > index)
      index = i;
  pclose(pf);

  return index;
}

void move_to_saved(ulong phase, int mode)
{
  slong last_saved_index, last_index, i, t;

  last_saved_index = last_stored_index(-1, mode);
  last_index = last_stored_index(phase, mode);

  printf("last saved index: %ld, last_index: %ld\n", last_saved_index, last_index);

  for(i = 0; i <= last_index; i ++)
  {
    for(t = 0; t < 2; t ++)
    {
      // backup -> saved
      if (mode == 1)
      {
        exec("mv -vn %s/n-e%lu-t%lu %s/s-e%lu-t%lu",
             folders.backup, i, t,
             folders.saved, i + last_saved_index + 1, t);
      }
      // tree -> saved
      else
      {
        exec("mv -vn %s/c%lu-l%lu-e%lu-t%lu %s/s-e%lu-t%lu",
             folders.tree, phase, state.last_h[phase], i, t,
             folders.saved, i + last_saved_index + 1, t);
      }
    }
  }
}
