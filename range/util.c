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

void exec(const char* format, ...)
{
  char buf[256];

  va_list args;
  va_start(args, format);
  vsprintf(buf, format, args);
  va_end(args);

  if (system(buf))
    error("system error: %s", buf);
}

void error(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  fprintf(stderr, "ERROR :: ");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
  va_end(args);
}

FILE *open_file(const char *mode, const char* format, ...)
{
  FILE *f;
  char buf[256];

  va_list args;
  va_start(args, format);
  vsprintf(buf, format, args);
  va_end(args);

  f = fopen(buf, mode);
  if (!f)
    error("file open error: %s", buf);
  return f;
}

void start_time(struct timespec *t_start)
{
  clock_gettime(CLOCK_MONOTONIC, t_start);
}

void print_time(struct timespec *t_start, struct timespec *t_finish)
{
  clock_gettime(CLOCK_MONOTONIC, t_finish);
  printf("\n\n--- %.3lfs ---\n\n\n",
         (t_finish->tv_sec - t_start->tv_sec) +
         (t_finish->tv_nsec - t_start->tv_nsec) / 1e9f);
}
