/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


void* safe_malloc(size_t n)
{
  void* ptr = malloc(n);
  if (ptr == NULL)
    fatal_error("out of memory, n = %zu", n);
  return ptr;
}


void* safe_realloc(void* ptr, size_t n)
{
  ptr = realloc(ptr, n);
  if (ptr == NULL)
    fatal_error("out of memory, n = %zu", n);
  return ptr;
}


void safe_free(void* ptr)
{
  free(ptr);
}


void list_insert(list_node_t** list, uint64_t* ptr)
{
  list_node_t* node = safe_malloc(sizeof(list_node_t));
  node->ptr = ptr;
  node->next = *list;
  *list = node;
}


uint64_t* list_extract(list_node_t** list)
{
  if (*list == NULL)
    return NULL;
  list_node_t* node = *list;
  uint64_t* ptr = node->ptr;
  *list = node->next;
  safe_free(node);
  return ptr;
}


void pool_init(pool_t* pool, size_t block_size)
{
  INIT_LOCK(&pool->lock);
  pool->block_size = block_size;
  pool->available = NULL;
  pool->from_heap = NULL;
}


void pool_clear(pool_t* pool)
{
  while (pool->from_heap != NULL)
    safe_free(list_extract(&pool->from_heap));

  while (pool->available != NULL)
    list_extract(&pool->available);

  DESTROY_LOCK(&pool->lock);
}


void pool_insert(pool_t* pool, uint64_t* ptr)
{
  SET_LOCK(&pool->lock);
  list_insert(&pool->available, ptr);
  UNSET_LOCK(&pool->lock);
}


uint64_t* pool_request(pool_t* pool)
{
  SET_LOCK(&pool->lock);

  uint64_t* ptr = list_extract(&pool->available);
  if (ptr == NULL)
    {
      // no available blocks; need to get one from the heap
      ptr = safe_malloc(pool->block_size * sizeof(uint64_t));
      list_insert(&pool->from_heap, ptr);
    }

  UNSET_LOCK(&pool->lock);
  return ptr;
}


void vbuf_init(vbuf_t* vbuf, size_t n, pool_t* pool)
{
  vbuf->wrap = NULL;

  INIT_LOCK(&vbuf->lock);
  vbuf->block_size = pool->block_size;
  vbuf->pool = pool;
  vbuf->n = n;

  vbuf->regions_alloc = 1;
  vbuf->regions_size = 0;
  vbuf->regions = safe_malloc(2 * sizeof(size_t));

  vbuf->blocks = safe_malloc(n * sizeof(uint64_t*));
  for (size_t i = 0; i < n; i++)
    vbuf->blocks[i] = NULL;
}


void vbuf_init_wrap(vbuf_t* vbuf, size_t block_size, uint64_t* wrap)
{
  vbuf->wrap = wrap;
  vbuf->block_size = block_size;
}


void vbuf_clear(vbuf_t* vbuf)
{
  if (vbuf->wrap == NULL)
    {
      for (size_t i = 0; i < vbuf->n; i++)
	{
	  if (vbuf->blocks[i] != NULL)
	    pool_insert(vbuf->pool, vbuf->blocks[i]);
	}
      safe_free(vbuf->blocks);
      safe_free(vbuf->regions);

      DESTROY_LOCK(&vbuf->lock);
    }
}


uint64_t* vbuf_get_block(vbuf_t* vbuf, size_t i)
{
  if (vbuf->wrap != NULL)
    {
      assert(i == 0);
      return vbuf->wrap;
    }

  SET_LOCK(&vbuf->lock);
  
  assert(i < vbuf->n);
  uint64_t* ptr = vbuf->blocks[i];
  if (ptr == NULL)
    ptr = vbuf->blocks[i] = pool_request(vbuf->pool);

  UNSET_LOCK(&vbuf->lock);
  return ptr;
}


uint64_t* vbuf_get(vbuf_t* vbuf, size_t i)
{
  if (vbuf->wrap != NULL)
    {
      assert(i < vbuf->block_size);
      return vbuf->wrap + i;
    }

  size_t block = i / vbuf->block_size;
  size_t offset = i - block * vbuf->block_size;
  return vbuf_get_block(vbuf, block) + offset;
}


void vbuf_free_region(vbuf_t* vbuf, size_t start, size_t end)
{
  assert(end >= start);

  if (vbuf->wrap != NULL)
    return;

  SET_LOCK(&vbuf->lock);

  // make room for one more pair in regions
  if (vbuf->regions_alloc == vbuf->regions_size)
    {
      vbuf->regions_alloc *= 2;
      vbuf->regions =
	safe_realloc(vbuf->regions, 2 * vbuf->regions_alloc * sizeof(size_t));
    }

  // search for where to add new region
  ptrdiff_t j;
  for (j = 0; j < vbuf->regions_size; j++)
    if (end <= vbuf->regions[2*j])
      break;    // found insertion point

  assert(j == 0 || start >= vbuf->regions[2*j - 1]);

  int merge_prev = (j >= 1) && (start == vbuf->regions[2*j - 1]);
  int merge_next = (j < vbuf->regions_size) && (end == vbuf->regions[2*j]);

  if (merge_prev)
    {
      if (merge_next)
	{
	  // merge with previous and next regions; delete one region
	  vbuf->regions_size--;
	  for (size_t i = 2 * j - 1; i < 2 * vbuf->regions_size; i++)
	    vbuf->regions[i] = vbuf->regions[i + 2];
	}
      else
	// merge with previous region only
	vbuf->regions[2*j - 1] = end;

      j--;
    }
  else
    {
      if (merge_next)
	// merge with next region only
	vbuf->regions[2*j] = start;
      else
	{
	  // merge with neither previous nor next region; insert one region
	  for (ptrdiff_t i = 2 * vbuf->regions_size - 1; i >= 2 * j; i--)
	    vbuf->regions[i + 2] = vbuf->regions[i];
	  vbuf->regions_size++;
	  vbuf->regions[2 * j] = start;
	  vbuf->regions[2 * j + 1] = end;
	}
    }

  // newly created region is new_start, new_end
  size_t new_start = vbuf->regions[2 * j];
  size_t new_end = vbuf->regions[2 * j + 1];

  // scan through pool blocks that overlap (start, end), and return them to
  // the memory pool if possible
  for (size_t i = start / vbuf->block_size;
       i < vbuf->n && i * vbuf->block_size < end; i++)
    {
      if (vbuf->blocks[i] != NULL &&
	  i * vbuf->block_size >= new_start &&
	  (i + 1) * vbuf->block_size <= new_end)
	{
	  pool_insert(vbuf->pool, vbuf->blocks[i]);
	  vbuf->blocks[i] = NULL;
	}
    }
  
  UNSET_LOCK(&vbuf->lock);
}
