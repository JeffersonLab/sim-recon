/*
 * memcheck - a simple memory management checking tool
 *
 *	Typically the management of a memory structure is restricted to
 *	a limited segment of code.  General malloc/free memory leak
 *	tools can be found that will trap every call to malloc or free.
 *	Often it is simpler just to insert some checkpoint calls around
 *	the relevant calls, and just study the behavior in that region.
 *	This is the purpose of the memcheck routines.
 *
 * Richard Jones  - July 18, 2000
 * University of Connecticut
 *
 *
 * Instructions:
 * -------------
 * 1) After each relevant malloc, insert a call to checkin(pointer) as 
 *	p = malloc(n);   	// old code
 *	checkin(p,string);	// user string helps trace memory leaks 
 *    or the following more compact form will have the same effect
 *	p = checkin(malloc(size_t),string);
 *
 * 2) Before each relevant free, insert a call to checkout(pointer) as
 *	checkout(p);		// new insertion
 *	free(p);        	// old code
 *
 * 3) Any time you think the memory balance should be zero do checkpoint()
 *	checkpoint();		// look for leaks
 *
 * 
 * Programmer's Notes:
 * -------------------
 * 1) The "bintree" binary tree package is used to store the allocation
 *    tables.
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <bintree.h>

typedef struct {
   void* addr;
   int count;
   char* tag;
} memblock_t;

binTree_t* memcheckTree = 0;
int* addressRef = 0;
int nodeCount = 0;

void* checkin (void* p, char* tag)
{
   int mark = (int*)p - addressRef;
   void** twig = getTwig(&memcheckTree, mark);
   if (*twig == 0)
   {
      memblock_t *blk = *twig = malloc(sizeof(memblock_t));
      if (tag)
      {
         blk->tag = malloc(strlen(tag)+1);
         strcpy(blk->tag,tag);
      }
      else
      {
         blk->tag = malloc(7);
         strcpy(blk->tag,"(null)");
      }
      blk->count = 1;
      blk->addr = p;
      nodeCount++;
   }
   else if (((memblock_t*) *twig)->count == 0)
   {
      memblock_t *blk = *twig;
      if (blk->tag)
      {
         free(blk->tag);
      }
      if (tag)
      {
         blk->tag = malloc(strlen(tag)+1);
         strcpy(blk->tag,tag);
      }
      else
      {
         blk->tag = malloc(7);
         strcpy(blk->tag,"(null)");
      }
      blk->count = 1;
      blk->addr = p;
   }
   else
   {
      memblock_t *blk = *twig;
      fprintf(stderr,"memcheck report:");
      fprintf(stderr," reallocation of allocated memory block\n");
      fprintf(stderr," original tag was %s\n",blk->tag);
      fprintf(stderr," second tag was %s\n",tag);
      assert (1 == 0);
   }
   return p;
}

void* checkout (void* p)
{
   int mark = (int*)p - addressRef;
   void** twig = getTwig(&memcheckTree, mark);
   if (*twig == 0)
   {
      fprintf(stderr,"memcheck report:");
      fprintf(stderr," attempt to free unallocated memory block\n");
//      assert (1 == 0);
   }
   else if (((memblock_t*) *twig)->count < 1)
   {
      memblock_t *blk = *twig;
      fprintf(stderr,"memcheck report:");
      fprintf(stderr," attempt to refree freed memory block\n");
      fprintf(stderr," tag was %s\n",blk->tag);
      assert (1 == 0);
   }
   else
   {
      memblock_t *blk = *twig;
      blk->count = 0;
   }
   return p;
}

void checkpoint ()
{
   memblock_t* node;
   int abort = 0;
   while (node = pickTwig(&memcheckTree))
   {
      if (node->count > 0)
      {
         fprintf(stderr,"memcheck report:");
         fprintf(stderr," checkpoint found allocated memory block\n");
         fprintf(stderr," tag was %s\n",node->tag);
         ++abort;
      }
      nodeCount--;
      free(node->tag);
      free(node);
   }
   if (abort)
   {
      fprintf(stderr," quitting because of above error%s.\n",
              ((abort == 1) ? "" : "s"));
      exit(1);
   }
}
