/*
 * bintree.c - library for managing binary tree of hits pointers
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 */

#include <stdlib.h>
#include <assert.h>
#include <bintree.h>

void** getTwig(binTree_t** tree, int mark)
{
   binTree_t* node = *tree;
   if (node == 0)
   {
      node = *tree = malloc(sizeof(binTree_t));
      node->mark = mark;
      node->left = 0;
      node->right = 0;
      node->this = 0;
      return &node->this;
   }
   else if (mark == node->mark)
   {
      return &node->this;
   }
   else if (mark < node->mark)
   {
      return getTwig(&node->left, mark);
   }
   else
   {
      assert (node->mark >= 0);
      return getTwig(&node->right, mark);
   }
}

void* pickTwig(binTree_t** tree)
{
   binTree_t* node = *tree;
   if (node == 0)
   {
      return 0;
   }
   else if (node->left)
   {
      return pickTwig(&node->left);
   }
   else
   {
      void* twig = node->this;
      *tree = node->right;
      free(node);
      return twig;
   }
}
