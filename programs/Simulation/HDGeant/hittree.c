/*
 * hittree.c - library for managing binary tree of hits pointers
 *
 *	This is a part of the hits package for the
 *	HDGeant simulation program for Hall D.
 *
 *	version 1.0 	-Richard Jones July 16, 2001
 */

#include <stdlib.h>

#include <hittree.h>

void** getTwig(hitTree_t** tree, int mark)
{
   hitTree_t* node = *tree;
   if (node == 0)
   {
      node = *tree = malloc(sizeof(hitTree_t));
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
      return getTwig(&node->right, mark);
   }
}

void* pickTwig(hitTree_t** tree)
{
   hitTree_t* node = *tree;
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
