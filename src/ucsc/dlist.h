/* dlist.h - Headers for generic doubly-linked list routines. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#ifndef DLIST_H
#define DLIST_H

#ifndef COMMON_H
#include "common.h"
#endif

struct dlNode
/* An element on a doubly linked list. */
    {
    struct dlNode *next;
    struct dlNode *prev;
    void *val;
    };

struct dlList
/* A doubly linked list. */
    {
    struct dlNode *head;
    struct dlNode *nullMiddle;
    struct dlNode *tail;
    };

#define dlEnd(node) (node->next == NULL)
/* True if node past end. */

#define dlStart(node) (node->prev == NULL)
/* True if node before start. */

/* Iterate on a doubly linked list as so:
    for (el = list->head; !dlEnd(el); el = el->next)
        val = el->val;
   or
    for (el = list->tail; !dlStart(el); el = el->prev)
        val = el->val;
 */

struct dlList *newDlList();
/* Return a new doubly linked list. */

#define dlListNew newDlList
/* Add object-first synonym. */

void dlListInit(struct dlList *dl);
/* Initialize list to be empty */

void dlListReset(struct dlList *dl);
/* Reset a list to the empty state (does not free values)  */

void freeDlList(struct dlList **pList);
/* Free up a doubly linked list and it's nodes (but not the node values). */
#define dlListFree freeDlList


void dlAddBefore(struct dlNode *anchor, struct dlNode *newNode);
/* Add a node to list before anchor member. */

void dlAddAfter(struct dlNode *anchor, struct dlNode *newNode);
/* Add a node to list after anchor member. */

void dlAddHead(struct dlList *list, struct dlNode *newNode);
/* Add a node to head of list. */

void dlAddTail(struct dlList *list, struct dlNode *newNode);
/* Add a node to tail of list. */

void dlRemove(struct dlNode *node);
/* Removes a node from list. Node is not freed. */

struct dlNode *dlPopHead(struct dlList *list);
/* Remove first node from list and return it. */


int dlCount(struct dlList *list);
/* Return length of list. */

boolean dlEmpty(struct dlList *list);
/* Return TRUE if list is empty. */

#define dlIsEmpty(list) ((list)->head->next == NULL)
/* Return TRUE if list is empty.  Macro version of above. */


#endif /* DLIST_H */


