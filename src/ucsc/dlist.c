/* dlist.c - Doubly-linked list routines. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */
#include "common.h"
#include "dlist.h"


void dlListInit(struct dlList *dl)
/* Initialize list to be empty */
{
dl->head = (struct dlNode *)(&dl->nullMiddle);
dl->nullMiddle = NULL;
dl->tail = (struct dlNode *)(&dl->head);
}

struct dlList *newDlList()
/* Return a new doubly linked list. */
{
struct dlList *dl;
AllocVar(dl);
dl->head = (struct dlNode *)(&dl->nullMiddle);
dl->tail = (struct dlNode *)(&dl->head);
return dl;
}

void dlListReset(struct dlList *dl)
/* Reset a list to the empty state (does not free values)  */
{
struct dlNode *node, *next;
for (node = dl->head; node->next != NULL; node = next)
    {
    next = node->next;
    freeMem(node);
    }
dl->head = (struct dlNode *)(&dl->nullMiddle);
dl->nullMiddle = NULL;
dl->tail = (struct dlNode *)(&dl->head);
}

void freeDlList(struct dlList **pList)
/* Free up a doubly linked list and it's nodes (but not the node values). */
{
struct dlList *list = *pList;
if (list != NULL)
    {
    dlListReset(list);
    freez(pList);
    }
}


void dlInsertBetween(struct dlNode *before, struct dlNode *after, struct dlNode *newNode)
{
before->next = newNode; 
newNode->prev = before; 
newNode->next = after;  
after->prev = newNode; 
}

void dlAddBefore(struct dlNode *anchor, struct dlNode *newNode)
/* Add a node to list before anchor member. */
{
dlInsertBetween(anchor->prev, anchor, newNode);
}

void dlAddAfter(struct dlNode *anchor, struct dlNode *newNode)
/* Add a node to list after anchor member. */
{
dlInsertBetween(anchor, anchor->next, newNode);
}

void dlAddHead(struct dlList *list, struct dlNode *newNode)
/* Add a node to head of list. */
{
struct dlNode *head = list->head;
dlInsertBetween(head->prev, head, newNode);
}

void dlAddTail(struct dlList *list, struct dlNode *newNode)
/* Add a node to tail of list. */
{
struct dlNode *tail = list->tail;
dlInsertBetween(tail, tail->next, newNode);
}

void dlRemove(struct dlNode *node)
/* Removes a node from list. Node is not freed. */
{
struct dlNode *before = node->prev;
struct dlNode *after = node->next;
before->next = after;
after->prev = before;
node->prev = NULL;
node->next = NULL;
}

struct dlNode *dlPopHead(struct dlList *list)
/* Remove first node from list and return it. */
{
struct dlNode *node = list->head;
if (node->next == NULL)
    return NULL;
dlRemove(node);
return node;
}

int dlCount(struct dlList *list)
/* Return length of list. */
{
return slCount(list->head) - 1;
}


struct dlSorter 
/* Helper structure for sorting dlNodes preserving order */
    {
    struct dlNode *node;
    };

static int (*compareFunc)(const void *elem1, const void *elem2);
/* Node comparison pointer, just used by dlSortNodes and helpers. */

static int dlNodeCmp(const void *elem1, const void *elem2)
/* Compare two dlSorters indirectly, by calling compareFunc. */
{
struct dlSorter *a = (struct dlSorter *)elem1;
struct dlSorter *b = (struct dlSorter *)elem2;
return compareFunc(&a->node->val, &b->node->val);
}
    
boolean dlEmpty(struct dlList *list)
/* Return TRUE if list is empty. */
{
return dlIsEmpty(list);
}
