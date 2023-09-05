/*
 * queue.c
 *
 * Basis queue implementation, used in various places.
 */


#include "symplectic_kernel.h"

Queue *resize_queue(Queue *);

Queue *init_queue(int size) {
    Queue *q    = NEW_STRUCT( Queue );

    q->front    = 0;
    q->rear     = -1;
    q->len      = 0;
    q->size     = MAX(size, 256);
    q->array    = NEW_ARRAY(q->size, int);

    return q;
}

Queue *enqueue(Queue *q, int i) {
    // Queue is full
    if ( q->size == q->len ) {
        q = resize_queue(q);
        q = enqueue(q, i);
    } else {
        q->rear = (q->rear + 1) % q->size;
        q->array[q->rear] = i;
        q->len++;
    }

    return q;
}

int dequeue(Queue *q) {
    // User to verify queue is not empty
    int i = q->array[q->front];

    q->front = (q->front + 1) % q->size;
    q->len--;

    return i;
}

Boolean empty_queue(Queue *q) {
    if (q->len > 0)
        return FALSE;

    return TRUE;
}

Queue *resize_queue(Queue *q) {
    int i;
    Queue *p = init_queue(2 * q->size);

    // Copy elements to new array
    while (!empty_queue(q)) {
        i = dequeue(q);
        enqueue(p, i);
    }

    free_queue(q);
    return p;
}

void free_queue(Queue *q) {
    my_free(q->array);
    my_free(q);
}