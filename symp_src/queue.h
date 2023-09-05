//
// Created by joshu on 5/09/2023.
//

#include "kernel.h"
#include "SnapPea.h"

typedef struct Queue {
    int                         front;
    int                         rear;
    int                         len;
    int                         size;
    int                         *array;
} Queue ;

Queue                   *init_queue(int);
Queue                   *enqueue(Queue *, int);
int                     dequeue(Queue *);
Queue                   *resize_queue(Queue *);
Boolean                 empty_queue(Queue *);
void                    free_queue(Queue *);
