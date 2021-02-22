#include "Queue.h"


static void _queueTest();


void _QueueTestSet() {
    _queueTest();
}



/*
 * Tests for working functions.
 */

static void _queueTest() {
    printf("\n**************** _queueTest ****************\n");
    uint64_t returnValue = 0;

    Queue* queueInstance = initQueue();

    returnValue = isQueueEmpty(queueInstance);
    printf("Queue is empty or not: %"PRIu64"\n", returnValue);

    QueueCell* queueCell = newQueueCell(1);
    printf("New queue cell's data: %"PRIu64"\n", queueCell->data);

    // dequeue an empty queue
    deQueue(queueInstance, queueCell);

    // construct a queue
    printf("Construct a queue ...\n");
    for(uint64_t i = 0; i < 10; i++) {
        enQueue(queueInstance, newQueueCell(i));
    }
    printf("Queue:\n");
    printQueue(queueInstance);

    // dequeue the queue
    printf("Dequeue the queue ...\n");
    printf("queue length: %"PRIu64"\n", queueInstance->length);
    returnValue = queueInstance->length;
    for(uint64_t i = 0; i < returnValue / 2; i++) {
        deQueue(queueInstance, queueCell);
        /**< \todo */
        printf("dequeue the queue cell: (%"PRIu64", 0x%p)\n", queueCell->data, queueCell->next);
    }
    printf("Queue:\n");
    printQueue(queueInstance);

    // clear the queue
    printf("Clear the queue (0x%p) ...\n", queueInstance);
    clearQueue(queueInstance);
}




/*
 * Working functions.
 */





QueueCell* newQueueCell(uint64_t data) {
    /**< \note memory required by malloc will not be freed when the function ends */
    QueueCell* queueCell = (QueueCell*)malloc(sizeof(QueueCell));
    if(queueCell == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    queueCell->data = data;
    queueCell->next = NULL;
    return queueCell;
}


Queue* initQueue() {
    Queue* queueInstance = (Queue*)malloc(sizeof(Queue));
    if(queueInstance == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    queueInstance->head = NULL;
    queueInstance->tail = NULL;
    queueInstance->length = 0;
    return queueInstance;
}


void setQueue(Queue* queueInstance, QueueCell* head, QueueCell* tail, uint64_t length) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred when setting a queue. \n");
        exit(EXIT_FAILURE);
    }
    if((head == NULL) ^ (tail == NULL)) {
        printf("ERROR: cannot set a queue with only a head or a tail. \n");
        exit(EXIT_FAILURE);
    }
    queueInstance->head = head;
    queueInstance->tail = tail;
    queueInstance->length = length;
}



void clearQueue(Queue* queueInstance) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred when clearing a queue. \n");
        exit(EXIT_FAILURE);
    }
    QueueCell* thisCell = NULL;
    QueueCell* nextCell = NULL;

    thisCell = queueInstance->head;
//    printf("queue length: %"PRIu64"\n", queueInstance->length);
    while(thisCell != NULL) {
//        printf("0x%p: (%"PRIu64", 0x%p) ->\n", thisCell, thisCell->data, thisCell->next);
        nextCell = thisCell->next;
        free(thisCell);
        queueInstance->length = queueInstance->length - 1;
        thisCell = nextCell;
    }
//    printf("NULL\n\n");
    queueInstance->head = NULL;
    queueInstance->tail = NULL;
    free(queueInstance);
}


uint64_t isQueueEmpty(Queue* queueInstance) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred judging whether a queue is empty. \n");
        exit(EXIT_FAILURE);
    }
    if(queueInstance->length == 0) {
        return QUEUE_EMPTY;
    }
    return QUEUE_NOT_EMPTY;
}


uint64_t queueLength(Queue* queueInstance) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred when getting the length of a queue. \n");
        exit(EXIT_FAILURE);
    }
    return queueInstance->length;
}


QueueCell* getQueueHead(Queue* queueInstance) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred when getting the head of a queue. \n");
        exit(EXIT_FAILURE);
    }
    return queueInstance->head;
}


QueueCell* getQueueTail(Queue* queueInstance) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred when getting the tail of a queue. \n");
        exit(EXIT_FAILURE);
    }
    return queueInstance->tail;
}


void deQueue(Queue* queueInstance, QueueCell* queueCell) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred when dequeuing a queue. \n");
        exit(EXIT_FAILURE);
    }
    QueueCell* queueHead = queueInstance->head;

    if(queueHead == NULL) {
        /**< \note if queue is empty  */
        printf("WARNING: Queue head is NULL!\n");
        queueCell = NULL;
        return;
    }

    queueCell->data = queueHead->data;
    queueCell->next = queueHead->next;

    free(queueHead);
    queueInstance->head = queueCell->next;
    queueInstance->length = queueInstance->length - 1;
    if(queueInstance->length == 0) {    // make the head consists with tail
        queueInstance->head = NULL;
        queueInstance->tail = NULL;
    }
}


void enQueue(Queue* queueInstance, QueueCell* queueCell) {
    if(queueInstance == NULL || queueCell == NULL) {
        printf("ERROR: null pointer occurred when dequeuing a queue. \n");
        exit(EXIT_FAILURE);
    }
    QueueCell* queueTail = queueInstance->tail;

    if(queueInstance->length == 0) {
        queueInstance->head = queueCell;
        queueInstance->tail = queueCell;
        queueInstance->length = queueInstance->length + 1;
        return;
    }

    queueTail->next = queueCell;

    queueInstance->tail = queueCell;
    queueInstance->length = queueInstance->length + 1;
}


void printQueue(Queue* queueInstance) {
    if(queueInstance == NULL) {
        printf("ERROR: null pointer occurred when printing a queue. \n");
        exit(EXIT_FAILURE);
    }
    QueueCell* queueCell = NULL;

    queueCell = queueInstance->head;
    printf("queue length: %"PRIu64"\n", queueInstance->length);
//    if(queueCell != NULL) {
//        printf("queue head: (%#"PRIx64", 0x%p)\n",
//               queueInstance->head->data, queueInstance->head->next);
//        printf("queue tail: (%#"PRIx64", 0x%p)\n",
//               queueInstance->tail->data, queueInstance->tail->next);
//    }
    while(queueCell != NULL) {
        printf("0x%p: (%#"PRIx64", 0x%p) ->\n", queueCell, queueCell->data, queueCell->next);
        queueCell = queueCell->next;
    }
    printf("(queue ends)\n\n");
}






/*
 * Static functions. (file-localized functions)
 */
