#include "thread_study.h"

pthread_mutex_t mutex_statistics = PTHREAD_MUTEX_INITIALIZER;

typedef struct _define_Statistics_thread {
  double time_taken;
} Statistics_thread;

typedef struct _define_Data_thread {
  int64_t id;
  int64_t start;
  int64_t end;
  Statistics_thread *statistics;
} Data_thread;

void *thread_iterate(void *thread_data) {
  Data_thread *data = (Data_thread *)thread_data;

  // printf("thread (%d) - start: %" PRId64 ", end: %" PRId64 "\n", data->id,
  //        data->start, data->end);

  // Iterate
  clock_t time_start = clock();
  int64_t iteration = 0;
  for (int64_t i = data->start; i < data->end; i++) {
    iteration++;
  }
  clock_t time_end = clock();

  float time_taken = (float)(time_end - time_start) / CLOCKS_PER_SEC;
  // printf("time start: %lu, time end: %lu\n", time_start, time_end);
  printf("time taken: %fs for thread (%" PRId64 ")\n", time_taken, data->id);

  // Lock
  if (pthread_mutex_lock(&mutex_statistics) != 0) {
    perror("pthread_mutex_lock");
    exit(EXIT_FAILURE);
  }
  data->statistics->time_taken += time_taken;
  // Unlock
  if (pthread_mutex_unlock(&mutex_statistics) != 0) {
    perror("pthread_mutex_unlock");
    exit(EXIT_FAILURE);
  }
  return (void *)(data->id);
}

void threadTest(int cnt_threads, int64_t cnt_iteration) {
  static Statistics_thread statistics;
  pthread_t threads[cnt_threads];
  Data_thread datas[cnt_threads];

  // Assign groups for data
  for (int i = 0; i < cnt_threads; i++) {
    datas[i].id = i;
    datas[i].start = (cnt_iteration / cnt_threads) * i;
    datas[i].statistics = &statistics;
  }
  for (int i = 0; i < cnt_threads - 1; i++) {
    datas[i].end = datas[i + 1].start - 1;
  }
  datas[cnt_threads - 1].end = cnt_iteration - 1;

  // Create threads
  for (int i = 0; i < cnt_threads; i++) {
    if (pthread_create(&threads[i], NULL, thread_iterate, (void *)&datas[i]) !=
        0) {
      fprintf(stderr, "Error: thread creation failed for thread %d\n", i);
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < cnt_threads; i++) {
    void *thread_ret = 0;
    pthread_join(threads[i], &thread_ret);
    printf("thread (%" PRId64 ") ended\n", (int64_t)thread_ret);
  }
}