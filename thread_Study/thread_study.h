#ifndef THREAD_STUDY_H_INCLUDED
#define THREAD_STUDY_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

void threadTest(int cnt_threads, int64_t cnt_iteration);


#endif