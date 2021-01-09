#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct threadpool_t threadpool_t;

threadpool_t *threadpool_create(int thread_count, int queue_size, int flags);

int threadpool_add(threadpool_t *pool, void (*routine)(void *),
                   void *arg, int flags);


typedef struct {
    void (*function)(void *);
    void *argument;
} threadpool_task_t;

struct threadpool_t {
  pthread_mutex_t lock;
  pthread_cond_t notify;
  pthread_t *threads;
  threadpool_task_t *queue;
  int thread_count;
  int queue_size;
  int head;
  int tail;
  int count;
  int shutdown;
  int started;
};


static void *threadpool_thread(void *threadpool);

int threadpool_free(threadpool_t *pool);

threadpool_t *threadpool_create(int thread_count, int queue_size, int flags)
{
    threadpool_t *pool;

    if((pool = (threadpool_t *)malloc(sizeof(threadpool_t))) == NULL) {
        threadpool_free(pool);
            return NULL;
    }

    pool->thread_count =pool->head = pool->tail = pool->count =pool->shutdown = pool->started = 0;
    pool->queue_size = queue_size;

    pool->threads = (pthread_t *)malloc(sizeof(pthread_t) * thread_count);
    pool->queue = (threadpool_task_t *)malloc(sizeof(threadpool_task_t) * queue_size);
    // all kinds of problems fail so free pool
    if((pthread_mutex_init(&(pool->lock), NULL) != 0) || (pthread_cond_init(&(pool->notify), NULL) != 0) || (pool->threads == NULL) || (pool->queue == NULL)) {
            threadpool_free(pool);
            return NULL;
    }

    for(int i = 0; i < thread_count; i++) {
        // create thread fail
	if(pthread_create(&(pool->threads[i]), NULL,threadpool_thread, (void*)pool) != 0) {
            threadpool_free(pool);
            return NULL;
        }
        pool->thread_count++;
        pool->started++;
    }

    return pool;

}

int threadpool_add(threadpool_t *pool, void (*function)(void *), void *argument, int flags){
    //this pool will deal with the queue to add new tasks into it
    int next;
    if(pthread_mutex_lock(&(pool->lock)) != 0){ 
        return -1;}

    next = pool->tail + 1;
    next = (next == pool->queue_size) ? 0 : next;

    do{

        pool->queue[pool->tail].function = function;
        pool->queue[pool->tail].argument = argument;

        pool->tail = next;
        pool->count += 1;

        if(pthread_cond_signal(&(pool->notify)) != 0)
         return -1;

    } while(0);

    if(pthread_mutex_unlock(&pool->lock) != 0) 
        return -1;
    return 0;
}

int threadpool_free(threadpool_t *pool){
    // no pool or pool still with remaining tasks
    if(pool == NULL || pool->started > 0)
        return -1;

    if(pool->threads) {
        free(pool->threads);
        free(pool->queue);
        pthread_mutex_lock(&(pool->lock));
        pthread_mutex_destroy(&(pool->lock));
        pthread_cond_destroy(&(pool->notify));
    }
    free(pool);
    return 0;
}


static void *threadpool_thread(void *threadpool){
    threadpool_t *pool = (threadpool_t *)threadpool;
    threadpool_task_t task;
    // this thread will continuously take task from queue
    while(1) {

        pthread_mutex_lock(&(pool->lock));
        //wait new task
        while(pool->count == 0 ) {
            pthread_cond_wait(&(pool->notify), &(pool->lock));
        }

        task.function = pool->queue[pool->head].function;
        task.argument = pool->queue[pool->head].argument;

        pool->head += 1;   
	if (pool->head == pool->queue_size) 
	 pool->head = 0;
        pool->count -= 1;


        pthread_mutex_unlock(&(pool->lock));

        (*(task.function))(task.argument);

    }

    pool->started--;

    pthread_mutex_unlock(&(pool->lock));
    pthread_exit(NULL);
    return NULL;
}


int* image;

//double pi;
     int iters;
    double left;
    double right;
    double lower;
    double upper;
    int width;
    int height;

int* tasks;
int left_n;
pthread_mutex_t lock;

void thread2thread(void *ptr2)
{
    int j=*(int*)ptr2;
 //   printf("%d \n",j);
    double y0 = j * ((upper - lower) / height) + lower;
          for (int i = 0; i < width; ++i) {
            double x0 = i * ((right - left) / width) + left;
            int repeats = 0;
            double x = 0;
            double y = 0;
            double length_squared = 0;
            while (repeats < iters && length_squared < 4) {
                double temp = x * x - y * y + x0;
                y = 2 * x * y + y0;
                x = temp;
                length_squared = x * x + y * y;
                ++repeats;
            }
            image[j * width + i] = repeats;
          }
    pthread_mutex_lock(&lock);
    left_n--;
    pthread_mutex_unlock(&lock);
}



void write_png(const char* filename, int iters, int width, int height, const int* buffer) {
    FILE* fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
    png_write_info(png_ptr, info_ptr);
    png_set_compression_level(png_ptr, 1);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; ++x) {
            int p = buffer[(height - 1 - y) * width + x];
            png_bytep color = row + x * 3;
            if (p != iters) {
                if (p & 16) {
                    color[0] = 240;
                    color[1] = color[2] = p % 16 * 16;
                } else {
                    color[0] = p % 16 * 16;
                }
            }
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

int main(int argc, char** argv) {
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    printf("%d cpus available\n", CPU_COUNT(&cpu_set));
    int n_threads= CPU_COUNT(&cpu_set);
    /* argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    iters = strtol(argv[2], 0, 10);
    left = strtod(argv[3], 0);
    right = strtod(argv[4], 0);
    lower = strtod(argv[5], 0);
    upper = strtod(argv[6], 0);
    width = strtol(argv[7], 0, 10);
    height = strtol(argv[8], 0, 10);
       threadpool_t *pool;

    /* allocate memory for image */
     image = (int*)malloc(width * height * sizeof(int));
     tasks = (int*)malloc( height * sizeof(int));
    assert(image);
   int copy=1;
    left_n=height;

  
    pthread_mutex_init(&lock, NULL);

    pool = threadpool_create(n_threads, height, 0);
    usleep(10);

    for(int i = 0; i < height; i++) {
        tasks[i] = i;
    }
 
    for(int i = 0; i < height; i++) {	
        threadpool_add(pool, &thread2thread, &tasks[i], 0);
    }   

    while(left_n > 0) {
        usleep(10);
    }

        threadpool_free(pool);

    pthread_mutex_destroy(&lock);

    /* draw and cleanup */
    write_png(filename, iters, width, height, image);
    free(image);
}
