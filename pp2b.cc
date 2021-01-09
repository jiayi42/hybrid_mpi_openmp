#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <assert.h>
#include <png.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#define WORK_TAG    0           
#define DATA_TAG    1          
#define STOP_TAG    2
#include <unistd.h>
int* image;
int done=0;
/*
struct arg_struct {
    int arg1;
    int arg2;
    double arg3;
    double arg4;
    int arg5;
    int arg6;
    int arg7;
};
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  cond  = PTHREAD_COND_INITIALIZER;
void* mworker(void *ptr2)
{
    pthread_detach(pthread_self()); 
    struct arg_struct *args = (struct arg_struct *)ptr2;
    double x0,y0,x,y,temp,repeats,length_squared;
    int i;
    int* data_msg;
    int lower=args->arg1;
    int left=args->arg2;
    double v= args->arg3;
    double h= args->arg4;
    int row=args->arg5;
    int width=args->arg6;
    int iters=args->arg7;
    data_msg = new int[(width)];
    for(int k=0;k<=row;k++){
        #pragma omp  parallel  for  schedule(dynamic) 
        for ( i = 0; i < width; i++) {
               //y0 =((row) * v)+lower;
               y0 =((k) * v)+lower;
               x0 = ( i *h) + left;
               repeats=0;
               x=0;
               y=0;
               length_squared = 0;
               while (repeats <iters && length_squared <4) {
                temp= x*x- y*y+ x0;
                y =  (x * y)+(x*y) + y0;
                x = temp;
                length_squared = x * x + y * y;
                ++repeats;
                }
             data_msg[i]=repeats;
          //image[k*width+i]=repeats;
        }
      //pthread_mutex_lock(&mutex);  
        #pragma omp  parallel  for    
	for (int col = 0; col < width; ++col) {
          //image[row*width+col]=data_msg[col];
          image[k*width+col]=data_msg[col];
        }
	//memcpy(image+k*width+1, data_msg, width*sizeof(int));

      //pthread_mutex_unlock(&mutex);
    }
    delete[] data_msg;
    pthread_mutex_lock(&mutex);  
    done=0;
    pthread_mutex_unlock(&mutex);
    pthread_exit(NULL);
}
*/
int  next_row;
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
         //#pragma omp parallel for  schedule(dynamic) 
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
void master_pgm(int n_worker, int width, int height, double left, double right,double lower, double upper, int iters,int* image,int nthreads){

    //int  next_row;
    int *data_msg;
    int *mdata_msg;
    data_msg = new int[(width+1)];
    mdata_msg = new int[(width+1)];
    MPI_Status status;
    int tasks_not_done;
    int id;
    MPI_Request	rst;
    next_row = 0;          
    tasks_not_done = 0;    
/*    
    pthread_t some_thread;
    struct arg_struct args;
    args.arg1 = lower;
    args.arg2 = left;
    double v= ((upper-lower)/height);
    args.arg3 = v;
    double h= ((right - left) / width);
    args.arg4 = h;
    args.arg6 = width;
    args.arg7 = iters;
   


    //args.arg5 = next_row;
   // pthread_create(&some_thread, NULL, &mworker, (void *)&args); 
   // pthread_join(some_thread, NULL);
   // ++next_row;
    
    pthread_mutex_lock(&mutex);  
    done=1;
    pthread_mutex_unlock(&mutex);  
    args.arg5 = 0;//height;
    pthread_create(&some_thread, NULL, &mworker, (void *)&args); 
    //pthread_join(some_thread, NULL);
    //pthread_detach(some_thread);
    next_row=0;//height/20+1;
*/
  if(n_worker>0){
    for (int p = 0; p < n_worker; ++p) {
        MPI_Isend(&next_row, 1, MPI_INT, p+1, WORK_TAG, MPI_COMM_WORLD,&rst);
        //MPI_Send(&next_row, 1, MPI_INT, p+1, WORK_TAG, MPI_COMM_WORLD);//,&rst);
        ++next_row;
        ++tasks_not_done;
    }
    while (tasks_not_done > 0) {
    /*
    if ( done==0 && k==1){
    pthread_join(some_thread, NULL);
    k=0;
       if (next_row < height) {
     k=1;
    pthread_mutex_lock(&mutex);  
    done=1;
    pthread_mutex_unlock(&mutex);
    args.arg5 = next_row;
    pthread_create(&some_thread, NULL, &mworker, (void *)&args); 
    //pthread_join(some_thread, NULL);
     ++next_row;
       }
    }*/
        MPI_Recv(data_msg, width+1, MPI_INT, MPI_ANY_SOURCE,DATA_TAG, MPI_COMM_WORLD, &status);
	--tasks_not_done;
        id = status.MPI_SOURCE;

       if (next_row < height) {
            MPI_Isend(&next_row, 1, MPI_INT, id, WORK_TAG, MPI_COMM_WORLD,&rst);
            ++next_row;
            ++tasks_not_done;
 	}
        else {
            MPI_Isend(&next_row, 0, MPI_INT, id, STOP_TAG, MPI_COMM_WORLD,&rst);
        }

    //pthread_mutex_lock(&mutex);  
        #pragma omp  parallel  for  schedule(dynamic)  
	for (int col = 0; col < width; ++col) {
          image[data_msg[0]*width+col]=data_msg[col+1];
        }
    //pthread_mutex_unlock(&mutex);
    }
  }
 else{
   
   


        #pragma omp  parallel  for  schedule(dynamic)  
       for (int j = 0; j < height; ++j) {
          for (int i = 0; i < width; ++i) {
            double y0 = j * ((upper - lower) / height) + lower;
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
   
       } 
   }
  /*
  while(done==1)
  {
  usleep(10);
  }   
    pthread_join(some_thread, NULL);
 
    pthread_mutex_destroy(&mutex);  
    pthread_cond_destroy(&cond); 
  */  
    //pthread_join(some_thread, NULL);
    delete[] data_msg;
    delete[] mdata_msg;
}

void worker_pgm(int rank_id, int width, int height, double left, double right,double lower, double upper, int iters,int nthreads) {

    MPI_Status status;
    int the_row;
    int *data_msg;
    data_msg = new int[(width+1)];

    double x0,y0,x,y,temp,repeats,length_squared;
    int i;
    double v= ((upper-lower)/height);
    double h= ((right - left) / width);

   MPI_Request	rst;
    while ( ((MPI_Recv(&the_row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status)) == MPI_SUCCESS)){

           if (status.MPI_TAG != WORK_TAG)  break;
        data_msg[0] = the_row;
        #pragma omp  parallel  for  schedule(dynamic) num_threads(nthreads) 
        for ( i = 0; i < width; i++) {
            y0 =((the_row) * v)+lower;
            x0 = ( i *h) + left;
            repeats=0;
            x=0;
            y=0;
            length_squared = 0;
            while (repeats <iters && length_squared <4) {
                temp= x*x- y*y+ x0;
                y =  (x * y)+(x*y) + y0;
                x = temp;
                length_squared = x * x + y * y;
                ++repeats;
            }
             data_msg[i+1] =repeats;
        }
        
        MPI_Isend(data_msg, width+1, MPI_INT, 0, DATA_TAG,MPI_COMM_WORLD,&rst);

    }

    delete[] data_msg;
   
}


int main(int argc, char** argv) {
    /* detect how many CPUs are available */
      //MPI_Init(&argc,&argv);
    MPI_Init(NULL,NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rank_size;
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);

    /* argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    int iters = strtol(argv[2], 0, 10);
    double left = strtod(argv[3], 0);
    double right = strtod(argv[4], 0);
    double lower = strtod(argv[5], 0);
    double upper = strtod(argv[6], 0);
    int width = strtol(argv[7], 0, 10);
    int height = strtol(argv[8], 0, 10);
    
    //int* image;
    if (rank == 0){
      image = new int[width * height ];//seq version
      assert(image);
    }
      if (rank == 0) {
        master_pgm(rank_size-1, width, height,left,right,lower,upper, iters,image,CPU_COUNT(&cpu_set));
      }
      else {
        worker_pgm(rank, width, height,left,right,lower,upper, iters,CPU_COUNT(&cpu_set));
      }
     if(rank==0){
        write_png(filename, iters, width, height, image);
        delete[] image;
      }
      // MPI_Finalize();
}

