#ifdef __cplusplus
extern "C" {
#endif

#ifndef CHROMOPAINTERERROR_H
#define CHROMOPAINTERERROR_H

#include <stdlib.h>
#include <stdio.h>
#include <setjmp.h>
  
  //extern jmp_buf cp_error_buf;
  // extern static jmp_buf cp_error_buf;

  void stop_on_error(int val,int errormode,jmp_buf cp_error_buf);

#endif

#ifdef __cplusplus
}
#endif
