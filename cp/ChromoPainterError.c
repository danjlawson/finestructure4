#include "ChromoPainterError.h"

void stop_on_error(int val,int errormode,jmp_buf cp_error_buf) {
    if(errormode==0) {
        exit(val);
    }else if(errormode==1){
        printf("Unrecoverable error, cancel execution via the GUI and reconfigure!\n");
       //# Sleep(10000);
        getchar();
        exit(val);
    }else{
      if(cp_error_buf==NULL) exit(val);
      longjmp(cp_error_buf,1);
    }
}
