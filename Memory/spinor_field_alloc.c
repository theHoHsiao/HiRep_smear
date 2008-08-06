#include <stdlib.h>
#include "suN.h"
#include "linear_algebra.h"
#include "error.h"
#include "memory.h"
#include "global.h"

suNf_spinor* alloc_spinor_field_f(unsigned int n)
{
   suNf_spinor *p;
	 unsigned int len;

	 get_spinor_len(&len);
   p=amalloc(n*len*sizeof(suNf_spinor),ALIGN);
   error(p==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
         "Could not allocate memory space for the spinor field");
   
   return p;
}

suNf_spinor_flt* alloc_spinor_field_f_flt(unsigned int n)
{
   suNf_spinor_flt *p;
	 unsigned int len;

	 get_spinor_len(&len);
   p=amalloc(n*len*sizeof(suNf_spinor_flt),ALIGN);
   error(p==NULL,1,"alloc_spinor_field_f [field_alloc.c]",
         "Could not allocate memory space for the spinor field");
   
   return p;
}


