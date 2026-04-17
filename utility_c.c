#include <stdio.h>

int copy_file(char *old_filename, char *new_filename) {
    char    buff[BUFSIZ];
    FILE    *ptr_old, *ptr_new;
    size_t  n;

    ptr_old = fopen(old_filename, "rb");
    ptr_new = fopen(new_filename, "wb");
    while ( (n = fread(buff,1,BUFSIZ,ptr_old)) != 0 ) {
        fwrite(buff, 1, n, ptr_new);
    }
    fclose(ptr_old);
    fclose(ptr_new);
    return 0;
}

/*
int main(void){
   char filename_src[101], filename_dest[101];
   printf("\nSource file: ");
   gets(filename_src);

   printf("\nDestination filename: ");
   gets(filename_dest);

   if(c_copy_file(filename_src, filename_dest) == 0)
      printf("Copy Successful\n");
   else
      fprintf(stderr, "Error during copy!");
}
*/
