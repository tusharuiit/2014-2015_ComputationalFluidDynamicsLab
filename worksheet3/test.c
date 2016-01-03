#include "helper.c"
int main(int argc, char **argv){
int** picture;
char* name;
name = malloc(strlen(argv[0]));
picture = read_pgm(argv[1]);
printf("Value of picture: %d\n", picture[490][1]);
free(name);
return 0;
}
