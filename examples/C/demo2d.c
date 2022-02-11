/* simple 2d demo of baobzi. to compile & run:

# here we assume baobzi installed in $HOME/local
export LIBRARY_PATH=$LIBRARY_PATH:$HOME/local/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/lib
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HOME/local/include

gcc demo2d.c -o demo2d -lbaobzi -lm && time OMP_NUM_THREADS=1 ./demo2d

(here math libs needed for our user function)

 */
#include <baobzi.h>
#include <stdio.h>
#include <math.h>


double f(const double *x, const void *data)
/* the user func to approximate: must have arguments (x,data)
   
   x = 2-element ordinate array, ie x- then y-coords
   data = user-available parameters array (here unused)
*/
{
  //  return x[0] * x[1];
  return sin(3.0 + 4*x[0] + 7*x[1]);   // plane wave
}


int main(int argc, char *argv[]) {
  
  baobzi_input_t input = {
    .func = f,                 // ptr to func
    .data = NULL,
    .dim = 2,
    .order = 16,                // high order. (should auto-choose given tol?)
    .tol = 1E-10
  };

  // define the rectangle to approximate over
  const double hl[2] = {1.0, 2.0};  // half-length and half-width
  const double center[2] = {0.0, 0.0};
  baobzi_t b = baobzi_init(&input, center, hl);    // b is an object
  
  // target pt
  const double x[2] = {1./3, 1./5};   // generic point avoids tree box edges

  double fe = f(x,NULL);              // plain func eval (2nd arg unused)
  double fb = baobzi_eval(b, x);      // the approximant at same target pt
  
  printf("true val  \t%.14g\nbaobzi eval\t%.14g \tabs err  %.3g (compare tol=%.3g)\n", fe,fb,fabs(fe-fb),input.tol);

  if (1) {     // demo save and load capability
    baobzi_save(b, "b.baobzi");
    b = baobzi_free(b);            // why b=... here?
    b = baobzi_restore("b.baobzi");
    fb = baobzi_eval(b, x);
    printf("restored eval\t%.14g\n", fb);
  }
  
  return 0;
}
