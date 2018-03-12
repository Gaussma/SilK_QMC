/* rand.c: Provide a Fortran interface of some gsl-random number functions
           for the kink-code.
   author: Frank Löffler, Nov 2012
*/

#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

unsigned int nstreams = 0;
gsl_rng **streams;

/* Initialize 'ns' many streams
 * (actually, initialize one more to seed the others ) */
int rand_init(unsigned int ns)
{
    const gsl_rng_type * T;
    T = gsl_rng_default;

    printf("Initializing %d streams\n", ns);
    nstreams = ns+1;

    gsl_rng_env_setup();

    // Allocate memory for the streams
    streams = (gsl_rng**) malloc(nstreams*sizeof(gsl_rng*));
    if (!streams)
        return 0;

    // Allocate all streams
    for (unsigned int i=0; i<nstreams; i++)
    {
        streams[i] = gsl_rng_alloc(T);
        // If we can't allocate, free what we have so far and return
        if (!streams[i])
        {
            for (unsigned int j=0; j<i; j++)
                gsl_rng_free(streams[j]);
            return 0;
        }
    }

    // Warm-up seed generator.
    // Not sure whether this helps, but it shouldn't hurt either.
    for (unsigned int i=1; i<10000; i++)
        gsl_rng_get(streams[0]);

    // Initialize all streams, nr. 0 is only used to seed the others
    for (unsigned int i=1; i<nstreams; i++)
    {
        gsl_rng_set(streams[i], gsl_rng_get(streams[0]));
    }
    return 1;
}

void rand_free()
{
    for (unsigned int i=0; i<nstreams; i++)
        gsl_rng_free (streams[i]);
    nstreams = 0;
}

unsigned long int rand_i(int s)
{
    assert(s>=0);
    assert(s<=nstreams);
    return gsl_rng_get(streams[s]);
}
double rand_d(int s)
{
    assert(s>=0);
    assert(s<=nstreams);
    return gsl_rng_uniform(streams[s]);
}

/* Fortran interfaces */
void rand_init_(int *ns, int *ierr) { *ierr = rand_init(*ns);  return; }
void rand_free_()                   {         rand_free();     return; }
int     rand_i_(int *s)             {  return rand_i(*s); }
double  rand_d_(int *s)             {  return rand_d(*s); }

void rand_init__(int *ns, int *ierr) { *ierr = rand_init(*ns);  return; }
void rand_free__()                   {         rand_free();     return; }
int     rand_i__(int *s)             {  return rand_i(*s); }
double  rand_d__(int *s)             {  return rand_d(*s); }

