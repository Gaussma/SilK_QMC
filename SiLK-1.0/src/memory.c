#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>

size_t getRSS()
{
  long rss = 0L;
  FILE* fp;
  if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
    return rss;
  fscanf( fp, "%*s%ld", &rss );
  fclose( fp );
  return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);
}

long getrss_() { return (long)getRSS(); }
long getrss__() { return (long)getRSS(); }
