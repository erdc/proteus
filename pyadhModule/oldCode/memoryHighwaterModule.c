#if _LINUX
#include <malloc.h>

size_t hiwm (void) {
/* info.arena - number of bytes allocated
* info.hblkhd - size of the mmap'ed space
* info.uordblks - number of bytes used (?)
*/
struct mallinfo info = mallinfo();
size_t s = (size_t) info.arena + (size_t) info.hblkhd;
return (s);
}

#elif _MAXOSX || _UNIX
#include <unistd.h>

size_t hiwm (void) {
size_t s = (size_t) sbrk(0);
return (s);
}

#elif _WINDOWS
size_t hiwm (void) {
size_t s = (size_t) 0; /* ??? */
return (s);
}

#endif

static PyObject* me(PyObject* self,
                                                          PyObject* args)
