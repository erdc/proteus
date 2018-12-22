#ifndef PROTEUS_SUPERLU_PROTO_H
#define PROTEUS_SUPERLU_PROTO_H

#ifdef __cplusplus
extern "C"
{
#endif

extern struct SuperLUStat_t;
extern struct SuperMatrix;
extern enum trans_t;
extern void dgstrs(trans_t, SuperMatrix *, SuperMatrix *, int *, int*, SuperMatrix *, SuperLUStat_t*, int *);

#ifdef __cplusplus
}
#endif
#endif
