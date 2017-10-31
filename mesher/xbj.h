/*
XBJV0001
起始8字节必须满足

int32 number of meshes NM
(int32 length of mesh name LN
char * LN + null 补齐4字节对齐
byte nbytes of double
byte nbytes of int
null
null
int32 number of verticis NV
float64 * 3 * NV
int32 number of faces NF
int32 * 3 * NF, 索引从1开始
)*NM
*/

#define MAX_XBJ_SIZE (0x800000 - 12)