/*
XBJV0001
��ʼ8�ֽڱ�������

int32 number of meshes NM
(int32 length of mesh name LN
char * LN + null ����4�ֽڶ���
byte nbytes of double
byte nbytes of int
null
null
int32 number of verticis NV
float64 * 3 * NV
int32 number of faces NF
int32 * 3 * NF, ������1��ʼ
)*NM
*/

#define MAX_XBJ_SIZE (0x800000 - 12)