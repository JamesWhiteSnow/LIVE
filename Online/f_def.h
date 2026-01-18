/* this file contains necessary definitions for the file simulator */
#ifndef _F_DEF
#define _F_DEF
//----------------------------------------------------------
class BNode;

#define BFHEAD_LENGTH (sizeof(int) * 2) // file header size
#define TRUE 1
#define FALSE 0
#define SEEK_CUR 1
#define SEEK_SET 0
#define SEEK_END 2
#define FLOATZERO 1e-20
#define B_MAXINT 99999999

enum BINSRT
{
    B_NRML,
    B_OVRFLW
};
enum BDEL
{
    B_NONE,
    B_UNDRFLW,
    B_NOTFOUND
};
// #define min(a, b) (((a) < (b)) ? (a) : (b))
// #define max(a, b) (((a) > (b)) ? (a) : (b))

typedef char Block[];

//----------------------------------------------------------
void error(char *msg, bool _exit);
float traverse(BNode *_bn);
//----------------------------------------------------------
#endif