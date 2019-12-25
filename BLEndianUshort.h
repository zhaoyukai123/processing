#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#define M_PI 3.14159265358979323846
short BLEndianUshort(short value)
{
    return ((value & 0x00FF) << 8 ) | ((value & 0xFF00) >> 8);
}