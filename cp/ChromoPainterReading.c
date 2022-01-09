
#include "ChromoPainterReading.h"

int reading(char **st, char *format, void *res)
{
    int i;
    char *rs;
    rs = *st;
    for(i = 0; isspace(rs[i]); i++) ;
    if (!rs[i]) return 0;
    for(; !isspace(rs[i]); i++) ;
    if (rs[i]) rs[i++] = 0;
    if (!sscanf(*st, format, res)) return 0;
    *st += i;
    return 1;
}
