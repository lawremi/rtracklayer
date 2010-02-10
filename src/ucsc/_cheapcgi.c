/* rtracklayer took out this little utility from cheapcgi.c */

#include "common.h"

void cgiDecode(char *in, char *out, int inLength)
/* Decode from cgi pluses-for-spaces format to normal.
 * Out will be a little shorter than in typically, and
 * can be the same buffer. */
{
  char c;
  int i;
  for (i=0; i<inLength;++i)
    {
      c = *in++;
      if (c == '+')
	*out++ = ' ';
      else if (c == '%')
	{
          int code;
          if (sscanf(in, "%2x", &code) != 1)
	    code = '?';
          in += 2;
          i += 2;
          *out++ = code;
	}
      else
	*out++ = c;
    }
  *out++ = 0;
}
