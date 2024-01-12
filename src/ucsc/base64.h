/* Base64 encoding and decoding.
 * by Galt Barber */

#ifndef BASE64_H
#define BASE64_H

#define B64CHARS "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

char *base64Encode(char *input, size_t inplen);
/* Use base64 to encode a string.  Returns one long encoded
 * string which need to be freeMem'd. Note: big-endian algorithm.
 * For some applications you may need to break the base64 output
 * of this function into lines no longer than 76 chars.
 */

#endif /* BASE64_H */
