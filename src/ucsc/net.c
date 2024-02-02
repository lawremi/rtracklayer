
/* net.c some stuff to wrap around net communications. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */
#define CURL_STATICLIB
#include <curl/curl.h>
#include "common.h"
#include "errAbort.h"
#include "net.h"
#include "hash.h"
#include "cheapcgi.h"

time_t header_get_last_modified(CURL *curl) {
    curl_off_t last_modified;
    CURLcode status = curl_easy_getinfo(curl, CURLINFO_FILETIME_T, &last_modified);

    if ((CURLE_OK == status) && (last_modified >= 0)) {
        struct tm *utc_tm_info = gmtime(&last_modified);
        time_t utc_time_t = mktime(utc_tm_info);
        return utc_time_t;
    }

    // could not retrieve time
    if (last_modified == -1)
        return 0;

    errAbort("curl_easy_getinfo() failed: %s\n", curl_easy_strerror(status));
}

long long header_get_content_length(CURL *curl) {
    curl_off_t content_length;
    CURLcode status = curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T,
                                        &content_length);

    // could not retrieve length
    if (content_length == -1)
        return 0;

    if (CURLE_OK == status)
        return content_length;

    errAbort("curl_easy_getinfo() failed: %s\n", curl_easy_strerror(status));
}

CURL *wrapped_curl_init() {
    CURLcode status = curl_global_init(CURL_GLOBAL_DEFAULT);
    if (status != 0)
        errAbort("curl_global_init() failed: %s\n", curl_easy_strerror(status));

    CURL *curl = curl_easy_init();
    if (curl == NULL)
        errAbort("curl_easy_init() failed\n");

    return curl;
}

void wrapped_curl_cleanup(CURL *curl) {
    curl_easy_cleanup(curl);
    curl_global_cleanup();
}

CURLcode wrapped_curl_perform(CURL *curl) {
    CURLcode status = curl_easy_perform(curl);
    if (CURLE_OK != status)
        errAbort("curl_easy_perform() failed: %s\n", curl_easy_strerror(status));
    return status;
}

size_t write_callback(void *buffer, size_t size, size_t nitems, void *userdata) {
    return size * nitems;
}

enum HTTP_METHODS {
    GET,
    HEAD
};

CURLcode wrapped_curl_request(CURL *curl, enum HTTP_METHODS type) {
    if (type == HEAD) {
        curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
        curl_easy_setopt(curl, CURLOPT_HEADER, 1L);
    }
    curl_easy_setopt(curl, CURLOPT_FILETIME, 1L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    return wrapped_curl_perform(curl);
}

boolean netGetFtpInfo(char *url, long long *retSize, time_t *retTime) {
    struct netParsedUrl npu;
    netParseUrl(url, &npu);
    if (!sameString(npu.protocol, "ftp"))
        errAbort("netGetFtpInfo: url (%s) is not for ftp.", url);

    // TODO maybe remove this workaround where udc cache wants info on URL "/" ?
    if (sameString(npu.file,"/")) {
        *retSize = 0;
        *retTime = time(NULL);
        return TRUE;
    }

    CURL *curl = wrapped_curl_init();
    curl_easy_setopt(curl, CURLOPT_URL, url);
    wrapped_curl_request(curl, HEAD);

    *retTime = header_get_last_modified(curl);
    *retSize = header_get_content_length(curl);

    wrapped_curl_cleanup(curl);
    return TRUE;
}

size_t header_callback(void *buffer, size_t size, size_t nitems, void *userdata) {
    size_t realsize = size * nitems;
    struct hash **header = (struct hash**)userdata;

    char *line = strtok((char*) buffer, "\n");
    if (line != NULL) {
        char *colon = memchr(line, ':', strlen(line));
        if (colon != NULL) {
            *colon = '\0';
            hashAdd(*header, strUpper(line), cloneString(colon+1));
        }
    }
    return realsize;
}

int netUrlHead(char *url, struct hash *hash) {
    CURL* curl = wrapped_curl_init();
    curl_easy_setopt(curl, CURLOPT_URL, url);

    curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, header_callback);
    curl_easy_setopt(curl, CURLOPT_HEADERDATA, &hash);

    CURLcode status = wrapped_curl_request(curl, HEAD);
    wrapped_curl_cleanup(curl);
    return status;
}

int netUrlOpenSockets(char *url, int *retCtrlSocket) {
    if (stringIn("://", url) == NULL)
        return open(url, O_RDONLY);
    else {
        CURL* curl = wrapped_curl_init();
        curl_easy_setopt(curl, CURLOPT_URL, url);

        if (startsWith("http://", url) || startsWith("https://", url)) {
            int sockfd;
            curl_easy_setopt(curl, CURLOPT_OPENSOCKETDATA, &sockfd);
            wrapped_curl_request(curl, GET);
            wrapped_curl_cleanup(curl);
            return sockfd;
        } else if (startsWith("ftp://", url)) {
            curl_socket_t ctrlSocket;
            CURLcode status = wrapped_curl_request(curl, GET);
            curl_easy_getinfo(curl, CURLINFO_ACTIVESOCKET, &ctrlSocket);

            if (retCtrlSocket != NULL)
                *retCtrlSocket = ctrlSocket;

            wrapped_curl_cleanup(curl);
            return ctrlSocket;
        } else {
            errAbort("Sorry, can only netUrlOpen http, https and ftp currently, not '%s'", url);
        }
    }
}

int netUrlOpen(char *url) {
    return netUrlOpenSockets(url, NULL);
}

int netUrlFakeHeadByGet(char *url, struct hash *hash) {
    CURL* curl = wrapped_curl_init();
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_RANGE, "0-0");

    curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, header_callback);
    curl_easy_setopt(curl, CURLOPT_HEADERDATA, &hash);

    CURLcode status = wrapped_curl_request(curl, GET);
    wrapped_curl_cleanup(curl);
    return status;
}

boolean netSkipHttpHeaderLinesHandlingRedirect(int sd, char *url, int *redirectedSd, char **redirectedUrl) {
    char *effectiveUrl;
    curl_socket_t nsd;

    CURL* curl = wrapped_curl_init();
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

    wrapped_curl_request(curl, GET);

    curl_easy_getinfo(curl, CURLINFO_EFFECTIVE_URL, &effectiveUrl);
    curl_easy_getinfo(curl, CURLINFO_ACTIVESOCKET, &nsd);
    if (sd != nsd) {
        close(sd);
        *redirectedSd = nsd;
    }
    if (!sameString(effectiveUrl, *redirectedUrl))
        *redirectedUrl = cloneString(effectiveUrl);
    wrapped_curl_cleanup(curl);
    return TRUE;
}

// UCSC Existing Code

boolean hasProtocol(char *urlOrPath) {
    return stringIn("://", urlOrPath) != NULL;
}

static void parseByteRange(const char *url, ssize_t *rangeStart, ssize_t *rangeEnd, boolean terminateAtByteRange)

/* parse the byte range information from url */
{
char *x;
/* default to no byte range specified */
*rangeStart = -1;
*rangeEnd = -1;
x = strrchr(url, ';');
if (x)
    {
    if (startsWith(";byterange=", x))
	{
	char *y=strchr(x, '=');
	++y;
	char *z=strchr(y, '-');
	if (z)
	    {
	    ++z;
	    if (terminateAtByteRange)
		*x = 0;
	    // TODO: use something better than atoll() ?
	    *rangeStart = atoll(y); 
	    if (z[0] != '\0')
		*rangeEnd = atoll(z);
	    }    
	}
    }

}

void cgiDecodeFull(char *in, char *out, int inLength)
/* Out will be a cgi-decoded version of in (no space from plus!).
 * Out will be a little shorter than in typically, and
 * can be the same buffer. */
{
  char c;
  int i;
  for (i=0; i<inLength;++i)
    {
      c = *in++;
      if (c == '%')
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

void netParseUrl(char *url, struct netParsedUrl *parsed)
/* Parse a URL into components.   A full URL is made up as so:
 *   http://user:password@hostName:port/file;byterange=0-499
 * User and password may be cgi-encoded.
 * This is set up so that the http:// and the port are optional. 
 */
{
char *s, *t, *u, *v, *w;
char buf[2024];

/* Make local copy of URL. */
if (strlen(url) >= sizeof(buf))
    errAbort("Url too long: '%s'", url);
strcpy(buf, url);
url = buf;

/* Find out protocol - default to http. */
s = trimSpaces(url);
s = stringIn("://", url);
if (s == NULL)
    {
    strcpy(parsed->protocol, "http");
    s = url;
    }
else
    {
    *s = 0;
    tolowers(url);
    safecpy(parsed->protocol, sizeof(parsed->protocol), url);
    s += 3;
    }

/* Split off file part. */
parsed->byteRangeStart = -1;  /* default to no byte range specified */
parsed->byteRangeEnd = -1;
u = strchr(s, '/');
if (u == NULL)
    strcpy(parsed->file, "/");
else
    {

    parseByteRange(u, &parsed->byteRangeStart, &parsed->byteRangeEnd, TRUE);

    if (sameWord(parsed->protocol,"http") ||
        sameWord(parsed->protocol,"https"))
	{
	// http servers expect the URL request to be URL-encoded already.
	/* need to encode spaces, but not ! other characters */
	char *t=replaceChars(u," ","%20");
	safecpy(parsed->file, sizeof(parsed->file), t);
	freeMem(t);
	}

    *u = 0; // terminate the host:port string

    if (sameWord(parsed->protocol,"ftp"))
	{
	++u; // that first slash is not considered part of the ftp path 
	// decode now because the FTP server does NOT expect URL-encoding.
	cgiDecodeFull(u,parsed->file,strlen(u));  // decodes %FF but not +
	}

    }

/* Split off user part */
v = strchr(s, '@');
if (v == NULL)
    {
    if (sameWord(parsed->protocol,"http") ||
        sameWord(parsed->protocol,"https"))
	{
	strcpy(parsed->user, "");
	strcpy(parsed->password, "");
	}
    if (sameWord(parsed->protocol,"ftp"))
	{
	strcpy(parsed->user, "anonymous");
	strcpy(parsed->password, "x@genome.ucsc.edu");
	}
    }
else
    {
    *v = 0;
    /* split off password part */
    w = strchr(s, ':');
    if (w == NULL)
	{
	safecpy(parsed->user, sizeof(parsed->user), s);
	strcpy(parsed->password, "");
	}
    else
	{
	*w = 0;
	safecpy(parsed->user, sizeof(parsed->user), s);
	safecpy(parsed->password, sizeof(parsed->password), w+1);
	}
    
    cgiDecode(parsed->user,parsed->user,strlen(parsed->user));
    cgiDecode(parsed->password,parsed->password,strlen(parsed->password));
    s = v+1;
    }


/* Save port if it's there.  If not default to 80. */
t = strchr(s, ':');
if (t == NULL)
    {
    if (sameWord(parsed->protocol,"http"))
	strcpy(parsed->port, "80");
    if (sameWord(parsed->protocol,"https"))
	strcpy(parsed->port, "443");
    if (sameWord(parsed->protocol,"ftp"))
	strcpy(parsed->port, "21");
    }
else
    {
    *t++ = 0;
    if (!isdigit(t[0]))
	errAbort("Non-numeric port name %s", t);
    safecpy(parsed->port, sizeof(parsed->port), t);
    }

/* What's left is the host. */
safecpy(parsed->host, sizeof(parsed->host), s);
}


char *urlFromNetParsedUrl(struct netParsedUrl *npu)
/* Build URL from netParsedUrl structure */
{
struct dyString *dy = newDyString(512);

dyStringAppend(dy, npu->protocol);
dyStringAppend(dy, "://");
if (npu->user[0] != 0)
    {
    char *encUser = cgiEncode(npu->user);
    dyStringAppend(dy, encUser);
    freeMem(encUser);
    if (npu->password[0] != 0)
	{
	dyStringAppend(dy, ":");
	char *encPassword = cgiEncode(npu->password);
	dyStringAppend(dy, encPassword);
	freeMem(encPassword);
	}
    dyStringAppend(dy, "@");
    }
dyStringAppend(dy, npu->host);
/* do not include port if it is the default */
if (!(
 (sameString(npu->protocol, "ftp"  ) && sameString("21", npu->port)) ||
 (sameString(npu->protocol, "http" ) && sameString("80", npu->port)) ||
 (sameString(npu->protocol, "https") && sameString("443",npu->port))
    ))
    {
    dyStringAppend(dy, ":");
    dyStringAppend(dy, npu->port);
    }
dyStringAppend(dy, npu->file);
if (npu->byteRangeStart != -1)
    {
    dyStringPrintf(dy, ";byterange=%lld-", (long long)npu->byteRangeStart);
    if (npu->byteRangeEnd != -1)
	dyStringPrintf(dy, "%lld", (long long)npu->byteRangeEnd);
    }

/* Clean up and return handle. */
return dyStringCannibalize(&dy);
}

char *transferParamsToRedirectedUrl(char *url, char *newUrl)
/* Transfer password, byteRange, and any other parameters from url to newUrl and return result.
 * freeMem result. */
{
    struct netParsedUrl npu, newNpu;
/* Parse the old URL to make parts available for graft onto the redirected url. */
/* This makes redirection work with byterange urls and user:password@ */
    netParseUrl(url, &npu);
    netParseUrl(newUrl, &newNpu);
    if (npu.byteRangeStart != -1)
    {
        newNpu.byteRangeStart = npu.byteRangeStart;
        newNpu.byteRangeEnd = npu.byteRangeEnd;
    }
    if ((npu.user[0] != 0) && (newNpu.user[0] == 0))
    {
        safecpy(newNpu.user,     sizeof newNpu.user,     npu.user);
        safecpy(newNpu.password, sizeof newNpu.password, npu.password);
    }
    return urlFromNetParsedUrl(&newNpu);
}
