
/* net.c - network wrappers using libcurl */

#define CURL_STATICLIB
#include <curl/curl.h>
#include "common.h"
#include "errAbort.h"
#include "net.h"
#include "hash.h"
#include "cheapcgi.h"

// Struct to hold data for libcurl write callback
struct webReader
    {
    char *buffer;
    size_t size;
    size_t received;
    };

// libcurl write callback
static size_t webWriteCallback(void *ptr, size_t size, size_t nmemb, void *userdata)
{
    struct webReader *reader = (struct webReader *)userdata;
    size_t bytes = size * nmemb;

    // prevent buffer overflow
    size_t canWrite = 0;
    if (reader->received < reader->size)
        canWrite = reader->size - reader->received;
    
    if (bytes > canWrite)
        bytes = canWrite;

    if (bytes > 0)
        {
        memcpy(reader->buffer + reader->received, ptr, bytes);
        reader->received += bytes;
        }

    return size * nmemb; // Must tell curl we processed all its data
}

// Struct to hold header data for libcurl header callback
struct headerData
    {
    struct hash *hash;
    char lowercased[256]; // buffer for lowercasing header keys
    };

// libcurl header callback
static size_t webHeaderCallback(char *buffer, size_t size, size_t nitems, void *userdata)
{
    struct headerData *header = (struct headerData *)userdata;
    size_t len = size * nitems;

    // Find colon
    char *colon = strchr(buffer, ':');
    if (colon)
        {
        int keyLen = colon - buffer;
        if (keyLen < ArraySize(header->lowercased))
            {
            // key is before colon, value is after. Trim whitespace.
            memcpy(header->lowercased, buffer, keyLen);
            header->lowercased[keyLen] = 0;
            tolowers(header->lowercased);
            char *key = trimSpaces(header->lowercased);

            char *value = colon + 1;
            value = trimSpaces(value);
            
            hashAdd(header->hash, key, cloneString(value));
            }
        }
    return len;
}

static void netGlobalInit()
{
    static bool initialized = FALSE;
    if (!initialized)
    {
        if (curl_global_init(CURL_GLOBAL_DEFAULT) != 0)
            errAbort("curl_global_init() failed");
        initialized = TRUE;
    }
}

// Initialize a curl handle with common options
static CURL *netInitCurl(char *url)
{
    netGlobalInit();
    CURL *curl = curl_easy_init();
    if (curl == NULL)
        errAbort("curl_easy_init() failed");
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L); // Follow redirects
    curl_easy_setopt(curl, CURLOPT_MAXREDIRS, 5L);
    // It's good practice to set a user-agent
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "ucsc-rtracklayer-client/1.0");
    // enable TCP keep-alive
    curl_easy_setopt(curl, CURLOPT_TCP_KEEPALIVE, 1L);
    return curl;
}

// Get info about a URL using HEAD or a ranged GET.
// Returns HTTP status code. 0 on other curl errors.
int netUrlHead(char *url, struct hash *hash, char **effectiveUrl)
{
    CURL *curl = netInitCurl(url);
    long http_code = 0;
    
    struct headerData header = {hash, ""};

    // First try with HEAD
    curl_easy_setopt(curl, CURLOPT_NOBODY, 1L);
    curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, webHeaderCallback);
    curl_easy_setopt(curl, CURLOPT_HEADERDATA, &header);

    CURLcode res = curl_easy_perform(curl);
    if (res == CURLE_OK)
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
    
    // Some servers don't like HEAD. If it fails, or doesn't give content-length,
    // try a 0-byte GET range request.
    if (res != CURLE_OK || http_code == 403 || (http_code == 200 && hashFindVal(hash, "content-length") == NULL))
        {
        // Reset options for GET
        curl_easy_reset(curl);
        curl_easy_setopt(curl, CURLOPT_URL, url);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl, CURLOPT_USERAGENT, "ucsc-rtracklayer-client/1.0");
        curl_easy_setopt(curl, CURLOPT_TCP_KEEPALIVE, 1L);
        curl_easy_setopt(curl, CURLOPT_HTTPGET, 1L);
        curl_easy_setopt(curl, CURLOPT_RANGE, "0-0");
        curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, webHeaderCallback);
        curl_easy_setopt(curl, CURLOPT_HEADERDATA, &header);

        struct webReader reader = {NULL, 0, 0}; // discard body
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, webWriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &reader);

        res = curl_easy_perform(curl);
        if (res == CURLE_OK)
            curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
        else
            http_code = 0; // curl error
        }

    if (http_code > 0 && effectiveUrl)
    {
        char *effUrl = NULL;
        curl_easy_getinfo(curl, CURLINFO_EFFECTIVE_URL, &effUrl);
        if (effUrl && !sameString(url, effUrl))
            *effectiveUrl = cloneString(effUrl);
        else
            *effectiveUrl = NULL;
    }

    curl_easy_cleanup(curl);
    return http_code;
}

// Fetch a range of bytes from a URL.
// Returns bytes read, or -1 on error.
int netUrlRead(char *url, bits64 offset, int size, void *buffer)
{
    if (size == 0) return 0;

    CURL *curl = netInitCurl(url);
    long http_code = 0;
    
    struct webReader reader = {buffer, size, 0};
    
    char range[128];
    safef(range, sizeof(range), "%lld-%lld", (long long)offset, (long long)(offset + size - 1));
    
    curl_easy_setopt(curl, CURLOPT_RANGE, range);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, webWriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &reader);
    
    CURLcode res = curl_easy_perform(curl);
    
    if (res != CURLE_OK)
        {
        warn("curl_easy_perform() failed for %s: %s", url, curl_easy_strerror(res));
        curl_easy_cleanup(curl);
        return -1;
        }

    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
    // 200 can happen for servers that don't support range requests.
    // The webWriteCallback will stop reading past the requested size.
    if (http_code != 206 && http_code != 200)
        {
        warn("HTTP error for %s: status %ld", url, http_code);
        }

    curl_easy_cleanup(curl);
    return reader.received;
}

// Get info about an FTP file.
boolean netGetFtpInfo(char *url, long long *retSize, time_t *retTime)
{
    CURL *curl = netInitCurl(url);
    CURLcode res;
    
    // To get size and time for FTP, we must not use NOBODY
    curl_easy_setopt(curl, CURLOPT_NOBODY, 0L);
    curl_easy_setopt(curl, CURLOPT_HEADER, 1L);
    curl_easy_setopt(curl, CURLOPT_FILETIME, 1L);
    // We don't want the file content, so we ask curl to stop after the header.
    // It will return CURLE_FTP_COULDNT_RETR_FILE, which is what we expect.
    struct webReader reader = {NULL, 0, 0};
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, webWriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &reader);
    
    res = curl_easy_perform(curl);
    if (res != CURLE_OK && res != CURLE_FTP_COULDNT_RETR_FILE)
        {
        warn("curl_easy_perform() failed for ftp info on %s: %s", url, curl_easy_strerror(res));
        curl_easy_cleanup(curl);
        return FALSE;
        }

    curl_off_t fileSize = -1;
    res = curl_easy_getinfo(curl, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &fileSize);
    if (res != CURLE_OK || fileSize < 0)
    {
        warn("Couldn't get size for FTP URL %s", url);
        fileSize = 0;
    }
    *retSize = fileSize;

    time_t fileTime = -1;
    #if LIBCURL_VERSION_NUM >= 0x073b00
        #define CI_FILETIME CURLINFO_FILETIME_T
    #else
        #define CI_FILETIME CURLINFO_FILETIME
    #endif
    res = curl_easy_getinfo(curl, CI_FILETIME, &fileTime);
    if (res != CURLE_OK || fileTime < 0)
    {
        warn("Couldn't get modification time for FTP URL %s", url);
        fileTime = 0; // 0 is a reasonable default
    }
    *retTime = fileTime;
    
    curl_easy_cleanup(curl);
    return TRUE;
}

boolean hasProtocol(char *urlOrPath)
{
    return stringIn("://", urlOrPath) != NULL;
}

time_t netParseDate(char *dateString)
/* Parse a date string (like from an HTTP header) into a time_t.
 * Returns -1 on error. */
{
    return curl_getdate(dateString, NULL);
}
