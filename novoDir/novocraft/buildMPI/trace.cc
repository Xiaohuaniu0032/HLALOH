/* Obtain a backtrace and print it to stdout. */
#include <stdio.h>
#include <stdlib.h>
#if defined(__SUNPRO_CC)
    void
    print_trace (void)
    {
    }
#elif defined(CRAY_MPICH2)
    void
    print_trace (void)
    {
    }
#elif defined(TCMALLOC)
    #define UNW_LOCAL_ONLY
    #include <libunwind.h>
    void
    print_trace (void)
    {
        fprintf(stderr, " Stack  Dump ...\n");
        unw_cursor_t cursor; unw_context_t uc;
        unw_word_t ip, sp;

        unw_getcontext(&uc);
        unw_init_local(&cursor, &uc);
        do {
            if(unw_get_reg(&cursor, UNW_REG_IP, &ip) < 0) break;
            if(unw_get_reg(&cursor, UNW_REG_SP, &sp) < 0) break;
            fprintf (stderr, "ip = 0x%16lx, sp = 0x%16lx\n", (long) ip, (long) sp);
        } while (unw_step(&cursor) > 0);
    }
#else
    #include <execinfo.h>
    #ifndef __APPLE__
        #include <malloc.h>
    #endif
    void
    print_trace (void)
    {
        void *array[100];
        size_t size;
        char **strings;
        size_t i;

        size = backtrace (array, 100);
        strings = backtrace_symbols (array, size);

        fprintf (stderr, "Obtained %zd stack frames.\n", size);

        for(i = 0; i < size; i++)
            fprintf (stderr, "%s\n", strings[i]);

        free (strings);
    }
#endif

//#include <openssl/md5.h>
#include <openssl/evp.h>
#include <openssl/engine.h>

void cleanCrypto() {
    ENGINE_cleanup();
    ERR_free_strings();
    #ifndef NOEVP_cleanup
    EVP_cleanup();
    #endif
    CRYPTO_cleanup_all_ex_data();
}

#include <sys/utsname.h>
static struct utsname unamebuf;
char * nodeName()
{
        uname(&unamebuf);
        return unamebuf.nodename;
}

