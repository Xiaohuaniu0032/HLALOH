#include <stdio.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "dummy"
#endif
#include "main.h"

int main(int argc, char *argv[])
{
	return maq_novo2maq(argc, argv);
}
