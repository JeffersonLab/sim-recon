#include "XString.hpp"
#include "XParsers.hpp"

void usage();
char* plural(const char* str);
char* simpleStructType(const char* tag);
char* listStructType(const char* tag);
void checkConsistency(DOMElement* el, int t);
void writeHeader(DOMElement* el);
int constructGroup(DOMElement* el);
void constructMakeFuncs();
void constructUnpackers();
void constructReadFunc(DOMElement* topEl);
void constructPackers();
void constructFlushFunc(DOMElement* el);
void writeMatcher();
void constructOpenFunc(DOMElement* el);
void constructInitFunc(DOMElement* el);
void constructCloseFunc(DOMElement* el);
void constructDocument(DOMElement* el);
