// $Id$

#include <stdlib.h>
#include <iostream>
using namespace std;

#include "DFactory_base.h"

//-------------
// printheader
//-------------
void DFactory_base::printheader(const char *header)
{
	/// Print the string given, along with a separator(----)
	/// and the factory name to signify the start of this factory's data
	
	// Print the header with separator
	cout<<dataClassName()<<endl;
	cout<<"---------------------------------------"<<endl;
	cout<<header<<endl;
	cout<<endl;

	// Find and record the column positions (just look for colons)
	char *c = (char*)header;
	_icol = 0;
	while(c = strstr(c,":")){
		_columns[_icol++] = (int)((unsigned long)c - (unsigned long)header);
		if(_icol>=99)break;
		c++;
	}
	for(int i=_icol;i<100;i++)_columns[i] = 100;
	_table = "";
}

//-------------
// printnewrow
//-------------
void DFactory_base::printnewrow(void)
{
	/// Initialize internal buffer in preparation of printing a new row.
	/// Call this before calling printcol().
	_row = string(79,' ');
	_icol = 0;
}

//-------------
// printnewrow
//-------------
void DFactory_base::printrow(void)
{
	/// Print the row to the screen. Make a call to printcol() for every
	/// column before calling this.
	_table += _row + "\n";
}

//-------------
// printcol
//-------------
void DFactory_base::printcol(const char *str)
{
	/// Print a formatted value to "str". Used by Print()
	_row.replace(_columns[_icol++]-strlen(str), strlen(str), str);
}
