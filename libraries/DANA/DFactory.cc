// $Id$

#include <stdlib.h>

#include "DFactory.h"
#include "DEvent.h"


//-------------
// DFactory
//-------------
DFactory::DFactory(DEvent *my_devent, char *my_name, int rowsize)
{
	/// This is a base class that specific factories inherit from.
	/// my_devent will be kept and used to allow this factory to
	/// access other factories.

	event = my_devent;
	name = my_name;
	_data = new DContainer(NULL, rowsize, my_name);

	// Allow any factory to have its debug_level set via environment variable
	debug_level = 0;
	char envar[256];
	sprintf(envar, "DEBUG_%s", name);
	char *ptr = getenv(envar);
	if(ptr)debug_level = atoi(ptr);
}

//-------------
// ~DFactory
//-------------
DFactory::~DFactory()
{
	delete _data;
}

//-------------
// Get
//-------------
DContainer* DFactory::Get()
{
	/// Return the list of values for the type of data this factory
	/// generates. This is function checks first if the data already
	/// exists and then calls the GetData method if it doesn't.
	
	// If evnt_called is set, then just return the _data pointer
	if(evnt_called)return _data;
	
	// Make sure we're initialized
	if(!init_called){
		init();
		init_called = 1;
	}
	
	// Call brun routine if run number has changed or it's not been called
	if(event->runnumber!=brun_runnumber){
		if(brun_called && !erun_called){
			erun();
			erun_called = 1;
		}
		brun_called = 0;
	}
	if(!brun_called){
		brun(event->runnumber);
		brun_called = 1;
		erun_called = 0;
		brun_runnumber = event->runnumber;
	}
	
	// Call evnt routine to generate data
	evnt(event->eventnumber);
	evnt_called = 1;
	
	return _data;
}

//-------------
// printheader
//-------------
void DFactory::printheader(const char *header)
{
	/// Print the string given, along with a separator(----)
	/// and the factory name to signify the start of this factory's data
	
	// Print the header with separator
	cout<<name<<endl;
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
}

//-------------
// printnewrow
//-------------
void DFactory::printnewrow(void)
{
	/// Initialize internal buffer in preparation of printing a new row.
	/// Call this before calling printcol().
	memset(_str,' ',80);
	_str[79] = 0;
	_icol = 0;
}

//-------------
// printnewrow
//-------------
void DFactory::printrow(void)
{
	/// Print the row to the screen. Make a call to printcol() for every
	/// column before calling this.
	cout<<_str<<endl;
}

//-------------
// printcol
//-------------
void DFactory::printcol(const char *val)
{
	/// Print a formatted value to "str". Used by Print()
	strncpy(&_str[_columns[_icol++]-strlen(val)], val, strlen(val));
}

//-------------
// printcol
//-------------
void DFactory::printcol(const char *format, float val)
{
	/// Print a formatted value to "str". Used by Print()
	char num[32];
	sprintf(num, format, val);
	strncpy(&_str[_columns[_icol++]-strlen(num)], num, strlen(num));
}

//-------------
// printcol
//-------------
void DFactory::printcol(const char *format, double val)
{
	/// Print a formatted value to "str". Used by Print()
	printcol(format,(float)val);
}


//-------------
// printcol
//-------------
void DFactory::printcol(const char *format, int val)
{
	/// Print a formatted value to "str". Used by Print()
	char num[32];
	sprintf(num, format, val);
	strncpy(&_str[_columns[_icol++]-strlen(num)], num, strlen(num));
}

