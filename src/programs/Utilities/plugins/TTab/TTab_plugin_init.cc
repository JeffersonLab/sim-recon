// $Id$
//
//    File: TTab_plugin_init.cc
// Created: Tue Jan  6 09:57:14 EST 2015
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#include <stdlib.h>
#include <iostream>
using namespace std;

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(void *arg){

	cout << endl;
	cout << endl;
	cout << "=== NOTICE: The TTab plugin has been converted into the TTAB ===" << endl;
	cout << "===         statically linked library. You should no longer  ===" << endl;
	cout << "===         specify it as a plugin. Exiting now to make sure ===" << endl;
	cout << "===         that you see this message.                       ===" << endl;
	cout << endl;
	cout << endl;
	exit(0);
}
} // "C"

