/*
 *  hddm-cpp :	tool that reads in a HDDM document (Hall D Data Model)
 *		and writes a c++ header file that embodies the model in
 *		c++ classes.  It also generates input/output member
 *		functions to translate the model between the memory
 *		representation and a default binary representation that
 *		is suitable for passing over a pipe or storing on disk.
 *
 *  Original version - Richard Jones, May 25 2001.
 *
 *
 *  Programmer's Notes:
 *  -------------------
 * This translator has yet to be implemented, but it should be simple
 * to do.  Just take the hddm-c.cpp code and make the following changes.
 *
 * 1. In the header file replace the make_x_yyy function declarations
 *    with class definitions named x_yyy, using the signature of the
 *    make_x_yyy for the constructor, and the default destructor.
 *
 * 2. Implement the data members of the c structures as public data 
 *    members of the parent class, for efficient access.
 *
 * 3. The pointers become pointers to class instance.
 *
 * 4. The i/o functions are implemented as member functions of the
 *    by top-level class x_HDDM.
 *
 * 5. The pop-stack helper functions are implemented in an auxilliary
 *    class popstack.
 *
 */
