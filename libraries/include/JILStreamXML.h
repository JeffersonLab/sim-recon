
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#ifndef _JILSTREAMXML_H_
#define _JILSTREAMXML_H_

#include <JILStream.h>

class JILStreamXML: public JILStream {
   public:
	
	JILStreamXML(string filename="", string mode="w"){
		tabs[0] = 0;
		this->filename = filename;
		if(mode == "r")JILStreamInit(STREAM_INPUT);
		else if(mode == "w")JILStreamInit(STREAM_OUTPUT);
		else{
			cerr<<"Unknown mode for JILStreamXML \""<<mode<<"\"  should be \"r\" or \"w\""<<endl;
			exit(-1);
		}
		if(iotype == STREAM_OUTPUT){
			// Output to file
			if(filename != ""){
				xmlout = new ofstream(filename.c_str());
			}else{
				xmlout = &cout;
			}
			(*xmlout)<<"<JILStreamXML>"<<std::endl;
			(*xmlout)<<JILMyDictionary();
			strcat(tabs,"\t");
		}else{
			// Input from file
			if(filename != ""){
				xmlin = new ifstream(filename.c_str());
			}else{
				xmlin = &cin;
			}
			
			// Read in the XML dictionary
			GetDictionary(xmlin);
		}
	}

	~JILStreamXML(){
		ClearNamed();
		if(iotype == STREAM_OUTPUT){
			(*xmlout)<<"</JILStreamXML>"<<std::endl;
			delete xmlout;
		}else{
			delete xmlin;
		}
	};
	
	// XML doesn't like angle brackets <> so we convert them into
	// curly brackets here {}
	const char* JILtypeid2nameXML(const std::type_info *t){
		static string str;
		str = JILtypeid2name(t);
		for(unsigned int i=0;i<str.size();i++){
			if(str[i] == '>')str[i] = '}';
			if(str[i] == '<')str[i] = '{';
		}

		return str.c_str();
	}
	
	// Atomic types. All atomic types (plus strings) are
	// explicitly defined.
   JILStreamXML& operator<<(short i){(*xmlout)<<tabs<<"<short val=\""<<i<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(int i){(*xmlout)<<tabs<<"<int val=\""<<i<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(long i){(*xmlout)<<tabs<<"<long val=\""<<i<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(unsigned short i){(*xmlout)<<tabs<<"<ushort val=\""<<i<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(unsigned int i){(*xmlout)<<tabs<<"<uint val=\""<<i<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(unsigned long i){(*xmlout)<<tabs<<"<ulong val=\""<<i<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(float f){(*xmlout)<<tabs<<"<float val=\""<<f<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(double f){(*xmlout)<<tabs<<"<double val=\""<<f<<"\"/>"<<std::endl; return *this;}
   JILStreamXML& operator<<(std::string s){(*xmlout)<<tabs<<"<string val=\""<<s<<"\"/>"<<std::endl; return *this;}
	JILStreamXML& operator<<(const std::type_info *t){
		(*xmlout)<<tabs<<"<object type=\""<<JILtypeid2nameXML(t)<<"\">"<<std::endl;
		strcat(tabs,"\t");
		return *this;
	}

	// Catch stream manipulators
	JILStreamXML& operator<<(JILStreamManipulator_t m){
		switch(m){
			case END_NAMED:
				named_depth--;
				if(strlen(tabs)>0)tabs[strlen(tabs)-1] =0;
				if(names.size()>0){
					string s = *names.begin();
					names.pop_front();
					(*xmlout)<<tabs<<"</"<<s<<">"<<std::endl;
				}else{
					std::cerr<<"END_NAMED manipulator reached without matching name!"<<std::endl;
				}
				break;
			case END_OBJECT:
				type_depth--;
				if(strlen(tabs)>0)tabs[strlen(tabs)-1] =0;
				(*xmlout)<<tabs<<"</object>"<<std::endl;
				break;
			case END_VECTOR:
				vector_depth--;
				if(strlen(tabs)>0)tabs[strlen(tabs)-1] =0;
				(*xmlout)<<tabs<<"</vector>"<<std::endl;
				break;
			case END_LIST:
				list_depth--;
				if(strlen(tabs)>0)tabs[strlen(tabs)-1] =0;
				(*xmlout)<<tabs<<"</list>"<<std::endl;
				break;
			case END_ARRAY:
				array_depth--;
				if(strlen(tabs)>0)tabs[strlen(tabs)-1] =0;
				(*xmlout)<<tabs<<"</array>"<<std::endl;
				break;
			case END_POINTER:
				pointer_depth--;
				if(strlen(tabs)>0)tabs[strlen(tabs)-1] =0;
				(*xmlout)<<tabs<<"</pointer>"<<std::endl;
				break;
		}
		return *this;
	}
	
	JILStream& UnknownOut(const type_info* type){
		(*xmlout)<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
	}
	
	JILStream& UnknownIn(const type_info* type){
		std::cerr<<"Attempting to convert unknown object!! (type="<<type->name()<<")"<<endl;
		return *this;
	}

	// Allow user to insert object-like tags
	void StartNamedWrite(const char *name){
		names.push_front(string(name));
		(*xmlout)<<tabs<<"<"<<name<<">"<<std::endl;
		strcat(tabs, "\t");
	}

	// Called before all items of a vector are streamed
	void StartVectorWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*xmlout)<<tabs<<"<vector type=\""<<JILtypeid2nameXML(&t)<<"\">"<<std::endl;
		strcat(tabs, "\t");
	}

	// Called before all items of a list are streamed
	void StartListWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*xmlout)<<tabs<<"<list type=\""<<JILtypeid2nameXML(&t)<<"\">"<<std::endl;
		strcat(tabs, "\t");
	}

	// Called before all items of an array are streamed
	void StartArrayWrite(const type_info &t, unsigned int size, unsigned int bytes_per_item){
		(*xmlout)<<tabs<<"<array type=\""<<JILtypeid2nameXML(&t)<<"\">"<<std::endl;
		strcat(tabs, "\t");
	}

	// Called for all pointers. If it returns "true", then the thing being pointed
	// to is streamed. Otherwise, it is ignored and NO corresponding
	// END_POINTER manipulator will be streamed.
	bool StartPointerWrite(const type_info &t, void* ptr){
		(*xmlout)<<tabs<<"<pointer type=\""<<JILtypeid2name(&t)<<"\">"<<std::endl;
		strcat(tabs, "\t");
		return true;
	}
	
	bool GetNamed(const char *name){
		return true;
	}
	
	bool GetPointerFromStream(const type_info* type, void* &ptr){
		// returning false means we have handled  this object and
		// the deserializer routine should NOT try filling from the
		// stream.
		return false;
	}

	private:
		ostream *xmlout;
		istream *xmlin;
		char tabs[256];
		list<string> names;
};

#endif //_JILSTREAMXML_H_

