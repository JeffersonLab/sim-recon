/*
 * gpp: GEANT preprocessor for MCFAST database files
 *
 * Summary:
 * MCFAST is a MonteCarlo event generator used by American HEP laboratories.
 * GEANT is a MonteCarlo event generator from CERN used by HEP worldwide.
 * The purpose of gpp is to input a geometry database file (db) in MCFAST
 * format and generate the equivalent Fortran code for GEANT.  It has been
 * tested and used with GEANT v3.21.
 *
 * Author: 	Richard Jones
 * Institution: University of Connecticut
 * Language:	ANSI C++
 * Original:	v 1.0, Jan. 12 1999
 * Updated:
 *
 * Comments:
 * > Richard Jones, July 5 2001
 *	This package has been superseded by the hdds-geant translator
 *	that interfaces Geant3 to the HDDS xml geometry database.
 */

#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>

#define _REENTRANT 1
#include <string.h>

ifstream fdbase;
ifstream fsource;
ofstream fdest;

int thisLine;
char thisFile[100];
char sourceFile[] = "mcgeom.f";

class TemplateItem;
class TemplateList;

class TemplateAtom {
friend class TemplateItem;
private:
   TemplateAtom *fNext;
public:
   char *fType;
   char *fName;
   TemplateAtom();
   TemplateAtom(char *type, char *name);
   ~TemplateAtom();
   TemplateAtom *Next() { return fNext; }
};

class TemplateItem {
friend class TemplateList;
private:
   TemplateAtom *firstAtom;
   TemplateItem *fNext;
public:
   char *fName;
   char *fNumber;
   TemplateItem();
   TemplateItem(char *name, char *number);
   ~TemplateItem();
   void Append(TemplateAtom *atom);
   void Delete(char *name);
   TemplateAtom *Find(char *name);
   TemplateAtom *First() { return firstAtom; }
   TemplateItem *Next() { return fNext; }
};

class TemplateList {
private:
   TemplateItem *firstItem;
public:
   TemplateList();
   ~TemplateList();
   void Append(TemplateItem *item);
   void Delete(char *name);
   TemplateItem *Find(char *name);
   TemplateItem *First() { return firstItem; }
};

TemplateList table;


TemplateAtom::TemplateAtom()
{
   fType = fName = 0;
   fNext = 0;
}

TemplateAtom::TemplateAtom(char *type, char *name)
{
   fType = new char[strlen(type)+1];
   strcpy(fType,type);
   fName = new char[strlen(name)+1];
   strcpy(fName,name);
   fNext = 0;
}

TemplateAtom::~TemplateAtom()
{
   if (fNext) delete fNext;
   if (fType) delete []fType;
   if (fName) delete []fName;
   fType = fName = 0;
   fNext = 0;
}

TemplateItem::TemplateItem()
{
   fName = fNumber = 0;
   firstAtom = 0;
   fNext = 0;
}

TemplateItem::TemplateItem(char *name, char *number)
{
   fName = new char[strlen(name)+1];
   strcpy(fName,name);
   fNumber = new char[strlen(number)+1];
   strcpy(fNumber,number);
   firstAtom = 0;
   fNext = 0;
}

TemplateItem::~TemplateItem()
{
   if (firstAtom) delete firstAtom;
   if (fNext) delete fNext;
   if (fName) delete []fName;
   if (fNumber) delete []fNumber;
   fName = fNumber = 0;
   firstAtom = 0;
   fNext = 0;
}

void TemplateItem::Append(TemplateAtom *atom)
{
   TemplateAtom *p=firstAtom;
   if (p == NULL) {
      firstAtom = atom;
      return;
   }
   while (p->fNext != NULL) {
      p = p->fNext;
   }
   p->fNext = atom;
}

void TemplateItem::Delete(char *name)
{
   TemplateAtom **p = &firstAtom;
   while (*p != NULL) {
      if (strcasecmp((*p)->fName,name) == 0) break;
      p = &(*p)->fNext;
   }
   if (*p) {
      TemplateAtom *old = *p;
      *p = (*p)->fNext;
      old->fNext = 0;
      delete old;
   }
}

TemplateAtom *TemplateItem::Find(char *name)
{
   TemplateAtom *p=firstAtom;
   while (p != NULL) {
      if (strcasecmp(p->fName,name) == 0) break;
      p = p->fNext;
   }
   return p;
}

TemplateList::TemplateList()
{
   firstItem = 0;
}

TemplateList::~TemplateList()
{
   if (firstItem) delete firstItem;
}

void TemplateList::Append(TemplateItem *item)
{
   TemplateItem *p=firstItem;
   if (p == NULL) {
      firstItem = item;
      return;
   }
   while (p->fNext != NULL) {
      p = p->fNext;
   }
   p->fNext = item;
}

void TemplateList::Delete(char *name)
{
   TemplateItem **p = &firstItem;
   while (*p != NULL) {
      if (strcasecmp((*p)->fName,name) == 0) break;
      p = &(*p)->fNext;
   }
   if (*p) {
      TemplateItem *old = *p;
      *p = (*p)->fNext;
      old->fNext = 0;
      delete old;
   }
}

TemplateItem *TemplateList::Find(char *name)
{
   TemplateItem *p=firstItem;
   while (p != NULL) {
      if (strcasecmp(p->fName,name) == 0) break;
      p = p->fNext;
   }
   return p;
}

int preProcessFile(ifstream &fin);
int preProcessLine(char *line, ifstream &fin);

int preProcessFile(ifstream &fin)
{
   char line[250];
   char token[250];
   while (!fin.eof()) {
      fin.getline(line,250), ++thisLine;
      strcpy(token,line);
      strtok(token," ");
      if (strcasecmp(token,"end") == 0) break;
      preProcessLine(line,fin);
   }
   return 0;
}

int preProcessLine(char *line, ifstream &fin)
{
   if ((strlen(line) == 0) || (line[0] == '!'))	return 0;

   char* key=strtok(line," ");
   if (strcasecmp(key,"database") == 0) {
      char *fout=strtok(NULL," ");
      strcat(fout,".f");
      fdest.open(fout);
      if (!fdest) {
         cerr << "gpp error: unable to open output file ";
         cerr << fout << endl;
         cerr << "bombing out in line " << thisLine;
         cerr << " of file " << thisFile << endl;
         exit(1);
      }
      fdest << "      " << "program " << strtok(fout,".") << endl;
      fdest << "      " << "call makeGeometry" << endl;
      return 0;
   }

   else if (!fdest) {
      return 0;
   }

   else if (strcasecmp(key,"include") == 0) {
      char thisFilePushed[100];
      strncpy(thisFilePushed,thisFile,100);
      int thisLinePushed=thisLine;
      strncpy(thisFile,strtok(NULL," "),100), thisLine = 0;
      ifstream includeFile(thisFile);
      preProcessFile(includeFile);
      includeFile.close();
      strncpy(thisFile,thisFilePushed,100);
      thisLine = thisLinePushed;
   }

   else if (strcasecmp(key,"template") == 0) {
      char *name=strtok(NULL," ()");
      char *number=strtok(NULL," ()");
      TemplateItem *item = new TemplateItem(name,number);
      table.Append(item);
      fdest << "      " << "end" << endl << endl;
      fdest << "      " << "subroutine " << name << "Def" <<  endl;
      while (!fin.eof()) {
         fin.getline(line,250), ++thisLine;
         if ((strlen(line) > 0) && (line[0] != '!')) {
            char *type=strtok(line," ");
            char *name=strtok(NULL," ");
            if (strcasecmp(type,"end") == 0) break;
            TemplateAtom *atom = new TemplateAtom(type,name);
            item->Append(atom);
         }
      }
      TemplateAtom *a=item->First();
      while (a) {
         if (strcasecmp(a->fType,"int") == 0)
            fdest << "      " << "integer " << a->fName << endl;
         else if (strcasecmp(a->fType,"real") == 0)
            fdest << "      " << "real " << a->fName << endl;
         else if (strcasecmp(a->fType,"char") == 0)
            fdest << "      " << "character*40 " << a->fName << endl;
         else if (strcasecmp(a->fType,"material") == 0)
            fdest << "      " << "character*40 " << a->fName << endl;
         a = a->Next();
      }
      a = item->First();
      fdest << "      " << "common /" << item->fName << "/";
      int col=15+strlen(item->fName);
      int needsComma=0;
      while (a) {
         char vname[100];
         strncpy(vname,a->fName,100);
         strtok(vname," ()");
         if ((strcasecmp(a->fType,"int") == 0) ||
             (strcasecmp(a->fType,"real") == 0) ||
             (strcasecmp(a->fType,"char") == 0) ||
             (strcasecmp(a->fType,"material") == 0)) {
            if (needsComma++) fdest << ",";
            col += strlen(vname)+1;
            if (col > 71) {
               fdest << endl;
               fdest << "     + ";
               col = strlen(vname)+7;
            }
            fdest << vname;
         }
         a = a->Next();
      }
      fdest << endl;
   }

   else if (strcasecmp(key,"make") == 0) {
      char *name = strtok(NULL," ");
      TemplateItem *item=table.Find(name);
      if (item == NULL) {
         cerr << "gpp error: keyword " << name << " not defined!" << endl;
         cerr << "bombing out in line " << thisLine;
         cerr << " of file " << thisFile << endl;
         exit(1);
      }
      TemplateAtom *a=item->First();
      while (a != NULL) {
         char vname[100];
         strncpy(vname,a->fName,100);
         char *lasts[100];
         char *numb=strtok_r(vname," ()",lasts);
         numb = strtok_r(NULL," ()",lasts);
         int nvalu=1;
         if (numb) sscanf(numb,"%d",&nvalu);
         if ((strcasecmp(a->fType,"child") == 0) ||
             (strcasecmp(a->fType,"parent") == 0)) nvalu = 0;
         for (int i=1; i<=nvalu; i++) {
            char *valu=strtok(NULL," ");
            if (valu == NULL) {
               // cerr << "gpp warning: missing argument";
               // cerr << " in line " << thisLine;
               // cerr << " in file " << thisFile << endl;
               break;
            }
            fdest << "      " << vname;
            if (numb) fdest << "(" << i << ")";
            fdest << " =";
            if (valu[0] == '"') {
               while (valu[strlen(valu)-1] != '"') {
                  fdest << " " << valu;
                  valu = strtok(NULL," ");
                  if (valu == NULL) {
                     cerr << "gpp error: quoted argument unterminated!";
                     cerr << endl;
                     cerr << "bombing out in line " << thisLine;
                     cerr << " of file " << thisFile << endl;
                     exit(1);
                  }
               }
            }
            fdest << " " << valu << endl;
         }
         a = a->Next();
      }
      fdest << "      " << "call make" << name << endl;
   }

   else {
      cerr << "gpp error: keyword " << key << " not valid!" << endl;
      cerr << "bombing out in line " << thisFile << endl;
      cerr << " of file " << thisFile;
      exit(1);
   }
}

void pruneTemplates()
{
   fsource.open(sourceFile);
   if (!fsource) {
      cerr << "gpp error: unable to open input source file ";
      cerr << sourceFile << endl;
      exit(1);
   }
   char line[250];
   char token[250];
   TemplateItem *item = 0;
   while (!fsource.eof()) {
      fsource.getline(line,250);
      strcpy(token,line);
      strtok(token," ");
      if (strcasecmp(token,"      subroutine") == 0) {
         item = table.Find(strtok(NULL," ")+4);
      }
      else if (strcasecmp(token,"      character*40") == 0) {
         if (item) {
            item->Delete(strtok(NULL," "));
         }
      }
      else if (strcasecmp(token,"      real") == 0) {
         if (item) {
            item->Delete(strtok(NULL," "));
         }
      }
      else if (strcasecmp(token,"      integer") == 0) {
         if (item) {
            item->Delete(strtok(NULL," "));
         }
      }
      else if (strcasecmp(token,"      call") == 0) {
         char *maker = strtok(NULL," ");
         maker[strlen(maker)-3] = 0;
         item = table.Find(maker);
         if (item) {
            TemplateAtom *atom = item->First();
            while (atom != NULL) {
               if ((strcasecmp(atom->fType,"child") != 0) &&
                   (strcasecmp(atom->fType,"parent") != 0)) break;
               atom = atom->Next();
            }
            if (atom == 0) {
               table.Delete(item->fName);
            }
         }
      }
      else if (strcasecmp(line,"      end") == 0) {
         item = 0;
      }
   }
   fsource.close();
}

void postProcessFile()
{
   TemplateItem *item=table.First();
   while (item) {
      cerr << ">> Make function for object " << item->fName;
      cerr << " not found, dummy subroutine inserted." << endl;
      fdest << "      " << "end" << endl << endl;
      fdest << "      " << "subroutine make" << item->fName << endl;
      TemplateAtom *a=item->First();
      while (a) {
         if (strcasecmp(a->fType,"int") == 0)
            fdest << "      " << "integer " << a->fName << endl;
         else if (strcasecmp(a->fType,"real") == 0)
            fdest << "      " << "real " << a->fName << endl;
         else if (strcasecmp(a->fType,"char") == 0)
            fdest << "      " << "character*40 " << a->fName << endl;
         else if (strcasecmp(a->fType,"material") == 0)
            fdest << "      " << "character*40 " << a->fName << endl;
         a = a->Next();
      }
      a = item->First();
      fdest << "      " << "common /" << item->fName << "/";
      int col=15+strlen(item->fName);
      int needsComma=0;
      while (a) {
         char vname[100];
         strncpy(vname,a->fName,100);
         strtok(vname," ()");
         if ((strcasecmp(a->fType,"int") == 0) ||
             (strcasecmp(a->fType,"real") == 0) ||
             (strcasecmp(a->fType,"char") == 0) ||
             (strcasecmp(a->fType,"material") == 0)) {
            if (needsComma++) fdest << ",";
            col += strlen(vname)+1;
            if (col > 71) {
               fdest << endl;
               fdest << "     + ";
               col = strlen(vname)+7;
            }
            fdest << vname;
         }
         a = a->Next();
      }
      fdest << endl << endl;
      fdest << "CC---> add the appropriate GEANT calls here" << endl << endl;
      item = item->Next();
   }
   fdest << "      " << "end" << endl << endl;
   fdest << "      " << "subroutine makeGeometry" << endl;
   item=table.First();
   while (item) {
      fdest << "      " << "call " << item->fName << "Def" << endl;
      item = item->Next();
   }
   fdest << "      " << "end" << endl;
}

int main(int argc, char** argv)
{
   for (int arg=1; arg<argc; arg++) {
      strncpy(thisFile,argv[arg],100);
      fdbase.open(thisFile), thisLine=0;
      if (!fdbase) {
         cerr << "gpp error: unable to open input db file ";
         cerr << argv[arg] << endl;
         exit(1);
      }
      preProcessFile(fdbase);
      fdbase.close();
      pruneTemplates();
      postProcessFile();
      fdest.close();
   }
   cout << "gpp processed " << argc-1 << " files" << endl;
}
