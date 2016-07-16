#include <fstream>
#include "hddm_x.hpp"

int report(hddm_x::HDDM &xrec);

int main()
{
   hddm_x::HDDM xrec;
   std::ofstream ofs("exam2.hddm");
   hddm_x::ostream *ostr = new hddm_x::ostream(ofs);

   //ostr->setCompression(hddm_x::k_z_compression);

   for (int n=0; n<1000000; ++n) {
      hddm_x::StudentList student = xrec.addStudents();
      student().setName("Humphrey Gaston");
      hddm_x::EnrolledList enrolled = student().addEnrolleds();
      enrolled().setYear(2005);
      enrolled().setSemester(2);
      hddm_x::CourseList course = enrolled().addCourses(3);
      course(0).setCredits(3);
      course(0).setTitle("Beginning Russian");
      course(0).addResults();
      course(0).getResult().setGrade("A-");
      course(0).getResult().setPass(true);
      course(1).setCredits(1);
      course(1).setTitle("Bohemian Poetry");
      course(1).addResults();
      course(1).getResult().setGrade("C");
      course(1).getResult().setPass(1);
      course(2).setCredits(4);
      course(2).setTitle("Developmental Psychology");
      course(2).addResults();
      course(2).getResult().setGrade("B+");
      course(2).getResult().setPass(true);
      *ostr << xrec;
      xrec.clear();
   }
   delete ostr;
   ofs.close();

   // try reading it back in from the output file just written
   std::ifstream ifs("exam2.hddm");
   hddm_x::istream istr(ifs);
   int count=0;
   while (ifs.good()) {
      istr >> xrec;
      if (count/100000*100000 == count) {
         std::cout << "event " << count << std::endl;
         report(xrec);
      }
      ++count;
   }
   std::cout << "finished after " << count << " events read." << std::endl;
   return 0;
}

int report(hddm_x::HDDM &xrec)
{
   hddm_x::CourseList courses = xrec.getCourses();
   int total_courses = courses.size();
   int total_enrolled = 0;
   int total_credits = 0;
   int total_passed = 0;
   hddm_x::CourseList::iterator iter;
   for (iter = courses.begin(); iter != courses.end(); ++iter) {
      if (iter->getResult().getPass()) {
         if (iter->getYear() > 1992) {
            total_credits += iter->getCredits();
         }
         ++total_passed;
      }
   }
   std::cout << courses().getName() << " enrolled in " 
             << total_courses << " courses "
             << "and passed " << total_passed << " of them, " << std::endl
             << "earning a total of " << total_credits
             << " credits." << std::endl;
   return 0;
}
