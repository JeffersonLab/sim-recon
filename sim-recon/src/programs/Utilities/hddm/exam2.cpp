#include <fstream>
#include "hddm_x.hpp"

int report(hddm_x::HDDM &xrec);

int main()
{
   std::ofstream ofs("exam2.hddm");
   hddm_x::ostream ostr(ofs);
   hddm_x::HDDM xrec;

   // ostr.setCompression(hddm_x::k_z_compression);

   for (int n=0; n<1000000; ++n) {
      hddm_x::StudentList student = xrec.addStudents();
      student().name("Humphrey Gaston");
      hddm_x::EnrolledList enrolled = student().addEnrolleds();
      enrolled().year(2005);
      enrolled().semester(2);
      hddm_x::CourseList course = enrolled().addCourses(3);
      course(0).credits(3);
      course(0).title("Beginning Russian");
      course(0).addResults();
      course(0).result().grade("A-");
      course(0).result().pass(true);
      course(1).credits(1);
      course(1).title("Bohemian Poetry");
      course(1).addResults();
      course(1).result().grade("C");
      course(1).result().pass(1);
      course(2).credits(4);
      course(2).title("Developmental Psychology");
      course(2).addResults();
      course(2).result().grade("B+");
      course(2).result().pass(true);
      ostr << xrec;
      xrec.clear();
   }
   ofs.close();
   exit(0);

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
   hddm_x::CourseList course = xrec.courses();
   int total_courses = course.size();
   int total_enrolled = 0;
   int total_credits = 0;
   int total_passed = 0;
   hddm_x::CourseList::iterator iter;
   for (iter = course.begin(); iter != course.end(); ++iter) {
      if (iter->result().pass()) {
         if (iter->year() > 1992) {
            total_credits += iter->credits();
         }
         ++total_passed;
      }
   }
   std::cout << course().name() << " enrolled in " 
             << total_courses << " courses "
             << "and passed " << total_passed << " of them, " << std::endl
             << "earning a total of " << total_credits
             << " credits." << std::endl;
   return 0;
}
